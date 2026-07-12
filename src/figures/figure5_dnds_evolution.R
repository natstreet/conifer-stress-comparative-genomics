suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork); library(ggtext); library(data.table)
})

# Figure 5 — evolutionary and genomic features of co-expression conservation.
# Panel a tests the core claim (conserved co-expression -> stronger purifying
# selection) CONTINUOUSLY, against conservation breadth (the number of the four
# stress-tissue comparisons in which an orthogroup has a significant conserved
# co-expressolog, 0-4) rather than the discrete categories, and reports a partial
# correlation controlling for expression level (co-expression significance is
# expression-dependent). dN/dS is the YN00 (Yang & Nielsen 2000) estimator.
# Panel b (SD enrichment by category) is unchanged.

INTEG <- "results/integration"
RDAT  <- "results/ComPlEx/RData"
OUT   <- "results/final_figures"

# ── PANEL a: dN/dS vs conservation breadth (threshold-free) ──────────────────
yn <- as.data.table(read.table(file.path(INTEG,"cross_species_dnds_yn00.tsv"), sep="\t", header=TRUE))
wp <- fread(file.path(RDAT,"weighted_gene_pairs.tsv"))
wp[, breadth := cold_needle_present + cold_root_present + drought_needle_present + drought_root_present]
gb  <- wp[, .(breadth = max(breadth, na.rm=TRUE)), by=.(pa_gene=Species1)]   # gene's best-conserved co-expressolog
deg <- fread(file.path(INTEG,"network_degree.tsv")); setnames(deg, c("pa_gene","degree"))
sc  <- fread("data/expression/SC_expression.txt"); sd_ <- fread("data/expression/SD_expression.txt")
setnames(sc,1,"pa_gene"); setnames(sd_,1,"pa_gene")
expr <- merge(sc[, .(pa_gene, e1=rowMeans(.SD)), .SDcols=setdiff(names(sc),"pa_gene")],
              sd_[, .(pa_gene, e2=rowMeans(.SD)), .SDcols=setdiff(names(sd_),"pa_gene")], by="pa_gene", all=TRUE)
expr[, mean_expr := rowMeans(cbind(e1,e2), na.rm=TRUE)]

d <- merge(yn[, .(pa_gene, ps_gene, dN, dS, dNdS)], gb, by="pa_gene", all.x=TRUE)
d[is.na(breadth), breadth := 0L]
d <- merge(d, deg[, .(pa_gene, degree)], by="pa_gene", all.x=TRUE)
d <- merge(d, expr[, .(pa_gene, mean_expr)], by="pa_gene", all.x=TRUE)
d <- d[dS>0 & dS<5 & dNdS<10 & !is.na(dNdS)]
n_usable <- nrow(d)

# statistics: Spearman + expression-controlled partial (rank residuals)
sp <- cor.test(d$breadth, d$dNdS, method="spearman", exact=FALSE)
de <- d[!is.na(mean_expr)]
pr <- cor.test(resid(lm(rank(dNdS) ~ rank(mean_expr), de)),
               resid(lm(rank(breadth) ~ rank(mean_expr), de)), method="pearson")
write.table(data.frame(n=n_usable, spearman_rho=unname(sp$estimate), spearman_p=sp$p.value,
                       partial_rho_expr=unname(pr$estimate), partial_p=pr$p.value),
            file.path(INTEG,"dnds_conservation_breadth_stats.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

dsum <- d[, .(med=median(dNdS), q25=quantile(dNdS,.25), q75=quantile(dNdS,.75), n=.N), by=breadth][order(breadth)]
gsum <- d[!is.na(degree), .(deg=median(degree), dq25=quantile(degree,.25), dq75=quantile(degree,.75)), by=breadth][order(breadth)]
write.table(dsum, file.path(INTEG,"dnds_by_conservation_breadth.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
cat(sprintf("Panel a: n=%d, Spearman rho=%.3f (p=%.2g), partial(expr) rho=%.3f (p=%.2g)\n",
            n_usable, sp$estimate, sp$p.value, pr$estimate, pr$p.value)); print(dsum)

yl <- c(0.20, 0.55)
d_lo <- floor(min(gsum$dq25)/10)*10; d_hi <- ceiling(max(gsum$dq75)/10)*10   # fit the IQR ribbon
sc_fn <- function(v) yl[1] + (v - d_lo)/(d_hi - d_lo) * (yl[2]-yl[1])
ORANGE <- "#C05621"; BLUE <- "#2B6CB0"
# Derive the P-value annotation from the computed value.
# If the Spearman p underflows to 0, floor at the smallest representable double.
p_exp  <- if (sp$p.value <= 0) ceiling(log10(.Machine$double.xmin)) else floor(log10(sp$p.value))
lab1   <- sprintf("paste('Spearman ', rho, ' = %.2f, ', italic(P), ' < 10'^%d)", sp$estimate, p_exp)
lab2 <- sprintf("paste('partial ', rho, ' (expression-controlled) = %.2f')", pr$estimate)
pa <- ggplot() +
  geom_ribbon(data=gsum, aes(x=breadth, ymin=sc_fn(dq25), ymax=sc_fn(dq75)), fill=BLUE, alpha=0.12) +
  geom_line(data=gsum, aes(x=breadth, y=sc_fn(deg)), colour=BLUE, linewidth=0.9, linetype="22") +
  geom_point(data=gsum, aes(x=breadth, y=sc_fn(deg)), colour=BLUE, size=2.6, shape=15) +
  geom_linerange(data=dsum, aes(x=breadth, ymin=q25, ymax=q75), colour=ORANGE, linewidth=1.2, alpha=0.5) +
  geom_line(data=dsum, aes(x=breadth, y=med), colour=ORANGE, linewidth=0.7, linetype="22") +
  geom_point(data=dsum, aes(x=breadth, y=med), colour=ORANGE, size=3.3) +
  annotate("text", x=0, y=yl[2]-0.004, hjust=0, vjust=1, size=2.9, colour="grey25", parse=TRUE, label=lab1) +
  annotate("text", x=0, y=yl[2]-0.045, hjust=0, vjust=1, size=2.9, colour="grey25", parse=TRUE, label=lab2) +
  scale_x_continuous(breaks=0:4,
                     labels=c("0\n(not\nco-expr.)","1","2","3","4\n(all\ncomparisons)"),
                     name="Conservation breadth (conserved co-expressolog comparisons)",
                     expand=expansion(add=0.4)) +
  scale_y_continuous(name=expression(italic(d)[N]/italic(d)[S]~"(median, IQR)"), limits=yl,
                     breaks=c(0.20,0.30,0.40,0.50),
                     sec.axis=sec_axis(~ d_lo + (. - yl[1])/(yl[2]-yl[1])*(d_hi-d_lo),
                                       name="co-expression network degree")) +
  labs(x=NULL, tag="a",
       subtitle="dN/dS (orange) falls and network degree (blue) rises with conservation breadth") +
  theme_classic(base_size=10) +
  theme(panel.grid.major.y=element_line(colour="grey92", linewidth=0.3), panel.grid.minor=element_blank(),
        axis.title.y.left=element_text(colour=ORANGE), axis.text.y.left=element_text(colour=ORANGE),
        axis.title.y.right=element_text(colour=BLUE), axis.text.y.right=element_text(colour=BLUE),
        axis.text.x=element_text(size=8, colour="grey20"), axis.line=element_line(colour="grey40"),
        plot.subtitle=element_text(size=7.3, colour="grey50"), plot.tag=element_text(size=12, face="bold"))

# ── PANEL b: SD enrichment (unchanged design) ────────────────────────────────
CAT_ORD <- c("conserved","multi_tissue","cold_specific","drought_specific","not_coex")
CAT_LAB <- c("Conserved","Multi-\ntissue","Cold-\nspecific","Drought-\nspecific","Not\nco-expressed")
sd <- read.table(file.path(INTEG,"sd_category_fisher.tsv"), sep="\t", header=TRUE) |>
  filter(sd_class %in% c("shared_SD","spruce_only_SD")) |>
  mutate(log2OR = log2(pmax(OR,0.01)),
         sig = ifelse(padj<0.001,"***", ifelse(padj<0.01,"**", ifelse(padj<0.05,"*",""))),
         sd_lab = factor(sd_class, c("shared_SD","spruce_only_SD"), c("Shared SD","*P. abies*-only SD")),
         cat_lab = factor(coex_category, CAT_ORD, CAT_LAB),
         enrich = log2OR > 0,
         sig_y = ifelse(log2OR >= 0, log2OR + 0.18, 0.18))

pb <- ggplot(sd, aes(cat_lab, log2OR)) +
  geom_col(aes(fill=enrich), width=0.68, colour="white", linewidth=0.3, alpha=0.88) +
  geom_hline(yintercept=0, colour="grey30", linewidth=0.6) +
  geom_text(aes(y=sig_y, label=sig), size=3.8, hjust=0.5, vjust=0, colour="grey15") +
  scale_fill_manual(values=c("TRUE"="#D55E00","FALSE"="#0072B2"), guide="none") +
  scale_y_continuous(name=expression(log[2]~"(odds ratio vs non-SD genes)"),
                     breaks=seq(-3,3,1), limits=c(-3.5,3.8)) +
  facet_wrap(~sd_lab, nrow=1, labeller=label_value) +
  labs(x=NULL, tag="b",
       subtitle="Orange bars: enriched (OR > 1)  |  Blue bars: depleted (OR < 1)  |  All bars: padj < 0.05",
       caption="Fisher exact test, Benjamini-Hochberg correction") +
  theme_classic(base_size=10) +
  theme(strip.text=element_markdown(size=9.5,face="bold"),
        strip.background=element_rect(fill="white",colour=NA),
        panel.grid.major.y=element_line(colour="grey92",linewidth=0.3),
        axis.text.x=element_text(size=8.5,angle=30,hjust=1),
        axis.line=element_line(colour="grey40"),
        plot.tag=element_text(size=12,face="bold"),
        plot.caption=element_text(size=7.5,colour="grey50",hjust=0),
        plot.subtitle=element_text(size=8,colour="grey40",hjust=0),
        plot.margin=margin(5,5,5,5))

fig5 <- pa / pb + plot_layout(heights=c(1, 1.15))
pdf(file.path(OUT,"Figure5.pdf"), width=18/2.54, height=18/2.54); print(fig5); dev.off()
ggsave(file.path(OUT,"Figure5.png"), fig5, width=18/2.54, height=18/2.54, dpi=300)
cat("Figure5.pdf / Figure5.png saved (breadth-based panel a)\n")
