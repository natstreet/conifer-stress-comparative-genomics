suppressPackageStartupMessages({library(data.table); library(ggplot2); library(dplyr)})

PROJ   <- here::here()
INTEG  <- file.path(PROJ, "results/integration")
FIGDIR <- file.path(INTEG, "figures")
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)  # ensure output dir exists in a fresh checkout
# Population-genomic selection and envGWAS gene sets (Kalman et al.), vendored into
# the deposit; see SOURCES.tsv for origin.
POPGEN <- file.path(PROJ, "data/popgen")
THEME  <- theme_bw(base_size=10) +
  theme(panel.grid.minor=element_blank(), strip.background=element_rect(fill="grey92"))

cat("=== Stress co-expression x Population selection signal ===\n")

# ── Per-gene co-expression category (5-category backbone) ────────────────────
backbone <- fread(file.path(INTEG, "integration_backbone_1to1.tsv"),
                  select=c("pa_gene","coex_category"))
backbone <- backbone[!is.na(coex_category)]
# Best category per gene
CAT_ORD <- c("conserved","cold_specific","drought_specific","multi_tissue","not_coex")
gene_cat <- backbone[, .(coex_category=coex_category[which.min(
              match(coex_category, CAT_ORD))]), by=pa_gene]
UNIV <- nrow(gene_cat)
cat(sprintf("Gene universe: %d PA genes with co-expression categories\n", UNIV))
print(gene_cat[, .N, by=coex_category][order(match(coex_category, CAT_ORD))])

# ── Load all selection gene sets ──────────────────────────────────────────────

# XP-EHH / iHS gene sets (8995 genes with signal type)
xpehh_ihs <- fread(file.path(POPGEN,"selection_xpehh_ihs_genes.tsv"),
                   header=FALSE, col.names=c("pa_gene","signal_type"))
cat(sprintf("\n8995 XP-EHH/iHS genes loaded\n"))
print(xpehh_ihs[, .N, by=signal_type])

# Separate into contrasting (xpehh) and shared (ihs) sweep genes
xpehh_genes <- xpehh_ihs[grepl("xpehh", signal_type), pa_gene]  # 4202 + shared
ihs_genes    <- xpehh_ihs[grepl("ihs",   signal_type), pa_gene]  # 7769 + shared  
xpehh_only   <- xpehh_ihs[signal_type=="xpehh", pa_gene]         # contrasting only
ihs_only     <- xpehh_ihs[signal_type=="ihs",   pa_gene]         # shared sweep only
both_xpehh_ihs <- xpehh_ihs[signal_type=="xpehh,ihs", pa_gene]  # both signals
cat(sprintf("  xpehh (any): %d, ihs (any): %d, both: %d\n",
            length(xpehh_genes), length(ihs_genes), length(both_xpehh_ihs)))

# Note: the earlier PCAdapt combined-score gene set (679) and its PCAdapt∩XP-EHH/iHS intersection
# (319) are intentionally not used here. The final Kalman et al. deposit ships a revised per-PC
# PCAdapt analysis (PC1/PC2) rather than that combined-score set, and no reported number in this
# study relied on the PCAdapt sets. Selection enrichment therefore uses the XP-EHH/iHS scans above
# (deposited as PopGen/gene-id_selected-genes_set-8995.txt) plus the envGWAS sets below.

# envGWAS — per climate variable
gwas_files <- list.files(POPGEN, pattern="_genes\\.uniq\\.txt$", full.names=TRUE)
gwas_by_var <- lapply(gwas_files, function(f) fread(f, header=FALSE)$V1)
names(gwas_by_var) <- sub("_genes\\.uniq\\.txt$","", basename(gwas_files))
cat(sprintf("\nenvGWAS gene counts per variable:\n"))
sapply(gwas_by_var, length) |> print()

# ── Enrichment function ───────────────────────────────────────────────────────
fisher_enrich <- function(sig_genes, label) {
  rbindlist(lapply(CAT_ORD, function(cat) {
    a <- sum(gene_cat$coex_category==cat &  gene_cat$pa_gene %in% sig_genes)
    b <- sum(gene_cat$coex_category==cat & !gene_cat$pa_gene %in% sig_genes)
    c <- sum(gene_cat$coex_category!=cat &  gene_cat$pa_gene %in% sig_genes)
    d <- UNIV - sum(gene_cat$coex_category==cat) - c
    ft <- fisher.test(matrix(c(a,b,c,d),2))
    data.table(signal=label, category=cat, n_sig=a, n_cat=a+b,
               pct=round(100*a/(a+b),1),
               OR=round(ft$estimate,3), pvalue=ft$p.value,
               ci_lo=round(ft$conf.int[1],3), ci_hi=round(ft$conf.int[2],3))
  }))
}

# ── Run enrichment for all signal types ──────────────────────────────────────
all_res <- rbindlist(list(
  fisher_enrich(xpehh_only,    "XP-EHH\n(contrasting N/S)"),
  fisher_enrich(ihs_only,      "iHS\n(shared sweep)"),
  fisher_enrich(both_xpehh_ihs,"XP-EHH+iHS\n(both signals)")
))
all_res[, padj := p.adjust(pvalue, "BH"), by=signal]
all_res[, sig  := ifelse(padj<0.001,"***",ifelse(padj<0.01,"**",
                  ifelse(padj<0.05,"*","ns")))]

cat("\n\nSignificant enrichments (padj<0.05):\n")
print(all_res[padj<0.05, .(signal, category, n_sig, n_cat, pct, OR, pvalue, padj, sig)])
fwrite(all_res, file.path(INTEG,"selection_category_enrichment.tsv"), sep="\t")

# ── envGWAS by variable ───────────────────────────────────────────────────────
cat("\n\nenvGWAS enrichment by climate variable:\n")
gwas_res <- rbindlist(lapply(names(gwas_by_var), function(nm) {
  fisher_enrich(gwas_by_var[[nm]], nm)
}))
gwas_res[, padj := p.adjust(pvalue,"BH"), by=signal]
gwas_res[, sig  := ifelse(padj<0.001,"***",ifelse(padj<0.01,"**",
                   ifelse(padj<0.05,"*","ns")))]
print(gwas_res[padj<0.05, .(signal, category, n_sig, n_cat, OR, padj, sig)])
fwrite(gwas_res, file.path(INTEG,"gwas_category_enrichment.tsv"), sep="\t")

# ── Figure: heatmap of OR by signal × co-expression category ─────────────────
fig_df <- all_res %>%
  mutate(category = factor(category, CAT_ORD),
         log2OR   = log2(pmax(OR, 0.01)),
         sig_text = ifelse(padj<0.001,"***",ifelse(padj<0.01,"**",
                    ifelse(padj<0.05,"*",""))))

p1 <- ggplot(fig_df, aes(x=category, y=signal, fill=log2OR)) +
  geom_tile(colour="white", linewidth=0.4) +
  geom_text(aes(label=sig_text), size=3.5, colour="grey20") +
  scale_fill_gradient2(low="#4393C3", mid="white", high="#D6604D",
                       midpoint=0, name="log2 OR") +
  scale_x_discrete(labels=c(conserved="Conserved",cold_specific="Cold",
                             drought_specific="Drought",multi_tissue="Multi",
                             not_coex="Not\nco-expr.")) +
  labs(x=NULL, y=NULL,
       title="Population selection signals x stress co-expression categories",
       subtitle="Fisher's exact OR (BH-adjusted); * p<0.05 ** p<0.01 *** p<0.001") +
  THEME + theme(axis.text.x=element_text(angle=30, hjust=1))
ggsave(file.path(FIGDIR,"selection_category_heatmap.pdf"), p1,
       width=14, height=10, units="cm", device="pdf")

# envGWAS heatmap — label by climate variable meaning
gwas_fig <- gwas_res %>%
  mutate(category=factor(category, CAT_ORD),
         log2OR=log2(pmax(OR,0.01)),
         sig_text=ifelse(padj<0.001,"***",ifelse(padj<0.01,"**",
                  ifelse(padj<0.05,"*",""))),
         var_lab=recode(signal,
           bio3="Isothermality (bio3)",
           bio7="Annual temp. range (bio7)",
           bio10="Warmest qtr. temp (bio10)",
           bio11="Coldest qtr. temp (bio11)",
           bio18="Warmest qtr. precip (bio18)",
           bio19="Coldest qtr. precip (bio19)"))

p2 <- ggplot(gwas_fig, aes(x=category, y=var_lab, fill=log2OR)) +
  geom_tile(colour="white", linewidth=0.4) +
  geom_text(aes(label=sig_text), size=3.5, colour="grey20") +
  scale_fill_gradient2(low="#4393C3", mid="white", high="#D6604D",
                       midpoint=0, name="log2 OR") +
  scale_x_discrete(labels=c(conserved="Conserved",cold_specific="Cold",
                             drought_specific="Drought",multi_tissue="Multi",
                             not_coex="Not\nco-expr.")) +
  labs(x=NULL, y=NULL,
       title="envGWAS climate-associated genes x co-expression categories",
       subtitle="Fisher's exact OR (BH-adjusted)") +
  THEME + theme(axis.text.x=element_text(angle=30, hjust=1))
ggsave(file.path(FIGDIR,"gwas_category_heatmap.pdf"), p2,
       width=16, height=10, units="cm", device="pdf")

# ── Population-selection periphery: co-expressolog degree + expression level ──────
# Universe: genes with a within-species network position (network_degree.tsv). A network node with no
# conserved cross-species co-expressolog is a genuine degree-0 (periphery) member and is kept, 0-filled;
# genes with no network position are excluded. Background = no-selection genes within this universe.
suppressPackageStartupMessages(library(DESeq2))
sel_class_of <- function(g) fifelse(g %in% both_xpehh_ihs, "XP-EHH+iHS",
                          fifelse(g %in% ihs_only,   "iHS only",
                          fifelse(g %in% xpehh_only, "XP-EHH only", "background")))
CLASS_ORD <- c("background","iHS only","XP-EHH only","XP-EHH+iHS")

# co-expressolog partner-count (conserved cross-species co-expressolog pairs per PA gene), 0-filled
# over the within-species network-node universe.
wp <- fread(file.path(PROJ,"results/ComPlEx/RData/weighted_gene_pairs.tsv"))
pc <- wp[, .(deg=.N), by=.(gene=Species1)]
nd <- fread(file.path(INTEG,"network_degree.tsv")); setnames(nd, 1, "gene")
deg_u <- merge(data.table(gene=unique(nd$gene)), pc, by="gene", all.x=TRUE)
deg_u[is.na(deg), deg:=0]
deg_u[, sel_class := sel_class_of(gene)]
deg_tab <- deg_u[, .(n=.N, median=as.numeric(median(deg)), mean=round(mean(deg),1),
                     pct_zero=round(100*mean(deg==0),1)), by=sel_class][order(match(sel_class,CLASS_ORD))]
fwrite(deg_tab, file.path(INTEG,"selection_network_degree.tsv"), sep="\t")
cat("\n=== Selection periphery: co-expressolog degree (network-node universe, 0-filled) ===\n"); print(deg_tab)

# expression level: per-gene MEAN DESeq2 baseMean across the FOUR Norway spruce stress experiments
# (SCN,SCR,SDN,SDR = cold/drought x needle/root), defined analogously to co-expressolog degree (the
# mean across the four spruce stress networks) so text = code = number for both metrics in the same
# paragraph. Restricted to the SAME within-species network-node universe as degree.
# Missing-dds rule: a gene's value is the mean over the spruce conditions where it is expressed
# (mean over AVAILABLE conditions, not requiring all four); a network-node gene absent from every
# spruce dds has no expression value and is excluded from the expression comparison only.
SPRUCE_DDS <- c("SCN","SCR","SDN","SDR")
bm_long <- rbindlist(lapply(SPRUCE_DDS, function(cond){
  e <- new.env()
  load(file.path(PROJ, sprintf("data/dds/dds_%s.rda", cond)), envir=e)
  d <- get(grep("^dds", ls(e), value=TRUE)[1], envir=e)
  d <- estimateSizeFactors(d)
  data.table(gene=rownames(d), bm=rowMeans(counts(d, normalized=TRUE)), cond=cond)
}))
bm_gene <- bm_long[, .(bm=mean(bm), n_cond=.N), by=gene]          # per-gene mean over available spruce conditions
bm_u <- merge(data.table(gene=unique(nd$gene)), bm_gene, by="gene")  # network-node universe with an expression value
bm_u[, sel_class := sel_class_of(gene)]
bm_tab <- bm_u[, .(n=.N, median_bm=as.numeric(median(bm)), mean_bm=mean(bm),
                   pct_lowexp=round(100*mean(bm<10),1)), by=sel_class][order(match(sel_class,CLASS_ORD))]
fwrite(bm_tab, file.path(INTEG,"selection_expression_level.tsv"), sep="\t")
cat("\n=== Selection periphery: mean baseMean across spruce SCN/SCR/SDN/SDR (network-node universe) ===\n")
cat(sprintf("   genes with an expression value: %d of %d network nodes; n_cond distribution: %s\n",
            nrow(bm_u), length(unique(nd$gene)),
            paste(sprintf("%d->%d", sort(unique(bm_gene$n_cond)),
                          tabulate(factor(bm_gene$n_cond, levels=sort(unique(bm_gene$n_cond))))), collapse=" ")))
print(bm_tab)

# Wilcoxon rank-sum, both-signals vs no-selection background (committed, not inline-only)
wd <- suppressWarnings(wilcox.test(deg_u[sel_class=="XP-EHH+iHS",deg], deg_u[sel_class=="background",deg]))
wb <- suppressWarnings(wilcox.test(bm_u[sel_class=="XP-EHH+iHS",bm],   bm_u[sel_class=="background",bm]))
wilcox_tab <- data.table(metric=c("coexpressolog_degree","baseMean_expression"),
  comparison="XP-EHH+iHS vs no-selection background",
  W=c(wd$statistic, wb$statistic), p_value=c(wd$p.value, wb$p.value))
fwrite(wilcox_tab, file.path(INTEG,"selection_periphery_wilcox.tsv"), sep="\t")
cat("\n=== Selection periphery: Wilcoxon (both vs background) ===\n"); print(wilcox_tab)

cat("\n=== Done ===\n")
