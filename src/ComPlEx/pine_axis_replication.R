#!/usr/bin/env Rscript
# pine_axis_replication.R — reciprocal Pinus sylvestris replication of the
# connectivity-constraint-conservation axis (degree-dN/dS + integrative_conservation_model
# Model A/B), computed from pine network degree and pine conservation breadth. Population-genetic
# / PAV / promoter-TE terms are spruce-only and are omitted here.
# CAVEAT (Model A SD term): no pine shared/pine-specific SD *class* table exists — only a pine SD
# membership set (Orthogroups/pine_all_in_duplicated_blocks_PGset.tsv.gz) — so Model A uses a single
# binary pine_SD flag rather than the spruce shared_SD/spruce_SD pair. Dropping it leaves the
# degree + expression effects unchanged in direction. Only Model B (no SD term) is cited in the
# manuscript, so no reported number depends on this choice.
# Mirrors: build_network_degree.py, degree_dnds_correlation.R, integrative_conservation_model.R.
# Outputs (results/integration/):
#   network_degree_PINE.tsv          per-PS-gene mean ComPlEx network degree
#   degree_dnds_correlation_PINE.tsv Spearman rho/P + dN/dS by degree quartile
#   pine_axis_models.tsv             Model A/B ordinal odds ratios + P (+ range-restriction rows)
USERLIB <- file.path(Sys.getenv("HOME"), "Library/R/arm64/4.6/library"); if (dir.exists(USERLIB)) .libPaths(c(USERLIB, .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(MASS) })

args <- commandArgs(trailingOnly = FALSE)
sp0  <- sub("--file=", "", args[grep("--file=", args)])
PROJ <- normalizePath(file.path(if (length(sp0) > 0) dirname(normalizePath(sp0)) else getwd(), "../.."))
I    <- file.path(PROJ, "results/integration")
NETS <- c("cold_needle","cold_root","drought_needle","drought_root")

# (1) pine network degree = mean Degree across the pine networks a gene is present in.
# Source: the ComPlEx Python-port per-network centrality (results/integration/centrality/
# centrality_pine_<net>.tsv), the SAME basis build_network_degree.py uses for spruce
# (centrality_spruce_<net>.tsv) — NOT the raw R-port results/ComPlEx/.../centrality_pine.tsv,
# so pine and spruce degree are computed on identical methodology.
frames <- lapply(NETS, function(n){
  d <- fread(file.path(PROJ, sprintf("results/integration/centrality/centrality_pine_%s.tsv", n)))
  setnames(d, 1:2, c("Genes","Degree")); d[, .(Genes, Degree)] })
deg <- Reduce(function(a,b) merge(a,b,by="Genes",all=TRUE),
              lapply(seq_along(frames), function(i){ x<-copy(frames[[i]]); setnames(x,"Degree",paste0("d",i)); x }))
deg[, mean_degree := rowMeans(.SD, na.rm=TRUE), .SDcols=paste0("d",1:4)]
fwrite(deg[, .(Genes, mean_degree)], file.path(I, "network_degree_PINE.tsv"), sep="\t")
deg <- deg[, .(ps_gene=Genes, degree=mean_degree)]

# (2) pine conservation breadth = sum of 4 *_present, keyed on Species2 (PS gene), max per gene
wp <- fread(file.path(PROJ, "results/ComPlEx/RData/weighted_gene_pairs.tsv"))
wp[, breadth := cold_needle_present + cold_root_present + drought_needle_present + drought_root_present]
d <- wp[, .(breadth = max(breadth)), by = .(ps_gene = Species2)]

# (3) pine mean expression over PC + PD
pc <- fread(file.path(PROJ,"data/expression/PC_expression.txt")); setnames(pc,1,"ps_gene")
pd <- fread(file.path(PROJ,"data/expression/PD_expression.txt")); setnames(pd,1,"ps_gene")
ex <- merge(pc[, .(ps_gene, e1=rowMeans(.SD)), .SDcols=setdiff(names(pc),"ps_gene")],
            pd[, .(ps_gene, e2=rowMeans(.SD)), .SDcols=setdiff(names(pd),"ps_gene")], by="ps_gene", all=TRUE)
ex[, expr := rowMeans(cbind(e1,e2), na.rm=TRUE)]

psd <- fread(cmd=sprintf("gunzip -c %s", file.path(PROJ,"Orthogroups/pine_all_in_duplicated_blocks_PGset.tsv.gz")), header=FALSE)[[1]]
psd <- sub("\\.mRNA\\.\\d+$","", psd)
yn  <- unique(fread(file.path(I,"cross_species_dnds_yn00.tsv"))[dS>0 & dS<5 & dNdS<10, .(ps_gene, dNdS)], by="ps_gene")

d <- merge(d, deg, by="ps_gene", all.x=TRUE)
d <- merge(d, ex[, .(ps_gene, expr)], by="ps_gene", all.x=TRUE)
d <- merge(d, yn, by="ps_gene", all.x=TRUE)
d[, pine_SD := as.integer(ps_gene %in% psd)]
d <- d[!is.na(degree) & !is.na(expr)]
z <- function(x) as.numeric(scale(x))
d[, `:=`(zdeg=z(degree), zexpr=z(expr), bf=factor(breadth, levels=0:4, ordered=TRUE))]

orl <- function(po, preds, label){ ct<-coef(summary(po)); co<-ct[seq_along(preds),1]; se<-ct[seq_along(preds),2]
  data.table(model=label, term=rownames(ct)[seq_along(preds)], OR=round(exp(co),3),
             ci_lo=round(exp(co-1.96*se),3), ci_hi=round(exp(co+1.96*se),3),
             p=signif(pnorm(abs(ct[seq_along(preds),3]),lower.tail=FALSE)*2,3)) }

# (4a) pine degree vs dN/dS (dNdS>0), Spearman + quartile medians
m <- d[!is.na(dNdS) & dNdS>0]
sp <- suppressWarnings(cor.test(m$degree, m$dNdS, method="spearman"))
m[, q := cut(degree, quantile(degree,0:4/4), include.lowest=TRUE, labels=paste0("Q",1:4))]
qmed <- m[, .(n=.N, median_dNdS=round(median(dNdS),4)), by=q][order(q)]
fwrite(rbind(data.table(statistic=c("n_usable","spearman_rho","spearman_P"),
                        value=c(nrow(m), round(unname(sp$estimate),4), signif(sp$p.value,3))),
             data.table(statistic=paste0(qmed$q,"_median_dNdS"), value=qmed$median_dNdS)),
       file.path(I,"degree_dnds_correlation_PINE.tsv"), sep="\t")
cat(sprintf("(4a) PINE degree vs dN/dS: n=%d rho=%.3f P=%.2e ; quartiles %s\n",
            nrow(m), unname(sp$estimate), sp$p.value, paste(qmed$median_dNdS, collapse="/")))

# (4b) Model A ; (4c) Model B
A <- polr(bf ~ zdeg + zexpr + pine_SD, data=d, Hess=TRUE)
dB <- d[!is.na(dNdS)]; dB[, zdnds := z(dNdS)]
B <- polr(bf ~ zdnds + zdeg + zexpr, data=dB, Hess=TRUE)
models <- rbind(orl(A, c("zdeg","zexpr","pine_SD"), "A_pine_genome_wide"),
                orl(B, c("zdnds","zdeg","zexpr"),   "B_pine_1to1_subset"))
# range-restriction rows
rr <- data.table(model="range_restriction",
  term=c("expr_sd_full","expr_sd_1to1","degree_sd_full","degree_sd_1to1","expr_sd_ratio","degree_sd_ratio"),
  OR=c(round(sd(d$expr),3),round(sd(dB$expr),3),round(sd(d$degree),3),round(sd(dB$degree),3),
       round(sd(dB$expr)/sd(d$expr),3), round(sd(dB$degree)/sd(d$degree),3)), ci_lo=NA, ci_hi=NA, p=NA)
fwrite(rbind(models, rr), file.path(I,"pine_axis_models.tsv"), sep="\t")
cat(sprintf("(4b) Model A pine n=%d ; (4c) Model B pine 1:1 n=%d\n", nrow(d), nrow(dB)))
print(models, row.names=FALSE)
cat("Wrote network_degree_PINE.tsv, degree_dnds_correlation_PINE.tsv, pine_axis_models.tsv\n")
