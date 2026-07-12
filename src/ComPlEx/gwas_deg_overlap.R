suppressPackageStartupMessages({library(data.table); library(DESeq2); library(matrixStats)})
PROJ  <- here::here()
INTEG <- file.path(PROJ, "results/integration")
# envGWAS climate-associated gene lists (Kalman et al.), vendored to the deposit.
POPGEN <- file.path(PROJ, "data/popgen")

# Load all GWAS gene lists
gwas_vars <- c("bio3","bio7","bio10","bio11","bio18","bio19")
gwas      <- lapply(gwas_vars, function(v)
  fread(file.path(POPGEN,sprintf("%s_genes.uniq.txt",v)),header=FALSE)$V1)
names(gwas) <- gwas_vars
all_gwas <- unique(unlist(gwas))

# Climate groupings
cold_gwas    <- unique(c(gwas$bio10, gwas$bio11))
drought_gwas <- unique(c(gwas$bio18, gwas$bio19))
temp_gwas    <- unique(c(gwas$bio3,  gwas$bio7,  gwas$bio10, gwas$bio11))
cat(sprintf("GWAS: cold=%d, drought=%d, all-temp=%d, all=%d\n",
            length(cold_gwas), length(drought_gwas), length(temp_gwas), length(all_gwas)))

# DESeq2 results — collect all 8 PA stress contrasts
load_results <- function(dds_file, coef_name) {
  e <- new.env()
  load(dds_file, envir=e)
  dds <- DESeq(e[[ls(envir=e)[1]]], quiet=TRUE)
  res <- as.data.frame(results(dds, name=coef_name,
           filter=rowMedians(counts(dds,normalized=TRUE))))
  res$gene <- rownames(res)
  res
}
cat("Loading DESeq2 results...\n")
DATA_DIR <- file.path(PROJ,"data/dds")
deg_res <- list(
  SCN = load_results(file.path(DATA_DIR,"dds_SCN.rda"),
                     "Condition_neg5C_10d_vs_20C_0h"),
  SCR = load_results(file.path(DATA_DIR,"dds_SCR.rda"),
                     "Condition_neg5C_10d_vs_20C_0h"),
  SDN = load_results(file.path(DATA_DIR,"dds_SDN.rda"),
                     "Level_Collapse_vs_80."),
  SDR = load_results(file.path(DATA_DIR,"dds_SDR.rda"),
                     "Level_Collapse_vs_80.")
)
cat("  Done.\n")

# For each contrast, get DE gene set (padj<0.01, |lfc|>=2)
degs <- lapply(deg_res, function(r)
  r$gene[!is.na(r$padj) & r$padj<0.01 & abs(r$log2FoldChange)>=2])
sapply(degs, length) |> print()

# Overlap GWAS x DEGs: does each GWAS variable hit DE genes in relevant stress?
cat("\n=== GWAS gene x DEG overlap (any direction) ===\n")
combos <- list(
  list(gwas_name="bio10+bio11 (cold temp)",  gwas=cold_gwas,    stress=c("SCN","SCR")),
  list(gwas_name="bio18+bio19 (precip)",     gwas=drought_gwas, stress=c("SDN","SDR")),
  list(gwas_name="bio3+bio7 (temp range)",   gwas=temp_gwas,    stress=c("SCN","SCR")),
  list(gwas_name="all GWAS",                 gwas=all_gwas,     stress=c("SCN","SCR","SDN","SDR"))
)
rows <- rbindlist(lapply(combos, function(x) {
  rbindlist(lapply(x$stress, function(s) {
    ov <- intersect(x$gwas, degs[[s]])
    data.table(gwas_set=x$gwas_name, stress=s,
               n_gwas=length(x$gwas), n_degs=length(degs[[s]]),
               n_overlap=length(ov), genes=paste(ov,collapse=","))
  }))
}))
print(rows[, .(gwas_set, stress, n_gwas, n_degs, n_overlap)])

# Specifically bio3 and CHS3
cat("\n=== bio3 GWAS ===\n")
cat(sprintf("bio3 genes: %s\n", paste(gwas$bio3, collapse=", ")))
for (s in names(degs)) {
  ov <- intersect(gwas$bio3, degs[[s]])
  if (length(ov)>0) cat(sprintf("  %s: %s\n", s, paste(ov,collapse=",")))
}

fwrite(rows[, .(gwas_set,stress,n_gwas,n_degs,n_overlap)],
       file.path(INTEG,"gwas_deg_overlap_summary.tsv"), sep="\t")
cat("\nSaved gwas_deg_overlap_summary.tsv\n")

# ── Gene-level GWAS x DEG tables (the Results-text up/down split is read from these) ──
# Temperature-associated (cold) GWAS genes = cold_gwas (bio10+bio11) differentially expressed in the
# spruce cold-needle contrast (SCN); drought = drought_gwas (bio18+bio19) in the spruce drought-needle
# contrast (SDN). Each gene carries its SCN/SDN log2FC and padj plus its co-expression category, taken
# from the canonical per-gene category table (coex_category_universe.tsv, build_category_universe.py).
cat_of <- fread(file.path(INTEG,"coex_category_universe.tsv"))[, .(pa_gene, coex_category)]
gene_level <- function(res, gwas_set) {
  d <- as.data.table(res)[gene %in% gwas_set & !is.na(padj) & padj < 0.01 & abs(log2FoldChange) >= 2]
  d <- merge(d[, .(pa_gene = gene, log2FC = log2FoldChange, padj)], cat_of, by = "pa_gene", all.x = TRUE)
  setorder(d, pa_gene)   # genes with no category-universe entry keep an empty coex_category (NA -> blank)
  d
}
report <- function(tag, d) cat(sprintf("%s: %d genes (%d up, %d down)\n", tag,
                                        nrow(d), sum(d$log2FC > 0), sum(d$log2FC < 0)))
cold_tab <- gene_level(deg_res$SCN, cold_gwas)
fwrite(cold_tab,    file.path(INTEG, "gwas_cold_deg_overlap.tsv"),    sep = "\t"); report("gwas_cold_deg_overlap.tsv", cold_tab)
drought_tab <- gene_level(deg_res$SDN, drought_gwas)
fwrite(drought_tab, file.path(INTEG, "gwas_drought_deg_overlap.tsv"), sep = "\t"); report("gwas_drought_deg_overlap.tsv", drought_tab)
