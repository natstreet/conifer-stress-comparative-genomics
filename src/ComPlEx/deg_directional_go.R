#!/usr/bin/env Rscript
# deg_directional_go.R — Supplementary Table S2 (broad GO enrichment)
#
# GO enrichment of the up- and down-regulated DEGs in each species x stress x organ
# dataset. DE is recomputed from the DESeq2 objects (DESeq(); Wald per coefficient;
# filter = per-gene median of normalised counts; DE = padj<0.01 & |log2FC|>=2 — the
# same definition as the DE notebooks). GO uses the single shared method
# (src/lib/go_enrichment.R): topGO weight01/Fisher, nodeSize 5, no BH, BP/MF/CC, all
# terms with a plant_consistent flag. Gene->GO comes from the SAME source for both
# species (eggNOG + InterPro): spruce = data/annotation/gene_annotation.tsv.gz, pine = pine_GO_annotation.tsv.gz.
#
# Run from AbioticStressConifers/.  Output: results/integration/go/deg_directional_go.tsv
suppressPackageStartupMessages({ library(DESeq2); library(matrixStats); library(data.table) })
source("src/lib/go_enrichment.R")

DDS_BASE <- "data/dds"
OUT      <- "results/integration/go/deg_directional_go.tsv"

# Gene -> GO from the single shared annotation source (eggNOG+InterPro), both species
cat("Loading gene->GO (eggNOG+InterPro) for both species...\n")
pa_gene2GO <- build_gene2go("data/annotation/gene_annotation.tsv.gz")                          # Picea abies
ps_gene2GO <- build_gene2go("data/annotation/pine_GO_annotation.tsv.gz")    # Pinus sylvestris
cat(sprintf("  spruce GO genes: %d ; pine GO genes: %d\n", length(pa_gene2GO), length(ps_gene2GO)))

DATASETS <- list(
  list(rda="dds_SCN.rda", species="PA", label="P. abies cold needle",
       names=c("Condition_5C_6h_vs_20C_0h","Condition_5C_24h_vs_20C_0h","Condition_5C_3d_vs_20C_0h",
               "Condition_5C_10d_vs_20C_0h","Condition_neg5C_6h_vs_20C_0h","Condition_neg5C_24h_vs_20C_0h",
               "Condition_neg5C_3d_vs_20C_0h","Condition_neg5C_10d_vs_20C_0h")),
  list(rda="dds_SCR.rda", species="PA", label="P. abies cold root",
       names=c("Condition_5C_6h_vs_20C_0h","Condition_5C_24h_vs_20C_0h","Condition_5C_3d_vs_20C_0h",
               "Condition_5C_10d_vs_20C_0h","Condition_neg5C_6h_vs_20C_0h","Condition_neg5C_24h_vs_20C_0h",
               "Condition_neg5C_3d_vs_20C_0h","Condition_neg5C_10d_vs_20C_0h")),
  list(rda="dds_SDN.rda", species="PA", label="P. abies drought needle",
       names=c("Level_60._vs_80.","Level_40._vs_80.","Level_30._vs_80.","Level_30.7d_vs_80.",
               "Level_Collapse_vs_80.","Level_C2d_vs_80.","Level_Rehydrate_vs_80.")),
  list(rda="dds_SDR.rda", species="PA", label="P. abies drought root",
       names=c("Level_60._vs_80.","Level_40._vs_80.","Level_30._vs_80.","Level_30.7d_vs_80.",
               "Level_Collapse_vs_80.","Level_C2d_vs_80.","Level_Rehydrate_vs_80.")),
  list(rda="dds_PCN.rda", species="PS", label="P. sylvestris cold needle",
       names=c("Condition_5C_6h_vs_20C_0h","Condition_5C_24h_vs_20C_0h","Condition_5C_3d_vs_20C_0h",
               "Condition_5C_10d_vs_20C_0h","Condition_neg5C_6h_vs_20C_0h","Condition_neg5C_24h_vs_20C_0h",
               "Condition_neg5C_3d_vs_20C_0h","Condition_neg5C_10d_vs_20C_0h")),
  list(rda="dds_PCR.rda", species="PS", label="P. sylvestris cold root",
       names=c("Condition_5C_6h_vs_20C_0h","Condition_5C_24h_vs_20C_0h","Condition_5C_3d_vs_20C_0h",
               "Condition_5C_10d_vs_20C_0h","Condition_neg5C_6h_vs_20C_0h","Condition_neg5C_24h_vs_20C_0h",
               "Condition_neg5C_3d_vs_20C_0h","Condition_neg5C_10d_vs_20C_0h")),
  list(rda="dds_PDN.rda", species="PS", label="P. sylvestris drought needle",
       names=c("Condition_FC60_vs_FC80","Condition_FC40_vs_FC80","Condition_FC30_vs_FC80",
               "Condition_FC30d7_vs_FC80","Condition_Collapsed_vs_FC80","Condition_Collapsed2d_vs_FC80",
               "Condition_Rehydrated_vs_FC80")),
  list(rda="dds_PDR.rda", species="PS", label="P. sylvestris drought root",
       names=c("Condition_FC60_vs_FC80","Condition_FC40_vs_FC80","Condition_FC30_vs_FC80",
               "Condition_FC30d7_vs_FC80","Condition_Collapsed_vs_FC80","Condition_Collapsed2d_vs_FC80",
               "Condition_Rehydrated_vs_FC80"))
)

all_results <- list()
for (ds in DATASETS) {
  cat("\n=== ", ds$label, " ===\n", sep="")
  e <- new.env(); load(file.path(DDS_BASE, ds$rda), envir=e)
  dds  <- DESeq(e[[ls(envir=e)[1]]], quiet=TRUE)
  filt <- rowMedians(counts(dds, normalized=TRUE))                 # matches the DE notebooks
  bg_df <- as.data.frame(results(dds, name=ds$names[1], filter=filt, independentFiltering=TRUE))
  bg_genes <- rownames(bg_df[!is.na(bg_df$padj), ])                # expressed & tested = background
  up_genes <- character(0); down_genes <- character(0)
  for (nm in ds$names) {
    r <- as.data.frame(results(dds, name=nm, filter=filt, independentFiltering=TRUE))
    r <- r[!is.na(r$padj), ]
    up_genes   <- union(up_genes,   rownames(r[r$padj<0.01 & r$log2FoldChange>=2,  ]))
    down_genes <- union(down_genes, rownames(r[r$padj<0.01 & r$log2FoldChange<=-2, ]))
  }
  cat(sprintf("  tested %d ; DEGs up %d / down %d\n", length(bg_genes), length(up_genes), length(down_genes)))
  g2go <- if (ds$species=="PA") pa_gene2GO else ps_gene2GO
  for (dir in c("UP","DOWN")) {
    sig <- if (dir=="UP") up_genes else down_genes
    out <- run_go(sig, bg_genes, g2go, ontologies=c("BP","MF","CC"))
    if (nrow(out) == 0) next
    out[, Label := paste0(ds$label, " [", dir, "]")]
    all_results[[paste0(ds$label, "_", dir)]] <- out
  }
}

final <- rbindlist(all_results, fill=TRUE)
setcolorder(final, c("Label","ont","GO.ID","Term","Annotated","Significant","Expected","p",
                     "plant_consistent","plant_consistency"))
fwrite(final, OUT, sep="\t")
cat(sprintf("\nSaved %d rows to %s (%d plant-consistent)\n", nrow(final), OUT, sum(final$plant_consistent)))
