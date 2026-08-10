#!/usr/bin/env Rscript
## no1to1_broad_species_class.R
## Broad-species classification of the no_1to1 stress-DE set using the 27-species N10 hierarchical
## orthogroups (Phylogenetic_Hierarchical_Orthogroups_N10.tsv.gz), which include angiosperm and other
## gymnosperm outgroups. Splits the "no 1:1 orthologue in the sister conifer" set into:
##   novel_gain       = the gene's HOG contains NO species other than the focal conifer (or the gene is in
##                      no HOG at all) -> genuine lineage-specific novelty; no homolog to compute dN/dS against.
##   lost_or_diverged = no member in the SISTER conifer, but >=1 member in some outgroup (angiosperm and/or
##                      other gymnosperm) -> ancient gene whose 1:1 sister partner was lost / strongly diverged.
##   expanded_family  = the sister conifer HAS >=1 member, but the relationship is non-1:1.
## NOTE: this is an orthogroup-level (N10 hierarchical) view — a DIFFERENT method from the reciprocal-best
## 1:1 backbone that set the no_1to1 flag; frame the text as "broad-species orthogroup analysis".
suppressPackageStartupMessages(library(data.table))
I   <- Sys.getenv("INTEG_DIR", "results/integration")
## N0 = root node of the 27-species OrthoFinder run: ALL species (angiosperms, other gymnosperms, moss)
## with our current Picea/Pinus gene IDs. (The N10 node holds only the 3 Pinaceae and cannot see angiosperms.)
HOG <- Sys.getenv("HOG_FILE", "Orthogroups/Phylogenetic_Hierarchical_Orthogroups_N0.tsv.gz")
strip_tx <- function(x) sub("\\.mRNA\\.[0-9]+$", "", x)

og <- fread(cmd = paste("gzip -dc", shQuote(HOG)), sep = "\t", header = TRUE, quote = "")
sp_cols       <- setdiff(names(og), c("HOG", "OG", "Gene Tree Parent Clade"))
stopifnot("Picea_abies" %in% sp_cols, "Pinus_sylvestris" %in% sp_cols)
## fread can type sparse species columns oddly / read empty cells as NA -> force character, NA => ""
og[, (sp_cols) := lapply(.SD, function(x) { x <- as.character(x); x[is.na(x)] <- ""; x }), .SDcols = sp_cols]
outgroup_cols <- setdiff(sp_cols, c("Picea_abies", "Pinus_sylvestris"))  # 25 outgroups (incl. Pinus_tabuliformis)
## angiosperm + moss columns: an orthologue here means the family predates the seed-plant split (~300 My),
## a taxon-sampling-ROBUST "ancient" signal (unlike the conifer congener Pinus_tabuliformis).
ancient_cols <- intersect(c("Amborella_trichopoda","Arabidopsis_thaliana","Beta_vulgaris","Betula_pendula",
  "Camellia_sinensis","Cinnamomum_kanehirae","Cucumis_sativus","Eucalyptus_grandis","Malania_oleifera",
  "Malus_domestica","Medicago_truncatula","Mimulus_guttatus","Nicotiana_tabacum","Populus_tremula",
  "Prunus_avium","Salix_purpurea","Vitis_vinifera","Zostera_marina","Physcomitrella_patens"), sp_cols)

## per-HOG presence flags (one row per HOG in N10)
flags <- data.table(
  HOG          = og$HOG,
  has_picea    = og$Picea_abies      != "",
  has_psylv    = og$Pinus_sylvestris != "",
  has_outgroup = og[, Reduce(`|`, lapply(.SD, function(x) x != "")), .SDcols = outgroup_cols],
  has_ancient  = og[, Reduce(`|`, lapply(.SD, function(x) x != "")), .SDcols = ancient_cols])

## gene -> HOG (strip transcript suffix)
melt_map <- function(col) {
  d <- og[get(col) != "", .(HOG, cell = get(col))]
  unique(d[, .(gene = strip_tx(trimws(unlist(strsplit(cell, ",", fixed = TRUE))))), by = HOG][gene != ""],
         by = "gene")
}
sp_map <- melt_map("Picea_abies")
pi_map <- melt_map("Pinus_sylvestris")

classify <- function(genes, gmap, sister_flag) {   # sister_flag: "has_psylv" (spruce) / "has_picea" (pine)
  d <- merge(data.table(gene = genes), gmap,  by = "gene", all.x = TRUE)
  d <- merge(d, flags, by = "HOG", all.x = TRUE)
  d[, class := fifelse(is.na(HOG), "novel_gain",
                fifelse(get(sister_flag) == TRUE, "expanded_family",
                fifelse(has_outgroup == TRUE, "lost_or_diverged", "novel_gain")))]
  ## depth axis (taxon-sampling-robust): does the family reach angiosperms/moss?
  d[, deepest_evidence := fifelse(class == "novel_gain", "conifer_lineage_specific",
                           fifelse(has_ancient == TRUE, "angiosperm_shared", "gymnosperm_only"))]
  d[]
}
sp_nc <- fread(file.path(I, "not_coex_de_genes.tsv"))[mechanism == "no_1to1_ortholog", pa_gene]
pi_nc <- fread(file.path(I, "not_coex_de_genes_PINE.tsv"))[mechanism == "no_1to1_ortholog", ps_gene]
stopifnot(length(sp_nc) == 1859L, length(pi_nc) == 2676L)   # sanity gate

sp_cls <- classify(sp_nc, sp_map, "has_psylv")[, species := "Picea_abies"]
pi_cls <- classify(pi_nc, pi_map, "has_picea")[, species := "Pinus_sylvestris"]

## per-gene output; outgroup_homolog_present (dN/dS computable in principle) = class != novel_gain
mk <- function(d) d[, .(species, gene, HOG = fcoalesce(HOG, ""), class, deepest_evidence,
                        outgroup_homolog_present = fifelse(class == "novel_gain", "no", "yes"))]
by_gene <- rbind(mk(sp_cls), mk(pi_cls))
fwrite(by_gene, file.path(I, "no1to1_broad_species_class_by_gene.tsv"), sep = "\t")

## primary summary: the three sister-conifer classes
summ <- by_gene[, .(n = .N), by = .(species, class)]
summ[, total := sum(n), by = species][, pct := round(100 * n / total, 1)]
setorder(summ, species, -n)
fwrite(summ, file.path(I, "no1to1_broad_species_class.tsv"), sep = "\t")

## depth summary (taxon-sampling-robust): conifer-lineage-specific / gymnosperm-only / angiosperm-shared
depth <- by_gene[, .(n = .N), by = .(species, deepest_evidence)]
depth[, total := sum(n), by = species][, pct := round(100 * n / total, 1)]
setorder(depth, species, -n)
fwrite(depth, file.path(I, "no1to1_broad_species_depth.tsv"), sep = "\t")

cat("Broad-species (N0, 27 species) no_1to1 sister-conifer classification:\n"); print(summ)
cat("\nDepth axis (robust to conifer taxon sampling):\n"); print(depth)
cat(sprintf("\nunassigned-in-N10 (-> novel_gain): spruce %d, pine %d\n",
            sum(is.na(sp_cls$HOG)), sum(is.na(pi_cls$HOG))))
