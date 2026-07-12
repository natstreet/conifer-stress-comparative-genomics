#!/usr/bin/env Rscript
# build_deg_overlap_counts.R
#
# Computes the Figure 2 / Results DEG-overlap counts quoted in the text ("834 genes in Norway
# spruce and 1,303 in Scots pine", etc.) as gene-level set operations on the deposited DE lists
# (data/DEG_lists/DE_all_<cond>_01_2L2FC.RData), per species:
#
#   cold-specific overlap    = (cold-needle INT cold-root)   \ (drought-needle UNION drought-root)
#   drought-specific overlap = (drought-needle INT drought-root) \ (cold-needle UNION cold-root)
#   shared core              = cold-needle INT cold-root INT drought-needle INT drought-root
#
# Note: these are GENE-level counts (the text's "orthogroup" phrasing is loose; the
# orthogroup-level counts differ). Run from the AbioticStressConifers directory.
# Output: results/deg_overlap_counts.tsv  (columns: species, overlap, count)

DEGL <- "data/DEG_lists"

ld <- function(cond) {
  f <- file.path(DEGL, paste0("DE_all_", cond, "_01_2L2FC.RData"))
  e <- new.env(); load(f, envir = e)
  get(ls(e)[grep("DE_all", ls(e))], envir = e)
}

overlaps <- function(CN, CR, DN, DR) {
  c(cold_specific    = length(setdiff(intersect(CN, CR), union(DN, DR))),
    drought_specific = length(setdiff(intersect(DN, DR), union(CN, CR))),
    shared_core      = length(Reduce(intersect, list(CN, CR, DN, DR))))
}

sp <- overlaps(ld("SCN"), ld("SCR"), ld("SDN"), ld("SDR"))
pi <- overlaps(ld("PCN"), ld("PCR"), ld("PDN"), ld("PDR"))

res <- data.frame(
  species = rep(c("Picea_abies", "Pinus_sylvestris"), each = 3),
  overlap = rep(names(sp), 2),
  count   = c(sp, pi),
  stringsAsFactors = FALSE)

dir.create("results", showWarnings = FALSE, recursive = TRUE)
write.table(res, "results/deg_overlap_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote results/deg_overlap_counts.tsv\n")
print(res)
