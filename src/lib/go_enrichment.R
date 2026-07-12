## go_enrichment.R — the shared GO-enrichment method used across the paper. Every GO script sources
## this so the framework is identical everywhere, with GO enrichment restricted to terms in the
## Arabidopsis GO universe.
##
## Settings:
##   * topGO, algorithm = "weight01", statistic = "fisher"
##   * nodeSize = 5  (drop GO terms with <5 annotated genes — unstable otherwise)
##   * NO multiple-testing correction: weight01 already accounts for GO-graph dependence
##     (Alexa et al. 2006); its p-values are not independent, so BH is not applied.
##   * ALL terms at p < 0.05 are returned, each carrying a `plant_consistent` flag
##     (present in the Arabidopsis GO universe, org.At.tair.db, and not curated-implausible).
##     Callers decide scope: SUPPLEMENTARY tables keep all ontologies + all terms + the flag;
##     MAIN-TEXT tables filter to ontology "BP" and plant_consistent == TRUE.
suppressPackageStartupMessages({ library(topGO); library(data.table); library(org.At.tair.db) })

## Arabidopsis GO universe + curated implausible terms (single definition, used everywhere)
.AT_GO <- unique(suppressWarnings(AnnotationDbi::keys(org.At.tair.db, keytype = "GO")))
.GO_IMPLAUSIBLE <- c("GO:0009294")  # DNA-mediated transformation: in At GO but implausible in conifer stress

## Build gene -> GO(list) from an eggNOG+InterPro annotation table (spruce tsv or gzipped pine tsv).
build_gene2go <- function(annot_file) {
  a <- if (grepl("\\.gz$", annot_file))
         fread(cmd = sprintf("gunzip -c %s", shQuote(annot_file)), quote = "")
       else fread(annot_file, quote = "")
  a[, gene_id := sub("\\.mRNA\\.\\d+$", "", id)]
  parse_go <- function(x) { x <- x[!is.na(x) & x != "NA" & x != ""]
    t <- unlist(strsplit(x, ",")); unique(t[grepl("^GO:", t)]) }
  g <- a[, .(go = list(parse_go(c(eggnog_go, interpro_go)))), by = gene_id][lengths(go) > 0]
  setNames(g$go, g$gene_id)
}

## Run GO enrichment for one significant gene set against a background.
## Returns a data.table: GO.ID, Term, Annotated, Significant, Expected, p, ont, plant_consistent,
## plant_consistency — every term at p<0.05 across the requested ontologies. No BH column (by design).
run_go <- function(sig_genes, background, gene2go,
                   ontologies = c("BP", "MF", "CC"), node_size = 5, p_cutoff = 0.05) {
  sig <- intersect(sig_genes, background)
  if (length(sig) < 5) return(data.table())
  gvec <- factor(as.integer(background %in% sig), levels = c(0, 1)); names(gvec) <- background
  one <- function(ont) {
    gd <- new("topGOdata", ontology = ont, allGenes = gvec, geneSel = function(x) x == 1,
              annot = annFUN.gene2GO, gene2GO = gene2go, nodeSize = node_size)
    res <- runTest(gd, algorithm = "weight01", statistic = "fisher")
    tb  <- GenTable(gd, p = res, orderBy = "p", topNodes = length(usedGO(gd)), numChar = 200)
    tb$ont <- ont; tb$p <- suppressWarnings(as.numeric(sub("< ", "", tb$p)))
    as.data.table(tb)[!is.na(p) & p <= p_cutoff]
  }
  out <- rbindlist(lapply(ontologies, one), fill = TRUE)
  if (nrow(out) == 0) return(out)
  out[, plant_consistent  := (GO.ID %in% .AT_GO) & !(GO.ID %in% .GO_IMPLAUSIBLE)]
  out[, plant_consistency := fifelse(GO.ID %in% .GO_IMPLAUSIBLE, "implausible (curated)",
                             fifelse(GO.ID %in% .AT_GO,           "plant-consistent",
                                                                  "cross-kingdom (absent from Arabidopsis GO)"))]
  out[]
}
