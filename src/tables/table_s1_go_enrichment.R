infile  <- "results/integration/go/deg_directional_go.tsv"
outfile <- "manuscript/supplementary/Supplementary_Table_S1_GO_enrichment.xlsx"

go <- read.delim(infile, sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, quote="", colClasses="character")
stopifnot(all(c("Label","GO.ID","Term","Annotated","Significant","Expected") %in% names(go)))

# The plant-consistency classification is computed in the source producer (deg_directional_go.R) and
# carried in deg_directional_go.tsv, so it is read here (single source of truth), not recomputed.
stopifnot("plant_consistency" %in% names(go))
cat(sprintf("Plant-consistency (from source) of %d rows: %s\n", nrow(go),
            paste(names(table(go$plant_consistency)), table(go$plant_consistency), sep="=", collapse=", ")))

pcol <- intersect(c("pval","weight01","classicFisher","p.value","p"), names(go))[1]
if (is.na(pcol)) pcol <- tail(names(go), 1)
go$Annotated <- as.integer(go$Annotated); go$Significant <- as.integer(go$Significant)
go$Expected  <- round(suppressWarnings(as.numeric(go$Expected)), 2)
go <- go[order(go$Label, suppressWarnings(as.numeric(sub("^<\\s*","",go[[pcol]])))), ]

s1 <- data.frame(Condition=go$Label, `GO ID`=go$`GO.ID`, Term=go$Term,
  Annotated=go$Annotated, Significant=go$Significant, Expected=go$Expected,
  `p (topGO weight01, Fisher's exact)`=go[[pcol]],
  `Plant consistency`=go$plant_consistency, check.names=FALSE)

caption <- paste0(
"Supplementary Table S1. Gene Ontology (biological process) enrichment of differentially expressed genes, by species, stress, ",
"tissue and direction of regulation. For each of the sixteen species x stress x tissue x direction groups, all significantly ",
"enriched GO terms are listed. Enrichment was tested with topGO using the weight01 algorithm and Fisher's exact test; because ",
"the weight01 algorithm accounts for the dependency structure of the GO directed acyclic graph (Alexa et al., 2006), the ",
"p-values are reported without further multiple-testing correction and interpreted as enrichment scores rather than unadjusted ",
"Fisher p-values. Annotated, genes annotated to the term in the analysed universe; Significant, DEGs annotated to the term; ",
"Expected, number expected under the null. All enriched terms are listed; the Plant consistency column ",
"classifies each term as a reproducibility aid. plant-consistent: the GO term is present in the ",
"Arabidopsis thaliana (Col-0) GO annotation (Bioconductor org.At.tair.db). cross-kingdom (absent from ",
"Arabidopsis GO): the term is not annotated in Arabidopsis and most likely reflects eggNOG/InterPro ",
"annotation transfer from non-plant orthologues. implausible (curated): the term is present in the ",
"Arabidopsis GO annotation but was manually flagged as biologically implausible in a conifer stress ",
"context (GO:0009294, DNA-mediated transformation). Only plant-consistent terms are discussed in the main text.")

if (requireNamespace("openxlsx", quietly=TRUE)) {
  library(openxlsx); wb <- createWorkbook(); addWorksheet(wb,"S1")
  writeData(wb,"S1",caption,startRow=1,startCol=1); mergeCells(wb,"S1",cols=1:8,rows=1)
  addStyle(wb,"S1",createStyle(wrapText=TRUE,valign="top",textDecoration="italic"),rows=1,cols=1)
  setRowHeights(wb,"S1",rows=1,heights=90)
  writeData(wb,"S1",s1,startRow=3,headerStyle=createStyle(textDecoration="bold"))
  setColWidths(wb,"S1",cols=1:8,widths=c(28,14,48,11,12,11,28,20)); freezePane(wb,"S1",firstActiveRow=4)
  saveWorkbook(wb,outfile,overwrite=TRUE); cat("Wrote (openxlsx, caption row 1):",outfile,"\n")
} else {
  library(writexl)
  write_xlsx(list(`Supplementary Table S1`=s1, Caption=data.frame(Caption=caption, check.names=FALSE)), outfile)
  cat("Wrote (writexl: table + Caption sheet):",outfile,"\n")
}
invisible(system2("python3", c("manuscript/supplementary/strip_xlsx_drawings.py", outfile)))  # submission-clean xlsx
cat(nrow(s1),"rows,",length(unique(s1$Condition)),"conditions\n")
