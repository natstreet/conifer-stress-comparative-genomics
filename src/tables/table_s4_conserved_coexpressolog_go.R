infile  <- "results/integration/fig4/table2_conserved_dynamics_go.tsv"
outfile <- "manuscript/supplementary/Supplementary_Table_S4_conserved_coexpressolog_GO.xlsx"

go <- read.delim(infile, sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, quote="", colClasses="character")
stopifnot(all(c("stress","direction","ont","GO.ID","Term","Annotated","Significant","Expected","p",
                "plant_consistency","n_genes") %in% names(go)))

go$Annotated  <- as.integer(go$Annotated)
go$Significant <- as.integer(go$Significant)
go$Expected   <- round(suppressWarnings(as.numeric(go$Expected)), 2)
go$n_genes    <- as.integer(go$n_genes)
go <- go[order(go$stress, go$direction, suppressWarnings(as.numeric(sub("^<\\s*","",go$p)))), ]
cat(sprintf("Rows: %d  (%s)\n", nrow(go),
            paste(names(table(paste(go$stress, go$direction))),
                  table(paste(go$stress, go$direction)), sep="=", collapse=", ")))

s4 <- data.frame(Stress=go$stress, Direction=go$direction, Ontology=go$ont,
  `GO ID`=go$`GO.ID`, Term=go$Term, Annotated=go$Annotated, Significant=go$Significant,
  Expected=go$Expected, `p (topGO weight01, Fisher's exact)`=go$p,
  `Plant consistency`=go$plant_consistency, Genes=go$n_genes, check.names=FALSE)

caption <- paste0(
"Supplementary Table S4. Gene Ontology enrichment of the genes in the conserved cross-species ",
"co-expressolog sets (co-expressologs conserved between Norway spruce and Scots pine), split by stress ",
"(cold, drought) and direction of regulation (up, down). Enrichment was tested with topGO using the ",
"weight01 algorithm and Fisher's exact test; because the weight01 algorithm accounts for the dependency ",
"structure of the GO directed acyclic graph (Alexa et al., 2006), the p-values are reported without ",
"further multiple-testing correction and interpreted as enrichment scores rather than unadjusted Fisher ",
"p-values. Annotated, genes annotated to the term in the analysed universe; Significant, conserved-set ",
"genes annotated to the term; Expected, number expected under the null; Genes, conserved-set genes in the ",
"term. The Plant consistency column classifies each term as a reproducibility aid. plant-consistent: the ",
"GO term is present in the Arabidopsis thaliana (Col-0) GO annotation (Bioconductor org.At.tair.db). ",
"cross-kingdom (absent from Arabidopsis GO): the term is not annotated in Arabidopsis and most likely ",
"reflects eggNOG/InterPro annotation transfer from non-plant orthologues. implausible (curated): the term ",
"is present in the Arabidopsis GO annotation but was manually flagged as biologically implausible in a ",
"conifer stress context (GO:0009294, DNA-mediated transformation). Only plant-consistent terms are ",
"discussed in the main text.")

if (requireNamespace("openxlsx", quietly=TRUE)) {
  library(openxlsx); wb <- createWorkbook(); addWorksheet(wb,"S4")
  writeData(wb,"S4",caption,startRow=1,startCol=1); mergeCells(wb,"S4",cols=1:11,rows=1)
  addStyle(wb,"S4",createStyle(wrapText=TRUE,valign="top",textDecoration="italic"),rows=1,cols=1)
  setRowHeights(wb,"S4",rows=1,heights=110)
  writeData(wb,"S4",s4,startRow=3,headerStyle=createStyle(textDecoration="bold"))
  setColWidths(wb,"S4",cols=1:11,widths=c(10,11,10,14,48,11,12,11,28,20,9)); freezePane(wb,"S4",firstActiveRow=4)
  saveWorkbook(wb,outfile,overwrite=TRUE); cat("Wrote:",outfile,"\n")
} else {
  library(writexl)
  write_xlsx(list(`Supplementary Table S4`=s4, Caption=data.frame(Caption=caption, check.names=FALSE)), outfile)
  cat("Wrote (writexl):",outfile,"\n")
}
invisible(system2("python3", c("manuscript/supplementary/strip_xlsx_drawings.py", outfile)))  # submission-clean xlsx
cat(nrow(s4),"rows\n")
