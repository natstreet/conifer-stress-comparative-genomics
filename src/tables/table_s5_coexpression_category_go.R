infile  <- "results/integration/go/coex_category_go_vs_bg.tsv"
outfile <- "manuscript/supplementary/Supplementary_Table_S5_coexpression_category_GO.xlsx"

go <- read.delim(infile, sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, quote="", colClasses="character")
stopifnot(all(c("GO.ID","Term","Annotated","Significant","Expected","p","ont",
                "plant_consistency","category") %in% names(go)))

go$Annotated   <- as.integer(go$Annotated)
go$Significant <- as.integer(go$Significant)
go$Expected    <- round(suppressWarnings(as.numeric(go$Expected)), 2)
cat_order <- c("conserved","cold_specific","drought_specific","multi_tissue","not_coex")
go$category <- factor(go$category, levels = intersect(cat_order, unique(go$category)))
go <- go[order(go$category, suppressWarnings(as.numeric(sub("^<\\s*","",go$p)))), ]
cat(sprintf("Rows: %d  (%s)\n", nrow(go),
            paste(names(table(as.character(go$category))), table(as.character(go$category)),
                  sep="=", collapse=", ")))

s5 <- data.frame(Category=as.character(go$category), Ontology=go$ont,
  `GO ID`=go$`GO.ID`, Term=go$Term, Annotated=go$Annotated, Significant=go$Significant,
  Expected=go$Expected, `p (topGO weight01, Fisher's exact)`=go$p,
  `Plant consistency`=go$plant_consistency, check.names=FALSE)

caption <- paste0(
"Supplementary Table S5. Gene Ontology enrichment of each co-expression conservation category ",
"(conserved, cold-specific, drought-specific, multi-tissue, and not_coex) tested against the ComPlEx ",
"expression universe (the background of all genes entering the co-expression analysis). Enrichment was ",
"tested with topGO using the weight01 algorithm and Fisher's exact test; because the weight01 algorithm ",
"accounts for the dependency structure of the GO directed acyclic graph (Alexa et al., 2006), the ",
"p-values are reported without further multiple-testing correction and interpreted as enrichment scores ",
"rather than unadjusted Fisher p-values. Annotated, genes annotated to the term in the background ",
"universe; Significant, category genes annotated to the term; Expected, number expected under the null. ",
"The Plant consistency column classifies each term as a reproducibility aid. plant-consistent: the GO ",
"term is present in the Arabidopsis thaliana (Col-0) GO annotation (Bioconductor org.At.tair.db). ",
"cross-kingdom (absent from Arabidopsis GO): the term is not annotated in Arabidopsis and most likely ",
"reflects eggNOG/InterPro annotation transfer from non-plant orthologues. implausible (curated): the term ",
"is present in the Arabidopsis GO annotation but was manually flagged as biologically implausible in a ",
"conifer stress context (GO:0009294, DNA-mediated transformation). Only plant-consistent terms are ",
"discussed in the main text.")

if (requireNamespace("openxlsx", quietly=TRUE)) {
  library(openxlsx); wb <- createWorkbook(); addWorksheet(wb,"S5")
  writeData(wb,"S5",caption,startRow=1,startCol=1); mergeCells(wb,"S5",cols=1:9,rows=1)
  addStyle(wb,"S5",createStyle(wrapText=TRUE,valign="top",textDecoration="italic"),rows=1,cols=1)
  setRowHeights(wb,"S5",rows=1,heights=110)
  writeData(wb,"S5",s5,startRow=3,headerStyle=createStyle(textDecoration="bold"))
  setColWidths(wb,"S5",cols=1:9,widths=c(16,10,14,48,11,12,11,28,20)); freezePane(wb,"S5",firstActiveRow=4)
  saveWorkbook(wb,outfile,overwrite=TRUE); cat("Wrote:",outfile,"\n")
} else {
  library(writexl)
  write_xlsx(list(`Supplementary Table S5`=s5, Caption=data.frame(Caption=caption, check.names=FALSE)), outfile)
  cat("Wrote (writexl):",outfile,"\n")
}
invisible(system2("python3", c("manuscript/supplementary/strip_xlsx_drawings.py", outfile)))  # submission-clean xlsx
cat(nrow(s5),"rows\n")
