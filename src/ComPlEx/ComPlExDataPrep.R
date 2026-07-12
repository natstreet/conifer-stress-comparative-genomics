#' ---
#' title: "ComPlEx data setup"
#' author: "Elena van Zalen"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' # Setup

#' 
#' Loading libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(here)
  library(ggplot2)
})
theme_set(theme_classic())
theme_update(plot.title = element_text(face="bold"))

#' Loading the data
load(here("data/dds/dds_SCN.rda"))
load(here("data/dds/dds_SCR.rda"))
load(here("data/dds/dds_PCN.rda"))
load(here("data/dds/dds_PCR.rda"))

load(here("data/dds/dds_SC.rda"))
load(here("data/dds/dds_PC.rda"))

load(here("data/dds/dds_SDN.rda"))
load(here("data/dds/dds_SDR.rda"))
load(here("data/dds/dds_PDN.rda"))
load(here("data/dds/dds_PDR.rda"))

load(here("data/dds/dds_SD.rda"))
load(here("data/dds/dds_PD.rda"))

#load(here("data/dds/dds_SAll.rda"))
#load(here("data/dds/dds_PAll.rda"))


# Generate expression files
#' # Variance Stabilising Transformation
#' Drought
#' SDN
vsd_SDN <- varianceStabilizingTransformation(dds_SDN, blind=TRUE)
vst_SDN <- assay(vsd_SDN)
vst_SDN <- vst_SDN - min(vst_SDN)
exp_SDN <- data.frame(Genes = row.names(vst_SDN), vst_SDN)
colnames(exp_SDN) <- gsub("^X", "SDN", colnames(exp_SDN))
write.table(exp_SDN, file = here("data/expression/SDN_expression.txt"), 
            sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)

#' SDR
vsd_SDR <- varianceStabilizingTransformation(dds_SDR, blind=TRUE)
vst_SDR <- assay(vsd_SDR)
vst_SDR <- vst_SDR - min(vst_SDR)
exp_SDR <- data.frame(Genes = row.names(vst_SDR), vst_SDR)
colnames(exp_SDR) <- gsub("^X", "SDR", colnames(exp_SDR))
write.table(exp_SDR, file = here("data/expression/SDR_expression.txt"), 
            sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)

#' PDN
vsd_PDN <- varianceStabilizingTransformation(dds_PDN, blind=TRUE)
vst_PDN <- assay(vsd_PDN)
vst_PDN <- vst_PDN - min(vst_PDN)
exp_PDN <- data.frame(Genes = row.names(vst_PDN), vst_PDN)
write.table(exp_PDN, file = here("data/expression/PDN_expression.txt"), 
            sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)

#' PDR
vsd_PDR <- varianceStabilizingTransformation(dds_PDR, blind=TRUE)
vst_PDR <- assay(vsd_PDR)
vst_PDR <- vst_PDR - min(vst_PDR)
exp_PDR <- data.frame(Genes = row.names(vst_PDR), vst_PDR)
write.table(exp_PDR, file = here("data/expression/PDR_expression.txt"), 
            sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)

#' SD
vsd_SD <- varianceStabilizingTransformation(dds_SD, blind=TRUE)
vst_SD <- assay(vsd_SD)
vst_SD <- vst_SD - min(vst_SD)
exp_SD <- data.frame(Genes = row.names(vst_SD), vst_SD)
colnames(exp_SD) <- gsub("^X", "SD", colnames(exp_SD))
write.table(exp_SD, file = here("data/expression/SD_expression.txt"), 
            sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)

#' PD
vsd_PD <- varianceStabilizingTransformation(dds_PD, blind=TRUE)
vst_PD <- assay(vsd_PD)
vst_PD <- vst_PD - min(vst_PD)
exp_PD <- data.frame(Genes = row.names(vst_PD), vst_PD)
write.table(exp_PD, file = here("data/expression/PD_expression.txt"), 
            sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)

#' Cold
#' SCN
vsd_SCN <- varianceStabilizingTransformation(dds_SCN, blind=TRUE)
vst_SCN <- assay(vsd_SCN)
vst_SCN <- vst_SCN - min(vst_SCN)
exp_SCN <- data.frame(Genes = row.names(vst_SCN), vst_SCN)
colnames(exp_SCN) <- gsub("^X", "SCN", colnames(exp_SCN))
write.table(exp_SCN, file = here("data/expression/SCN_expression.txt"), 
            sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)

#' SCR
vsd_SCR <- varianceStabilizingTransformation(dds_SCR, blind=TRUE)
vst_SCR <- assay(vsd_SCR)
vst_SCR <- vst_SCR - min(vst_SCR)
exp_SCR <- data.frame(Genes = row.names(vst_SCR), vst_SCR) 
colnames(exp_SCR) <- gsub("^X", "SCR", colnames(exp_SCR))
write.table(exp_SCR, file = here("data/expression/SCR_expression.txt"), 
            sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)

#' PCN
vsd_PCN <- varianceStabilizingTransformation(dds_PCN, blind=TRUE)
vst_PCN <- assay(vsd_PCN)
vst_PCN <- vst_PCN - min(vst_PCN)
exp_PCN <- data.frame(Genes = row.names(vst_PCN), vst_PCN)
colnames(exp_PCN) <- gsub("^X", "PCN", colnames(exp_PCN))
write.table(exp_PCN, file = here("data/expression/PCN_expression.txt"), 
            sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)

#' PCR
vsd_PCR <- varianceStabilizingTransformation(dds_PCR, blind=TRUE)
vst_PCR <- assay(vsd_PCR)
vst_PCR <- vst_PCR - min(vst_PCR)
exp_PCR <- data.frame(Genes = row.names(vst_PCR), vst_PCR)
colnames(exp_PCR) <- gsub("^X", "PCR", colnames(exp_PCR))
write.table(exp_PCR, file = here("data/expression/PCR_expression.txt"), 
            sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)

#' SC
vsd_SC <- varianceStabilizingTransformation(dds_SC, blind=TRUE)
vst_SC <- assay(vsd_SC)
vst_SC <- vst_SC - min(vst_SC)
exp_SC <- data.frame(Genes = row.names(vst_SC), vst_SC)
colnames(exp_SC) <- gsub("^X", "SC", colnames(exp_SC))
write.table(exp_SC, file = here("data/expression/SC_expression.txt"), 
            sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)
#' PC
vsd_PC <- varianceStabilizingTransformation(dds_PC, blind=TRUE)
vst_PC <- assay(vsd_PC)
vst_PC <- vst_PC - min(vst_PC)
exp_PC <- data.frame(Genes = row.names(vst_PC), vst_PC)
colnames(exp_PC) <- gsub("^X", "PC", colnames(exp_PC))
write.table(exp_PC, file = here("data/expression/PC_expression.txt"), 
            sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)


#' #' Spruce All
#' vsd_SAll <- varianceStabilizingTransformation(dds_SAll, blind=TRUE)
#' vst_SAll <- assay(vsd_SAll)
#' vst_SAll <- vst_SAll - min(vst_SAll)
#' exp_SAll <- data.frame(Genes = row.names(vst_SAll), vst_SAll)
#' write.table(exp_SAll, file = here("data/expression/SAll_expression.txt"), 
#'             sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)
#' 
#' #' Pine All
#' vsd_PAll <- varianceStabilizingTransformation(dds_PAll, blind=TRUE)
#' vst_PAll <- assay(vsd_PAll)
#' vst_PAll <- vst_PAll - min(vst_PAll)
#' exp_PAll <- data.frame(Genes = row.names(vst_PAll), vst_PAll)
#' write.table(exp_PAll, file = here("data/expression/PAll_expression.txt"), 
#'             sep = "\t", dec = ".", col.names = TRUE, row.names = FALSE)







