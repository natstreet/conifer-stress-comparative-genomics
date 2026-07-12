#!/usr/bin/env Rscript
# figure1c_deg_counts.R — Figure 1c: per-timepoint DEG counts from the DESeq2 contrasts
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
})
args <- commandArgs(trailingOnly=FALSE)
sp   <- sub("--file=","",args[grep("--file=",args)])
PD   <- normalizePath(file.path(if(length(sp)>0) dirname(normalizePath(sp)) else getwd(),"../.."))
DEGL <- file.path(PD,"data/DEG_lists")
FIGS <- file.path(PD,"results")

ld_bar <- function(cond) {
  f <- file.path(DEGL,paste0("DEG_BARPLOTcount_",cond,"_01-2.RData"))
  if(!file.exists(f)){warning("Missing: ",f); return(NULL)}
  load(f); get(ls()[grep("deg_counts",ls())])
}

# Condition order for each stress
cold_order   <- c("5c_6h","5c_24h","5c_3d","5c_10d","neg5c_6h","neg5c_24h","neg5c_3d","neg5c_10d")
drought_order<- c("fc60","fc40","fc30","fc30d7","c1d","c2d","reh")

cold_labs    <- c("5\u00b0C\n6h","5\u00b0C\n24h","5\u00b0C\n3d","5\u00b0C\n10d",
                  "-5\u00b0C\n6h","-5\u00b0C\n24h","-5\u00b0C\n3d","-5\u00b0C\n10d")
drought_labs <- c("FC60%","FC40%","FC30%","FC30%\n7d","Collapsed","Col.\n2d","Rehydrated")

prep_df <- function(cond, order_vec, labels) {
  df <- ld_bar(cond)
  if(is.null(df)) return(NULL)
  df |>
    mutate(condition=tolower(as.character(condition)),
           regulation=factor(regulation, c("Up","Down"))) |>
    filter(condition %in% order_vec) |>
    mutate(condition=factor(condition, order_vec, labels),
           n_signed=ifelse(regulation=="Down",-count,count))
}

theme_bar <- theme_classic(base_size=8) +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=6),
        strip.text=element_text(size=8,face="bold"),
        strip.background=element_blank(),
        legend.position="none",
        axis.line=element_line(colour="grey40"))

make_panel <- function(conds, stress_lab, order_vec, labels) {
  df <- bind_rows(lapply(conds, function(c) {
    d <- prep_df(c, order_vec, labels)
    if(is.null(d)) return(NULL)
    d |> mutate(dataset=c)
  }))
  if(nrow(df)==0) return(NULL)

  # dataset labels
  ds_labs <- c(SCN="P.a. needle",SCR="P.a. root",SDN="P.a. needle",SDR="P.a. root",
               PCN="P.s. needle",PCR="P.s. root",PDN="P.s. needle",PDR="P.s. root")
  df <- df |> mutate(tissue=ds_labs[dataset])

  ggplot(df, aes(condition, n_signed, fill=regulation)) +
    geom_col(width=0.75) +
    geom_hline(yintercept=0, colour="grey40", linewidth=0.3) +
    facet_wrap(~tissue, nrow=1, scales="free_x") +
    scale_fill_manual(values=c(Up="#D55E00",Down="#0072B2"), guide="none") +
    scale_y_continuous(labels=function(x) abs(x), name="DEGs") +
    labs(subtitle=stress_lab, x=NULL) +
    theme_bar
}

p_sc <- make_panel(c("SCN","SCR"), "Spruce cold",   cold_order,    cold_labs)
p_pc <- make_panel(c("PCN","PCR"), "Pine cold",     cold_order,    cold_labs)
p_sd <- make_panel(c("SDN","SDR"), "Spruce drought",drought_order, drought_labs)
p_pd <- make_panel(c("PDN","PDR"), "Pine drought",  drought_order, drought_labs)

panels <- Filter(Negate(is.null), list(p_sc, p_pc, p_sd, p_pd))
if (length(panels) == 4) {
  p_all <- (p_sc | p_pc) / (p_sd | p_pd)
  ggsave(file.path(FIGS,"fig_deg_barplots.pdf"), p_all,
         width=22, height=12, units="cm", device="pdf")
  message("Saved fig_deg_barplots.pdf — verified per-timepoint DESeq2 counts")
} else message("Some panels missing")
