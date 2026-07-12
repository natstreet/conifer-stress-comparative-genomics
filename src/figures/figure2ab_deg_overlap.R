#!/usr/bin/env Rscript
# figure2ab_deg_overlap.R — Figure 2a/b: within-species DEG overlap (2-panel patchwork)
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork); library(ComplexUpset)
})
source("src/lib/fig_palette.R")
DEGL <- "data/DEG_lists"
FIGS <- "results"

load_deg <- function(cond) {
  f <- file.path(DEGL, paste0("DE_all_", cond, "_01_2L2FC.RData"))
  if (!file.exists(f)) return(character(0))
  load(f); get(ls()[grep("DE_all", ls())])
}
make_binary_matrix <- function(gene_lists) {
  gene_lists <- gene_lists[sapply(gene_lists, length) > 0]
  all_genes  <- unique(unlist(gene_lists))
  as.data.frame(lapply(gene_lists, function(g) as.logical(all_genes %in% g)),
                row.names=all_genes, check.names=FALSE)
}
make_upset <- function(mat, title, barfill) {
  set_names <- names(mat)
  upset(mat, set_names,
    base_annotations=list(
      "Intersection size"=intersection_size(text=list(size=2.8), fill=barfill) +
        scale_y_continuous(expand=expansion(mult=c(0,0.22)),
                           labels=scales::comma)
    ),
    set_sizes=(upset_set_size(geom=geom_bar(fill=barfill, width=0.6)) +
      scale_y_continuous(labels=scales::comma)),
    stripes="white",                                  # house style: white matrix background
    themes=upset_modify_themes(list(
      "intersections_matrix"=theme(text=element_text(size=8.5),
                                   panel.background=element_rect(fill="white", colour=NA)),
      "overall_sizes"=theme(text=element_text(size=8.5),
                            panel.background=element_rect(fill="white", colour=NA))
    )),
    sort_sets=FALSE, width_ratio=0.25
  ) + labs(title=title) +
    theme(plot.title=element_text(size=9, face="bold"))
}

sp_mat <- make_binary_matrix(list(
  "Cold\nneedle"=load_deg("SCN"), "Cold\nroot"=load_deg("SCR"),
  "Drought\nneedle"=load_deg("SDN"), "Drought\nroot"=load_deg("SDR")
))
pi_mat <- make_binary_matrix(list(
  "Cold\nneedle"=load_deg("PCN"), "Cold\nroot"=load_deg("PCR"),
  "Drought\nneedle"=load_deg("PDN"), "Drought\nroot"=load_deg("PDR")
))

p_sp <- make_upset(sp_mat, "Spruce (P. abies)", PAL$spruce)
p_pi <- make_upset(pi_mat, "Pine (P. sylvestris)", PAL$pine)

# One panel tag per UpSet plot (a = spruce, b = pine) to match the caption "(a) P. abies and
# (b) P. sylvestris"; wrap_elements makes each UpSet a single taggable unit so the three internal
# UpSet components are NOT tagged.
p_combined <- wrap_elements(full = p_sp) / wrap_elements(full = p_pi) +
  plot_annotation(title="Differentially expressed gene overlap across stress conditions",
                  tag_levels = list(c("a", "b")),
                  theme=theme(plot.title=element_text(size=10, face="bold"))) &
  theme(plot.tag = element_text(size=12, face="bold"))

ggsave(file.path(FIGS,"fig_deg_overlap.pdf"), p_combined,
       width=16, height=18, units="cm", device="pdf")
message("Saved fig_deg_overlap.pdf (2-panel, no clipping)")
