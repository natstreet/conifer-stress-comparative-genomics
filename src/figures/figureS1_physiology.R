#!/usr/bin/env Rscript
# Figure S1 — Scots pine drought experiment physiology.
# Reproduces the published 6-panel figure from the raw measurement data.
# Source data: D. Castro (data/physiology/). Run from AbioticStressConifers/.
#
# Fixes over the original Plots.R:
#   * correct, vectorised mean +/- 95% CI (the original conint() indexed a tibble
#     by rownames and referenced a global, giving unreliable intervals);
#   * root-RWC SD was computed from RWC_shoot (copy-paste bug) — corrected;
#   * panels combined into the single relabelled 6-panel figure used in the paper;
#   * shared Control/Drought legend + stress-phase colour key.

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })

PHYS  <- "data/physiology"
meta  <- read.csv(file.path(PHYS, "metadata.csv"), sep = ";", dec = ".")
photo <- read.csv(file.path(PHYS, "scots_pine_drought_physiology.csv"), sep = ";", dec = ".")
meta$Treatment  <- factor(meta$Treatment,  c("Control", "Drought"))
photo$Treatment <- factor(photo$Treatment, c("Control", "Drought"))

# mean +/- 95% CI per Day x Treatment
summ <- function(df, col) {
  df %>% group_by(Day, Treatment) %>%
    summarise(m = mean(.data[[col]], na.rm = TRUE),
              s = sd(.data[[col]],   na.rm = TRUE),
              n = sum(!is.na(.data[[col]])), .groups = "drop") %>%
    mutate(err = qt(0.975, pmax(n - 1, 1)) * s / sqrt(n), lo = m - err, hi = m + err)
}

phases <- data.frame(
  phase = factor(c("Pre-treatment","Mild stress","Moderate stress","Severe stress","Recovery"),
                 levels = c("Pre-treatment","Mild stress","Moderate stress","Severe stress","Recovery")),
  x = c(-1.5, 1, 14, 21, 29), xend = c(1, 14, 21, 29, 34))
phase_cols <- c("Pre-treatment"="steelblue3","Mild stress"="palegreen",
                "Moderate stress"="yellow2","Severe stress"="tomato","Recovery"="skyblue")
phase_bar <- function(y) geom_segment(data = phases, inherit.aes = FALSE,
  aes(x = x, xend = xend, y = y, yend = y, colour = phase), linewidth = 2.2, show.legend = FALSE)

source("src/lib/fig_palette.R")
base_theme <- theme_paper(base_size = 9, major_y = TRUE) +
  theme(plot.tag = element_text(face = "bold", size = 12),
        legend.position = "none", axis.text = element_text(colour = "black"))

pl_panel <- function(d, ylab, ylim, ybreaks, bar_y) {
  ggplot(d, aes(Day, m, shape = Treatment, fill = Treatment, ymin = lo, ymax = hi)) +
    phase_bar(bar_y) +
    geom_line(aes(group = Treatment), linewidth = 0.6, colour = "black") +
    geom_errorbar(width = 0.7, linewidth = 0.35, colour = "black") +
    geom_point(size = 2.4, stroke = 0.6, colour = "black") +
    scale_shape_manual(values = c(Control = 24, Drought = 21)) +
    scale_fill_manual(values  = c(Control = "white", Drought = "black")) +
    scale_colour_manual(values = phase_cols) +
    scale_x_continuous("Days of experiment", breaks = seq(0,35,5), limits = c(-1.5,35)) +
    scale_y_continuous(ylab, limits = ylim, breaks = ybreaks) + base_theme
}

fc  <- summ(meta,  "FC_percent")
rws <- summ(meta,  "RWC_shoot")
rwr <- summ(meta,  "RWC_root")
pho <- summ(photo, "Photo_percent")
con <- summ(photo, "cond_percent")
wp  <- summ(meta,  "mPa") %>% filter(Treatment == "Drought")

pa <- pl_panel(fc,  "Soil water content [% FC]",        c(0,110),   seq(0,100,20),  105) + labs(tag="a")
pc <- pl_panel(rws, "Shoot RWC [%]",                    c(35,115),  seq(40,100,20), 112) + labs(tag="c")
pd <- pl_panel(rwr, "Root RWC [%]",                     c(0,90),    seq(0,80,20),    88) + labs(tag="d")
pe <- pl_panel(pho, expression(italic(A)[net]~"[%]"),  c(-60,240), seq(0,200,100), 232) + labs(tag="e")
pf <- pl_panel(con, expression(italic(g)[s]~"[%]"),    c(-30,240), seq(0,200,100), 232) + labs(tag="f")

pb <- ggplot(wp, aes(Day, m, ymin = lo, ymax = hi)) +
  phase_bar(0) +
  geom_col(width = 1.4, fill = "grey45", colour = "black", linewidth = 0.2) +
  geom_errorbar(width = 0.7, linewidth = 0.35, colour = "black") +
  scale_colour_manual(values = phase_cols) +
  scale_x_continuous("Days of experiment", breaks = seq(0,35,5), limits = c(-1.5,35)) +
  scale_y_continuous(expression(Psi[shoot]~"[MPa]"), limits = c(-1.7, 0.05)) +
  base_theme + labs(tag = "b")

# legend strip: stress-phase colour key + Control/Drought symbols
leg_phase <- ggplot(phases, aes(x, 1, colour = phase)) + geom_point(size = 3, alpha = 0) +
  scale_colour_manual(values = phase_cols, name = NULL) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(shape = 15, size = 5, alpha = 1))) +
  theme_void() + theme(legend.position = "top", legend.text = element_text(size = 8))
leg_trt <- ggplot(data.frame(Treatment = factor(c("Control","Drought"), c("Control","Drought"))),
                  aes(1, 1, shape = Treatment, fill = Treatment)) +
  geom_point(size = 3, colour = "black", alpha = 0) +
  scale_shape_manual(values = c(Control = 24, Drought = 21), name = NULL) +
  scale_fill_manual(values = c(Control = "white", Drought = "black"), name = NULL) +
  guides(shape = guide_legend(override.aes = list(alpha = 1)),
         fill  = guide_legend(override.aes = list(alpha = 1))) +
  theme_void() + theme(legend.position = "top", legend.text = element_text(size = 8))

body    <- (pa | pb) / (pc | pd) / (pe | pf)
legends <- (leg_phase | leg_trt) + plot_layout(widths = c(3, 1))
final   <- legends / body + plot_layout(heights = c(0.06, 1))

OUT <- "results/final_figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)  # ensure output dir exists in a fresh checkout
ggsave(file.path(OUT, "FigureS1.pdf"), final, width = 21/2.54, height = 26/2.54)
ggsave(file.path(OUT, "FigureS1.png"), final, width = 21/2.54, height = 26/2.54, dpi = 200)
cat("Figure S1 written to", OUT, "\n")
