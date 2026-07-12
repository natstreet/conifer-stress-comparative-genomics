#!/usr/bin/env Rscript
# figure6_chs3_case_study.R — Figure 6: CHS3 segmental-duplication gene-pair case study
# PA_chr09_G004115 and PA_chr09_G004116
# Both on minus strand, ~145 kb apart on chr9

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(readr)
})

# Run from the AbioticStressConifers/ project root (standalone or sourced by assemble_figures.R);
# getwd() is the project root in both cases. (commandArgs("--file=") is unreliable under source().)
PROJ <- normalizePath(getwd())
EXPR    <- file.path(PROJ, "data/expression")               # individual-sample VST matrices
PAVFILE <- file.path(PROJ, "data/chs3/chs3_pav_north_south.tsv")  # see file header for provenance
OUTDIR  <- file.path(PROJ, "results/integration/chs3")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

theme_fig <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        plot.title = element_text(size = 9, face = "bold"))

# ── Panel A: Genomic positions (verified from GFF3 only) ─────────────────────
# G004115: 1,096,275,812–1,096,278,064 (minus, 2 exons, 2.25 kb)
# G004116: 1,096,423,140–1,096,425,179 (mRNA.3, minus, 2.04 kb)
# Both on PA_chr09, same strand, ~145 kb apart

locus_df <- tibble(
  gene   = c("G004115", "G004116"),
  pos_mb = c(1096277, 1096424),   # midpoint in bp → label position
  label  = c("G004115\n(1,096,276-1,096,278 kb)", "G004116\n(1,096,423-1,096,425 kb)"),
  fill   = c("#D6604D", "#4393C3"),
  size_kb= c(2.3, 2.0)
)

p_a <- ggplot() +
  # chromosome line
  annotate("segment", x=1096270, xend=1096430, y=1, yend=1,
           colour="grey70", linewidth=1) +
  # gene rectangles
  annotate("rect", xmin=1096275.8, xmax=1096278.1, ymin=0.7, ymax=1.3,
           fill="#D6604D", colour="white") +
  annotate("rect", xmin=1096423.1, xmax=1096425.2, ymin=0.7, ymax=1.3,
           fill="#4393C3", colour="white") +
  # labels
  annotate("text", x=1096277, y=1.5, label="G004115\n2.3 kb",
           size=2.8, fontface="bold", colour="#D6604D") +
  annotate("text", x=1096424, y=1.5, label="G004116\n2.0 kb",
           size=2.8, fontface="bold", colour="#4393C3") +
  # distance bracket
  annotate("segment", x=1096278.1, xend=1096423.1, y=0.5, yend=0.5,
           colour="grey50", linewidth=0.5,
           arrow=arrow(ends="both", length=unit(0.15,"cm"), type="open")) +
  annotate("text", x=1096350, y=0.38, label="~145 kb", size=2.8, colour="grey40") +
  # GWAS label
  annotate("text", x=1096350, y=1.7,
           label="GWAS: isothermality (bio3)  |  both minus strand",
           size=2.5, colour="grey40", fontface="italic") +
  scale_x_continuous(name="PA_chr09 position (kb)", labels=scales::comma) +
  scale_y_continuous(limits=c(0.2, 2.0), breaks=NULL, name=NULL) +
  labs(title="a   Genomic locus (chr9)") +
  theme_fig

# ── Panel B: PAV north/south ──────────────────────────────────────────────────
# PAV counts read from the deposit data file
pav_raw <- readr::read_tsv(PAVFILE, comment = "#", show_col_types = FALSE)
pav_raw$Gene <- sub("^PA_chr09_", "", pav_raw$gene)
pav_df <- pav_raw %>%
  transmute(Gene,
            Population = sprintf("%s\n(n=%d)", population, n_total),
            Present = n_present, Absent = n_absent) %>%
  pivot_longer(c(Present,Absent), names_to="Status", values_to="n") %>%
  mutate(Status = factor(Status, levels=c("Absent","Present")))

# PAV summary string for the figure subtitle, derived from the same data file
pav_str <- function(g) {
  n <- pav_raw[pav_raw$Gene == g & pav_raw$population == "Northern", ]
  s <- pav_raw[pav_raw$Gene == g & pav_raw$population == "Southern", ]
  sprintf("%s present in %d/%d north, %d/%d south",
          g, n$n_present, n$n_total, s$n_present, s$n_total)
}
pav_text <- paste0(pav_str("G004115"), "  |  ", pav_str("G004116"))

p_b <- ggplot(pav_df, aes(x=Population, y=n, fill=Status)) +
  geom_col(width=0.65) +
  facet_wrap(~Gene, nrow=1) +
  scale_fill_manual(values=c("Present"="#4393C3","Absent"="#E0E0E0"), name=NULL) +
  scale_y_continuous(limits=c(0,28), breaks=c(0,13,26)) +
  labs(title="b   Presence-absence variation (PAV)",
       x=NULL, y="Individuals") +
  theme_fig +
  theme(legend.position="bottom", legend.key.size=unit(0.35,"cm"))

# ── Panel C: Expression ───────────────────────────────────────────────────────
# Derived from the deposit's combined cold (SC) and drought (SD) sample matrices,
# split into needle/root by the _N / _R column-name suffix.
genes_oi <- c("PA_chr09_G004115","PA_chr09_G004116")
ctx_defs <- list(
  "Cold\nneedles"    = list("SC_expression.txt", "_N$"),
  "Cold\nroots"      = list("SC_expression.txt", "_R$"),
  "Drought\nneedles" = list("SD_expression.txt", "_N$"),
  "Drought\nroots"   = list("SD_expression.txt", "_R$"))

expr_df <- lapply(names(ctx_defs), function(ctx) {
  m <- read_tsv(file.path(EXPR, ctx_defs[[ctx]][[1]]), show_col_types=FALSE) %>%
    filter(Genes %in% genes_oi) %>%
    column_to_rownames("Genes")
  m <- m[, grepl(ctx_defs[[ctx]][[2]], colnames(m)), drop=FALSE]
  tibble(
    Gene    = rownames(m),
    Context = ctx,
    mean_vst= rowMeans(m),
    se_vst  = apply(m, 1, sd) / sqrt(ncol(m))
  )
}) %>% bind_rows() %>%
  mutate(
    Context = factor(Context, levels=names(ctx_defs)),
    Gene    = recode(Gene,
                     "PA_chr09_G004115" = "G004115",
                     "PA_chr09_G004116" = "G004116")
  )

p_c <- ggplot(expr_df, aes(x=Context, y=mean_vst, fill=Gene)) +
  geom_col(position=position_dodge(0.7), width=0.65) +
  geom_errorbar(aes(ymin=mean_vst-se_vst, ymax=mean_vst+se_vst),
                position=position_dodge(0.7), width=0.25, linewidth=0.4) +
  scale_fill_manual(values=c("G004115"="#D6604D","G004116"="#4393C3"), name=NULL) +
  labs(title="c   Mean expression across stress conditions\n    (VST, P. abies seedlings)",
       x=NULL, y=expression("Mean VST" %+-% "SE")) +
  theme_fig +
  theme(legend.position="bottom", legend.key.size=unit(0.35,"cm"))

# ── CHS3 promoter-motif summary (one-line subtitle only; per-motif panel removed) ──
# The FIMO/PlantTFDB proximal-2 kb scan (src/ComPlEx/chs3_promoter_motifs.py, run FIRST)
# still supplies the compressed one-line CHS3 motif statement (the copies share most motif
# families). The former per-motif bar chart (old panel d) has been removed — Figure 6 is
# now three panels a/b/c.
motif_summ <- readr::read_tsv(file.path(OUTDIR, "chs3_promoter_motif_summary.tsv"),
                              show_col_types = FALSE)

# ── Assemble (a/b/c) ──────────────────────────────────────────────────────────
p_final <- (p_a | p_b) / p_c +
  plot_annotation(
    title    = "CHS3 SD gene pair: PA_chr09_G004115 / PA_chr09_G004116",
    subtitle = paste0(
      "dN/dS = ", motif_summ$kaks, " (purifying selection)  |  GWAS: isothermality (bio3)  |  ",
      "PAV: ", pav_text, "  |  ",
      "Promoter TF motifs (2 kb): ", motif_summ$jaccard_pct, "% shared (Jaccard)"
    ),
    theme = theme(
      plot.title    = element_text(face="bold", size=11),
      plot.subtitle = element_text(size=8, colour="grey30")
    )
  )

out <- file.path(OUTDIR, "chs3_sd_pair_figure")
ggsave(paste0(out,".pdf"), p_final, width=14, height=7.5, device="pdf")
ggsave(paste0(out,".png"), p_final, width=14, height=7.5, dpi=150)
cat("Saved", out, ".pdf/.png\n")
