## fig_palette.R — paper-wide colour palette + house theme for van Zalen et al.
## Source this from every figure script so colours and styling are identical
## across all manuscript figures:  source("src/lib/fig_palette.R")
##
## Palette rationale:
##   Species colours are a NEW convention — spruce and pine had no prior colour
##   because they were always shown in separate panels; Figure 2c and the timing
##   supplementary are the first figures showing both together.
suppressPackageStartupMessages({ library(ggplot2) })

PAL <- list(
  ## species (new convention)
  spruce        = "#762A83",  # purple
  pine          = "#35978F",  # teal
  ## stress (from Figure 4)
  cold          = "#4D90D6",  # blue
  drought       = "#D6A44D",  # gold
  ## organ (from the heatmap Tissue annotation)
  needle        = "#2CA25F",  # green
  root          = "#8C2D04",  # brown
  ## "shared across both stresses" (from Figure 4)
  stress_shared = "#6BAF6B",
  ## neutral greys
  shared_dark   = "#404040",  # Fig 2c shared bar
  specific_grey = "#dcdcdc"   # species-specific portion (Fig 2c / timing stacks)
)

## Convenience named vectors for scale_*_manual()
PAL_SPECIES <- c(spruce = PAL$spruce, pine = PAL$pine)
PAL_STRESS  <- c(cold = PAL$cold, drought = PAL$drought)
PAL_ORGAN   <- c(needle = PAL$needle, root = PAL$root)

## Figure 2c bar fills: spruce-specific / shared / pine-specific
PAL_FIG2C <- c(`spruce-specific` = PAL$spruce,
               shared            = PAL$shared_dark,
               `pine-specific`   = PAL$pine)

## Timing stacked-bar fills: in-shared-OG (species colour) vs species-specific (grey)
PAL_TIMING <- c(`spruce shared` = PAL$spruce,
                `pine shared`   = PAL$pine,
                `species-specific` = PAL$specific_grey)

## House theme — white background, journal style (Addendum 2 of the task spec):
##   white panel, thin dark axes, NO gridlines by default. Set major_y=TRUE for a
##   single faint major-y gridline where a numeric read-off genuinely helps
##   (e.g. overlap-count bars, dN/dS breadth). No minor gridlines, ever.
theme_paper <- function(base_size = 11, base_family = "", major_y = FALSE) {
  th <- theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      axis.line        = element_line(colour = "grey20", linewidth = 0.4),
      axis.ticks       = element_line(colour = "grey20", linewidth = 0.4),
      axis.text        = element_text(colour = "grey20"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white", colour = NA),
      strip.text       = element_text(face = "bold", size = base_size),
      plot.title       = element_text(face = "bold", size = base_size + 1),
      legend.key       = element_rect(fill = "white", colour = NA)
    )
  if (major_y) {
    th <- th + theme(panel.grid.major.y =
                       element_line(colour = "grey88", linewidth = 0.3))
  }
  th
}
