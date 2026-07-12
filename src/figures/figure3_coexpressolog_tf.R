suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork)
  library(ComplexUpset); library(ggplotify)
})

WDIR <- getwd()  # Run from AbioticStressConifers/ project root
RDAT <- "results/ComPlEx/RData"
OUT  <- file.path(WDIR,"results/final_figures")
load(file.path(RDAT,"co_expressologs.RData"))
load(file.path(RDAT,"orthogroup_coexpressolog_presence.RData"))
load(file.path(RDAT,"weighted_gene_pairs.RData"))

ld<-function(c){f<-file.path(WDIR,paste0("data/DEG_lists/DE_all_",c,"_01_2L2FC.RData"));if(!file.exists(f))return(character(0));load(f);get(ls()[grep("DE_all",ls())])}

# ── PANEL a: UpSet ─────────────────────────────────────────────────────────────
oz <- as.data.frame(orthogroup_coexpressolog_presence)
SET_NAMES <- c("Cold\nNeedle","Cold\nRoot","Drought\nNeedle","Drought\nRoot")
names(oz) <- SET_NAMES
oz[] <- lapply(oz, as.logical)

# Keep only intersections with >= 1 member to reduce clutter
oz_nz <- oz[rowSums(oz) > 0, ]

raw_upset <- upset(
  oz_nz, SET_NAMES,
  stripes = "white",                 # house style: white matrix background
  # Intersection bar: rotate counts to avoid overlap
  base_annotations = list(
    "Intersection\nsize" = (
      intersection_size(
        text = list(size=2.3, angle=90, vjust=0.5, hjust=-0.1),
        bar_number_threshold=1
      ) +
      scale_y_continuous(expand=expansion(mult=c(0,0.35)),
                         labels=scales::comma) +
      theme(
        panel.grid    = element_blank(),
        axis.line.y   = element_line(colour="grey40"),
        axis.ticks.y  = element_line(colour="grey40"),
        axis.line.x   = element_blank(),
        axis.ticks.x  = element_blank()
      )
    )
  ),
  # Matrix dots: keep UpSet built-in contrast (dark = member, light grey = non-member).
  # (A previous fixed fill="grey25" painted every dot the same, so members and non-members
  # were indistinguishable and every column read as fully connected.)
  matrix = intersection_matrix(
    geom = geom_point(size=2.6)
  ) + theme(panel.grid=element_blank()),
  # Set size bars: no coord_flip so names stay on y-axis (no overlap)
  set_sizes = (
    upset_set_size(
      geom=geom_bar(fill="#555555", width=0.65)
    ) +
    scale_y_continuous(
      labels=scales::comma,
      expand=expansion(mult=c(0.05,0.2))
    ) +
    labs(y="Set size") +
    theme(
      panel.grid   = element_blank(),
      axis.line    = element_line(colour="grey40"),
      axis.ticks   = element_line(colour="grey40"),
      axis.text.y  = element_text(size=8.5),
      axis.text.x  = element_text(size=5.5, angle=45, hjust=1)
    )
  ),
  themes = upset_modify_themes(list(
    "intersections_matrix" = theme(
      text        = element_text(size=8.5),
      panel.grid  = element_blank(),
      axis.line.x = element_line(colour="grey40"),
      axis.ticks.x = element_line(colour="grey40")
    ),
    "overall_sizes" = theme(
      text        = element_text(size=8.5),
      panel.grid  = element_blank()
    )
  )),
  sort_sets        = FALSE,
  sort_intersections_by = "cardinality",
  height_ratio     = 0.35,
  width_ratio      = 0.20,
  min_size         = 10
)

p3a <- wrap_elements(full = raw_upset)
cat("Panel a built\n")

# ── PANEL b: TF families in all-condition conserved set ───────────────────────
c_og <- rownames(as.data.frame(orthogroup_coexpressolog_presence))[rowSums(as.data.frame(orthogroup_coexpressolog_presence))==4]
tf_pa <- read.table(
  "data/annotation/Picab02_230926_at01_longest_no_TE_aa_TF_predictions.tsv",
  sep="\t", header=TRUE, quote="", fill=TRUE
) |> mutate(gene_id=sub("\\.mRNA.*$","",Gene_ID)) |> select(gene_id,Family)

cons_tf <- co_expressologs |>
  filter(OrthoGroup %in% c_og,
         suppressWarnings(as.numeric(MaxpVal)) < 0.01) |>
  group_by(OrthoGroup) |>
  slice_min(suppressWarnings(as.numeric(MaxpVal)), n=1, with_ties=FALSE) |>
  ungroup() |>
  distinct(OrthoGroup, Species1) |>
  left_join(tf_pa, by=c("Species1"="gene_id")) |>
  filter(!is.na(Family))

fam_top <- c("WRKY","NAC","MYB","ERF","bHLH","bZIP","MYB_related",
             "MIKC_MADS","C2H2","Trihelix","G2-like","LBD","HSF","SBP","TCP")
fam_pal <- c(WRKY="#0072B2",NAC="#009E73",MYB="#E69F00",ERF="#D55E00",
             bHLH="#CC79A7",bZIP="#56B4E9",MYB_related="#F0E442",
             MIKC_MADS="#7B2D8B",C2H2="#999999",Trihelix="#E97DBD",
             "G2-like"="#A6611A",LBD="#543005",HSF="#BF812D",
             SBP="#80CDC1",TCP="#DFC27D",Other="#CCCCCC")

fam_cnt <- cons_tf |>
  mutate(Fam2=ifelse(Family %in% fam_top, Family, "Other")) |>
  count(Fam2) |> arrange(desc(n)) |>
  mutate(Fam2=factor(Fam2, rev(unique(Fam2))))

p3b <- ggplot(fam_cnt, aes(Fam2, n, fill=Fam2)) +
  geom_col(width=0.72, alpha=0.9) +
  geom_text(aes(label=n), hjust=-0.15, size=2.8) +
  scale_fill_manual(values=fam_pal, guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,0.28))) +
  coord_flip() +
  labs(x=NULL, y="Orthogroups (conserved\nacross all 4 comparisons)") +
  theme_classic(base_size=9) +
  theme(panel.grid=element_blank(),
        axis.text.y=element_text(size=8.5),
        axis.text.x=element_text(size=7.5))
cat("Panel b built\n")

# ── PANEL c: TF up/down — fixed scale, legend, labels ────────────────────────
pa_d<-unique(c(ld("SDN"),ld("SDR"))); ps_d<-unique(c(ld("PDN"),ld("PDR")))
pa_c<-unique(c(ld("SCN"),ld("SCR"))); ps_c<-unique(c(ld("PCN"),ld("PCR")))
drought_og <- weighted_gene_pairs|>filter(drought_specific==TRUE,Species1%in%pa_d,Species2%in%ps_d)|>group_by(OrthoGroup)|>slice_max(DroughtSum,n=1,with_ties=FALSE)|>ungroup()|>filter(DroughtSum>=60)|>pull(OrthoGroup)
cold_og    <- weighted_gene_pairs|>filter(cold_specific==TRUE,Species1%in%pa_c,Species2%in%ps_c)|>group_by(OrthoGroup)|>slice_max(ColdSum,n=1,with_ties=FALSE)|>ungroup()|>filter(ColdSum>=40)|>pull(OrthoGroup)

get_tf<-function(og,pa_de,ps_de,tissues,label){
  co_expressologs|>mutate(pv=suppressWarnings(as.numeric(MaxpVal)))|>
    filter(OrthoGroup%in%og,!is.na(pv),pv<0.01,Tissue%in%tissues,
           Species1%in%pa_de,Species2%in%ps_de)|>
    group_by(OrthoGroup)|>slice_min(pv,n=1,with_ties=FALSE)|>ungroup()|>
    distinct(OrthoGroup,Species1)|>left_join(tf_pa,by=c("Species1"="gene_id"))|>
    filter(!is.na(Family))|>mutate(stress=label)
}
add_dir<-function(df,mat,ctrl_pat,trt_pat){
  g<-intersect(df$Species1,rownames(mat))
  ctrl<-grep(ctrl_pat,names(mat),value=TRUE); trt<-grep(trt_pat,names(mat),value=TRUE)
  if(!length(g)||!length(ctrl)||!length(trt))return(df|>mutate(direction=NA))
  dir<-data.frame(Species1=g,direction=ifelse(rowMeans(mat[g,trt,drop=FALSE],na.rm=TRUE)>
    rowMeans(mat[g,ctrl,drop=FALSE],na.rm=TRUE),"Up","Down"))
  left_join(df,dir,by="Species1")
}
med_sd<-read.table(file.path(WDIR,"data/expression/expr_median/SDmed_expression.txt"),header=TRUE,sep="\t",check.names=FALSE); rownames(med_sd)<-med_sd[,1]; med_sd<-med_sd[,-c(1,2)]
med_sc<-read.table(file.path(WDIR,"data/expression/expr_median/SCmed_expression.txt"),header=TRUE,sep="\t",check.names=FALSE); rownames(med_sc)<-med_sc[,1]; med_sc<-med_sc[,-c(1,2)]

d_tf<-get_tf(drought_og,pa_d,ps_d,c("drought_needle","drought_root"),"Drought")|>add_dir(med_sd,"FC80","FC[0-9]|C1d|C2d|Reh")|>filter(!is.na(direction))|>mutate(Fam2=ifelse(Family%in%fam_top,Family,"Other"))
c_tf<-get_tf(cold_og,pa_c,ps_c,c("cold_needle","cold_root"),"Cold")|>add_dir(med_sc,"20C_0h","5C|neg5C")|>filter(!is.na(direction))|>mutate(Fam2=ifelse(Family%in%fam_top,Family,"Other"))
all_tf <- bind_rows(d_tf,c_tf)|>mutate(stress=factor(stress,c("Cold","Drought")))
cat("TF totals: Cold =",sum(all_tf$stress=="Cold"),"Drought =",sum(all_tf$stress=="Drought"),"\n")

# Only include TF families actually present in the data
pres_fams <- unique(all_tf$Fam2)
fam_pal_sub <- fam_pal[names(fam_pal) %in% pres_fams]
fam_sum <- all_tf|>count(stress,direction,Fam2)|>
  mutate(n_sign=ifelse(direction=="Down",-n,n),
         Fam2=factor(Fam2, names(fam_pal_sub)))

# Correct limits: based on TOTAL stacked height per stress×direction
tot_h <- all_tf|>count(stress,direction)
max_up   <- max(tot_h$n[tot_h$direction=="Up"])
max_down <- max(tot_h$n[tot_h$direction=="Down"])
y_max <- max_up   * 1.25
y_min <- -max_down * 1.25
cat(sprintf("Scale: [%.0f, %.0f]\n", y_min, y_max))

p3c <- ggplot(fam_sum, aes(stress, n_sign, fill=Fam2)) +
  geom_col(position="stack", width=0.68,
           colour="white", linewidth=0.2) +
  geom_hline(yintercept=0, colour="grey25", linewidth=0.6) +
  # Up / Down axis side labels — placed on y-axis as annotations outside clip
  scale_fill_manual(values=fam_pal_sub, name="TF family",
                    guide=guide_legend(ncol=1, reverse=FALSE)) +
  scale_y_continuous(
    labels=abs,
    name="Number of co-expressologs
← Down  |  Up →",
    limits=c(y_min, y_max),
    breaks=seq(-floor(max_down), floor(max_up), by=ceiling(max_up/4))
  ) +
  # Secondary axis to label Up / Down
  scale_x_discrete(name=NULL) +

  coord_flip(clip="off") +
  theme_classic(base_size=9) +
  theme(
    panel.grid   = element_blank(),
    axis.line.x  = element_line(colour="grey40"),
    axis.ticks.x = element_line(colour="grey40"),
    axis.line.y  = element_line(colour="grey40"),
    legend.key.size = unit(0.28,"cm"),
    legend.text  = element_text(size=7.5),
    legend.title = element_text(size=8, face="bold"),
    plot.margin  = margin(5, 5, 5, 8)
  )
cat("Panel c built\n")

# ── ASSEMBLE with correct tags ─────────────────────────────────────────────────
right_col <- p3b / p3c + plot_layout(heights=c(1,1))

fig3 <- (p3a | right_col) +
  plot_layout(widths=c(1.5, 1)) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag=element_text(size=11, face="bold"))
  )

# Key numbers, printed for traceability (top UpSet intersection + panel-b families)
{
  ozc <- oz_nz; combo <- apply(ozc, 1, function(r) paste(names(ozc)[r], collapse="+"))
  cat("Fig3 panel a: top UpSet intersections:\n"); print(head(sort(table(combo), decreasing=TRUE), 6))
  cat("Fig3 panel b: conserved (all-4) TF-family counts:\n"); print(fam_cnt)
}

pdf(file.path(OUT,"Figure3.pdf"), width=22/2.54, height=14/2.54)
print(fig3)
dev.off()
cat("Figure3.pdf saved\n")
# Also write a matching PNG from the same figure object so the two never drift.
ggsave(file.path(OUT,"Figure3.png"), fig3, width=22/2.54, height=14/2.54, dpi=300)
cat("Figure3.png saved\n")
