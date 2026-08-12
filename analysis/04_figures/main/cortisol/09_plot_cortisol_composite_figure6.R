# Manuscript section: Cortisol analyses
# Analysis family: manuscript figure assembly
# Original source path: legacy/cortisol_reports/markdown_source_bundle/4_cortisol_maxmin.Rmd
# Primary input dataset(s): cortisol manuscript analysis tables under data/derived/analysis/
# Primary output(s): outputs/figures/main/Figure_6.pdf
# Known TODOs: this follows the wide publication layout retained in the markdown workflow; alternate legacy variants remain archived
# Scientific logic note: figure assembly sources the topical cortisol plot scripts so panel edits stay local to those files

library(patchwork)
library(here)
source(here("01_analysis", "04_figures", "shared", "custom_theme_settings_publication.R"))
source(here("01_analysis", "04_figures", "main", "cortisol", "05_plot_raw_cortisol.R"))
source(here("01_analysis", "04_figures", "main", "cortisol", "07_plot_auc_predictor_effects.R"))
source(here("01_analysis", "04_figures", "main", "cortisol", "08_plot_t2t3_effects.R"))

mm_to_in <- function(mm) {
  mm / 25.4
}

style_publication_panel <- function(p, drop_guides = FALSE) {
  p <- p +
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
  if (drop_guides) {
    p <- p + theme(legend.position = "none")
  }
  p
}

strip_all_guides <- function(p, keep_guides = character()) {
  if ("fill" %in% keep_guides && !"colour" %in% keep_guides) {
    return(p + guides(colour = "none"))
  }
  if ("colour" %in% keep_guides && !"fill" %in% keep_guides) {
    return(p + guides(fill = "none"))
  }
  if (length(keep_guides) == 0) {
    return(p + theme(legend.position = "none"))
  }
  p
}

raw_plot_clean <- style_publication_panel(raw_plot, drop_guides = TRUE)
reactivity_barplot_clean <- style_publication_panel(reactivity_barplot, drop_guides = TRUE)
interaction_t2t3_clean <- style_publication_panel(interaction_t2t3_plot, drop_guides = TRUE)

kill_layer_legends <- function(p) {
  p$layers <- lapply(p$layers, function(layer) {
    layer$show.legend <- FALSE
    layer
  })
  p
}

interaction_auci_keep <- interaction_auci_plot + labs(title = "AUCi and Eye Gaze")
interaction_auci_keep <- kill_layer_legends(interaction_auci_keep)
interaction_auci_keep <- style_publication_panel(interaction_auci_keep)
interaction_auci_keep <- interaction_auci_keep +
  aes(fill = Group) +
  publication_palette(aesthetics = "fill")

interaction_auci_keep <- interaction_auci_keep +
  geom_point(
    data = data.frame(Group = c("ASD", "CTRL"), x = 0, y = 0),
    mapping = aes(x = x, y = y, fill = Group),
    inherit.aes = FALSE,
    shape = 22,
    size = 0,
    alpha = 0,
    show.legend = TRUE
  )

interaction_auci_keep <- strip_all_guides(interaction_auci_keep, keep_guides = "fill")

figure_6 <- (raw_plot_clean | reactivity_barplot_clean | interaction_auci_keep | interaction_t2t3_clean) +
  plot_layout(nrow = 1, guides = "collect") +
  plot_annotation(tag_levels = "A") +
  plot_annotation(theme = publication_legend_theme())

figure_6

ggsave(
  filename = here("02_outputs", "figures", "main", "Figure_6.pdf"),
  plot = figure_6,
  width = mm_to_in(190),
  height = mm_to_in(70),
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)
