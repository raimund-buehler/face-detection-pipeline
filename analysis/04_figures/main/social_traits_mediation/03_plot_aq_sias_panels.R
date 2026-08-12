# Manuscript section: Dimensional analyses
# Analysis family: figure generation
# Original source path: scripts/publication/analysis/AQ_SIAS.r
# Primary input dataset(s): data prepared by social_traits_mediation/shared_aq_sias_plot_utils.R from manuscript analysis tables
# Primary output(s): outputs/figures/main/Figure_4.png and the AQ/SIAS panel figures
# Known TODOs: helper still includes plotting-specific logic and mediation grob utilities
# Scientific logic note: figure assembly was separated from the original mixed script without changing figure behavior

library(here)
source(here("01_analysis", "04_figures", "shared", "custom_plot_settings.R"))
source(here("01_analysis", "02_main_manuscript", "social_traits_mediation", "shared_aq_sias_plot_utils.R"))

aq_sias_components <- build_aq_sias_panels(
  base_size = 8,
  axis_text_size = 7,
  axis_title_size = 8,
  title_size = 9,
  legend_text_size = 7,
  strip_text_size = 7,
  annotation_size = 3.2
)

compact_figure <- assemble_compact_aq_sias_mediation(
  panels = aq_sias_components$plots,
  base_size = 8
)

compact_figure

ggsave(here("02_outputs", "figures", "main", "Figure_4.png"), compact_figure, units = "mm", width = 190, height = 118)
