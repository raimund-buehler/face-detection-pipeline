# Manuscript section: Supplementary / control analyses
# Analysis family: menstrual cycle analysis
# Original source path: scripts/publication/analysis/menst_analysis.r
# Primary input dataset(s): data/derived/analysis/df_analysis_pub_prep.csv
# Primary output(s): outputs/figures/supplementary/menst_phases.svg; outputs/figures/supplementary/Supp_Figure_1.png
# Known TODOs: phase provenance depends on the external menstrual-prep step
# Scientific logic note: plotting logic was separated from the original mixed script without changing figure behavior

library(tidyverse)
library(here)
source(here("01_analysis", "04_figures", "shared", "custom_plot_settings.R"))

df <- read_csv(here("00_data", "derived", "analysis", "df_analysis_pub_prep.csv"))

df$fix <- factor(
  df$fix,
  levels = c("Eyes", "Mouth", "Face", "Background"),
  labels = c("Eyes", "Mouth", "Face", "Background"),
  ordered = TRUE
)

df_menst <- df %>%
  ungroup() %>%
  filter(phase != "Missing")

summary_df <- df_menst %>%
  group_by(Group, fix, phase) %>%
  summarise(
    mean_var = mean(sqrt_rate_all, na.rm = TRUE),
    se = sd(sqrt_rate_all, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

menst_phase_plot <- ggplot(summary_df, aes(x = phase, y = mean_var, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(
    aes(ymin = mean_var - se, ymax = mean_var + se),
    position = position_dodge(0.9),
    width = 0.2
  ) +
  facet_wrap(~fix, scales = "fixed") +
  theme_bw() +
  labs(y = "Fixation Rate (Hz)", x = "Phase") +
  theme(legend.position = "top") +
  apply_custom_settings()

menst_phase_plot

ggsave(
  here("02_outputs", "figures", "supplementary", "menstrual_cycle", "svg", "menst_phases.svg"),
  menst_phase_plot,
  units = "mm",
  width = 89,
  height = 89
)
ggsave(
  here("02_outputs", "figures", "supplementary", "menstrual_cycle", "png", "Supp_Figure_1.png"),
  menst_phase_plot,
  width = 8.27,
  height = 8.27 / 1.5
)
