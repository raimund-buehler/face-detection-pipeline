# Manuscript section: Timecourse and question-type analysis
# Analysis family: figure generation
# Original source path: scripts/publication/analysis/within_questions.R
# Primary input dataset(s): data/derived/analysis/df_analysis_pub.csv
# Primary output(s): outputs/figures/main/timecourse_plot.svg; outputs/figures/main/Figure_3.png
# Known TODOs: significance stars are recomputed from the label model summary object exposed by the model script
# Scientific logic note: plotting logic was separated from the original mixed script without changing figure behavior

library(tidyverse)
library(ggh4x)
library(here)
source(here("01_analysis", "04_figures", "shared", "custom_plot_settings.R"))
source(here("01_analysis", "02_main_manuscript", "eye_tracking", "02_fit_eye_question_timecourse_model.R"))

medication_levels <- c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA")
medication_labels <- c("BOTH", "NAL", "OXT", "PLA")

timecourse_group_stats <- timecourse_group_stats %>%
  mutate(medication = factor(medication, levels = medication_levels, labels = medication_labels))

timecourse_group_stats_by_sex <- df %>%
  group_by(Group, Sex, medication, label) %>%
  summarise(
    mean_sqrt_rate_q = mean(sqrt_rate_q, na.rm = TRUE),
    sem = sd(sqrt_rate_q, na.rm = TRUE) / sqrt(n()),
    lower_ci = mean_sqrt_rate_q - qt(0.975, df = n() - 1) * sem,
    upper_ci = mean_sqrt_rate_q + qt(0.975, df = n() - 1) * sem,
    .groups = "drop"
  ) %>%
  mutate(medication = factor(medication, levels = medication_levels, labels = medication_labels))

timecourse_plot <- ggplot(
  timecourse_group_stats,
  aes(x = label, y = mean_sqrt_rate_q, group = Group, color = Group)
) +
  geom_line(size = 0.25) +
  geom_point(size = 1) +
  geom_ribbon(
    aes(
      ymin = mean_sqrt_rate_q - sem,
      ymax = mean_sqrt_rate_q + sem,
      fill = Group
    ),
    alpha = 0.4,
    size = 0.25
  ) +
  geom_text(
    aes(x = label, y = 1.1, label = stars),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 4.5
  ) +
  theme_minimal() +
  labs(x = "Question", y = "Fixations on Eyes / Second") +
  theme(legend.position = "top") +
  facet_wrap2(~medication, axes = "all", remove_labels = "all") +
  apply_custom_settings()

timecourse_plot

ggsave(here("02_outputs", "figures", "main", "timecourse_plot.svg"), timecourse_plot, units = "mm", width = 183)
ggsave(here("02_outputs", "figures", "main", "Figure_3.png"), timecourse_plot, units = "mm", width = 183, height = 120)

timecourse_plot_by_sex <- ggplot(
  timecourse_group_stats_by_sex,
  aes(x = label, y = mean_sqrt_rate_q, group = Group, color = Group)
) +
  geom_line(size = 0.25) +
  geom_point(size = 0.9) +
  geom_ribbon(
    aes(
      ymin = mean_sqrt_rate_q - sem,
      ymax = mean_sqrt_rate_q + sem,
      fill = Group
    ),
    alpha = 0.35,
    size = 0.2
  ) +
  labs(x = "Question", y = "Fixations on Eyes / Second") +
  theme_minimal() +
  theme(legend.position = "top") +
  facet_grid2(Sex ~ medication, axes = "all", remove_labels = "all") +
  apply_custom_settings()

timecourse_plot_by_sex

ggsave(
  here("02_outputs", "figures", "supplementary", "timecourse", "svg", "timecourse_plot_by_sex.svg"),
  timecourse_plot_by_sex,
  units = "mm",
  width = 183,
  height = 150
)
ggsave(
  here("02_outputs", "figures", "supplementary", "timecourse", "png", "Supp_Figure_timecourse_by_sex.png"),
  timecourse_plot_by_sex,
  units = "mm",
  width = 183,
  height = 150
)
