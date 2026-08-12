# Manuscript section: Exploratory duration-based longitudinal analysis
# Analysis family: figure generation
# Original source path: adapted from the sliding-window eyes plot
# Primary input dataset(s): 00_data/derived/preprocessing/duration_sliding_window_data.csv
# Primary output(s): supplementary normalized sliding-window background figures
# Known TODOs: compare visually against the face/background collapsed version if this becomes manuscript-relevant
# Scientific logic note: plots bounded within-window fixation-duration proportions over overlapping normalized interview progress for the background AOI

library(tidyverse)
library(ggh4x)
library(here)
source(here("01_analysis", "04_figures", "shared", "custom_plot_settings.R"))
source(here("01_analysis", "03_supplementary", "10_fit_sliding_window_background_models.R"))

medication_levels <- c("BOTH", "NAL", "OXT", "PLA")
interval_colors <- c("ASD_gt_CTRL" = "#d95f02", "ASD_lt_CTRL" = "#1b9e77")

background_group_stats <- background_group_stats %>%
  mutate(medication = factor(medication, levels = medication_levels))

background_group_stats_by_sex <- background_group_stats_by_sex %>%
  mutate(
    medication = factor(medication, levels = medication_levels),
    Sex = factor(Sex, levels = c("f", "m"), labels = c("Female", "Male"))
  )

background_intervals <- background_intervals %>%
  mutate(
    medication = factor(medication, levels = medication_levels),
    band_color = unname(interval_colors[direction])
  )

background_intervals_by_sex <- background_intervals_by_sex %>%
  mutate(
    medication = factor(medication, levels = medication_levels),
    Sex = factor(Sex, levels = c("f", "m"), labels = c("Female", "Male")),
    band_color = unname(interval_colors[direction])
  )

interval_positions <- background_group_stats %>%
  group_by(medication) %>%
  summarise(interval_y = max(upper_ci, na.rm = TRUE) + 0.03, .groups = "drop")

interval_positions_by_sex <- background_group_stats_by_sex %>%
  group_by(Sex, medication) %>%
  summarise(interval_y = max(upper_ci, na.rm = TRUE) + 0.03, .groups = "drop")

background_intervals <- background_intervals %>%
  left_join(interval_positions, by = "medication")

background_intervals_by_sex <- background_intervals_by_sex %>%
  left_join(interval_positions_by_sex, by = c("Sex", "medication"))

background_group_stats_combined <- bind_rows(
  background_group_stats %>%
    mutate(Panel = factor("Overall", levels = c("Overall", "Female", "Male"))),
  background_group_stats_by_sex %>%
    mutate(Panel = factor(as.character(Sex), levels = c("Overall", "Female", "Male")))
)

background_intervals_combined <- bind_rows(
  background_intervals %>%
    mutate(Panel = factor("Overall", levels = c("Overall", "Female", "Male"))),
  background_intervals_by_sex %>%
    mutate(Panel = factor(as.character(Sex), levels = c("Overall", "Female", "Male")))
)

combined_interval_positions <- background_group_stats_combined %>%
  group_by(Panel, medication) %>%
  summarise(interval_y = max(upper_ci, na.rm = TRUE) + 0.03, .groups = "drop")

background_intervals_combined <- background_intervals_combined %>%
  select(-interval_y) %>%
  left_join(combined_interval_positions, by = c("Panel", "medication"))

background_plot_combined <- ggplot(
  background_group_stats_combined,
  aes(x = progress_midpoint, y = mean_window_fix_duration_prop, group = Group, color = Group)
) +
  geom_line(linewidth = 0.35) +
  geom_ribbon(
    aes(
      ymin = mean_window_fix_duration_prop - sem,
      ymax = mean_window_fix_duration_prop + sem,
      fill = Group
    ),
    alpha = 0.5,
    linewidth = 0.2
  ) +
  geom_segment(
    data = background_intervals_combined,
    aes(x = interval_start, xend = interval_end, y = interval_y, yend = interval_y),
    inherit.aes = FALSE,
    color = background_intervals_combined$band_color,
    linewidth = 2.2,
    lineend = "round"
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = c("0%", "25%", "50%", "75%", "100%")
  ) +
  labs(x = "Normalized Interview Progress", y = "Fixation Duration Proportion on Background") +
  theme_minimal() +
  theme(legend.position = "top") +
  facet_grid2(Panel ~ medication, axes = "all", remove_labels = "all") +
  apply_custom_settings()

ggsave(
  here("02_outputs", "figures", "supplementary", "timecourse", "svg", "duration_sliding_window_background_plot_combined.svg"),
  background_plot_combined,
  units = "mm",
  width = 183,
  height = 190
)
ggsave(
  here("02_outputs", "figures", "supplementary", "timecourse", "png", "Supp_Figure_duration_sliding_window_background_combined.png"),
  background_plot_combined,
  units = "mm",
  width = 183,
  height = 190
)
