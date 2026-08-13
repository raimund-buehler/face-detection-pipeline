# Manuscript section: Main Figure 5
# Analysis family: cortisol descriptive figure
# Primary input dataset(s): 00_data/derived/analysis/df_cortisol_merged.csv
# Primary output(s): raw cortisol timecourse by group and medication

library(tidyverse)
library(here)

source(here("01_analysis", "04_figures", "shared", "custom_theme_settings_publication.R"))

output_dir <- here("02_outputs", "figures", "main", "current_story")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_png <- file.path(output_dir, "Figure_5_cortisol_raw_timecourse.png")
output_svg <- file.path(output_dir, "Figure_5_cortisol_raw_timecourse.svg")
output_alias_png <- file.path(output_dir, "Figure_5_cortisol.png")
output_alias_svg <- file.path(output_dir, "Figure_5_cortisol.svg")

medication_levels_raw <- c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA")
medication_levels_plot <- c("PLA", "OXT", "NAL", "BOTH")
medication_values <- c(
  "PLA" = "#0B7A6B",
  "OXT" = "#93AD7C",
  "NAL" = "#C7B655",
  "BOTH" = "#FDC010"
)

raw_cortisol <- read_csv(
  here("00_data", "derived", "analysis", "df_cortisol_merged.csv"),
  show_col_types = FALSE
) %>%
  filter(fix == "Eyes", !is.na(cortisol_1), !is.na(medication), !is.na(Group)) %>%
  mutate(
    medication = factor(
      medication,
      levels = medication_levels_raw,
      labels = c("BOTH", "NAL", "OXT", "PLA")
    ),
    medication = factor(medication, levels = medication_levels_plot),
    Group = factor(Group, levels = c("ASD", "CTRL")),
    timepoint = factor(
      timepoint,
      levels = c(1, 2, 3),
      labels = c("T1", "T2", "T3")
    )
  )

raw_summary <- raw_cortisol %>%
  group_by(Group, medication, timepoint) %>%
  summarise(
    mean_cortisol = mean(cortisol_1, na.rm = TRUE),
    sem = sd(cortisol_1, na.rm = TRUE) / sqrt(sum(!is.na(cortisol_1))),
    n = sum(!is.na(cortisol_1)),
    .groups = "drop"
  )

figure_5 <- ggplot(
  raw_summary,
  aes(x = timepoint, y = mean_cortisol, color = medication, group = medication)
) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.7) +
  geom_errorbar(
    aes(ymin = mean_cortisol - sem, ymax = mean_cortisol + sem),
    width = 0.08,
    linewidth = 0.35
  ) +
  facet_wrap(~ Group, nrow = 1) +
  scale_color_manual(values = medication_values, breaks = medication_levels_plot) +
  labs(
    title = "Cortisol",
    x = "Timepoint",
    y = "Cortisol (nmol/l)",
    color = "Medication"
  ) +
  publication_panel_theme(base_size = 8) +
  publication_legend_theme(position = "top", text_size = 8) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(b = 2)),
    legend.justification = "center",
    legend.key.width = grid::unit(9, "pt"),
    legend.spacing.x = grid::unit(2, "pt"),
    strip.text = element_text(size = 8),
    panel.spacing.x = grid::unit(10, "pt"),
    axis.title.x = element_text(margin = margin(t = 4)),
    axis.title.y = element_text(margin = margin(r = 4))
  )

ggsave(output_png, figure_5, units = "mm", width = 120, height = 70, dpi = 300)
ggsave(output_svg, figure_5, units = "mm", width = 120, height = 70)
ggsave(output_alias_png, figure_5, units = "mm", width = 120, height = 70, dpi = 300)
ggsave(output_alias_svg, figure_5, units = "mm", width = 120, height = 70)
