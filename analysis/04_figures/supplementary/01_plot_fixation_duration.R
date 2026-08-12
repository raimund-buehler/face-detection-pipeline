# Manuscript section: Supplementary / control analyses
# Analysis family: fixation duration alternative DV
# Original source path: scripts/publication/analysis/duration_analysis.r
# Primary input dataset(s): data/derived/preprocessing/duration_data.csv
# Primary output(s): outputs/figures/supplementary/Figure_Supp1_Fix_Duration.png
# Known TODOs: significance annotations remain hard-coded from the source figure script
# Scientific logic note: plotting logic was separated from the original mixed script without changing figure behavior

library(tidyverse)
library(ggdist)
library(ggh4x)
library(ggsignif)
library(here)
source(here("01_analysis", "04_figures", "shared", "custom_plot_settings.R"))

df <- read_csv(here("00_data", "derived", "preprocessing", "duration_data.csv"))

df$fix <- factor(
  df$fix,
  levels = c("fix_on_eyes", "fix_on_mouth", "fix_on_face", "fix_on_background"),
  labels = c("Eyes", "Mouth", "Face", "Background"),
  ordered = TRUE
)

df$medication <- factor(
  df$medication,
  levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA"),
  labels = c("BOTH", "NAL", "OXT", "PLA")
)

bar_plot <- function(data, var) {
  var <- rlang::sym(var)

  df_summary <- data %>%
    group_by(Group, fix, medication) %>%
    summarise(
      mean_rate_all = mean(!!var),
      sem = sd(!!var) / sqrt(n()),
      .groups = "drop"
    )

  ggplot(df_summary, aes(x = Group, y = mean_rate_all, fill = medication)) +
    geom_bar(stat = "identity", position = position_dodge(0.85), width = 0.4) +
    geom_errorbar(
      aes(ymin = mean_rate_all - sem, ymax = mean_rate_all + sem),
      width = 0.1,
      position = position_dodge(0.85)
    ) +
    stat_halfeye(
      data = data,
      aes(x = Group, y = !!var, fill = medication),
      adjust = 0.6,
      position = position_dodge(0.85),
      justification = -1.1,
      .width = 0,
      point_color = NA,
      alpha = 0.4,
      normalize = "groups",
      scale = 0.25
    ) +
    labs(y = "% Fixation Duration") +
    facet_wrap2(~fix, scales = "fixed", axes = "all", remove_labels = "all")
}

df_summary <- df %>%
  group_by(Group, fix, medication) %>%
  summarise(
    mean_rate_all = mean(percentage_fix_duration),
    sem = sd(percentage_fix_duration) / sqrt(n()),
    .groups = "drop"
  )

y_position_main <- max(df_summary %>% filter(fix == "Eyes") %>% pull(mean_rate_all)) + 0.08

add_significance_lines <- function(summary_df, y_position_main, aoi, annotation) {
  y_position_branch <- y_position_main - 0.03

  list(
    geom_signif(
      data = summary_df %>% filter(fix == aoi),
      aes(xmin = "ASD", xmax = "CTRL", annotations = annotation),
      y_position = y_position_main,
      tip_length = 0.05,
      textsize = 5
    ),
    geom_signif(
      data = summary_df %>% filter(fix == aoi, Group == "ASD"),
      aes(xmin = 0.6, xmax = 1.4, annotations = ""),
      y_position = y_position_branch,
      tip_length = 0.03,
      textsize = 4,
      manual = TRUE
    ),
    geom_signif(
      data = summary_df %>% filter(fix == aoi, Group == "CTRL"),
      aes(xmin = 1.6, xmax = 2.4, annotations = ""),
      y_position = y_position_branch,
      tip_length = 0.03,
      textsize = 4,
      manual = TRUE
    )
  )
}

group_by_drug <- bar_plot(df, "percentage_fix_duration") +
  apply_custom_settings(values = c("#FDC010", "#C7B655", "#93AD7C", "#0B7A6B")) +
  add_significance_lines(df_summary, y_position_main, "Eyes", "***") +
  add_significance_lines(df_summary, y_position_main, "Background", "***")

group_by_drug

ggsave(
  here("02_outputs", "figures", "supplementary", "fixation_duration", "png", "Figure_Supp1_Fix_Duration.png"),
  group_by_drug,
  units = "mm",
  width = 183,
  height = 120
)

df_summary_by_sex <- df %>%
  filter(Sex %in% c("f", "m")) %>%
  group_by(Sex, Group, fix, medication) %>%
  summarise(
    mean_rate_all = mean(percentage_fix_duration),
    sem = sd(percentage_fix_duration) / sqrt(n()),
    .groups = "drop"
  )

group_by_drug_by_sex <- ggplot(
  df_summary_by_sex,
  aes(x = Group, y = mean_rate_all, fill = medication)
) +
  geom_bar(stat = "identity", position = position_dodge(0.85), width = 0.4) +
  geom_errorbar(
    aes(ymin = mean_rate_all - sem, ymax = mean_rate_all + sem),
    width = 0.1,
    position = position_dodge(0.85)
  ) +
  stat_halfeye(
    data = df %>% filter(Sex %in% c("f", "m")),
    aes(x = Group, y = percentage_fix_duration, fill = medication),
    adjust = 0.6,
    position = position_dodge(0.85),
    justification = -1.1,
    .width = 0,
    point_color = NA,
    alpha = 0.4,
    normalize = "groups",
    scale = 0.25
  ) +
  labs(y = "% Fixation Duration") +
  facet_grid(Sex ~ fix) +
  apply_custom_settings(values = c("#FDC010", "#C7B655", "#93AD7C", "#0B7A6B"))

group_by_drug_by_sex

ggsave(
  here("02_outputs", "figures", "supplementary", "fixation_duration", "svg", "fixation_duration_proportion_barplot_by_sex.svg"),
  group_by_drug_by_sex,
  units = "mm",
  width = 183,
  height = 165
)
ggsave(
  here("02_outputs", "figures", "supplementary", "fixation_duration", "png", "Supp_Figure_fixation_duration_by_sex.png"),
  group_by_drug_by_sex,
  units = "mm",
  width = 183,
  height = 165
)
ggsave(
  here("02_outputs", "figures", "supplementary", "fixation_duration", "svg", "fixation_duration_proportion_barplot.svg"),
  group_by_drug,
  units = "mm",
  width = 183,
  height = 120
)
