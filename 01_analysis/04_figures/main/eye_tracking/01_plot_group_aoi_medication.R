# Manuscript section: Main eye-tracking model
# Analysis family: figure generation
# Original source path: scripts/publication/analysis/overall_model_pub.R
# Primary input dataset(s): data/derived/analysis/df_analysis_pub_prep.csv
# Primary output(s): outputs/figures/main/group_by_drug.svg; outputs/figures/main/Figure_2.png; outputs/figures/main/submission/Figure_2.pdf; outputs/figures/main/submission/Figure_2.png
# Known TODOs: significance annotations remain hard-coded from the source figure script
# Scientific logic note: plotting logic was separated from the original mixed script without changing figure behavior

library(tidyverse)
library(ggdist)
library(ggh4x)
library(ggsignif)
library(here)
source(here("01_analysis", "04_figures", "shared", "custom_plot_settings.R"))

df <- read_csv(here("00_data", "derived", "analysis", "df_analysis_pub_prep.csv"))

df$fix <- factor(
  df$fix,
  levels = c("Eyes", "Mouth", "Face", "Background"),
  labels = c("Eyes", "Mouth", "Face", "Background"),
  ordered = TRUE
)

df$medication_simple <- factor(
  df$medication,
  levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA"),
  labels = c("BOTH", "NAL", "OXT", "PLA")
)

bar_plot <- function(data, var) {
  var <- rlang::sym(var)

  df_summary <- data %>%
    group_by(Group, fix, medication_simple) %>%
    summarise(
      mean_rate_all = mean(!!var),
      sem = sd(!!var) / sqrt(n()),
      .groups = "drop"
    )

  ggplot(df_summary, aes(x = Group, y = mean_rate_all, fill = medication_simple)) +
    geom_bar(stat = "identity", position = position_dodge(0.85), width = 0.4) +
    geom_errorbar(
      aes(ymin = mean_rate_all - sem, ymax = mean_rate_all + sem),
      width = 0.1,
      position = position_dodge(0.85)
    ) +
    stat_halfeye(
      data = data,
      aes(x = Group, y = !!var, fill = medication_simple),
      adjust = 0.6,
      position = position_dodge(0.85),
      justification = -1.1,
      .width = 0,
      point_color = NA,
      alpha = 0.4,
      normalize = "groups",
      scale = 0.25
    ) +
    labs(x = NULL, y = "Fixations / Second") +
    facet_wrap2(~fix, scales = "fixed", axes = "all", remove_labels = "all")
}

df_summary <- df %>%
  group_by(Group, fix, medication_simple) %>%
  summarise(
    mean_rate_all = mean(sqrt_rate_all),
    sem = sd(sqrt_rate_all) / sqrt(n()),
    .groups = "drop"
  )

y_position_main <- max(df_summary %>% filter(fix == "Eyes") %>% pull(mean_rate_all)) + 0.75

add_significance_lines <- function(summary_df, y_position_main, aoi) {
  y_position_branch <- y_position_main - 0.1

  list(
    geom_signif(
      data = summary_df %>% filter(fix == aoi),
      aes(xmin = "ASD", xmax = "CTRL", annotations = "*"),
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

group_by_drug <- bar_plot(df, "sqrt_rate_all") +
  apply_custom_settings(
    values = c("#FDC010", "#C7B655", "#93AD7C", "#0B7A6B"),
    base_size = 8,
    axis_text_size = 8,
    axis_title_size = 8,
    title_size = 9,
    legend_text_size = 7,
    strip_text_size = 9
  ) +
  theme(
    legend.position = "right",
    legend.key.spacing.y = unit(4, "mm")
  ) +
  add_significance_lines(df_summary, y_position_main, "Eyes") +
  add_significance_lines(df_summary, y_position_main, "Background") +
  geom_signif(
    data = df_summary %>% filter(fix == "Eyes"),
    aes(xmin = 1.3, xmax = 2.3, annotations = "*"),
    y_position = 1.3,
    tip_length = 0.03,
    textsize = 5
  )

group_by_drug

ggsave(here("02_outputs", "figures", "main", "group_by_drug.svg"), group_by_drug, units = "mm", width = 183)
ggsave(here("02_outputs", "figures", "main", "Figure_2.png"), group_by_drug, units = "mm", width = 183, height = 120)
ggsave(here("02_outputs", "figures", "main", "submission", "Figure_2.pdf"), group_by_drug, units = "mm", width = 190, height = 110)
ggsave(here("02_outputs", "figures", "main", "submission", "Figure_2.png"), group_by_drug, units = "mm", width = 190, height = 110, dpi = 600)
