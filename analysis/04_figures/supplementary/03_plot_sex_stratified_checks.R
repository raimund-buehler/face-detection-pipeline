# Manuscript section: Supplementary / control analyses
# Analysis family: sex-stratified robustness checks
# Original source path: scripts/publication/analysis/sex_diffs.r
# Primary input dataset(s): data/derived/analysis/df_analysis_pub_prep.csv
# Primary output(s): exploratory sex-stratified descriptive plots
# Known TODOs: no canonical manuscript export path was defined for these plots in the source workflow
# Scientific logic note: descriptive plotting logic was separated from the original mixed script without changing plot construction

library(tidyverse)
library(ggdist)
library(here)

df <- read_csv(here("00_data", "derived", "analysis", "df_analysis_pub_prep.csv"))

df$fix <- factor(
  df$fix,
  levels = c("Eyes", "Mouth", "Face", "Background"),
  labels = c("Eyes", "Mouth", "Face", "Background"),
  ordered = TRUE
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
    labs(y = deparse(var)) +
    facet_wrap(~fix, scales = "fixed") +
    theme_bw() +
    theme(legend.position = "top")
}

barplot_m <- bar_plot(df %>% filter(Sex == "m"), "rate_all")
barplot_f <- bar_plot(df %>% filter(Sex == "f"), "rate_all")

barplot_m
barplot_f
