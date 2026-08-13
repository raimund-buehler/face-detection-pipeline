# Manuscript section: Supplementary / control analyses
# Analysis family: experimenter effects robustness checks
# Original source path: scripts/publication/analysis/experimenter_diffs.r
# Primary input dataset(s): data/derived/analysis/df_analysis_pub_prep.csv
# Primary output(s): exploratory experimenter-comparison plot
# Known TODOs: no canonical manuscript export path was defined for this plot in the source workflow
# Scientific logic note: descriptive plotting logic was separated from the original mixed script without changing plot construction

library(tidyverse)
library(here)

df <- read_csv(here("00_data", "derived", "analysis", "df_analysis_pub_prep.csv"))

df$fix <- factor(
  df$fix,
  levels = c("Eyes", "Mouth", "Face", "Background"),
  labels = c("Eyes", "Mouth", "Face", "Background"),
  ordered = TRUE
)

experimenter_plot <- ggplot(df, aes(x = Experimenter, y = rate_all, fill = Group)) +
  geom_boxplot() +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 3,
    color = "red",
    position = position_dodge(width = 0.75)
  ) +
  facet_wrap(~fix, scales = "free") +
  theme_bw()

experimenter_plot
