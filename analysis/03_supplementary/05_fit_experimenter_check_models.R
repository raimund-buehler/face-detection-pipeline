# Manuscript section: Supplementary / control analyses
# Analysis family: experimenter effects robustness checks
# Original source path: scripts/publication/analysis/experimenter_diffs.r
# Primary input dataset(s): data/derived/analysis/df_analysis_pub_prep.csv
# Primary output(s): experimenter-effects model summaries
# Known TODOs: figure generation moved to analysis/04_figures/supplementary/04_plot_experimenter_checks.R
# Scientific logic note: scientific logic and model formulas are unchanged from source; plotting was separated for structural cleanup

library(tidyverse)
library(lmerTest)
library(emmeans)
library(here)
source(here("01_analysis", "shared", "model_output_utils.R"))

df <- read_csv(here("00_data", "derived", "analysis", "df_analysis_pub_prep.csv"))

df$fix <- factor(
  df$fix,
  levels = c("Eyes", "Mouth", "Face", "Background"),
  labels = c("Eyes", "Mouth", "Face", "Background"),
  ordered = TRUE
)

model <- lmer(
  sqrt_rate_all ~ Group * fix * Experimenter + (fix | ID),
  data = df
)

output_dir <- here("02_outputs", "model_outputs", "supplementary", "experimenter")

summary(model)
anova(model)
write_model_anova(model, output_dir, "01_experimenter_model_anova.csv", "experimenter_model")

emmeans(model, pairwise ~ fix, by = "Experimenter")
emmeans(model, pairwise ~ Group, by = "Experimenter", at = list(fix = "Eyes"))
