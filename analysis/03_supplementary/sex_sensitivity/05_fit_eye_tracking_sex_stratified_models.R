# Manuscript section: Supplementary / control analyses
# Analysis family: sex-stratified robustness checks
# Original source path: scripts/publication/analysis/sex_diffs.r
# Primary input dataset(s): data/derived/analysis/df_analysis_pub_prep.csv
# Primary output(s): sex-stratified model summaries under 02_outputs/model_outputs/sensitivity_analyses/sex_sensitivity/eye_tracking_stratified/
# Known TODOs: figure generation moved to analysis/04_figures/supplementary/03_plot_sex_stratified_checks.R
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

df_m <- df %>% filter(Sex == "m")

model_m <- lmer(
  rate_all ~ Group * fix * medication + session + (1 | ID),
  data = df_m
)

output_dir <- here("02_outputs", "model_outputs", "sensitivity_analyses", "sex_sensitivity", "eye_tracking_stratified")

summary(model_m)
anova(model_m)
write_model_anova(model_m, output_dir, "01_male_stratified_model_anova.csv", "male_stratified_model")

df_f <- df %>% filter(Sex == "f")

model_f <- lmer(
  rate_all ~ Group * fix * medication + session + (1 | ID),
  data = df_f
)

summary(model_f)
anova(model_f)
write_model_anova(model_f, output_dir, "02_female_stratified_model_anova.csv", "female_stratified_model")

emmeans(model_f, pairwise ~ Group, by = c("fix"))
emmeans(model_f, pairwise ~ Group, by = c("fix", "medication"), at = list(fix = "Eyes"))
emmeans(model_f, pairwise ~ medication, by = c("fix", "Group"), at = list(fix = "Eyes"))
