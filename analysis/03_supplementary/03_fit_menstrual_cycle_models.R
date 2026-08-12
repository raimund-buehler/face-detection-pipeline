# Manuscript section: Supplementary / control analyses
# Analysis family: menstrual cycle analysis
# Original source path: scripts/publication/analysis/menst_analysis.r
# Primary input dataset(s): data/derived/analysis/df_analysis_pub_prep.csv
# Primary output(s): menstrual-cycle model summaries
# Known TODOs: verify dependency on phase labels from external prep step; figure generation moved to analysis/04_figures/supplementary/02_plot_menstrual_cycle.R
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

df_menst <- df %>%
  ungroup() %>%
  filter(phase != "Missing")

var <- sym("sqrt_rate_all")

df_menst_eyes <- df_menst %>% filter(fix == "Eyes")

model_phase <- lmer(
  sqrt_rate_all ~ Group * medication * phase + session + (1 | ID),
  data = df_menst_eyes
)

output_dir <- here("02_outputs", "model_outputs", "supplementary", "menstrual_cycle")

anova(model_phase)
write_model_anova(model_phase, output_dir, "01_eyes_phase_model_anova.csv", "eyes_phase_model")
emmeans(model_phase, pairwise ~ Group)
emmeans(model_phase, pairwise ~ phase, by = c("Group"))

model_phase_med <- lmer(
  sqrt_rate_all ~ Group * fix * medication * phase + (1 | ID),
  data = df_menst
)

anova(model_phase_med)
write_model_anova(model_phase_med, output_dir, "02_full_phase_model_anova.csv", "full_phase_model")
emmeans(model_phase_med, pairwise ~ Group, by = c("fix", "phase"))
