# Manuscript section: Cortisol analyses
# Analysis family: raw cortisol over timepoints
# Original source path: legacy/cortisol_reports/markdown_source_bundle/1_cortisol_raw.Rmd
# Primary input dataset(s): data/derived/analysis/df_cortisol_merged.csv
# Primary output(s): raw cortisol model summaries and emmeans contrasts
# Known TODOs: plotting moved to analysis/04_figures/main/cortisol/05_plot_raw_cortisol.R
# Scientific logic note: extracted from the manuscript markdown workflow without changing the model formula

library(tidyverse)
library(here)
library(lmerTest)
library(emmeans)
source(here("01_analysis", "shared", "model_output_utils.R"))

df <- read_csv(here("00_data", "derived", "analysis", "df_cortisol_merged.csv"))

df <- df %>%
  filter(fix == "Eyes") %>%
  mutate(timepoint = as.factor(timepoint))

raw_cortisol_model <- lmer(
  log(cortisol_1) ~ Group * medication * timepoint + (1 | ID/session),
  data = df
)

output_dir <- here("02_outputs", "model_outputs", "main_manuscript", "cortisol", "raw_cortisol")

summary(raw_cortisol_model)
anova(raw_cortisol_model)
write_model_anova(raw_cortisol_model, output_dir, "01_raw_cortisol_model_anova.csv", "raw_cortisol_model")

emm_med <- emmeans(raw_cortisol_model, pairwise ~ medication)
emm_med_x_group <- emmeans(raw_cortisol_model, pairwise ~ medication, by = "Group")
emm_time <- emmeans(raw_cortisol_model, pairwise ~ timepoint)
emm_time_x_med <- emmeans(raw_cortisol_model, pairwise ~ timepoint, by = "medication")
emm_med_x_time <- emmeans(raw_cortisol_model, pairwise ~ medication, by = "timepoint")
emm_group_x_med <- emmeans(raw_cortisol_model, pairwise ~ Group, by = "medication")
emm_group_x_med_x_time <- emmeans(raw_cortisol_model, pairwise ~ Group, by = c("medication", "timepoint"))

emm_med
emm_med_x_group
emm_time
emm_time_x_med
emm_med_x_time
emm_group_x_med
emm_group_x_med_x_time
