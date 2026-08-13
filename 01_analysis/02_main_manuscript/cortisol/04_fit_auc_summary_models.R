# Manuscript section: Cortisol analyses
# Analysis family: AUCg and AUCi summary models
# Original source path: legacy/cortisol_reports/markdown_source_bundle/2_cortisol_auc.Rmd
# Primary input dataset(s): data/derived/analysis/df_cortisol_merged_with_auc.csv
# Primary output(s): AUCg and AUCi model summaries and emmeans contrasts
# Known TODOs: plotting moved to analysis/04_figures/main/cortisol/06_plot_auc_summary.R
# Scientific logic note: extracted from the manuscript markdown workflow without changing the model formulas

library(tidyverse)
library(here)
library(lmerTest)
library(emmeans)
library(performance)
source(here("01_analysis", "shared", "model_output_utils.R"))

df <- read_csv(here("00_data", "derived", "analysis", "df_cortisol_merged_with_auc.csv"))

aucg_model <- lmer(AUCg ~ Group * medication + (1 | ID), data = df)
output_dir <- here("02_outputs", "model_outputs", "main_manuscript", "cortisol", "auc_summary")

summary(aucg_model)
anova(aucg_model)
write_model_anova(aucg_model, output_dir, "01_aucg_model_anova.csv", "aucg_model")

emm_aucg_med <- emmeans(aucg_model, pairwise ~ medication)
emm_aucg_med_x_group <- emmeans(aucg_model, pairwise ~ medication | Group)
emm_aucg_group_x_med <- emmeans(aucg_model, pairwise ~ Group | medication)

emm_aucg_med
emm_aucg_med_x_group
emm_aucg_group_x_med

auci_model <- lmer(AUCi ~ medication * Group + (1 | ID), data = df)

summary(auci_model)
anova(auci_model)
write_model_anova(auci_model, output_dir, "02_auci_model_anova.csv", "auci_model")

emm_auci_med <- emmeans(auci_model, pairwise ~ medication)
emm_auci_group <- emmeans(auci_model, pairwise ~ Group)
emm_auci_med_x_group <- emmeans(auci_model, pairwise ~ medication | Group)
emm_auci_group_x_med <- emmeans(auci_model, pairwise ~ Group | medication)

emm_auci_med
emm_auci_group
emm_auci_med_x_group
emm_auci_group_x_med

check_model(auci_model)
