# Manuscript section: Cortisol analyses
# Analysis family: cortisol reactivity predicting fixations on eyes
# Original source path: legacy/cortisol_reports/markdown_source_bundle/3_cortisol_auc_pred.Rmd
# Primary input dataset(s): data/derived/analysis/df_cortisol_merged_with_auc.csv
# Primary output(s): AUCg and AUCi predictor model summaries and emtrends results
# Known TODOs: plotting moved to analysis/04_figures/main/cortisol/07_plot_auc_predictor_effects.R
# Scientific logic note: extracted from the manuscript markdown workflow without changing the model formulas

library(tidyverse)
library(here)
library(lmerTest)
library(emmeans)
library(performance)
library(effectsize)
source(here("01_analysis", "shared", "model_output_utils.R"))

df_session <- read_csv(here("00_data", "derived", "analysis", "df_cortisol_merged_with_auc.csv"))

aucg_predictor_model <- lmer(
  sqrt_rate_all ~ medication * Group * AUCg + (1 | ID),
  data = df_session
)

output_dir <- here("02_outputs", "model_outputs", "main_manuscript", "cortisol", "auc_predictors")

summary(aucg_predictor_model)
anova(aucg_predictor_model)
write_model_anova(aucg_predictor_model, output_dir, "01_aucg_predictor_model_anova.csv", "aucg_predictor_model")

emm_aucg_group <- emmeans(aucg_predictor_model, pairwise ~ Group)
emm_aucg_group_x_med <- emmeans(aucg_predictor_model, pairwise ~ Group, by = c("medication"))
emm_aucg_trend_x_group <- emtrends(aucg_predictor_model, ~Group, var = "AUCg", infer = c(TRUE, TRUE))
emm_aucg_trend_x_group_med <- emtrends(aucg_predictor_model, ~medication * Group, var = "AUCg", infer = c(TRUE, TRUE))

emm_aucg_group
emm_aucg_group_x_med
emm_aucg_trend_x_group
emm_aucg_trend_x_group_med

auci_predictor_model <- lmer(
  sqrt_rate_all ~ medication * Group * AUCi + (1 | ID),
  data = df_session
)

summary(auci_predictor_model)
anova(auci_predictor_model)
write_model_anova(auci_predictor_model, output_dir, "02_auci_predictor_model_anova.csv", "auci_predictor_model")

emm_auci_group_x_med <- emmeans(auci_predictor_model, pairwise ~ Group, by = c("medication"))
emm_auci_trend_x_group <- emtrends(auci_predictor_model, ~Group, var = "AUCi", infer = c(TRUE, TRUE))
emm_auci_trend_x_group_med <- emtrends(auci_predictor_model, ~ medication * Group, var = "AUCi", infer = c(TRUE, TRUE))

emm_auci_group_x_med
emm_auci_trend_x_group
emm_auci_trend_x_group_med

emtrends_df <- as.data.frame(emm_auci_trend_x_group)
trend_asd <- emtrends_df %>% filter(Group == "ASD")
t_to_r(trend_asd$t.ratio, df = trend_asd$df)

check_model(auci_predictor_model)
