# Manuscript section: Cortisol analyses
# Analysis family: T2T3 / MinMax summary and predictor models
# Original source path: legacy/cortisol_reports/markdown_source_bundle/4_cortisol_maxmin.Rmd
# Primary input dataset(s): data/derived/analysis/df_cortisol_min_max.csv
# Primary output(s): T2T3 summary and predictor model summaries
# Known TODOs: plotting moved to analysis/04_figures/main/cortisol/08_plot_t2t3_effects.R
# Scientific logic note: extracted from the manuscript markdown workflow without changing the model formulas

library(tidyverse)
library(here)
library(lmerTest)
library(emmeans)
library(performance)
library(effectsize)
source(here("01_analysis", "shared", "model_output_utils.R"))

df <- read_csv(here("00_data", "derived", "analysis", "df_cortisol_min_max.csv")) %>%
  rename(T2T3 = MinMax, T2T3_scaled = MinMax_scaled)

t2t3_summary_model <- lmer(T2T3 ~ Group * medication + (1 | ID), data = df)
output_dir <- here("02_outputs", "model_outputs", "main_manuscript", "cortisol", "t2t3")

summary(t2t3_summary_model)
anova(t2t3_summary_model)
write_model_anova(t2t3_summary_model, output_dir, "01_t2t3_summary_model_anova.csv", "t2t3_summary_model")

emm_t2t3_med <- emmeans(t2t3_summary_model, pairwise ~ medication)
emm_t2t3_med_x_group <- emmeans(t2t3_summary_model, pairwise ~ medication | Group)
emm_t2t3_group_x_med <- emmeans(t2t3_summary_model, pairwise ~ Group | medication)

emm_t2t3_med
emm_t2t3_med_x_group
emm_t2t3_group_x_med

check_model(t2t3_summary_model)

t2t3_predictor_model <- lmer(
  sqrt_rate_all ~ medication * Group * T2T3 + (1 | ID),
  data = df
)

summary(t2t3_predictor_model)
anova(t2t3_predictor_model)
write_model_anova(t2t3_predictor_model, output_dir, "02_t2t3_predictor_model_anova.csv", "t2t3_predictor_model")

emm_t2t3_group_x_med_pred <- emmeans(t2t3_predictor_model, pairwise ~ Group, by = c("medication"))
emm_t2t3_trend_x_group <- emtrends(t2t3_predictor_model, ~Group, var = "T2T3", infer = c(TRUE, TRUE))
emm_t2t3_trend_x_med <- emtrends(t2t3_predictor_model, ~medication, var = "T2T3", infer = c(TRUE, TRUE))
emm_t2t3_trend_x_group_med <- emtrends(t2t3_predictor_model, ~ medication * Group, var = "T2T3", infer = c(TRUE, TRUE))

emm_t2t3_group_x_med_pred
emm_t2t3_trend_x_group
emm_t2t3_trend_x_med
emm_t2t3_trend_x_group_med

contrast(emm_t2t3_trend_x_group_med, method = "pairwise", adjust = "tukey")
contrast(emm_t2t3_trend_x_med, method = "pairwise", adjust = "tukey")

emtrends_df <- as.data.frame(emm_t2t3_trend_x_group_med)
pla_pla_trend <- emtrends_df %>%
  filter(medication == "PLA/PLA", Group == "CTRL")
t_to_r(pla_pla_trend$t.ratio, df = pla_pla_trend$df)

check_model(t2t3_predictor_model)
