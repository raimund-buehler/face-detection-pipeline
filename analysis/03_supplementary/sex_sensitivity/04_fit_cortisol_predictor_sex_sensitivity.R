# Manuscript section: Supplementary / control analyses
# Analysis family: sex sensitivity for cortisol predictor models of eye fixations
# Original source path: derived from the main cortisol predictor scripts
# Primary input dataset(s): 00_data/derived/analysis/df_cortisol_merged_with_auc.csv; 00_data/derived/analysis/df_cortisol_min_max.csv
# Primary output(s): omnibus ANOVA and selected group/sex slope-contrast tables under 02_outputs/model_outputs/sensitivity_analyses/sex_sensitivity/cortisol_predictors/
# Known TODOs: one non-binary participant is excluded because full sex interactions are only estimable on the binary subset
# Scientific logic note: preserves the original AUCi and T2T3 predictor models while adding sex interactions on the binary f/m subset

library(tidyverse)
library(here)
library(lmerTest)
library(emmeans)
source(here("01_analysis", "shared", "model_output_utils.R"))

sex_lookup <- read_csv(here("00_data", "derived", "analysis", "df_analysis_pub_prep.csv")) %>%
  distinct(ID, session, Sex)

output_root <- here("02_outputs", "model_outputs", "sensitivity_analyses", "sex_sensitivity", "cortisol_predictors")
auci_dir <- file.path(output_root, "auci_predictor")
t2t3_dir <- file.path(output_root, "t2t3_predictor")
for (dir_path in c(auci_dir, t2t3_dir)) {
  dir.create(file.path(dir_path, "01_omnibus"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dir_path, "02_group_trend_contrasts"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dir_path, "03_sex_trend_contrasts"), recursive = TRUE, showWarnings = FALSE)
}

write_trend_contrasts_only <- function(emm_result, output_dir, prefix) {
  readr::write_csv(
    summary(emm_result$contrasts, infer = c(TRUE, TRUE)) %>% as.data.frame() %>% as_tibble(),
    file.path(output_dir, paste0(prefix, "_contrasts.csv"))
  )
}

df_auc <- read_csv(here("00_data", "derived", "analysis", "df_cortisol_merged_with_auc.csv")) %>%
  left_join(sex_lookup, by = c("ID", "session")) %>%
  filter(Sex %in% c("f", "m"))

auci_predictor_sex_model <- lmer(
  sqrt_rate_all ~ medication * Group * AUCi * Sex + (1 | ID),
  data = df_auc
)
write_model_anova(
  auci_predictor_sex_model,
  file.path(auci_dir, "01_omnibus"),
  "01_auci_predictor_sex_interaction_anova.csv",
  "auci_predictor_sex_interaction"
)
write_trend_contrasts_only(
  emtrends(auci_predictor_sex_model, pairwise ~ Group | Sex, var = "AUCi"),
  file.path(auci_dir, "02_group_trend_contrasts"),
  "01_group_trend_within_sex"
)
write_trend_contrasts_only(
  emtrends(auci_predictor_sex_model, pairwise ~ Group | medication * Sex, var = "AUCi"),
  file.path(auci_dir, "02_group_trend_contrasts"),
  "02_group_trend_within_medication_by_sex"
)
write_trend_contrasts_only(
  emtrends(auci_predictor_sex_model, pairwise ~ Sex | Group, var = "AUCi"),
  file.path(auci_dir, "03_sex_trend_contrasts"),
  "01_sex_trend_within_group"
)
write_trend_contrasts_only(
  emtrends(auci_predictor_sex_model, pairwise ~ Sex | medication * Group, var = "AUCi"),
  file.path(auci_dir, "03_sex_trend_contrasts"),
  "02_sex_trend_within_medication_by_group"
)

df_t2t3 <- read_csv(here("00_data", "derived", "analysis", "df_cortisol_min_max.csv")) %>%
  rename(T2T3 = MinMax, T2T3_scaled = MinMax_scaled) %>%
  left_join(sex_lookup, by = c("ID", "session")) %>%
  filter(Sex %in% c("f", "m"))

t2t3_predictor_sex_model <- lmer(
  sqrt_rate_all ~ medication * Group * T2T3 * Sex + (1 | ID),
  data = df_t2t3
)
write_model_anova(
  t2t3_predictor_sex_model,
  file.path(t2t3_dir, "01_omnibus"),
  "01_t2t3_predictor_sex_interaction_anova.csv",
  "t2t3_predictor_sex_interaction"
)
write_trend_contrasts_only(
  emtrends(t2t3_predictor_sex_model, pairwise ~ Group | Sex, var = "T2T3"),
  file.path(t2t3_dir, "02_group_trend_contrasts"),
  "01_group_trend_within_sex"
)
write_trend_contrasts_only(
  emtrends(t2t3_predictor_sex_model, pairwise ~ Group | medication * Sex, var = "T2T3"),
  file.path(t2t3_dir, "02_group_trend_contrasts"),
  "02_group_trend_within_medication_by_sex"
)
write_trend_contrasts_only(
  emtrends(t2t3_predictor_sex_model, pairwise ~ Sex | Group, var = "T2T3"),
  file.path(t2t3_dir, "03_sex_trend_contrasts"),
  "01_sex_trend_within_group"
)
write_trend_contrasts_only(
  emtrends(t2t3_predictor_sex_model, pairwise ~ Sex | medication * Group, var = "T2T3"),
  file.path(t2t3_dir, "03_sex_trend_contrasts"),
  "02_sex_trend_within_medication_by_group"
)
