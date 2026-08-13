# Manuscript section: Supplementary / control analyses
# Analysis family: sex sensitivity for raw cortisol and cortisol summary levels
# Original source path: derived from the main cortisol scripts
# Primary input dataset(s): 00_data/derived/analysis/df_cortisol_merged.csv; 00_data/derived/analysis/df_cortisol_merged_with_auc.csv; 00_data/derived/analysis/df_cortisol_min_max.csv
# Primary output(s): omnibus ANOVA and selected group/sex contrast tables under 02_outputs/model_outputs/sensitivity_analyses/sex_sensitivity/cortisol_levels/
# Known TODOs: one non-binary participant is excluded because full sex interactions are only estimable on the binary subset
# Scientific logic note: preserves the original cortisol level models while adding sex interactions on the binary f/m subset

library(tidyverse)
library(here)
library(lmerTest)
library(emmeans)
source(here("01_analysis", "shared", "model_output_utils.R"))

sex_lookup <- read_csv(here("00_data", "derived", "analysis", "df_analysis_pub_prep.csv")) %>%
  distinct(ID, session, Sex)

output_root <- here("02_outputs", "model_outputs", "sensitivity_analyses", "sex_sensitivity", "cortisol_levels")
raw_dir <- file.path(output_root, "raw_cortisol")
auci_dir <- file.path(output_root, "auci_summary")
t2t3_dir <- file.path(output_root, "t2t3_summary")
for (dir_path in c(raw_dir, auci_dir, t2t3_dir)) {
  dir.create(file.path(dir_path, "01_omnibus"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dir_path, "02_group_contrasts"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dir_path, "03_sex_contrasts"), recursive = TRUE, showWarnings = FALSE)
}

write_contrasts_only <- function(emm_result, output_dir, prefix) {
  readr::write_csv(
    summary(emm_result$contrasts, infer = c(TRUE, TRUE)) %>% as.data.frame() %>% as_tibble(),
    file.path(output_dir, paste0(prefix, "_contrasts.csv"))
  )
}

df_raw <- read_csv(here("00_data", "derived", "analysis", "df_cortisol_merged.csv")) %>%
  filter(fix == "Eyes", Sex %in% c("f", "m")) %>%
  mutate(timepoint = factor(timepoint))

raw_cortisol_sex_model <- lmer(
  log(cortisol_1) ~ Group * medication * timepoint * Sex + (1 | ID/session),
  data = df_raw
)
write_model_anova(
  raw_cortisol_sex_model,
  file.path(raw_dir, "01_omnibus"),
  "01_raw_cortisol_sex_interaction_anova.csv",
  "raw_cortisol_sex_interaction"
)
write_contrasts_only(
  emmeans(raw_cortisol_sex_model, pairwise ~ Group | timepoint * Sex),
  file.path(raw_dir, "02_group_contrasts"),
  "01_group_within_timepoint_by_sex"
)
write_contrasts_only(
  emmeans(raw_cortisol_sex_model, pairwise ~ Sex | Group * timepoint),
  file.path(raw_dir, "03_sex_contrasts"),
  "01_sex_within_group_by_timepoint"
)

df_auc <- read_csv(here("00_data", "derived", "analysis", "df_cortisol_merged_with_auc.csv")) %>%
  left_join(sex_lookup, by = c("ID", "session")) %>%
  filter(Sex %in% c("f", "m"))

auci_summary_sex_model <- lmer(AUCi ~ Group * medication * Sex + (1 | ID), data = df_auc)
write_model_anova(
  auci_summary_sex_model,
  file.path(auci_dir, "01_omnibus"),
  "01_auci_summary_sex_interaction_anova.csv",
  "auci_summary_sex_interaction"
)
write_contrasts_only(
  emmeans(auci_summary_sex_model, pairwise ~ Group | Sex),
  file.path(auci_dir, "02_group_contrasts"),
  "01_group_within_sex"
)
write_contrasts_only(
  emmeans(auci_summary_sex_model, pairwise ~ Sex | Group),
  file.path(auci_dir, "03_sex_contrasts"),
  "01_sex_within_group"
)

df_t2t3 <- read_csv(here("00_data", "derived", "analysis", "df_cortisol_min_max.csv")) %>%
  rename(T2T3 = MinMax, T2T3_scaled = MinMax_scaled) %>%
  left_join(sex_lookup, by = c("ID", "session")) %>%
  filter(Sex %in% c("f", "m"))

t2t3_summary_sex_model <- lmer(T2T3 ~ Group * medication * Sex + (1 | ID), data = df_t2t3)
write_model_anova(
  t2t3_summary_sex_model,
  file.path(t2t3_dir, "01_omnibus"),
  "01_t2t3_summary_sex_interaction_anova.csv",
  "t2t3_summary_sex_interaction"
)
write_contrasts_only(
  emmeans(t2t3_summary_sex_model, pairwise ~ Group | Sex),
  file.path(t2t3_dir, "02_group_contrasts"),
  "01_group_within_sex"
)
write_contrasts_only(
  emmeans(t2t3_summary_sex_model, pairwise ~ Sex | Group),
  file.path(t2t3_dir, "03_sex_contrasts"),
  "01_sex_within_group"
)
