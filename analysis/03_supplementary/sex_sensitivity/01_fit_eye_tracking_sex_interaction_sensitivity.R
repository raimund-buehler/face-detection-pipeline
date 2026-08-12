# Manuscript section: Supplementary / control analyses
# Analysis family: sex-interaction sensitivity analysis
# Original source path: new structured robustness script
# Primary input dataset(s): 00_data/derived/analysis/df_analysis_pub_prep.csv
# Primary output(s): binary interaction omnibus ANOVA and exploratory contrast tables under 02_outputs/model_outputs/sensitivity_analyses/sex_sensitivity/eye_tracking_interaction/
# Known TODOs: interaction model is fit on the binary sex subset only because a single non-binary participant makes the full interaction rank-deficient
# Scientific logic note: this script focuses on the binary `Group * fix * medication * Sex` sensitivity model and its exploratory follow-ups

library(tidyverse)
library(lmerTest)
library(emmeans)
library(here)
source(here("01_analysis", "shared", "model_output_utils.R"))

output_dir <- here("02_outputs", "model_outputs", "sensitivity_analyses", "sex_sensitivity", "eye_tracking_interaction")
anova_dir <- file.path(output_dir, "01_omnibus")
group_dir <- file.path(output_dir, "02_group_contrasts")
sex_dir <- file.path(output_dir, "03_sex_contrasts")

dir.create(anova_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sex_dir, recursive = TRUE, showWarnings = FALSE)

df <- read_csv(here("00_data", "derived", "analysis", "df_analysis_pub_prep.csv")) %>%
  mutate(
    fix = factor(
      fix,
      levels = c("Eyes", "Mouth", "Face", "Background"),
      labels = c("Eyes", "Mouth", "Face", "Background"),
      ordered = TRUE
    ),
    Group = factor(Group),
    medication = factor(medication),
    Sex = factor(Sex)
  )

df_binary <- df %>% filter(Sex %in% c("f", "m"))

write_contrasts_only <- function(emm_result, output_dir, prefix) {
  readr::write_csv(
    summary(emm_result$contrasts, infer = c(TRUE, TRUE)) %>% as.data.frame() %>% as_tibble(),
    file.path(output_dir, paste0(prefix, "_contrasts.csv"))
  )
}

model_interaction_binary <- lmer(
  sqrt_rate_all ~ Group * fix * medication * Sex + session + (1 | ID),
  data = df_binary
)

results_interaction_binary <- write_model_anova(
  model_interaction_binary,
  anova_dir,
  "01_full_interaction_binary_subset_anova.csv",
  "full_interaction_binary_subset"
)

print(results_interaction_binary, n = nrow(results_interaction_binary))

# Current omnibus snapshot as of 2026-03-31:
# `full_interaction_binary_subset`: `fix` p < .001, `Group:fix` p = .024,
# `Group:Sex` p = .049, `fix:Sex` p < .001, `Group:fix:Sex` p = .018

# Exploratory follow-ups for the significant omnibus terms in the full
# `Group * fix * medication * Sex` sensitivity model:
# - `Group:Sex`: Group contrasts within Sex, collapsed over AOI/medication
# - `Group:fix`: Group contrasts within AOI, collapsed over medication/Sex
# - `Group:fix:Sex`: Group contrasts within each AOI x Sex cell and
#   complementary Sex contrasts within each Group x AOI cell

emm_interaction_group_by_sex <- emmeans(model_interaction_binary, pairwise ~ Group | Sex)
write_contrasts_only(
  emm_interaction_group_by_sex,
  group_dir,
  "01_group_within_sex"
)

emm_interaction_group_by_fix <- emmeans(model_interaction_binary, pairwise ~ Group | fix)
write_contrasts_only(
  emm_interaction_group_by_fix,
  group_dir,
  "02_group_within_fix"
)

emm_interaction_group_by_fix_sex <- emmeans(model_interaction_binary, pairwise ~ Group | fix * Sex)
write_contrasts_only(
  emm_interaction_group_by_fix_sex,
  group_dir,
  "03_group_within_fix_by_sex"
)

emm_interaction_sex_by_group_fix <- emmeans(model_interaction_binary, pairwise ~ Sex | Group * fix)
write_contrasts_only(
  emm_interaction_sex_by_group_fix,
  sex_dir,
  "01_sex_within_group_by_fix"
)
