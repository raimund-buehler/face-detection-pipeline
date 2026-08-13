# Manuscript section: Supplementary / control analyses
# Analysis family: sex sensitivity for timecourse and question-type models
# Original source path: derived from 01_analysis/02_main_manuscript/eye_tracking/02_fit_eye_question_timecourse_model.R
# Primary input dataset(s): 00_data/derived/analysis/df_analysis_pub.csv
# Primary output(s): omnibus ANOVA and selected group/sex contrast tables under 02_outputs/model_outputs/sensitivity_analyses/sex_sensitivity/timecourse/
# Known TODOs: confirm whether df_analysis_pub.csv is the intended canonical input for the label model
# Scientific logic note: preserves the original timecourse models while adding sex interactions on the binary f/m subset

library(tidyverse)
library(lmerTest)
library(emmeans)
library(here)
source(here("01_analysis", "shared", "model_output_utils.R"))

output_root <- here("02_outputs", "model_outputs", "sensitivity_analyses", "sex_sensitivity", "timecourse")
q_type_dir <- file.path(output_root, "q_type")
label_dir <- file.path(output_root, "label")
dir.create(file.path(q_type_dir, "01_omnibus"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(q_type_dir, "02_group_contrasts"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(q_type_dir, "03_sex_contrasts"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(label_dir, "01_omnibus"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(label_dir, "02_group_contrasts"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(label_dir, "03_sex_contrasts"), recursive = TRUE, showWarnings = FALSE)

write_contrasts_only <- function(emm_result, output_dir, prefix) {
  readr::write_csv(
    summary(emm_result$contrasts, infer = c(TRUE, TRUE)) %>% as.data.frame() %>% as_tibble(),
    file.path(output_dir, paste0(prefix, "_contrasts.csv"))
  )
}

df <- read_csv(here("00_data", "derived", "analysis", "df_analysis_pub.csv")) %>%
  filter(!is.na(rate_q), !is.na(medication), Sex %in% c("f", "m")) %>%
  select(ID, session, Sex, Group, medication, label, fix, rate_q, accuracy_degrees) %>%
  distinct(ID, session, fix, label, .keep_all = TRUE) %>%
  ungroup() %>%
  mutate(
    sqrt_rate_q = sqrt(rate_q),
    fix = factor(
      fix,
      levels = c("fix_on_eyes", "fix_on_mouth", "fix_on_face", "fix_on_background"),
      labels = c("Eyes", "Mouth", "Face", "Background")
    )
  )

ordered_labels <- c("Intro", paste("Question", 1:13))
df$label <- factor(df$label, levels = ordered_labels, labels = c("Intro", paste(1:13)), ordered = TRUE)
df <- df %>% filter(fix == "Eyes", label != "Intro")
df <- df %>% mutate(q_type = ifelse(label %in% paste(1:6), "closed", "open"))

q_type_model_sex <- lmer(
  sqrt_rate_q ~ Group * q_type * medication * Sex + session + (1 | ID/q_type),
  data = df
)

q_type_anova <- write_model_anova(
  q_type_model_sex,
  file.path(q_type_dir, "01_omnibus"),
  "01_q_type_sex_interaction_anova.csv",
  "q_type_sex_interaction"
)
print(q_type_anova, n = nrow(q_type_anova))

# Current omnibus snapshot should be checked on rerun; targeted follow-ups below
# focus on Group:Sex, Group:q_type, and Group:q_type:Sex if present.
write_contrasts_only(
  emmeans(q_type_model_sex, pairwise ~ Group | Sex),
  file.path(q_type_dir, "02_group_contrasts"),
  "01_group_within_sex"
)
write_contrasts_only(
  emmeans(q_type_model_sex, pairwise ~ Group | q_type),
  file.path(q_type_dir, "02_group_contrasts"),
  "02_group_within_question_type"
)
write_contrasts_only(
  emmeans(q_type_model_sex, pairwise ~ Group | q_type * Sex),
  file.path(q_type_dir, "02_group_contrasts"),
  "03_group_within_question_type_by_sex"
)
write_contrasts_only(
  emmeans(q_type_model_sex, pairwise ~ Sex | Group * q_type),
  file.path(q_type_dir, "03_sex_contrasts"),
  "01_sex_within_group_by_question_type"
)

label_model_sex <- lmer(
  sqrt_rate_q ~ Group * label * medication * Sex + (1 | ID),
  data = df
)

label_anova <- write_model_anova(
  label_model_sex,
  file.path(label_dir, "01_omnibus"),
  "01_label_sex_interaction_anova.csv",
  "label_sex_interaction"
)
print(label_anova, n = nrow(label_anova))

write_contrasts_only(
  emmeans(label_model_sex, pairwise ~ Group | label),
  file.path(label_dir, "02_group_contrasts"),
  "01_group_within_label"
)
write_contrasts_only(
  emmeans(label_model_sex, pairwise ~ Group | label * Sex),
  file.path(label_dir, "02_group_contrasts"),
  "02_group_within_label_by_sex"
)
write_contrasts_only(
  emmeans(label_model_sex, pairwise ~ Sex | Group * label),
  file.path(label_dir, "03_sex_contrasts"),
  "01_sex_within_group_by_label"
)
