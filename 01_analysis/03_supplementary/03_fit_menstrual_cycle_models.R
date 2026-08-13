# Manuscript section: Supplementary / control analyses
# Analysis family: menstrual cycle analysis
# Original source path: scripts/publication/analysis/menst_analysis.r
# Primary input dataset(s): data/derived/analysis/df_analysis_pub.csv
# Primary output(s): menstrual-cycle model summaries
# Known TODOs: verify dependency on phase labels from external prep step; figure generation moved to analysis/04_figures/supplementary/02_plot_menstrual_cycle.R
# Scientific logic note: scientific logic and model formulas are unchanged from source; plotting was separated for structural cleanup

library(tidyverse)
library(glmmTMB)
library(car)
library(emmeans)
library(here)
source(here("01_analysis", "shared", "model_output_utils.R"))

normalize_session <- function(x) {
  x %>%
    str_trim() %>%
    str_to_lower() %>%
    str_replace_all(" ", "_")
}

write_glmmtmb_coefficients <- function(model, output_dir, filename, model_name) {
  coef_table <- summary(model)$coefficients$cond %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    as_tibble()

  coef_table <- coef_table %>%
    rename(p_value = `Pr(>|z|)`) %>%
    mutate(
      model = model_name,
      p_display = format_p_value(p_value),
      sig = significance_stars(p_value),
      significant_0_05 = !is.na(p_value) & p_value < 0.05,
      highlight = case_when(
        significant_0_05 ~ "significant",
        !is.na(p_value) & p_value < 0.1 ~ "trend",
        TRUE ~ ""
      )
    ) %>%
    relocate(model, term, p_value, p_display, sig, significant_0_05, highlight)

  write_csv(coef_table, file.path(output_dir, filename))
  coef_table
}

session_metadata <- read_csv(
  here("00_data", "derived", "analysis", "df_analysis_pub.csv"),
  show_col_types = FALSE
) %>%
  distinct(ID, session, phase, Sex, Group, medication) %>%
  mutate(session_norm = normalize_session(session))

df <- read_csv(
  here("00_data", "derived", "preprocessing", "duration_data.csv"),
  show_col_types = FALSE
) %>%
  mutate(ID = Sub_ID, session_norm = normalize_session(session)) %>%
  left_join(session_metadata, by = c("ID", "session_norm")) %>%
  mutate(session = coalesce(session.y, session.x)) %>%
  filter(!is.na(phase), Sex == "f") %>%
  mutate(
    fix = factor(
      fix,
      levels = c("fix_on_eyes", "fix_on_mouth", "fix_on_face", "fix_on_background"),
      labels = c("Eyes", "Mouth", "Face", "Background"),
      ordered = TRUE
    ),
    phase = factor(phase, levels = c("Follicular", "Luteal")),
    medication = factor(
      medication,
      levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA")
    ),
    Group = factor(Group),
    session = factor(session)
  ) %>%
  select(ID, session, Group, medication, phase, fix, percentage_fix_duration)

df_menst <- df %>%
  ungroup() %>%
  filter(!is.na(fix))

df_menst_eyes <- df_menst %>% filter(fix == "Eyes")

model_phase <- glmmTMB(
  percentage_fix_duration ~ Group * medication * phase + session + (1 | ID),
  family = beta_family(),
  data = df_menst_eyes
)

output_dir <- here("02_outputs", "model_outputs", "supplementary", "menstrual_cycle")

phase_anova <- car::Anova(model_phase, type = 3) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term") %>%
  as_tibble()
write_anova_table_from_df(phase_anova, output_dir, "01_eyes_phase_model_anova.csv", "eyes_phase_duration_beta")
write_glmmtmb_coefficients(model_phase, output_dir, "02_eyes_phase_model_coefficients.csv", "eyes_phase_duration_beta")
write_csv(broom::tidy(emmeans(model_phase, pairwise ~ Group | phase * medication)$contrasts),
          file.path(output_dir, "03_eyes_group_by_phase_medication.csv"))
write_csv(broom::tidy(emmeans(model_phase, pairwise ~ phase | Group * medication)$contrasts),
          file.path(output_dir, "04_eyes_phase_by_group_medication.csv"))

model_phase_med <- glmmTMB(
  percentage_fix_duration ~ Group * fix * medication * phase + session + (1 | ID),
  family = beta_family(),
  data = df_menst
)

phase_med_anova <- car::Anova(model_phase_med, type = 3) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term") %>%
  as_tibble()
write_anova_table_from_df(phase_med_anova, output_dir, "05_full_phase_model_anova.csv", "full_phase_duration_beta")
write_glmmtmb_coefficients(model_phase_med, output_dir, "06_full_phase_model_coefficients.csv", "full_phase_duration_beta")
write_csv(broom::tidy(emmeans(model_phase_med, pairwise ~ Group | fix * phase * medication)$contrasts),
          file.path(output_dir, "07_group_by_fix_phase_medication.csv"))
