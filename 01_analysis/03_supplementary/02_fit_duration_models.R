# Manuscript section: Supplementary / control analyses
# Analysis family: fixation duration alternative DV
# Original source path: scripts/publication/analysis/duration_analysis.r
# Primary input dataset(s): data/derived/preprocessing/duration_data.csv
# Primary output(s): duration model summaries and diagnostics for the alternative fixation-duration DV
# Known TODOs: clarify whether duration_data.csv or a prep variant is the intended canonical input; figure generation moved to analysis/04_figures/supplementary/01_plot_fixation_duration.R
# Scientific logic note: scientific logic and model formulas are unchanged from source; plotting was separated for structural cleanup

library(tidyverse)
library(lmerTest)
library(emmeans)
library(glmmTMB)
library(car)
library(DHARMa)
library(performance)
library(here)
source(here("01_analysis", "shared", "model_output_utils.R"))

write_emmeans_contrasts <- function(contrast_obj, output_dir, filename, model_name) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  contrast_df <- summary(contrast_obj, infer = TRUE) %>%
    as_tibble()

  p_col <- intersect(c("p.value", "p_value"), names(contrast_df))
  if (length(p_col) == 0) {
    contrast_df <- contrast_df %>% mutate(p_value = NA_real_)
  } else {
    contrast_df <- contrast_df %>% rename(p_value = all_of(p_col[[1]]))
  }

  contrast_df <- contrast_df %>%
    mutate(
      model = model_name,
      p_display = format_p_value(p_value),
      sig = significance_stars(p_value),
      significant_0_05 = !is.na(p_value) & p_value < 0.05,
      highlight = case_when(
        significant_0_05 ~ "significant",
        !is.na(p_value) & p_value < 0.1 ~ "trend",
        TRUE ~ ""
      ),
      p_display = if_else(sig == "", p_display, paste0(p_display, " ", sig))
    ) %>%
    relocate(model, .before = 1)

  write_csv(contrast_df, file.path(output_dir, filename))
  contrast_df
}

df <- read_csv(here("00_data", "derived", "preprocessing", "duration_data.csv"))

df$fix <- factor(
  df$fix,
  levels = c("fix_on_eyes", "fix_on_mouth", "fix_on_face", "fix_on_background"),
  labels = c("Eyes", "Mouth", "Face", "Background"),
  ordered = TRUE
)

df$medication <- factor(
  df$medication,
  levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA"),
  labels = c("BOTH", "NAL", "OXT", "PLA")
)

model <- lmer(
  percentage_fix_duration ~ Group * fix * medication + Sex + session + (1 | Sub_ID),
  data = df
)

output_dir <- here("02_outputs", "model_outputs", "supplementary", "duration")
omnibus_dir <- file.path(output_dir, "01_omnibus")
contrast_dir <- file.path(output_dir, "02_group_contrasts")
sex_dir <- file.path(output_dir, "03_sex_sensitivity")
diagnostic_dir <- file.path(output_dir, "04_diagnostics")

dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)

summary(model)
anova(model)
write_model_anova(model, omnibus_dir, "01_duration_lmer_anova.csv", "duration_lmer")

emmeans(model, pairwise ~ fix)
emmeans(model, pairwise ~ Group, by = c("fix"))

pairwise_results <- emmeans(model, pairwise ~ Group, by = c("fix", "medication"))
summary(pairwise_results$contrasts, infer = TRUE) %>%
  filter(p.value < 0.05)

emmeans(model, pairwise ~ medication, by = c("fix"))
emmeans(model, pairwise ~ medication, by = c("fix", "Group"))

df_ctrl <- df %>%
  filter(Group == "CTRL")

model_ctrl <- lmer(
  percentage_fix_duration ~ fix * medication + Sex + session + (1 | Sub_ID),
  data = df_ctrl
)

summary(model_ctrl)
anova(model_ctrl)
write_model_anova(model_ctrl, omnibus_dir, "02_duration_ctrl_lmer_anova.csv", "duration_ctrl_lmer")
emmeans(model_ctrl, pairwise ~ medication, by = c("fix"))

model_beta <- glmmTMB(
  percentage_fix_duration ~ Group * fix * medication + Sex + session + (1 | Sub_ID),
  family = beta_family(),
  data = df
)

summary(model_beta)
beta_anova <- Anova(model_beta, type = 2) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term") %>%
  as_tibble()
write_anova_table_from_df(beta_anova, omnibus_dir, "03_duration_beta_glmmtmb_anova.csv", "duration_beta_glmmtmb")

# Current beta-model omnibus pattern on 2026-03-31:
# Group:fix is significant (p = 1.41e-07), indicating a robust AOI allocation
# difference between groups when duration proportions are used.
# Beta-model group contrasts by AOI:
# Eyes ASD - CTRL p = 0.0005
# Background ASD - CTRL p < 0.0001
# Face ASD - CTRL p = 0.0607
beta_group_by_fix <- emmeans(model_beta, pairwise ~ Group, by = "fix")
write_emmeans_contrasts(
  beta_group_by_fix$contrasts,
  contrast_dir,
  "01_duration_beta_group_by_fix_contrasts.csv",
  "duration_beta_glmmtmb"
)

beta_group_by_fix_medication <- emmeans(model_beta, pairwise ~ Group, by = c("fix", "medication"))
write_emmeans_contrasts(
  beta_group_by_fix_medication$contrasts,
  contrast_dir,
  "02_duration_beta_group_by_fix_medication_contrasts.csv",
  "duration_beta_glmmtmb"
)

sim_res <- simulateResiduals(fittedModel = model_beta)
plot(sim_res)
testDispersion(sim_res)
testZeroInflation(sim_res)
testOutliers(sim_res)
plotResiduals(sim_res, df$Group)
plotResiduals(sim_res, df$fix)

outlier_indices <- outliers(sim_res)
print(outlier_indices)
write_csv(tibble(row_index = outlier_indices), file.path(diagnostic_dir, "01_duration_beta_outlier_indices.csv"))

df_filtered <- df[-outlier_indices, ]

model_beta_filtered <- glmmTMB(
  percentage_fix_duration ~ Group * fix * medication + Sex + session + (1 | Sub_ID),
  family = beta_family(),
  data = df_filtered
)

summary(model_beta_filtered)
beta_filtered_anova <- Anova(model_beta_filtered, type = 2) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term") %>%
  as_tibble()
write_anova_table_from_df(beta_filtered_anova, omnibus_dir, "04_duration_beta_filtered_glmmtmb_anova.csv", "duration_beta_filtered_glmmtmb")

df_sex <- df %>%
  filter(Sex %in% c("f", "m")) %>%
  droplevels()

model_beta_sex <- glmmTMB(
  percentage_fix_duration ~ Group * fix * medication * Sex + session + (1 | Sub_ID),
  family = beta_family(),
  data = df_sex
)

beta_sex_anova <- Anova(model_beta_sex, type = 2) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term") %>%
  as_tibble()

# Current binary-sex beta-model omnibus pattern on 2026-03-31:
# fix:Sex is significant (p = 7.25e-12)
# Group:fix:Sex is significant (p = 0.000296)
write_anova_table_from_df(
  beta_sex_anova,
  sex_dir,
  "01_duration_beta_sex_interaction_anova.csv",
  "duration_beta_glmmtmb_sex_interaction"
)

beta_group_by_fix_sex <- emmeans(model_beta_sex, pairwise ~ Group, by = c("fix", "Sex"))
write_emmeans_contrasts(
  beta_group_by_fix_sex$contrasts,
  sex_dir,
  "02_duration_beta_group_by_fix_sex_contrasts.csv",
  "duration_beta_glmmtmb_sex_interaction"
)

beta_sex_by_group_fix <- emmeans(model_beta_sex, pairwise ~ Sex, by = c("Group", "fix"))
write_emmeans_contrasts(
  beta_sex_by_group_fix$contrasts,
  sex_dir,
  "03_duration_beta_sex_by_group_fix_contrasts.csv",
  "duration_beta_glmmtmb_sex_interaction"
)
