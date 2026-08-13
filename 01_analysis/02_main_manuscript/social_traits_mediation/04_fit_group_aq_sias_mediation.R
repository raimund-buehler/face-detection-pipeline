# Manuscript section: Mediation analyses
# Analysis family: group -> SIAS -> eye fixations; AQ -> SIAS -> eye fixations
# Original source path: scripts/publication/analysis/Mediation.r
# Primary input dataset(s): data prepared by social_traits_mediation/shared_aq_sias_plot_utils.R from manuscript analysis tables
# Primary output(s): mediation model summaries for Group -> SIAS -> Eyes and AQ -> SIAS -> Eyes
# Known TODOs: mediation figure assembly moved to analysis/04_figures/main/social_traits_mediation/04_plot_mediation_figure.R; mediation panel PNGs remain external figure assets
# Scientific logic note: scientific logic is unchanged from source; plotting was separated for structural cleanup

library(here)
library(tidyverse)
library(lmerTest)
library(emmeans)
library(car)
library(effectsize)
source(here("01_analysis", "02_main_manuscript", "social_traits_mediation", "shared_aq_sias_plot_utils.R"))
source(here("01_analysis", "shared", "model_output_utils.R"))

df <- prepare_aq_sias_data()
df_cor <- df %>% distinct(ID, .keep_all = TRUE) %>% filter(!is.na(AQ))

df_mod <- df %>% filter(fix == "Eyes")
df_mod$Group <- factor(df_mod$Group, levels = c("CTRL", "ASD"))

df_mediator_group <- df_mod[complete.cases(df_mod$Group, df_mod$SIAS), ]
df_outcome_group <- df_mod[complete.cases(df_mod$Group, df_mod$SIAS, df_mod$sqrt_rate_all), ]

mediator_model_group <- lm(SIAS ~ Group, df_mediator_group)
outcome_model_group <- lm(sqrt_rate_all ~ Group + SIAS, df_outcome_group)
output_dir <- here("02_outputs", "model_outputs", "main_manuscript", "social_traits_mediation", "mediation")

write_model_anova(mediator_model_group, output_dir, "01_group_mediator_model_anova.csv", "group_mediator_model")
write_lm_coefficients(mediator_model_group, output_dir, "01_group_mediator_model_coefficients.csv", "group_mediator_model")
write_model_anova(outcome_model_group, output_dir, "02_group_outcome_model_anova.csv", "group_outcome_model")
write_lm_coefficients(outcome_model_group, output_dir, "02_group_outcome_model_coefficients.csv", "group_outcome_model")

med_fit_group <- tryCatch(
  mediation::mediate(
    mediator_model_group,
    outcome_model_group,
    treat = "Group",
    mediator = "SIAS",
    boot = FALSE,
    sims = 1000,
    cluster = df_outcome_group$ID
  ),
  error = function(e) e
)

summary(mediator_model_group)
t_to_d(15.92, 237)

summary(outcome_model_group)
if (inherits(med_fit_group, "error")) {
  write_text_output(
    c("mediation::mediate failed in this environment", conditionMessage(med_fit_group)),
    output_dir,
    "03_group_mediation_summary.txt"
  )
} else {
  summary(med_fit_group)
  write_text_output(capture.output(summary(med_fit_group)), output_dir, "03_group_mediation_summary.txt")
}

df_mean <- df_mod %>%
  group_by(ID) %>%
  mutate(sqrt_rate_all = mean(sqrt_rate_all)) %>%
  distinct(ID, .keep_all = TRUE)

effect_model_group <- lm(SIAS ~ Group, df_mean)
write_model_anova(effect_model_group, output_dir, "04_group_effect_model_anova.csv", "group_effect_model")
write_lm_coefficients(effect_model_group, output_dir, "04_group_effect_model_coefficients.csv", "group_effect_model")

summary(effect_model_group)
t_to_d(8.249, 62)

df_mediator_aq <- df_mod[complete.cases(df_mod$AQ, df_mod$SIAS), ]
df_outcome_aq <- df_mod[complete.cases(df_mod$AQ, df_mod$SIAS, df_mod$sqrt_rate_all), ]

mediator_model_aq <- lm(SIAS ~ AQ, data = df_mediator_aq)
outcome_model_aq <- lm(sqrt_rate_all ~ AQ + SIAS, data = df_outcome_aq)
write_model_anova(mediator_model_aq, output_dir, "05_aq_mediator_model_anova.csv", "aq_mediator_model")
write_lm_coefficients(mediator_model_aq, output_dir, "05_aq_mediator_model_coefficients.csv", "aq_mediator_model")
write_model_anova(outcome_model_aq, output_dir, "06_aq_outcome_model_anova.csv", "aq_outcome_model")
write_lm_coefficients(outcome_model_aq, output_dir, "06_aq_outcome_model_coefficients.csv", "aq_outcome_model")

med_fit_aq <- tryCatch(
  mediation::mediate(
    mediator_model_aq,
    outcome_model_aq,
    treat = "AQ",
    mediator = "SIAS",
    boot = FALSE,
    sims = 1000,
    cluster = df_outcome_aq$ID
  ),
  error = function(e) e
)

summary(mediator_model_aq)
summary(outcome_model_aq)
if (inherits(med_fit_aq, "error")) {
  write_text_output(
    c("mediation::mediate failed in this environment", conditionMessage(med_fit_aq)),
    output_dir,
    "07_aq_mediation_summary.txt"
  )
} else {
  summary(med_fit_aq)
  write_text_output(capture.output(summary(med_fit_aq)), output_dir, "07_aq_mediation_summary.txt")
}

t_to_r(25.08, 238)

cor(df_cor$AQ, df_cor$SIAS, method = "pearson")
cor.test(df_cor$AQ, df_cor$SIAS, method = "pearson")
