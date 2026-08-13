# Manuscript section: Main eye-tracking model (fixation duration)
# Analysis family: group x AOI x medication x sex
# Primary input dataset(s): 00_data/derived/preprocessing/duration_data.csv; 00_data/derived/analysis/df_analysis_pub_prep.csv
# Primary output(s): final duration-model summaries, post-hocs, and diagnostics for manuscript reporting

library(tidyverse)
library(glmmTMB)
library(car)
library(emmeans)
library(DHARMa)
library(performance)
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

output_dir <- here("02_outputs", "model_outputs", "main_manuscript", "eye_tracking", "group_aoi_medication_duration")
diag_dir <- file.path(output_dir, "diagnostics")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

session_metadata <- read_csv(
  here("00_data", "derived", "analysis", "df_analysis_pub_prep.csv"),
  show_col_types = FALSE
) %>%
  distinct(ID, session, Sex, Group, medication) %>%
  mutate(session_norm = normalize_session(session))

df <- read_csv(
  here("00_data", "derived", "preprocessing", "duration_data.csv"),
  show_col_types = FALSE
) %>%
  mutate(ID = Sub_ID, session_norm = normalize_session(session)) %>%
  left_join(session_metadata, by = c("ID", "session_norm")) %>%
  mutate(session = coalesce(session.y, session.x)) %>%
  filter(!is.na(Group), Sex %in% c("f", "m")) %>%
  mutate(
    fix = factor(
      fix,
      levels = c("fix_on_eyes", "fix_on_mouth", "fix_on_face", "fix_on_background"),
      labels = c("Eyes", "Mouth", "Face", "Background"),
      ordered = TRUE
    ),
    medication = factor(
      medication,
      levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA")
    ),
    Group = factor(Group),
    Sex = factor(Sex),
    session = factor(session)
  ) %>%
  select(ID, session, Sex, Group, medication, fix, percentage_fix_duration)

model <- glmmTMB(
  percentage_fix_duration ~ Group * fix * medication + Group * fix * Sex + session,
  family = beta_family(),
  data = df
)

anova_df <- car::Anova(model, type = 3) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term") %>%
  as_tibble()

write_anova_table_from_df(
  anova_df,
  output_dir,
  "01_main_duration_model_anova.csv",
  "main_duration_no_random"
)

write_glmmtmb_coefficients(
  model,
  output_dir,
  "02_main_duration_model_coefficients.csv",
  "main_duration_no_random"
)

emm_fix <- emmeans(model, pairwise ~ Group | fix)
emm_fix_sex <- emmeans(model, pairwise ~ Group | fix * Sex)
emm_fix_med <- emmeans(model, pairwise ~ Group | fix * medication)
emm_fix_med_sex <- emmeans(model, pairwise ~ Group | fix * medication * Sex)
emm_aoi <- pairs(regrid(emmeans(model, ~ fix), transform = "response"), adjust = "tukey")

write_csv(broom::tidy(emm_fix$contrasts), file.path(output_dir, "03_group_by_fix.csv"))
write_csv(broom::tidy(emm_fix_sex$contrasts), file.path(output_dir, "04_group_by_fix_sex.csv"))
write_csv(broom::tidy(emm_fix_med$contrasts), file.path(output_dir, "05_group_by_fix_medication.csv"))
write_csv(broom::tidy(emm_fix_med_sex$contrasts), file.path(output_dir, "06_group_by_fix_medication_sex.csv"))
write_csv(broom::tidy(emm_aoi), file.path(output_dir, "07_aoi_pairwise_response_scale.csv"))

descriptives <- df %>%
  group_by(fix) %>%
  summarise(
    mean = mean(percentage_fix_duration, na.rm = TRUE),
    sd = sd(percentage_fix_duration, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
write_csv(descriptives, file.path(output_dir, "08_aoi_descriptives.csv"))

sim_res <- DHARMa::simulateResiduals(model, plot = FALSE)

png(file.path(diag_dir, "01_dharma_overview.png"), width = 1600, height = 1200, res = 150)
plot(sim_res)
dev.off()

png(file.path(diag_dir, "02_dharma_residuals_by_predictor.png"), width = 1800, height = 1200, res = 150)
par(mfrow = c(2, 2))
plotResiduals(sim_res, df$Group, main = "Residuals by Group")
plotResiduals(sim_res, df$fix, main = "Residuals by AOI")
plotResiduals(sim_res, df$medication, main = "Residuals by Medication")
plotResiduals(sim_res, df$Sex, main = "Residuals by Sex")
dev.off()

png(file.path(diag_dir, "03_fitted_vs_observed.png"), width = 1600, height = 1200, res = 150)
plot(
  fitted(model),
  df$percentage_fix_duration,
  xlab = "Fitted values",
  ylab = "Observed fixation-duration proportion",
  pch = 16,
  col = scales::alpha("#0B7A6B", 0.25)
)
abline(0, 1, col = "red", lwd = 2)
dev.off()

diag_lines <- c(
  "Main fixation-duration model diagnostics (no random effects)",
  "",
  "Model:",
  "percentage_fix_duration ~ Group * fix * medication + Group * fix * Sex + session",
  "",
  paste("N observations:", nrow(df)),
  paste("N participants:", n_distinct(df$ID)),
  paste("AIC:", round(AIC(model), 3)),
  paste("DHARMa uniformity test p =", round(DHARMa::testUniformity(sim_res)$p.value, 4)),
  paste("DHARMa dispersion test p =", round(DHARMa::testDispersion(sim_res)$p.value, 4)),
  paste("DHARMa outlier test p =", round(DHARMa::testOutliers(sim_res)$p.value, 4)),
  paste("DHARMa zero-inflation test p =", round(DHARMa::testZeroInflation(sim_res)$p.value, 4)),
  paste("check_singularity() =", as.logical(performance::check_singularity(model))),
  paste("check_overdispersion() dispersion ratio:", round(performance::check_overdispersion(model)$dispersion_ratio, 4)),
  paste("check_overdispersion() p =", round(performance::check_overdispersion(model)$p_value, 4))
)

write_text_output(diag_lines, diag_dir, "00_diagnostic_summary.txt")
