# Manuscript section: Main eye-tracking model sensitivity
# Analysis family: AQ-adjusted group x AOI x sex sensitivity
# Primary input dataset(s): 00_data/derived/preprocessing/duration_data.csv; 00_data/derived/analysis/df_analysis_pub.csv
# Primary output(s): AQ imbalance checks and AQ-adjusted duration-model sensitivity summaries

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

beta_squeeze <- function(x) {
  n <- sum(!is.na(x))
  ((x * (n - 1)) + 0.5) / n
}

tidy_anova <- function(model, model_name) {
  car::Anova(model, type = 3) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    as_tibble() %>%
    mutate(
      model = model_name,
      p_display = format_p_value(`Pr(>Chisq)`),
      sig = significance_stars(`Pr(>Chisq)`),
      significant_0_05 = !is.na(`Pr(>Chisq)`) & `Pr(>Chisq)` < 0.05
    ) %>%
    relocate(model, term)
}

tidy_lm_anova <- function(model, model_name) {
  car::Anova(model, type = 3) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    as_tibble() %>%
    mutate(
      model = model_name,
      p_display = format_p_value(`Pr(>F)`),
      sig = significance_stars(`Pr(>F)`),
      significant_0_05 = !is.na(`Pr(>F)`) & `Pr(>F)` < 0.05
    ) %>%
    relocate(model, term)
}

extract_group_contrasts <- function(model, model_name) {
  bind_rows(
    broom::tidy(emmeans(model, pairwise ~ Group | fix)$contrasts) %>%
      mutate(contrast_family = "group_by_aoi"),
    broom::tidy(emmeans(model, pairwise ~ Group | fix * Sex)$contrasts) %>%
      mutate(contrast_family = "group_by_aoi_sex")
  ) %>%
    mutate(
      model = model_name,
      p_display = format_p_value(p.value),
      sig = significance_stars(p.value),
      significant_0_05 = !is.na(p.value) & p.value < 0.05
    ) %>%
    relocate(model, contrast_family)
}

extract_aq_slopes <- function(model, model_name) {
  bind_rows(
    broom::tidy(emtrends(model, ~ fix, var = "AQ_z")) %>%
      mutate(contrast_family = "aq_slope_by_aoi"),
    broom::tidy(emtrends(model, ~ fix * Sex, var = "AQ_z")) %>%
      mutate(contrast_family = "aq_slope_by_aoi_sex")
  ) %>%
    mutate(
      model = model_name,
      p_display = format_p_value(p.value),
      sig = significance_stars(p.value),
      significant_0_05 = !is.na(p.value) & p.value < 0.05
    ) %>%
    relocate(model, contrast_family)
}

output_dir <- here(
  "02_outputs", "model_outputs", "main_manuscript", "eye_tracking",
  "group_aoi_medication_duration", "sensitivity", "aq_adjusted_group_sex"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

session_metadata <- read_csv(
  here("00_data", "derived", "analysis", "df_analysis_pub.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    ID,
    session,
    session_norm = normalize_session(session),
    Sex,
    Group,
    medication,
    AQ = `AQ-K Score`,
    SIAS = `SIAS Score`
  ) %>%
  distinct()

df <- read_csv(
  here("00_data", "derived", "preprocessing", "duration_data.csv"),
  show_col_types = FALSE
) %>%
  mutate(ID = Sub_ID, session_norm = normalize_session(session)) %>%
  left_join(session_metadata, by = c("ID", "session_norm")) %>%
  mutate(session = coalesce(session.y, session.x)) %>%
  filter(!is.na(Group), Sex %in% c("f", "m"), !is.na(AQ)) %>%
  mutate(
    fix = factor(
      fix,
      levels = c("fix_on_eyes", "fix_on_mouth", "fix_on_face", "fix_on_background"),
      labels = c("Eyes", "Mouth", "Face", "Background"),
      ordered = TRUE
    ),
    medication = factor(
      medication,
      levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA"),
      labels = c("BOTH", "NAL", "OXT", "PLA")
    ),
    Group = factor(Group),
    Sex = factor(Sex),
    session = factor(session),
    AQ_z = as.numeric(scale(AQ)),
    percentage_fix_duration_beta = beta_squeeze(percentage_fix_duration)
  ) %>%
  select(ID, session, Sex, Group, medication, fix, percentage_fix_duration, percentage_fix_duration_beta, AQ, AQ_z, SIAS)

participant_aq <- df %>%
  distinct(ID, Group, Sex, AQ, AQ_z, SIAS)

write_csv(
  participant_aq %>%
    group_by(Group, Sex) %>%
    summarise(
      n = n(),
      mean_AQ = mean(AQ, na.rm = TRUE),
      sd_AQ = sd(AQ, na.rm = TRUE),
      median_AQ = median(AQ, na.rm = TRUE),
      min_AQ = min(AQ, na.rm = TRUE),
      max_AQ = max(AQ, na.rm = TRUE),
      .groups = "drop"
    ),
  file.path(output_dir, "01_aq_descriptives_by_group_sex.csv")
)

aq_imbalance_model <- lm(AQ ~ Group * Sex, data = participant_aq)
write_csv(
  tidy_lm_anova(aq_imbalance_model, "aq_imbalance_group_sex"),
  file.path(output_dir, "02_aq_imbalance_anova.csv")
)
write_csv(
  broom::tidy(emmeans(aq_imbalance_model, pairwise ~ Group | Sex)$contrasts) %>%
    mutate(p_display = format_p_value(p.value), sig = significance_stars(p.value)),
  file.path(output_dir, "03_aq_group_contrasts_by_sex.csv")
)
write_csv(
  broom::tidy(emmeans(aq_imbalance_model, pairwise ~ Sex | Group)$contrasts) %>%
    mutate(p_display = format_p_value(p.value), sig = significance_stars(p.value)),
  file.path(output_dir, "04_aq_sex_contrasts_by_group.csv")
)

models <- list(
  aq_complete_main = glmmTMB(
    percentage_fix_duration_beta ~ Group * fix * medication + Group * fix * Sex + session,
    family = beta_family(),
    data = df,
    control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
  ),
  aq_adjusted = glmmTMB(
    percentage_fix_duration_beta ~ Group * fix * medication + Group * fix * Sex + AQ_z * fix + session,
    family = beta_family(),
    data = df,
    control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
  ),
  aq_adjusted_by_sex = glmmTMB(
    percentage_fix_duration_beta ~ Group * fix * medication + Group * fix * Sex + AQ_z * fix * Sex + session,
    family = beta_family(),
    data = df,
    control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
  )
)

model_comparison <- tibble(
  model = names(models),
  formula = map_chr(models, ~ paste(deparse(formula(.x)), collapse = " ")),
  n_observations = map_int(models, nobs),
  AIC = map_dbl(models, AIC),
  BIC = map_dbl(models, BIC),
  logLik = map_dbl(models, ~ as.numeric(logLik(.x)))
) %>%
  arrange(AIC) %>%
  mutate(
    delta_AIC = AIC - min(AIC),
    delta_BIC = BIC - min(BIC),
    AIC_rank = row_number()
  )
write_csv(model_comparison, file.path(output_dir, "05_model_comparison.csv"))

write_csv(
  imap_dfr(models, ~ tidy_anova(.x, .y)),
  file.path(output_dir, "06_model_anova.csv")
)

write_csv(
  imap_dfr(models, ~ extract_group_contrasts(.x, .y)),
  file.path(output_dir, "07_group_contrasts.csv")
)

write_csv(
  imap_dfr(models[c("aq_adjusted", "aq_adjusted_by_sex")], ~ extract_aq_slopes(.x, .y)),
  file.path(output_dir, "08_aq_slopes.csv")
)

summary_lines <- c(
  "AQ-adjusted Group x AOI x Sex sensitivity",
  "",
  "AQ imbalance model:",
  "AQ ~ Group * Sex",
  "",
  "Duration sensitivity models:",
  "aq_complete_main: percentage_fix_duration_beta ~ Group * fix * medication + Group * fix * Sex + session",
  "aq_adjusted: percentage_fix_duration_beta ~ Group * fix * medication + Group * fix * Sex + AQ_z * fix + session",
  "aq_adjusted_by_sex: percentage_fix_duration_beta ~ Group * fix * medication + Group * fix * Sex + AQ_z * fix * Sex + session",
  "",
  "AQ_z is standardized across the AQ-complete analysis sample.",
  "Outputs 06-08 contain the omnibus model tests, group contrasts, and AQ simple slopes."
)
write_text_output(summary_lines, output_dir, "00_readme.txt")
