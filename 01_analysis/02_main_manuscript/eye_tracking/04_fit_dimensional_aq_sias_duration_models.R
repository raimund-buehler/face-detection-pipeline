# Manuscript section: Dimensional eye-tracking extensions
# Analysis family: AQ/SIAS x AOI x medication x sex beta-regression models
# Primary input dataset(s): 00_data/derived/preprocessing/duration_data.csv; 00_data/derived/analysis/df_analysis_pub.csv
# Primary output(s): dimensional AQ/SIAS duration-model comparisons, omnibus tests, and AOI-specific slopes

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

tidy_coefficients <- function(model, model_name) {
  summary(model)$coefficients$cond %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    as_tibble() %>%
    rename(p_value = `Pr(>|z|)`) %>%
    mutate(
      model = model_name,
      p_display = format_p_value(p_value),
      sig = significance_stars(p_value),
      significant_0_05 = !is.na(p_value) & p_value < 0.05
    ) %>%
    relocate(model, term)
}

fit_candidate_models <- function(data, trait_col) {
  candidate_formulas <- list(
    no_trait = "percentage_fix_duration ~ fix * medication + fix * Sex + session",
    trait_main = paste0("percentage_fix_duration ~ ", trait_col, " + fix * medication + fix * Sex + session"),
    trait_aoi = paste0("percentage_fix_duration ~ ", trait_col, " * fix + fix * medication + fix * Sex + session"),
    trait_aoi_medication = paste0("percentage_fix_duration ~ ", trait_col, " * fix * medication + fix * Sex + session"),
    trait_aoi_sex = paste0("percentage_fix_duration ~ ", trait_col, " * fix * Sex + fix * medication + session"),
    main_analogue = paste0("percentage_fix_duration ~ ", trait_col, " * fix * medication + ", trait_col, " * fix * Sex + session")
  )

  map(
    candidate_formulas,
    ~ glmmTMB(
      as.formula(.x),
      family = beta_family(),
      data = data,
      control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
    )
  )
}

extract_trait_slopes <- function(model, trait_col, model_name) {
  bind_rows(
    broom::tidy(emtrends(model, ~ fix, var = trait_col)) %>%
      mutate(contrast_family = "trait_slope_by_aoi"),
    broom::tidy(emtrends(model, ~ fix * Sex, var = trait_col)) %>%
      mutate(contrast_family = "trait_slope_by_aoi_sex"),
    broom::tidy(emtrends(model, ~ fix * medication, var = trait_col)) %>%
      mutate(contrast_family = "trait_slope_by_aoi_medication")
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
  "dimensional_aq_sias_duration"
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
      levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA"),
      labels = c("BOTH", "NAL", "OXT", "PLA")
    ),
    Sex = factor(Sex),
    session = factor(session),
    AQ_z = as.numeric(scale(AQ)),
    SIAS_z = as.numeric(scale(SIAS))
  ) %>%
  select(ID, session, Sex, Group, medication, fix, percentage_fix_duration, AQ, SIAS, AQ_z, SIAS_z)

write_csv(
  df %>%
    summarise(
      n_observations = n(),
      n_participants = n_distinct(ID),
      n_participants_aq = n_distinct(ID[!is.na(AQ_z)]),
      n_participants_sias = n_distinct(ID[!is.na(SIAS_z)]),
      n_observations_aq = sum(!is.na(AQ_z)),
      n_observations_sias = sum(!is.na(SIAS_z))
    ),
  file.path(output_dir, "00_dimensional_model_sample_sizes.csv")
)

traits <- tibble(
  trait = c("AQ", "SIAS"),
  trait_col = c("AQ_z", "SIAS_z")
)

all_comparisons <- list()
all_anova <- list()
all_coefficients <- list()
all_slopes <- list()

for (i in seq_len(nrow(traits))) {
  trait_name <- traits$trait[[i]]
  trait_col <- traits$trait_col[[i]]
  model_data <- df %>% filter(!is.na(.data[[trait_col]]))
  models <- fit_candidate_models(model_data, trait_col)

  comparison <- tibble(
    trait = trait_name,
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

  write_csv(comparison, file.path(output_dir, paste0("01_", tolower(trait_name), "_model_comparison.csv")))
  all_comparisons[[trait_name]] <- comparison

  preferred_model <- models[["main_analogue"]]
  best_aic_model_name <- comparison$model[[1]]
  best_aic_model <- models[[best_aic_model_name]]

  all_anova[[paste0(trait_name, "_main_analogue")]] <- tidy_anova(preferred_model, paste0(trait_name, "_main_analogue")) %>%
    mutate(trait = trait_name, selected_model = "main_analogue")
  all_anova[[paste0(trait_name, "_best_aic")]] <- tidy_anova(best_aic_model, paste0(trait_name, "_best_aic")) %>%
    mutate(trait = trait_name, selected_model = best_aic_model_name)

  all_coefficients[[paste0(trait_name, "_main_analogue")]] <- tidy_coefficients(preferred_model, paste0(trait_name, "_main_analogue")) %>%
    mutate(trait = trait_name, selected_model = "main_analogue")
  all_coefficients[[paste0(trait_name, "_best_aic")]] <- tidy_coefficients(best_aic_model, paste0(trait_name, "_best_aic")) %>%
    mutate(trait = trait_name, selected_model = best_aic_model_name)

  all_slopes[[paste0(trait_name, "_main_analogue")]] <- extract_trait_slopes(
    preferred_model,
    trait_col,
    paste0(trait_name, "_main_analogue")
  ) %>%
    mutate(trait = trait_name, selected_model = "main_analogue")
  all_slopes[[paste0(trait_name, "_best_aic")]] <- extract_trait_slopes(
    best_aic_model,
    trait_col,
    paste0(trait_name, "_best_aic")
  ) %>%
    mutate(trait = trait_name, selected_model = best_aic_model_name)
}

write_csv(bind_rows(all_comparisons), file.path(output_dir, "01_all_model_comparisons.csv"))
write_csv(bind_rows(all_anova), file.path(output_dir, "02_dimensional_model_anova.csv"))
write_csv(bind_rows(all_coefficients), file.path(output_dir, "03_dimensional_model_coefficients.csv"))
write_csv(bind_rows(all_slopes), file.path(output_dir, "04_dimensional_trait_slopes.csv"))

summary_lines <- c(
  "Dimensional AQ/SIAS duration models",
  "",
  "Preferred confirmatory analogue:",
  "percentage_fix_duration ~ trait_z * fix * medication + trait_z * fix * Sex + session",
  "",
  "Candidate models are ranked by AIC in 01_all_model_comparisons.csv.",
  "Omnibus type-III tests for the main analogue and AIC-best model are in 02_dimensional_model_anova.csv.",
  "AOI-, sex-, and medication-specific trait slopes are in 04_dimensional_trait_slopes.csv."
)

write_text_output(summary_lines, output_dir, "00_readme.txt")
