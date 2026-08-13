# Manuscript section: Main eye-tracking model sensitivity
# Analysis family: medication/comorbidity sensitivity using anonymized inputs
# Primary input dataset(s): 00_data/derived/sensitivity/medication_disorder_anonymized/*.csv
# Primary output(s): model summaries testing ASD medication/comorbidity robustness

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

extract_group_contrasts <- function(model, model_name) {
  bind_rows(
    broom::tidy(emmeans(model, pairwise ~ Group | fix)$contrasts) %>%
      mutate(contrast_family = "group_by_aoi"),
    broom::tidy(emmeans(model, pairwise ~ Group | fix * medication)$contrasts) %>%
      mutate(contrast_family = "group_by_aoi_medication"),
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

extract_group_contrasts_safe <- function(model, model_name) {
  if (!"Group" %in% all.vars(stats::formula(model))) {
    return(tibble())
  }

  extract_group_contrasts(model, model_name)
}

read_anon_csv <- function(filename) {
  read_csv(
    here("00_data", "derived", "sensitivity", "medication_disorder_anonymized", filename),
    show_col_types = FALSE
  ) %>%
    rename_with(str_trim)
}

fit_glmmtmb_beta <- function(formula, data) {
  glmmTMB(
    formula,
    family = beta_family(),
    data = data,
    control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
  )
}

fit_sensitivity_model <- function(data, formula, model_name) {
  model <- fit_glmmtmb_beta(formula, data)
  model_formula_text <- paste(deparse(stats::formula(model)), collapse = " ")
  n_observations <- stats::nobs(model)
  model_aic <- AIC(model)
  model_bic <- BIC(model)
  model_loglik <- as.numeric(logLik(model))
  model_converged <- model$fit$convergence == 0
  model_pdhess <- isTRUE(model$sdr$pdHess)

  list(
    model_name = model_name,
    model = model,
    comparison = tibble(
      model = model_name,
      formula = model_formula_text,
      n_observations = n_observations,
      n_participants = n_distinct(data$anon_id),
      AIC = model_aic,
      BIC = model_bic,
      logLik = model_loglik,
      converged = model_converged,
      pdHess = model_pdhess
    ),
    anova = tidy_anova(model, model_name),
    coefficients = tidy_coefficients(model, model_name),
    contrasts = extract_group_contrasts_safe(model, model_name)
  )
}

extract_lrt_row <- function(reduced_model, full_model, test_label, predictor_label) {
  comparison <- anova(reduced_model, full_model, test = "Chisq") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("model_step") %>%
    as_tibble()

  comparison[2, ] %>%
    transmute(
      predictor = predictor_label,
      test = test_label,
      reduced_formula = paste(deparse(stats::formula(reduced_model)), collapse = " "),
      full_formula = paste(deparse(stats::formula(full_model)), collapse = " "),
      Df = Df,
      AIC = AIC,
      BIC = BIC,
      logLik = logLik,
      deviance = deviance,
      Chisq = Chisq,
      Chi_Df = `Chi Df`,
      p_value = `Pr(>Chisq)`,
      p_display = format_p_value(p_value),
      sig = significance_stars(p_value),
      significant_0_05 = !is.na(p_value) & p_value < 0.05
    )
}

fit_type_influence_tests <- function(data, predictor, predictor_label) {
  base_formula <- as.formula(
    "percentage_fix_duration ~ fix * medication + fix * Sex + session"
  )
  overall_formula <- as.formula(paste(
    "percentage_fix_duration ~ fix * medication + fix * Sex + session +",
    predictor
  ))
  aoi_formula <- as.formula(paste(
    "percentage_fix_duration ~ fix * medication + fix * Sex + session +",
    predictor, "* fix"
  ))
  medication_formula <- as.formula(paste(
    "percentage_fix_duration ~ fix * medication + fix * Sex + session +",
    predictor, "* fix +", predictor, "* medication"
  ))
  responsiveness_formula <- as.formula(paste(
    "percentage_fix_duration ~ fix * medication + fix * Sex + session +",
    predictor, "* fix * medication"
  ))

  models <- list(
    base = fit_glmmtmb_beta(base_formula, data),
    overall = fit_glmmtmb_beta(overall_formula, data),
    aoi = fit_glmmtmb_beta(aoi_formula, data),
    medication = fit_glmmtmb_beta(medication_formula, data),
    responsiveness = fit_glmmtmb_beta(responsiveness_formula, data)
  )

  list(
    models = models,
    model_summary = imap_dfr(models, ~ tibble(
      predictor = predictor_label,
      model_step = .y,
      formula = paste(deparse(stats::formula(.x)), collapse = " "),
      n_observations = stats::nobs(.x),
      n_participants = n_distinct(data$anon_id),
      AIC = AIC(.x),
      BIC = BIC(.x),
      logLik = as.numeric(logLik(.x)),
      converged = .x$fit$convergence == 0,
      pdHess = isTRUE(.x$sdr$pdHess)
    )),
    lrt = bind_rows(
      extract_lrt_row(models$base, models$overall, "overall gaze profile", predictor_label),
      extract_lrt_row(models$overall, models$aoi, "AOI profile", predictor_label),
      extract_lrt_row(models$aoi, models$medication, "study-medication response", predictor_label),
      extract_lrt_row(models$medication, models$responsiveness, "AOI-specific study-medication response", predictor_label)
    ),
    anova = imap_dfr(models, ~ tidy_anova(.x, paste(predictor_label, .y, sep = "_"))),
    coefficients = imap_dfr(models, ~ tidy_coefficients(.x, paste(predictor_label, .y, sep = "_")))
  )
}

output_dir <- here(
  "02_outputs", "model_outputs", "main_manuscript", "eye_tracking",
  "group_aoi_medication_duration", "sensitivity", "medication_disorder"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

session_metadata <- read_anon_csv("df_analysis_pub_prep_anonymized.csv") %>%
  transmute(
    anon_id,
    session,
    session_norm = normalize_session(session),
    Sex,
    Group,
    medication
  ) %>%
  distinct()

medication_flags <- read_anon_csv("df_meds_anonymized.csv") %>%
  distinct(anon_id, meds_name, type) %>%
  group_by(anon_id) %>%
  summarise(
    medication_record_count = n(),
    medication_type_count = n_distinct(type, na.rm = TRUE),
    any_medication_record = TRUE,
    .groups = "drop"
  )

medication_type_groups <- read_anon_csv("df_meds_anonymized.csv") %>%
  distinct(anon_id, type) %>%
  mutate(
    medication_type_class = case_when(
      type == "SSRI/SNRI" ~ "SSRI/SNRI",
      TRUE ~ "Other/multiple medication"
    )
  ) %>%
  distinct(anon_id, medication_type_class) %>%
  group_by(anon_id) %>%
  summarise(
    medication_type_group = case_when(
      n_distinct(medication_type_class) == 1 & first(medication_type_class) == "SSRI/SNRI" ~ "SSRI/SNRI only",
      TRUE ~ "Other/multiple medication"
    ),
    .groups = "drop"
  )

disorder_flags <- read_anon_csv("df_disorders_anonymized.csv") %>%
  distinct(anon_id, name_clean) %>%
  group_by(anon_id) %>%
  summarise(
    disorder_record_count = n(),
    any_disorder_record = TRUE,
    .groups = "drop"
  )

disorder_type_groups <- read_anon_csv("df_disorders_anonymized.csv") %>%
  distinct(anon_id, name_clean) %>%
  mutate(
    disorder_type_class = case_when(
      name_clean == "ADHS" ~ "ADHD",
      name_clean == "Depression" ~ "Depression",
      TRUE ~ "Other/multiple disorder"
    )
  ) %>%
  distinct(anon_id, disorder_type_class) %>%
  group_by(anon_id) %>%
  summarise(
    disorder_type_group = case_when(
      n_distinct(disorder_type_class) == 1 & first(disorder_type_class) == "ADHD" ~ "ADHD only",
      n_distinct(disorder_type_class) == 1 & first(disorder_type_class) == "Depression" ~ "Depression only",
      TRUE ~ "Other/multiple disorder"
    ),
    .groups = "drop"
  )

participant_flags <- session_metadata %>%
  distinct(anon_id, Group, Sex) %>%
  left_join(medication_flags, by = "anon_id") %>%
  left_join(medication_type_groups, by = "anon_id") %>%
  left_join(disorder_flags, by = "anon_id") %>%
  left_join(disorder_type_groups, by = "anon_id") %>%
  mutate(
    across(
      c(medication_record_count, medication_type_count, disorder_record_count),
      ~ replace_na(.x, 0L)
    ),
    any_medication_record = replace_na(any_medication_record, FALSE),
    any_disorder_record = replace_na(any_disorder_record, FALSE),
    medication_type_group = replace_na(medication_type_group, "No medication record"),
    disorder_type_group = replace_na(disorder_type_group, "No disorder record"),
    asd_medication_record = Group == "ASD" & any_medication_record,
    asd_disorder_record = Group == "ASD" & any_disorder_record,
    asd_medication_or_disorder_record = asd_medication_record | asd_disorder_record
  )

df <- read_anon_csv("duration_data_anonymized.csv") %>%
  mutate(session_norm = normalize_session(session)) %>%
  left_join(session_metadata, by = c("anon_id", "session_norm")) %>%
  mutate(session = coalesce(session.y, session.x)) %>%
  left_join(
    participant_flags %>%
      select(
        anon_id,
        any_medication_record,
        any_disorder_record,
        asd_medication_record,
        asd_disorder_record,
        asd_medication_or_disorder_record,
        medication_type_group,
        disorder_type_group
      ),
    by = "anon_id"
  ) %>%
  filter(!is.na(Group), Sex %in% c("f", "m"), !is.na(percentage_fix_duration)) %>%
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
    session = factor(session),
    asd_medication_record_factor = factor(asd_medication_record, levels = c(FALSE, TRUE)),
    asd_disorder_record_factor = factor(asd_disorder_record, levels = c(FALSE, TRUE)),
    any_medication_record_factor = factor(any_medication_record, levels = c(FALSE, TRUE)),
    any_disorder_record_factor = factor(any_disorder_record, levels = c(FALSE, TRUE)),
    medication_type_group = factor(
      medication_type_group,
      levels = c("No medication record", "SSRI/SNRI only", "Other/multiple medication")
    ),
    disorder_type_group = factor(
      disorder_type_group,
      levels = c("No disorder record", "ADHD only", "Depression only", "Other/multiple disorder")
    ),
    asd_med_ssri = as.integer(Group == "ASD" & medication_type_group == "SSRI/SNRI only"),
    asd_med_other = as.integer(Group == "ASD" & medication_type_group == "Other/multiple medication"),
    asd_dis_adhd = as.integer(Group == "ASD" & disorder_type_group == "ADHD only"),
    asd_dis_dep = as.integer(Group == "ASD" & disorder_type_group == "Depression only"),
    asd_dis_other = as.integer(Group == "ASD" & disorder_type_group == "Other/multiple disorder")
  ) %>%
  select(
    anon_id,
    session,
    Sex,
    Group,
    medication,
    fix,
    percentage_fix_duration,
    any_medication_record,
    any_disorder_record,
    asd_medication_record,
    asd_disorder_record,
    asd_medication_or_disorder_record,
    asd_medication_record_factor,
    asd_disorder_record_factor,
    any_medication_record_factor,
    any_disorder_record_factor,
    medication_type_group,
    disorder_type_group,
    asd_med_ssri,
    asd_med_other,
    asd_dis_adhd,
    asd_dis_dep,
    asd_dis_other
  )

scenario_definitions <- tribble(
  ~model, ~description,
  "full_anonymized", "Main duration model fitted to the anonymized sensitivity dataset.",
  "exclude_asd_medication_records", "Main duration model excluding ASD participants with any medication record.",
  "exclude_asd_disorder_records", "Main duration model excluding ASD participants with any comorbidity/disorder record.",
  "exclude_asd_medication_or_disorder_records", "Main duration model excluding ASD participants with any medication or comorbidity/disorder record.",
  "adjust_asd_medication_disorder_records", "Main duration model with ASD medication and comorbidity/disorder flags added as covariates."
)

write_csv(
  participant_flags %>%
    group_by(Group) %>%
    summarise(
      n_participants = n_distinct(anon_id),
      n_any_medication_record = sum(any_medication_record),
      n_any_disorder_record = sum(any_disorder_record),
      n_any_medication_or_disorder_record = sum(any_medication_record | any_disorder_record),
      .groups = "drop"
    ),
  file.path(output_dir, "01_participant_flag_summary.csv")
)

write_csv(
  participant_flags %>%
    filter(Group == "ASD") %>%
    count(medication_type_group, name = "n_participants"),
  file.path(output_dir, "10_asd_medication_type_group_counts.csv")
)

write_csv(
  participant_flags %>%
    filter(Group == "ASD") %>%
    count(disorder_type_group, name = "n_participants"),
  file.path(output_dir, "11_asd_disorder_type_group_counts.csv")
)

sex_balance <- bind_rows(
  participant_flags %>%
    filter(Group == "ASD", Sex %in% c("f", "m")) %>%
    count(predictor = "medication_type", group = medication_type_group, Sex, name = "n"),
  participant_flags %>%
    filter(Group == "ASD", Sex %in% c("f", "m")) %>%
    count(predictor = "disorder_type", group = disorder_type_group, Sex, name = "n")
)
write_csv(sex_balance, file.path(output_dir, "16_type_group_sex_counts.csv"))

sex_balance_tests <- bind_rows(
  tibble(
    predictor = "medication_type",
    test = "Fisher exact test of type group by sex in ASD participants",
    p_value = fisher.test(table(
      participant_flags %>% filter(Group == "ASD", Sex %in% c("f", "m")) %>% pull(medication_type_group),
      participant_flags %>% filter(Group == "ASD", Sex %in% c("f", "m")) %>% pull(Sex)
    ))$p.value
  ),
  tibble(
    predictor = "disorder_type",
    test = "Fisher exact test of type group by sex in ASD participants",
    p_value = fisher.test(table(
      participant_flags %>% filter(Group == "ASD", Sex %in% c("f", "m")) %>% pull(disorder_type_group),
      participant_flags %>% filter(Group == "ASD", Sex %in% c("f", "m")) %>% pull(Sex)
    ))$p.value
  )
) %>%
  mutate(
    p_display = format_p_value(p_value),
    sig = significance_stars(p_value),
    significant_0_05 = !is.na(p_value) & p_value < 0.05
  )
write_csv(sex_balance_tests, file.path(output_dir, "17_type_group_sex_balance_tests.csv"))

scenario_samples <- list(
  full_anonymized = df,
  exclude_asd_medication_records = df %>% filter(!(Group == "ASD" & any_medication_record)),
  exclude_asd_disorder_records = df %>% filter(!(Group == "ASD" & any_disorder_record)),
  exclude_asd_medication_or_disorder_records = df %>%
    filter(!(Group == "ASD" & asd_medication_or_disorder_record))
)

sample_summary <- imap_dfr(
  scenario_samples,
  ~ .x %>%
    distinct(anon_id, Group, Sex) %>%
    count(Group, Sex, name = "n_participants") %>%
    mutate(model = .y, .before = 1)
)
write_csv(sample_summary, file.path(output_dir, "02_model_sample_summary.csv"))

main_formula <- percentage_fix_duration ~ Group * fix * medication + Group * fix * Sex + session
adjusted_formula <- percentage_fix_duration ~ Group * fix * medication + Group * fix * Sex +
  asd_medication_record_factor + asd_disorder_record_factor + session
full_sample_influence_formula <- percentage_fix_duration ~ Group * fix * medication + Group * fix * Sex +
  asd_medication_record_factor * fix + asd_disorder_record_factor * fix + session
asd_medication_influence_formula <- percentage_fix_duration ~ medication * fix + any_medication_record_factor * fix + Sex + session
asd_disorder_influence_formula <- percentage_fix_duration ~ medication * fix + any_disorder_record_factor * fix + Sex + session
asd_joint_influence_formula <- percentage_fix_duration ~ medication * fix +
  any_medication_record_factor * fix + any_disorder_record_factor * fix + Sex + session
asd_medication_responsiveness_formula <- percentage_fix_duration ~ any_medication_record_factor * fix * medication + Sex + session
asd_disorder_responsiveness_formula <- percentage_fix_duration ~ any_disorder_record_factor * fix * medication + Sex + session

asd_df <- df %>%
  filter(Group == "ASD") %>%
  mutate(
    any_medication_record_factor = droplevels(any_medication_record_factor),
    any_disorder_record_factor = droplevels(any_disorder_record_factor),
    medication_type_group = droplevels(medication_type_group),
    disorder_type_group = droplevels(disorder_type_group)
  )

type_tests <- list(
  medication_type = fit_type_influence_tests(
    asd_df,
    "medication_type_group",
    "medication_type"
  ),
  disorder_type = fit_type_influence_tests(
    asd_df,
    "disorder_type_group",
    "disorder_type"
  )
)

full_sample_type_adjustment_models <- list(
  main = fit_glmmtmb_beta(
    percentage_fix_duration ~ Group * fix * medication + Group * fix * Sex + session,
    df
  ),
  medication_type_aoi_adjusted = fit_glmmtmb_beta(
    percentage_fix_duration ~ Group * fix * medication + Group * fix * Sex + session +
      (asd_med_ssri + asd_med_other) * fix,
    df
  ),
  disorder_type_aoi_adjusted = fit_glmmtmb_beta(
    percentage_fix_duration ~ Group * fix * medication + Group * fix * Sex + session +
      (asd_dis_adhd + asd_dis_dep + asd_dis_other) * fix,
    df
  ),
  medication_and_disorder_type_aoi_adjusted = fit_glmmtmb_beta(
    percentage_fix_duration ~ Group * fix * medication + Group * fix * Sex + session +
      (asd_med_ssri + asd_med_other + asd_dis_adhd + asd_dis_dep + asd_dis_other) * fix,
    df
  )
)

full_sample_type_adjustment_key_terms <- imap_dfr(
  full_sample_type_adjustment_models,
  ~ tidy_anova(.x, .y)
) %>%
  filter(term %in% c(
    "Group:fix",
    "Group:fix:medication",
    "Group:fix:Sex",
    "asd_med_ssri",
    "asd_med_other",
    "fix:asd_med_ssri",
    "fix:asd_med_other",
    "asd_dis_adhd",
    "asd_dis_dep",
    "asd_dis_other",
    "fix:asd_dis_adhd",
    "fix:asd_dis_dep",
    "fix:asd_dis_other"
  )) %>%
  select(model, term, Chisq, Df, `Pr(>Chisq)`, p_display, sig, significant_0_05)

results <- list(
  full_anonymized = fit_sensitivity_model(
    scenario_samples$full_anonymized,
    main_formula,
    "full_anonymized"
  ),
  exclude_asd_medication_records = fit_sensitivity_model(
    scenario_samples$exclude_asd_medication_records,
    main_formula,
    "exclude_asd_medication_records"
  ),
  exclude_asd_disorder_records = fit_sensitivity_model(
    scenario_samples$exclude_asd_disorder_records,
    main_formula,
    "exclude_asd_disorder_records"
  ),
  exclude_asd_medication_or_disorder_records = fit_sensitivity_model(
    scenario_samples$exclude_asd_medication_or_disorder_records,
    main_formula,
    "exclude_asd_medication_or_disorder_records"
  ),
  adjust_asd_medication_disorder_records = fit_sensitivity_model(
    df,
    adjusted_formula,
    "adjust_asd_medication_disorder_records"
  ),
  full_sample_influence = fit_sensitivity_model(
    df,
    full_sample_influence_formula,
    "full_sample_influence"
  ),
  asd_medication_influence = fit_sensitivity_model(
    asd_df,
    asd_medication_influence_formula,
    "asd_medication_influence"
  ),
  asd_disorder_influence = fit_sensitivity_model(
    asd_df,
    asd_disorder_influence_formula,
    "asd_disorder_influence"
  ),
  asd_joint_influence = fit_sensitivity_model(
    asd_df,
    asd_joint_influence_formula,
    "asd_joint_influence"
  ),
  asd_medication_responsiveness = fit_sensitivity_model(
    asd_df,
    asd_medication_responsiveness_formula,
    "asd_medication_responsiveness"
  ),
  asd_disorder_responsiveness = fit_sensitivity_model(
    asd_df,
    asd_disorder_responsiveness_formula,
    "asd_disorder_responsiveness"
  )
)

model_comparison <- map_dfr(results, "comparison") %>%
  arrange(AIC) %>%
  mutate(
    delta_AIC = AIC - min(AIC),
    delta_BIC = BIC - min(BIC)
  )
anova_results <- map_dfr(results, "anova")
coefficient_results <- map_dfr(results, "coefficients")
contrast_results <- map_dfr(results, "contrasts")

key_terms <- anova_results %>%
  filter(term %in% c(
    "Group",
    "fix",
    "medication",
    "Sex",
    "Group:fix",
    "Group:medication",
    "fix:medication",
    "Group:fix:medication",
    "Group:fix:Sex",
    "asd_medication_record_factor",
    "asd_disorder_record_factor",
    "asd_medication_record_factor:fix",
    "asd_disorder_record_factor:fix",
    "any_medication_record_factor",
    "any_disorder_record_factor",
    "any_medication_record_factor:fix",
    "any_disorder_record_factor:fix",
    "any_medication_record_factor:medication",
    "any_disorder_record_factor:medication",
    "any_medication_record_factor:fix:medication",
    "any_disorder_record_factor:fix:medication"
  )) %>%
  select(model, term, Chisq, Df, `Pr(>Chisq)`, p_display, sig, significant_0_05) %>%
  arrange(term, model)

influence_terms <- anova_results %>%
  filter(str_detect(term, "asd_medication_record|asd_disorder_record|any_medication_record|any_disorder_record")) %>%
  select(model, term, Chisq, Df, `Pr(>Chisq)`, p_display, sig, significant_0_05) %>%
  arrange(model, term)

write_csv(scenario_definitions, file.path(output_dir, "00_readme_scenarios.csv"))
write_csv(model_comparison, file.path(output_dir, "03_model_comparison.csv"))
write_csv(anova_results, file.path(output_dir, "04_model_anova.csv"))
write_csv(key_terms, file.path(output_dir, "05_key_terms.csv"))
write_csv(coefficient_results, file.path(output_dir, "06_model_coefficients.csv"))
write_csv(contrast_results, file.path(output_dir, "07_group_contrasts.csv"))
write_csv(influence_terms, file.path(output_dir, "09_influence_terms.csv"))
write_csv(map_dfr(type_tests, "model_summary"), file.path(output_dir, "12_type_model_summary.csv"))
write_csv(map_dfr(type_tests, "lrt"), file.path(output_dir, "13_type_lrt_tests.csv"))
write_csv(map_dfr(type_tests, "anova"), file.path(output_dir, "14_type_model_anova.csv"))
write_csv(map_dfr(type_tests, "coefficients"), file.path(output_dir, "15_type_model_coefficients.csv"))
write_csv(
  full_sample_type_adjustment_key_terms,
  file.path(output_dir, "19_full_sample_type_adjusted_key_terms.csv")
)

summary_lines <- c(
  "Medication/comorbidity sensitivity analysis for the main fixation-duration model",
  "",
  "Input data are anonymized and linked only through anon_id.",
  "",
  "Primary model:",
  "percentage_fix_duration ~ Group * fix * medication + Group * fix * Sex + session",
  "",
  "Sensitivity checks:",
  paste(scenario_definitions$model, scenario_definitions$description, sep = ": "),
  "",
  "Participant flags by group:",
  capture.output(print(read_csv(file.path(output_dir, "01_participant_flag_summary.csv"), show_col_types = FALSE))),
  "",
  "Key omnibus terms:",
  capture.output(print(key_terms %>% filter(term %in% c("Group:fix", "Group:fix:medication", "Group:fix:Sex")))),
  "",
  "Medication/comorbidity influence terms:",
  capture.output(print(influence_terms))
)
write_lines(summary_lines, file.path(output_dir, "08_interpretation_note.txt"))

print(model_comparison)
print(key_terms)
