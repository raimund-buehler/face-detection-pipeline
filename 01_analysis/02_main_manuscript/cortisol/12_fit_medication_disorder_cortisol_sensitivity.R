# Manuscript section: Cortisol sensitivity analysis
# Analysis family: medication/comorbidity type sensitivity for cortisol reactivity
# Primary input dataset(s): df_cortisol_min_max.csv; anonymized medication/disorder files; private ID lookup
# Primary output(s): simple tests of whether ASD heterogeneity moderates medication effects on T2-T3

library(tidyverse)
library(lmerTest)
library(car)
library(here)

source(here("01_analysis", "shared", "model_output_utils.R"))

read_anon_csv <- function(filename) {
  read_csv(
    here("00_data", "derived", "sensitivity", "medication_disorder_anonymized", filename),
    show_col_types = FALSE
  ) %>%
    rename_with(str_trim)
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

extract_lrt_row <- function(reduced_model, full_model, test_label, predictor_label) {
  comparison <- anova(reduced_model, full_model, test = "Chisq") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("model_step") %>%
    as_tibble()
  comparison_row <- comparison[2, ]
  deviance_value <- if ("deviance" %in% names(comparison_row)) {
    comparison_row$deviance
  } else if ("deviance()" %in% names(comparison_row)) {
    comparison_row[["deviance()"]]
  } else if ("-2*log(L)" %in% names(comparison_row)) {
    comparison_row[["-2*log(L)"]]
  } else {
    NA_real_
  }
  chi_df_value <- if ("Chi Df" %in% names(comparison_row)) {
    comparison_row[["Chi Df"]]
  } else if ("Chi.Df" %in% names(comparison_row)) {
    comparison_row[["Chi.Df"]]
  } else if ("Df" %in% names(comparison_row)) {
    comparison_row$Df
  } else if ("Df" %in% names(comparison)) {
    comparison$Df[[2]] - comparison$Df[[1]]
  } else {
    NA_real_
  }

  comparison_row %>%
    transmute(
      predictor = predictor_label,
      test = test_label,
      reduced_formula = paste(deparse(stats::formula(reduced_model)), collapse = " "),
      full_formula = paste(deparse(stats::formula(full_model)), collapse = " "),
      n_observations = nobs(full_model),
      AIC = AIC,
      BIC = BIC,
      logLik = logLik,
      deviance = deviance_value,
      Chisq = Chisq,
      Chi_Df = chi_df_value,
      p_value = `Pr(>Chisq)`,
      p_display = format_p_value(p_value),
      sig = significance_stars(p_value),
      significant_0_05 = !is.na(p_value) & p_value < 0.05
    )
}

fit_lmer_ml <- function(formula, data) {
  lmer(formula, data = data, REML = FALSE)
}

output_dir <- here(
  "02_outputs", "model_outputs", "main_manuscript", "cortisol",
  "sensitivity", "medication_disorder"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

lookup_path <- here("00_data", "private", "id_lookup.csv")
if (!file.exists(lookup_path)) {
  stop(
    "Private lookup table not found at 00_data/private/id_lookup.csv. ",
    "Run the anonymization helper first, or place the lookup there.",
    call. = FALSE
  )
}

id_lookup <- read_csv(lookup_path, show_col_types = FALSE) %>%
  transmute(
    anon_id = as.integer(anon_id),
    ID = as.character(original_id)
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

df <- read_csv(
  here("00_data", "derived", "analysis", "df_cortisol_min_max.csv"),
  show_col_types = FALSE
) %>%
  rename(T2T3 = MinMax) %>%
  mutate(ID = as.character(ID)) %>%
  left_join(id_lookup, by = "ID") %>%
  filter(!is.na(anon_id), !is.na(T2T3), !is.na(medication), Group %in% c("ASD", "CTRL")) %>%
  left_join(medication_type_groups, by = "anon_id") %>%
  left_join(disorder_type_groups, by = "anon_id") %>%
  mutate(
    medication_type_group = replace_na(medication_type_group, "No medication record"),
    disorder_type_group = replace_na(disorder_type_group, "No disorder record"),
    Group = factor(Group, levels = c("ASD", "CTRL")),
    medication = factor(
      medication,
      levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA"),
      labels = c("BOTH", "NAL", "OXT", "PLA")
    ),
    anon_id = factor(anon_id),
    asd_med_ssri = as.integer(Group == "ASD" & medication_type_group == "SSRI/SNRI only"),
    asd_med_other = as.integer(Group == "ASD" & medication_type_group == "Other/multiple medication"),
    asd_dis_adhd = as.integer(Group == "ASD" & disorder_type_group == "ADHD only"),
    asd_dis_dep = as.integer(Group == "ASD" & disorder_type_group == "Depression only"),
    asd_dis_other = as.integer(Group == "ASD" & disorder_type_group == "Other/multiple disorder")
  )

write_csv(
  df %>%
    distinct(anon_id, Group, medication_type_group, disorder_type_group) %>%
    group_by(Group) %>%
    summarise(
      n_participants = n_distinct(anon_id),
      n_ssri_snri_only = sum(medication_type_group == "SSRI/SNRI only"),
      n_other_multiple_medication = sum(medication_type_group == "Other/multiple medication"),
      n_adhd_only = sum(disorder_type_group == "ADHD only"),
      n_depression_only = sum(disorder_type_group == "Depression only"),
      n_other_multiple_disorder = sum(disorder_type_group == "Other/multiple disorder"),
      .groups = "drop"
    ),
  file.path(output_dir, "01_participant_type_counts.csv")
)

main_model <- fit_lmer_ml(
  T2T3 ~ Group * medication + (1 | anon_id),
  df
)

medication_adjusted_model <- fit_lmer_ml(
  T2T3 ~ Group * medication + asd_med_ssri + asd_med_other + (1 | anon_id),
  df
)
medication_moderation_model <- fit_lmer_ml(
  T2T3 ~ Group * medication + (asd_med_ssri + asd_med_other) * medication + (1 | anon_id),
  df
)

disorder_adjusted_model <- fit_lmer_ml(
  T2T3 ~ Group * medication + asd_dis_adhd + asd_dis_dep + asd_dis_other + (1 | anon_id),
  df
)
disorder_moderation_model <- fit_lmer_ml(
  T2T3 ~ Group * medication + (asd_dis_adhd + asd_dis_dep + asd_dis_other) * medication + (1 | anon_id),
  df
)

both_adjusted_model <- fit_lmer_ml(
  T2T3 ~ Group * medication +
    asd_med_ssri + asd_med_other + asd_dis_adhd + asd_dis_dep + asd_dis_other +
    (1 | anon_id),
  df
)
both_moderation_model <- fit_lmer_ml(
  T2T3 ~ Group * medication +
    (asd_med_ssri + asd_med_other + asd_dis_adhd + asd_dis_dep + asd_dis_other) * medication +
    (1 | anon_id),
  df
)

models <- list(
  main = main_model,
  medication_type_adjusted = medication_adjusted_model,
  medication_type_moderation = medication_moderation_model,
  disorder_type_adjusted = disorder_adjusted_model,
  disorder_type_moderation = disorder_moderation_model,
  medication_and_disorder_type_adjusted = both_adjusted_model,
  medication_and_disorder_type_moderation = both_moderation_model
)

model_summary <- imap_dfr(models, ~ tibble(
  model = .y,
  formula = paste(deparse(stats::formula(.x)), collapse = " "),
  n_observations = nobs(.x),
  n_participants = n_distinct(df$anon_id),
  AIC = AIC(.x),
  BIC = BIC(.x),
  logLik = as.numeric(logLik(.x)),
  singular = lme4::isSingular(.x)
))

anova_results <- imap_dfr(models, ~ tidy_anova(.x, .y))

key_terms <- anova_results %>%
  filter(term %in% c(
    "medication",
    "Group:medication",
    "asd_med_ssri",
    "asd_med_other",
    "asd_dis_adhd",
    "asd_dis_dep",
    "asd_dis_other",
    "medication:asd_med_ssri",
    "medication:asd_med_other",
    "medication:asd_dis_adhd",
    "medication:asd_dis_dep",
    "medication:asd_dis_other"
  )) %>%
  select(model, term, Chisq, Df, `Pr(>Chisq)`, p_display, sig, significant_0_05)

moderation_tests <- bind_rows(
  extract_lrt_row(
    medication_adjusted_model,
    medication_moderation_model,
    "Does prescribed medication type moderate the study-medication effect on T2-T3?",
    "prescribed_medication_type"
  ),
  extract_lrt_row(
    disorder_adjusted_model,
    disorder_moderation_model,
    "Does disorder type moderate the study-medication effect on T2-T3?",
    "disorder_type"
  ),
  extract_lrt_row(
    both_adjusted_model,
    both_moderation_model,
    "Do medication and disorder types jointly moderate the study-medication effect on T2-T3?",
    "medication_and_disorder_type"
  )
)

write_csv(model_summary, file.path(output_dir, "02_model_summary.csv"))
write_csv(anova_results, file.path(output_dir, "03_model_anova.csv"))
write_csv(key_terms, file.path(output_dir, "04_key_terms.csv"))
write_csv(moderation_tests, file.path(output_dir, "05_type_moderation_lrt_tests.csv"))

write_lines(
  c(
    "Cortisol sensitivity analysis: ASD medication/disorder type and T2-T3 cortisol reactivity",
    "",
    "Primary model:",
    "T2T3 ~ Group * medication + (1 | anon_id)",
    "",
    "Sensitivity question:",
    "Do ASD medication/disorder type indicators moderate the study-medication effect on T2-T3?",
    "",
    "Type indicators are ASD-specific dummy variables. Controls are coded 0 on all indicators.",
    "This avoids treating absence of ASD medication/disorder records in controls as a comparable clinical category.",
    "",
    "Primary reviewer-facing output:",
    "05_type_moderation_lrt_tests.csv"
  ),
  file.path(output_dir, "00_readme.txt")
)

print(model_summary)
print(moderation_tests)
