# Manuscript section: Main eye-tracking model sensitivity analysis
# Analysis family: crossover carryover/order sensitivity
# Primary input dataset(s): 00_data/derived/preprocessing/duration_data.csv; 00_data/derived/analysis/df_analysis_pub_prep.csv
# Primary output(s): carryover sensitivity tables for fixation-duration model

library(tidyverse)
library(glmmTMB)
library(car)
library(here)
source(here("01_analysis", "shared", "model_output_utils.R"))

normalize_session <- function(x) {
  x %>%
    str_trim() %>%
    str_to_lower() %>%
    str_replace_all(" ", "_")
}

extract_session_number <- function(x) {
  as.integer(str_extract(as.character(x), "\\d+"))
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

output_dir <- here(
  "02_outputs", "model_outputs", "main_manuscript", "eye_tracking",
  "group_aoi_medication_duration", "sensitivity", "carryover"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

session_metadata <- read_csv(
  here("00_data", "derived", "analysis", "df_analysis_pub_prep.csv"),
  show_col_types = FALSE
) %>%
  distinct(ID, session, Sex, Group, medication) %>%
  mutate(
    session_norm = normalize_session(session),
    session_number = extract_session_number(session_norm)
  )

previous_session_metadata <- session_metadata %>%
  arrange(ID, session_number) %>%
  group_by(ID) %>%
  mutate(
    previous_medication = lag(medication),
    previous_nal = case_when(
      is.na(previous_medication) ~ NA_character_,
      str_detect(previous_medication, "NAL") ~ "previous_NAL",
      TRUE ~ "no_previous_NAL"
    ),
    previous_oxt = case_when(
      is.na(previous_medication) ~ NA_character_,
      str_detect(previous_medication, "OXT") ~ "previous_OXT",
      TRUE ~ "no_previous_OXT"
    )
  ) %>%
  ungroup() %>%
  select(ID, session_norm, previous_medication, previous_nal, previous_oxt)

df <- read_csv(
  here("00_data", "derived", "preprocessing", "duration_data.csv"),
  show_col_types = FALSE
) %>%
  mutate(ID = Sub_ID, session_norm = normalize_session(session)) %>%
  left_join(session_metadata, by = c("ID", "session_norm")) %>%
  left_join(previous_session_metadata, by = c("ID", "session_norm")) %>%
  mutate(session = coalesce(session.y, session.x)) %>%
  filter(!is.na(Group), Sex %in% c("f", "m"), !is.na(previous_medication)) %>%
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
    previous_medication = factor(
      previous_medication,
      levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA")
    ),
    previous_nal = factor(previous_nal, levels = c("no_previous_NAL", "previous_NAL")),
    previous_oxt = factor(previous_oxt, levels = c("no_previous_OXT", "previous_OXT")),
    Group = factor(Group),
    Sex = factor(Sex),
    session = factor(session)
  ) %>%
  select(
    ID, session, Sex, Group, medication, previous_medication,
    previous_nal, previous_oxt, fix, percentage_fix_duration
  )

base_formula <- percentage_fix_duration ~ Group * fix * medication + Group * fix * Sex + session
previous_medication_formula <- update(base_formula, . ~ . + previous_medication)
previous_medication_aoi_formula <- update(base_formula, . ~ . + previous_medication * fix)
previous_components_formula <- update(base_formula, . ~ . + previous_nal + previous_oxt)
previous_components_aoi_formula <- update(base_formula, . ~ . + previous_nal * fix + previous_oxt * fix)

fit_beta <- function(formula) {
  glmmTMB(formula, family = beta_family(), data = df)
}

models <- list(
  base_sessions_2_to_4 = fit_beta(base_formula),
  previous_medication_main = fit_beta(previous_medication_formula),
  previous_medication_by_aoi = fit_beta(previous_medication_aoi_formula),
  previous_components_main = fit_beta(previous_components_formula),
  previous_components_by_aoi = fit_beta(previous_components_aoi_formula)
)

anova_tables <- imap_dfr(
  models,
  ~ car::Anova(.x, type = 3) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    as_tibble() %>%
    prepare_anova_table_from_df(.y)
)
write_csv(anova_tables, file.path(output_dir, "01_carryover_model_anova.csv"))

walk2(
  models,
  names(models),
  ~ write_glmmtmb_coefficients(
    .x,
    output_dir,
    paste0("02_coefficients_", .y, ".csv"),
    .y
  )
)

model_comparison <- tibble(
  model = names(models),
  formula = map_chr(models, ~ paste(deparse(formula(.x)), collapse = " ")),
  n_observations = map_int(models, nobs),
  n_participants = n_distinct(df$ID),
  aic = map_dbl(models, AIC),
  bic = map_dbl(models, BIC),
  logLik = map_dbl(models, ~ as.numeric(logLik(.x)))
)
write_csv(model_comparison, file.path(output_dir, "03_carryover_model_comparison.csv"))

prepare_lrt <- function(base_model, comparison_model, comparison_name) {
  anova(base_model, comparison_model) %>%
    as.data.frame(check.names = TRUE) %>%
    tibble::rownames_to_column("model_index") %>%
    as_tibble(.name_repair = "unique") %>%
    mutate(comparison = comparison_name) %>%
    relocate(comparison, model_index)
}

lrt_results <- bind_rows(
  prepare_lrt(
    models$base_sessions_2_to_4,
    models$previous_medication_main,
    "base_vs_previous_medication_main"
  ),
  prepare_lrt(
    models$base_sessions_2_to_4,
    models$previous_medication_by_aoi,
    "base_vs_previous_medication_by_aoi"
  ),
  prepare_lrt(
    models$base_sessions_2_to_4,
    models$previous_components_main,
    "base_vs_previous_components_main"
  ),
  prepare_lrt(
    models$base_sessions_2_to_4,
    models$previous_components_by_aoi,
    "base_vs_previous_components_by_aoi"
  )
)
write_csv(lrt_results, file.path(output_dir, "04_carryover_lrt.csv"))

carryover_terms <- anova_tables %>%
  filter(str_detect(term, "^previous_")) %>%
  select(model, term, Chisq, Df, p_value, p_display, sig, significant_0_05)
write_csv(carryover_terms, file.path(output_dir, "05_carryover_terms.csv"))

key_terms <- anova_tables %>%
  filter(term %in% c("medication", "Group:medication", "fix:medication", "Group:fix:medication")) %>%
  select(model, term, Chisq, Df, p_value, p_display, sig, significant_0_05)
write_csv(key_terms, file.path(output_dir, "06_key_medication_terms_across_models.csv"))

summary_lines <- c(
  "Carryover sensitivity analysis for the main fixation-duration model",
  "",
  "Current-session model on sessions 2-4:",
  paste(deparse(base_formula), collapse = " "),
  "",
  "Carryover predictors encode the treatment administered in the immediately preceding session.",
  "Session 1 is excluded because it has no preceding treatment.",
  "",
  paste("N observations:", nrow(df)),
  paste("N participants:", n_distinct(df$ID)),
  "",
  "Carryover term tests:",
  apply(
    carryover_terms,
    1,
    function(x) {
      paste0(
        x[["model"]], " / ", x[["term"]],
        ": chi-square(", x[["Df"]], ") = ", sprintf("%.2f", as.numeric(x[["Chisq"]])),
        ", p = ", x[["p_display"]]
      )
    }
  ),
  "",
  "Key current-medication terms across models:",
  apply(
    key_terms,
    1,
    function(x) {
      paste0(
        x[["model"]], " / ", x[["term"]],
        ": chi-square(", x[["Df"]], ") = ", sprintf("%.2f", as.numeric(x[["Chisq"]])),
        ", p = ", x[["p_display"]]
      )
    }
  )
)
write_text_output(summary_lines, output_dir, "00_carryover_sensitivity_summary.txt")
