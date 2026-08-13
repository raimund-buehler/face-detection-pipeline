# Manuscript section: Exploratory cortisol sensitivity
# Analysis family: AUCi eye-duration model selection
# Primary input dataset(s): 00_data/derived/analysis/df_cortisol_min_max.csv; 00_data/derived/preprocessing/duration_data.csv
# Primary output(s): fixed/random structure model selection for AUCi -> eye fixation-duration proportion

library(tidyverse)
library(glmmTMB)
library(car)
library(emmeans)
library(performance)
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

safe_fit <- function(formula, data) {
  tryCatch(
    glmmTMB(
      formula,
      family = beta_family(),
      data = data,
      control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
    ),
    error = function(e) e
  )
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

extract_group_slopes <- function(model, model_name) {
  broom::tidy(emtrends(model, ~ Group, var = "AUCi_z")) %>%
    mutate(
      model = model_name,
      p_display = format_p_value(p.value),
      sig = significance_stars(p.value),
      significant_0_05 = !is.na(p.value) & p.value < 0.05
    ) %>%
    relocate(model)
}

extract_group_contrast_of_slopes <- function(model, model_name) {
  contrast(emtrends(model, ~ Group, var = "AUCi_z"), method = "pairwise") %>%
    broom::tidy() %>%
    mutate(
      model = model_name,
      p_display = format_p_value(p.value),
      sig = significance_stars(p.value),
      significant_0_05 = !is.na(p.value) & p.value < 0.05
    ) %>%
    relocate(model)
}

output_dir <- here(
  "02_outputs", "model_outputs", "main_manuscript", "cortisol",
  "duration_exploratory", "model_selection"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

medication_levels_raw <- c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA")
medication_levels_plot <- c("PLA", "OXT", "NAL", "BOTH")

duration_eyes <- read_csv(
  here("00_data", "derived", "preprocessing", "duration_data.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    ID = Sub_ID,
    session_norm = normalize_session(session),
    fix,
    percentage_fix_duration
  ) %>%
  filter(fix == "fix_on_eyes")

df <- read_csv(
  here("00_data", "derived", "analysis", "df_cortisol_min_max.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    session_norm = normalize_session(session),
    medication = factor(
      medication,
      levels = medication_levels_raw,
      labels = c("BOTH", "NAL", "OXT", "PLA")
    ),
    medication = factor(medication, levels = medication_levels_plot),
    Group = factor(Group, levels = c("ASD", "CTRL")),
    ID = factor(ID),
    session = factor(session)
  ) %>%
  left_join(duration_eyes, by = c("ID", "session_norm")) %>%
  filter(
    !is.na(percentage_fix_duration),
    !is.na(AUCi),
    !is.na(Group),
    !is.na(medication)
  ) %>%
  mutate(
    AUCi_z = as.numeric(scale(AUCi)),
    percentage_fix_duration_beta = beta_squeeze(percentage_fix_duration)
  )

write_csv(
  df %>%
    summarise(
      n_observations = n(),
      n_participants = n_distinct(ID),
      min_sessions_per_participant = min(as.numeric(table(ID))),
      median_sessions_per_participant = median(as.numeric(table(ID))),
      max_sessions_per_participant = max(as.numeric(table(ID)))
    ),
  file.path(output_dir, "00_sample_summary.csv")
)

fixed_terms <- tribble(
  ~fixed_name, ~fixed_rhs,
  "group_auci", "Group * AUCi_z",
  "group_auci_med_add", "Group * AUCi_z + medication",
  "group_auci_session_add", "Group * AUCi_z + session",
  "group_auci_med_session_add", "Group * AUCi_z + medication + session",
  "group_auci_med_interaction", "Group * AUCi_z + Group * medication",
  "group_auci_by_medication", "Group * AUCi_z * medication",
  "group_auci_session_med", "Group * AUCi_z + medication + session"
)

random_terms <- tribble(
  ~random_name, ~random_rhs,
  "none", "",
  "id_intercept", " + (1 | ID)",
  "id_medication_slope_diag", " + diag(1 + medication | ID)",
  "id_session_slope_diag", " + diag(1 + session | ID)"
)

model_grid <- crossing(fixed_terms, random_terms) %>%
  mutate(
    model_name = paste(fixed_name, random_name, sep = "__"),
    formula_text = paste0("percentage_fix_duration_beta ~ ", fixed_rhs, random_rhs),
    formula = map(formula_text, as.formula)
  )

fits <- model_grid %>%
  mutate(
    fit = map(formula, safe_fit, data = df),
    fit_ok = map_lgl(fit, ~ inherits(.x, "glmmTMB")),
    error_message = map_chr(fit, ~ if (inherits(.x, "glmmTMB")) "" else .x$message),
    converged = map_lgl(fit, ~ if (inherits(.x, "glmmTMB")) isTRUE(.x$sdr$pdHess) else FALSE),
    singular = map_lgl(fit, ~ if (inherits(.x, "glmmTMB")) {
      out <- tryCatch(as.logical(performance::check_singularity(.x)), error = function(e) NA)
      out
    } else NA),
    n_observations = map_int(fit, ~ if (inherits(.x, "glmmTMB")) nobs(.x) else NA_integer_),
    AIC = map_dbl(fit, ~ if (inherits(.x, "glmmTMB")) AIC(.x) else NA_real_),
    BIC = map_dbl(fit, ~ if (inherits(.x, "glmmTMB")) BIC(.x) else NA_real_),
    logLik = map_dbl(fit, ~ if (inherits(.x, "glmmTMB")) as.numeric(logLik(.x)) else NA_real_)
  )

comparison <- fits %>%
  select(
    model_name, fixed_name, random_name, formula_text, fit_ok, converged, singular,
    n_observations, AIC, BIC, logLik, error_message
  ) %>%
  arrange(AIC) %>%
  mutate(
    delta_AIC = AIC - min(AIC, na.rm = TRUE),
    delta_BIC = BIC - min(BIC, na.rm = TRUE),
    AIC_rank = row_number()
  )
write_csv(comparison, file.path(output_dir, "01_model_selection_comparison.csv"))

valid_fits <- fits %>%
  filter(fit_ok, converged) %>%
  arrange(AIC)

anova_results <- valid_fits %>%
  mutate(anova = map2(fit, model_name, tidy_anova)) %>%
  select(model_name, anova) %>%
  unnest(anova)
write_csv(anova_results, file.path(output_dir, "02_model_selection_anova.csv"))

slopes <- valid_fits %>%
  mutate(slopes = map2(fit, model_name, extract_group_slopes)) %>%
  select(model_name, slopes) %>%
  unnest(slopes)
write_csv(slopes, file.path(output_dir, "03_group_auci_slopes_by_model.csv"))

slope_contrasts <- valid_fits %>%
  mutate(contrasts = map2(fit, model_name, extract_group_contrast_of_slopes)) %>%
  select(model_name, contrasts) %>%
  unnest(contrasts)
write_csv(slope_contrasts, file.path(output_dir, "04_group_auci_slope_contrasts_by_model.csv"))

interaction_summary <- anova_results %>%
  filter(term == "Group:AUCi_z") %>%
  left_join(
    comparison %>% select(model_name, fixed_name, random_name, AIC, BIC, delta_AIC, delta_BIC, AIC_rank, singular),
    by = "model_name"
  ) %>%
  arrange(AIC_rank)
write_csv(interaction_summary, file.path(output_dir, "05_group_auci_interaction_by_model.csv"))

best_non_singular_name <- comparison %>%
  filter(fit_ok, converged, !is.na(singular), !singular) %>%
  arrange(AIC) %>%
  slice_head(n = 1) %>%
  pull(model_name)

best_aic_name <- comparison %>%
  filter(fit_ok, converged) %>%
  arrange(AIC) %>%
  slice_head(n = 1) %>%
  pull(model_name)

write_text_output(
  c(
    "AUCi eye-duration model-selection grid",
    "",
    paste("Best converged AIC model:", best_aic_name),
    paste("Best converged non-singular AIC model:", ifelse(length(best_non_singular_name) == 0, "none", best_non_singular_name)),
    "",
    "Fixed-effect candidates vary medication/session adjustment and whether medication moderates Group x AUCi.",
    "Random-effect candidates: none, participant intercept, diagonal medication slopes by participant, diagonal session slopes by participant.",
    "",
    "Use 05_group_auci_interaction_by_model.csv to inspect whether Group x AUCi is stable across candidate structures."
  ),
  output_dir,
  "00_readme.txt"
)
