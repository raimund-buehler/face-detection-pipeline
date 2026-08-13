# Manuscript section: Exploratory cortisol sensitivity
# Analysis family: robustness of cortisol-reactivity to fixation-duration association
# Primary input dataset(s): 00_data/derived/analysis/df_cortisol_min_max.csv; 00_data/derived/preprocessing/duration_data.csv
# Primary output(s): robustness checks for T2T3/AUCi association with eye fixation-duration proportion

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

winsorize <- function(x, probs = c(0.05, 0.95)) {
  qs <- quantile(x, probs = probs, na.rm = TRUE)
  pmin(pmax(x, qs[[1]]), qs[[2]])
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

extract_group_slopes <- function(model, var_name, model_name) {
  broom::tidy(emtrends(model, ~ Group, var = var_name)) %>%
    mutate(
      model = model_name,
      predictor = var_name,
      p_display = format_p_value(p.value),
      sig = significance_stars(p.value),
      significant_0_05 = !is.na(p.value) & p.value < 0.05
    ) %>%
    relocate(model, predictor)
}

output_dir <- here(
  "02_outputs", "model_outputs", "main_manuscript", "cortisol",
  "duration_exploratory", "robustness"
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
    session = factor(session),
    T2T3 = MinMax
  ) %>%
  left_join(duration_eyes, by = c("ID", "session_norm")) %>%
  filter(
    !is.na(percentage_fix_duration),
    !is.na(T2T3),
    !is.na(AUCi),
    !is.na(Group),
    !is.na(medication)
  ) %>%
  mutate(
    T2T3_z = as.numeric(scale(T2T3)),
    AUCi_z = as.numeric(scale(AUCi)),
    T2T3_winsor = winsorize(T2T3),
    T2T3_winsor_z = as.numeric(scale(T2T3_winsor)),
    AUCi_winsor = winsorize(AUCi),
    AUCi_winsor_z = as.numeric(scale(AUCi_winsor)),
    T2T3_rank_z = as.numeric(scale(rank(T2T3, ties.method = "average"))),
    AUCi_rank_z = as.numeric(scale(rank(AUCi, ties.method = "average"))),
    T2T3_trim_95 = T2T3 >= quantile(T2T3, 0.025, na.rm = TRUE) &
      T2T3 <= quantile(T2T3, 0.975, na.rm = TRUE),
    T2T3_trim_90 = T2T3 >= quantile(T2T3, 0.05, na.rm = TRUE) &
      T2T3 <= quantile(T2T3, 0.95, na.rm = TRUE),
    AUCi_trim_95 = AUCi >= quantile(AUCi, 0.025, na.rm = TRUE) &
      AUCi <= quantile(AUCi, 0.975, na.rm = TRUE),
    AUCi_trim_90 = AUCi >= quantile(AUCi, 0.05, na.rm = TRUE) &
      AUCi <= quantile(AUCi, 0.95, na.rm = TRUE),
    percentage_fix_duration_beta = beta_squeeze(percentage_fix_duration)
  )

write_csv(
  df %>%
    summarise(
      n_observations = n(),
      n_participants = n_distinct(ID),
      T2T3_min = min(T2T3),
      T2T3_q025 = quantile(T2T3, 0.025),
      T2T3_q05 = quantile(T2T3, 0.05),
      T2T3_q95 = quantile(T2T3, 0.95),
      T2T3_q975 = quantile(T2T3, 0.975),
      T2T3_max = max(T2T3),
      AUCi_min = min(AUCi),
      AUCi_q025 = quantile(AUCi, 0.025),
      AUCi_q05 = quantile(AUCi, 0.05),
      AUCi_q95 = quantile(AUCi, 0.95),
      AUCi_q975 = quantile(AUCi, 0.975),
      AUCi_max = max(AUCi)
    ),
  file.path(output_dir, "00_predictor_distribution_summary.csv")
)

write_csv(
  df %>%
    mutate(abs_T2T3_z = abs(T2T3_z), abs_AUCi_z = abs(AUCi_z)) %>%
    arrange(desc(abs_T2T3_z)) %>%
    select(ID, Group, session, medication, T2T3, T2T3_z, AUCi, AUCi_z, percentage_fix_duration) %>%
    slice_head(n = 20),
  file.path(output_dir, "01_top_t2t3_extreme_sessions.csv")
)

write_csv(
  df %>%
    mutate(abs_T2T3_z = abs(T2T3_z), abs_AUCi_z = abs(AUCi_z)) %>%
    arrange(desc(abs_AUCi_z)) %>%
    select(ID, Group, session, medication, T2T3, T2T3_z, AUCi, AUCi_z, percentage_fix_duration) %>%
    slice_head(n = 20),
  file.path(output_dir, "02_top_auci_extreme_sessions.csv")
)

fit_beta_model <- function(data, predictor, model_name) {
  model <- glmmTMB(
    as.formula(paste0(
      "percentage_fix_duration_beta ~ Group * ", predictor,
      " + medication + session + (1 | ID)"
    )),
    family = beta_family(),
    data = data,
    control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
  )

  list(
    model_name = model_name,
    predictor = predictor,
    n_observations = nobs(model),
    n_participants = n_distinct(data$ID),
    model = model,
    anova = tidy_anova(model, model_name),
    slopes = extract_group_slopes(model, predictor, model_name)
  )
}

model_specs <- tribble(
  ~model_name, ~predictor, ~filter_expr,
  "t2t3_raw_full", "T2T3_z", "TRUE",
  "t2t3_winsor_5_95", "T2T3_winsor_z", "TRUE",
  "t2t3_rank", "T2T3_rank_z", "TRUE",
  "t2t3_trim_2_5_97_5", "T2T3_z", "T2T3_trim_95",
  "t2t3_trim_5_95", "T2T3_z", "T2T3_trim_90",
  "auci_raw_full", "AUCi_z", "TRUE",
  "auci_winsor_5_95", "AUCi_winsor_z", "TRUE",
  "auci_rank", "AUCi_rank_z", "TRUE",
  "auci_trim_2_5_97_5", "AUCi_z", "AUCi_trim_95",
  "auci_trim_5_95", "AUCi_z", "AUCi_trim_90"
)

results <- pmap(model_specs, function(model_name, predictor, filter_expr) {
  data_i <- df %>% filter(eval(parse(text = filter_expr)))
  fit_beta_model(data_i, predictor, model_name)
})

write_csv(
  map_dfr(results, ~ tibble(
    model = .x$model_name,
    predictor = .x$predictor,
    n_observations = .x$n_observations,
    n_participants = .x$n_participants,
    AIC = AIC(.x$model),
    BIC = BIC(.x$model)
  )),
  file.path(output_dir, "03_model_fit_summary.csv")
)

write_csv(
  map_dfr(results, "anova"),
  file.path(output_dir, "04_robustness_model_anova.csv")
)

write_csv(
  map_dfr(results, "slopes"),
  file.path(output_dir, "05_robustness_group_slopes.csv")
)

write_csv(
  map_dfr(results, function(x) {
    x$anova %>%
      filter(term == paste0("Group:", x$predictor)) %>%
      transmute(
        model,
        predictor = x$predictor,
        n_observations = x$n_observations,
        n_participants = x$n_participants,
        interaction_chisq = Chisq,
        interaction_df = Df,
        interaction_p = `Pr(>Chisq)`,
        interaction_p_display = p_display,
        interaction_sig = sig
      )
  }),
  file.path(output_dir, "06_interaction_robustness_summary.csv")
)

write_text_output(
  c(
    "Cortisol-reactivity to eye fixation-duration robustness checks",
    "",
    "Outcome: percentage_fix_duration_beta for fix_on_eyes.",
    "Base model: outcome ~ Group * predictor + medication + session + (1 | ID).",
    "",
    "Robustness variants:",
    "raw full sample",
    "winsorized 5th/95th percentile predictor",
    "rank-transformed predictor",
    "trimmed 2.5th/97.5th percentile predictor",
    "trimmed 5th/95th percentile predictor",
    "",
    "See 06_interaction_robustness_summary.csv for the Group x cortisol robustness summary."
  ),
  output_dir,
  "00_readme.txt"
)
