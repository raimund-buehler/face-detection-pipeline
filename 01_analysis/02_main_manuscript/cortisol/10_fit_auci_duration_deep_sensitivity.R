# Manuscript section: Exploratory cortisol sensitivity
# Analysis family: deep AUCi robustness for eye fixation-duration association
# Primary input dataset(s): 00_data/derived/analysis/df_cortisol_min_max.csv; 00_data/derived/preprocessing/duration_data.csv
# Primary output(s): leave-one-out, robust, and permutation checks for AUCi x Group

library(tidyverse)
library(glmmTMB)
library(car)
library(emmeans)
library(MASS)
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

get_interaction_row <- function(model, model_name, predictor = "AUCi_z") {
  tidy_anova(model, model_name) %>%
    filter(term == paste0("Group:", predictor)) %>%
    transmute(
      model,
      interaction_chisq = Chisq,
      interaction_df = Df,
      interaction_p = `Pr(>Chisq)`,
      interaction_p_display = p_display,
      interaction_sig = sig
    )
}

fit_beta <- function(data, predictor = "AUCi_z") {
  glmmTMB(
    as.formula(paste0(
      "percentage_fix_duration_beta ~ Group * ", predictor,
      " + medication + session + (1 | ID)"
    )),
    family = beta_family(),
    data = data,
    control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
  )
}

extract_slopes <- function(model, model_name, predictor = "AUCi_z") {
  broom::tidy(emtrends(model, ~ Group, var = predictor)) %>%
    mutate(
      model = model_name,
      predictor = predictor,
      p_display = format_p_value(p.value),
      sig = significance_stars(p.value)
    ) %>%
    relocate(model, predictor)
}

output_dir <- here(
  "02_outputs", "model_outputs", "main_manuscript", "cortisol",
  "duration_exploratory", "robustness", "auci_deep"
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
    AUCi_rank_z = as.numeric(scale(rank(AUCi, ties.method = "average"))),
    percentage_fix_duration_beta = beta_squeeze(percentage_fix_duration),
    eye_logit = qlogis(percentage_fix_duration_beta)
  )

base_model <- fit_beta(df, "AUCi_z")
base_summary <- get_interaction_row(base_model, "auci_raw_full", "AUCi_z")
base_slopes <- extract_slopes(base_model, "auci_raw_full", "AUCi_z")

write_csv(base_summary, file.path(output_dir, "00_base_interaction.csv"))
write_csv(base_slopes, file.path(output_dir, "01_base_group_slopes.csv"))

write_csv(
  df %>%
    mutate(abs_AUCi_z = abs(AUCi_z)) %>%
    arrange(desc(abs_AUCi_z)) %>%
    dplyr::select(ID, Group, session, medication, AUCi, AUCi_z, percentage_fix_duration) %>%
    slice_head(n = 40),
  file.path(output_dir, "02_top_auci_extreme_sessions.csv")
)

session_loso <- map_dfr(seq_len(nrow(df)), function(i) {
  data_i <- df[-i, ]
  model_i <- fit_beta(data_i, "AUCi_z")
  row <- get_interaction_row(model_i, paste0("drop_session_", i), "AUCi_z")
  row %>%
    mutate(
      dropped_row = i,
      dropped_ID = as.character(df$ID[[i]]),
      dropped_Group = as.character(df$Group[[i]]),
      dropped_session = as.character(df$session[[i]]),
      dropped_medication = as.character(df$medication[[i]]),
      dropped_AUCi = df$AUCi[[i]],
      dropped_AUCi_z = df$AUCi_z[[i]],
      dropped_eye_duration = df$percentage_fix_duration[[i]]
    )
})
write_csv(session_loso, file.path(output_dir, "03_leave_one_session_out_interactions.csv"))

subject_loso <- map_dfr(sort(unique(as.character(df$ID))), function(id_i) {
  data_i <- df %>% filter(ID != id_i)
  model_i <- fit_beta(data_i, "AUCi_z")
  row <- get_interaction_row(model_i, paste0("drop_subject_", id_i), "AUCi_z")
  subject_info <- df %>%
    filter(ID == id_i) %>%
    summarise(
      dropped_ID = first(as.character(ID)),
      dropped_Group = first(as.character(Group)),
      n_sessions = n(),
      mean_AUCi = mean(AUCi, na.rm = TRUE),
      max_abs_AUCi_z = max(abs(AUCi_z), na.rm = TRUE),
      mean_eye_duration = mean(percentage_fix_duration, na.rm = TRUE),
      .groups = "drop"
    )
  bind_cols(row, subject_info)
})
write_csv(subject_loso, file.path(output_dir, "04_leave_one_subject_out_interactions.csv"))

lm_model <- lm(
  eye_logit ~ Group * AUCi_z + medication + session,
  data = df
)
lm_cooks <- cooks.distance(lm_model)
lm_hat <- hatvalues(lm_model)
lm_aug <- df %>%
  mutate(
    cooks_distance = lm_cooks,
    leverage = lm_hat,
    abs_AUCi_z = abs(AUCi_z)
  ) %>%
  arrange(desc(cooks_distance)) %>%
  dplyr::select(ID, Group, session, medication, AUCi, AUCi_z, percentage_fix_duration, cooks_distance, leverage)
write_csv(lm_aug, file.path(output_dir, "05_linear_model_influence_sessions.csv"))

cook_cutoff <- 4 / nrow(df)
high_influence_df <- df[lm_cooks <= cook_cutoff, ]
lm_no_high_cook <- lm(
  eye_logit ~ Group * AUCi_z + medication + session,
  data = high_influence_df
)
robust_lm <- MASS::rlm(
  eye_logit ~ Group * AUCi_z + medication + session,
  data = df,
  maxit = 100
)
rank_lm <- lm(
  eye_logit ~ Group * AUCi_rank_z + medication + session,
  data = df
)

lm_summaries <- bind_rows(
  broom::tidy(lm_model) %>% mutate(model = "linear_full"),
  broom::tidy(lm_no_high_cook) %>% mutate(model = "linear_drop_cook_gt_4_over_n"),
  broom::tidy(robust_lm) %>% mutate(model = "robust_rlm"),
  broom::tidy(rank_lm) %>% mutate(model = "linear_rank_auci")
) %>%
  mutate(
    p_display = if_else(!is.na(p.value), format_p_value(p.value), NA_character_),
    sig = if_else(!is.na(p.value), significance_stars(p.value), NA_character_)
  ) %>%
  relocate(model)
write_csv(lm_summaries, file.path(output_dir, "06_linear_robust_model_coefficients.csv"))

set.seed(20260421)
n_perm <- 1000
observed_chisq <- base_summary$interaction_chisq[[1]]
perm_results <- map_dfr(seq_len(n_perm), function(i) {
  perm_df <- df %>%
    mutate(Group = sample(Group))
  model_i <- fit_beta(perm_df, "AUCi_z")
  row_i <- get_interaction_row(model_i, paste0("perm_", i), "AUCi_z")
  tibble(
    permutation = i,
    interaction_chisq = row_i$interaction_chisq[[1]],
    interaction_p = row_i$interaction_p[[1]]
  )
})
perm_p <- mean(perm_results$interaction_chisq >= observed_chisq, na.rm = TRUE)
write_csv(perm_results, file.path(output_dir, "07_permutation_interactions.csv"))
write_csv(
  tibble(
    observed_chisq = observed_chisq,
    n_permutations = n_perm,
    empirical_p = perm_p
  ),
  file.path(output_dir, "08_permutation_summary.csv")
)

write_csv(
  tibble(
    check = c(
      "base_glmmTMB_interaction_p",
      "leave_one_session_min_p",
      "leave_one_session_max_p",
      "leave_one_subject_min_p",
      "leave_one_subject_max_p",
      "n_high_cook_sessions",
      "cook_cutoff_4_over_n",
      "permutation_empirical_p"
    ),
    value = c(
      base_summary$interaction_p[[1]],
      min(session_loso$interaction_p, na.rm = TRUE),
      max(session_loso$interaction_p, na.rm = TRUE),
      min(subject_loso$interaction_p, na.rm = TRUE),
      max(subject_loso$interaction_p, na.rm = TRUE),
      sum(lm_cooks > cook_cutoff, na.rm = TRUE),
      cook_cutoff,
      perm_p
    )
  ),
  file.path(output_dir, "09_deep_sensitivity_summary.csv")
)

write_text_output(
  c(
    "AUCi deep sensitivity for eye fixation-duration association",
    "",
    "Base model: percentage_fix_duration_beta ~ Group * AUCi_z + medication + session + (1 | ID).",
    "Additional checks:",
    "leave-one-session-out GLMMs",
    "leave-one-subject-out GLMMs",
    "Cook's distance from comparable logit-linear model",
    "linear model after removing Cook's D > 4/n",
    "robust linear model using MASS::rlm",
    "rank-transformed AUCi model",
    "1000 permutations of Group labels for the Group x AUCi interaction"
  ),
  output_dir,
  "00_readme.txt"
)
