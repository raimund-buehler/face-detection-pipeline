# Manuscript section: Exploratory duration-based longitudinal analysis
# Analysis family: normalized interview-time overlapping sliding-window models for background AOI
# Original source path: adapted from the sliding-window eyes model
# Primary input dataset(s): 00_data/derived/preprocessing/duration_sliding_window_data.csv
# Primary output(s): omnibus and contrast outputs for normalized sliding-window background models
# Known TODOs: compare against collapsed face/background formulation if this becomes manuscript-relevant
# Scientific logic note: mirrors the exploratory eyes overlapping sliding-window model for the background region

library(tidyverse)
library(glmmTMB)
library(emmeans)
library(car)
library(here)
source(here("01_analysis", "shared", "model_output_utils.R"))

input_path <- here("00_data", "derived", "preprocessing", "duration_sliding_window_data.csv")

beta_squeeze <- function(x) {
  n <- sum(!is.na(x))
  ((x * (n - 1)) + 0.5) / n
}

df <- read_csv(input_path, show_col_types = FALSE) %>%
  filter(!is.na(window_fix_duration_prop), !is.na(medication)) %>%
  mutate(
    ID = Sub_ID,
    Group = factor(Group, levels = c("ASD", "CTRL")),
    fix = factor(
      fix,
      levels = c("fix_on_eyes", "fix_on_mouth", "fix_on_face", "fix_on_background"),
      labels = c("Eyes", "Mouth", "Face", "Background")
    ),
    medication = factor(
      medication,
      levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA"),
      labels = c("BOTH", "NAL", "OXT", "PLA")
    ),
    progress_window = factor(progress_window, levels = unique(progress_window[order(progress_window_index)]), ordered = TRUE)
  ) %>%
  filter(fix == "Background") %>%
  mutate(window_fix_duration_prop_beta = beta_squeeze(window_fix_duration_prop))

build_significance_intervals <- function(contrast_df, group_vars) {
  effect_col <- intersect(c("estimate", "odds.ratio", "ratio"), names(contrast_df))

  if (length(effect_col) == 0) {
    stop("No usable effect column found in contrast data.")
  }

  effect_col <- effect_col[[1]]

  contrast_df %>%
    filter(!is.na(p.value), p.value < 0.05) %>%
    arrange(across(all_of(c(group_vars, "progress_window_index")))) %>%
    mutate(
      direction = case_when(
        .data[[effect_col]] > 1 ~ "ASD_gt_CTRL",
        .data[[effect_col]] < 1 ~ "ASD_lt_CTRL",
        TRUE ~ "no_difference"
      )
    ) %>%
    filter(direction != "no_difference") %>%
    group_by(across(all_of(group_vars)), direction) %>%
    mutate(run_id = cumsum(c(TRUE, diff(progress_window_index) != 1))) %>%
    group_by(across(all_of(group_vars)), direction, run_id) %>%
    summarise(
      interval_start = min(progress_start, na.rm = TRUE),
      interval_end = max(progress_end, na.rm = TRUE),
      first_window = min(progress_window_index, na.rm = TRUE),
      last_window = max(progress_window_index, na.rm = TRUE),
      .groups = "drop"
    )
}

output_dir <- here("02_outputs", "model_outputs", "supplementary", "duration_sliding_window_background")
omnibus_dir <- file.path(output_dir, "01_omnibus")
group_contrast_dir <- file.path(output_dir, "02_group_contrasts")

dir.create(omnibus_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(group_contrast_dir, recursive = TRUE, showWarnings = FALSE)

main_model_background <- glmmTMB(
  window_fix_duration_prop_beta ~ Group * progress_window * medication + Sex + session + (1 | ID),
  family = beta_family(),
  data = df
)

main_anova_background <- Anova(main_model_background, type = 2) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term") %>%
  as_tibble()

write_anova_table_from_df(
  main_anova_background,
  omnibus_dir,
  "01_duration_sliding_window_background_beta_anova.csv",
  "duration_sliding_window_background_beta_model"
)

df_binary_sex <- df %>%
  filter(Sex %in% c("f", "m")) %>%
  droplevels()

sex_interaction_model_background <- glmmTMB(
  window_fix_duration_prop_beta ~ Group * progress_window * medication * Sex + session + (1 | ID),
  family = beta_family(),
  data = df_binary_sex
)

sex_interaction_anova_background <- Anova(sex_interaction_model_background, type = 2) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term") %>%
  as_tibble()

write_anova_table_from_df(
  sex_interaction_anova_background,
  omnibus_dir,
  "02_duration_sliding_window_background_beta_sex_interaction_anova.csv",
  "duration_sliding_window_background_beta_sex_interaction_model"
)

background_group_stats <- df %>%
  group_by(Group, medication, progress_window, progress_window_index, progress_midpoint, progress_start, progress_end) %>%
  summarise(
    mean_window_fix_duration_prop = mean(window_fix_duration_prop, na.rm = TRUE),
    sem = sd(window_fix_duration_prop, na.rm = TRUE) / sqrt(n()),
    lower_ci = mean_window_fix_duration_prop - qt(0.975, df = n() - 1) * sem,
    upper_ci = mean_window_fix_duration_prop + qt(0.975, df = n() - 1) * sem,
    .groups = "drop"
  )

background_contrasts <- emmeans(
  main_model_background,
  pairwise ~ Group | progress_window * medication,
  type = "response"
)

background_contrasts <- summary(background_contrasts$contrasts) %>%
  as_tibble() %>%
  mutate(
    progress_window = as.character(progress_window)
  )

readr::write_csv(
  background_contrasts,
  file.path(group_contrast_dir, "01_duration_sliding_window_background_group_by_window_medication_contrasts.csv")
)

background_intervals <- background_contrasts %>%
  left_join(
    background_group_stats %>%
      distinct(medication, progress_window, progress_window_index, progress_start, progress_end),
    by = c("medication", "progress_window")
  ) %>%
  build_significance_intervals(group_vars = "medication")

readr::write_csv(
  background_intervals,
  file.path(group_contrast_dir, "02_duration_sliding_window_background_group_by_window_medication_intervals.csv")
)

background_group_stats_by_sex <- df_binary_sex %>%
  group_by(Group, Sex, medication, progress_window, progress_window_index, progress_midpoint, progress_start, progress_end) %>%
  summarise(
    mean_window_fix_duration_prop = mean(window_fix_duration_prop, na.rm = TRUE),
    sem = sd(window_fix_duration_prop, na.rm = TRUE) / sqrt(n()),
    lower_ci = mean_window_fix_duration_prop - qt(0.975, df = n() - 1) * sem,
    upper_ci = mean_window_fix_duration_prop + qt(0.975, df = n() - 1) * sem,
    .groups = "drop"
  )

background_contrasts_by_sex <- emmeans(
  sex_interaction_model_background,
  pairwise ~ Group | progress_window * medication * Sex,
  type = "response"
)

background_contrasts_by_sex <- summary(background_contrasts_by_sex$contrasts) %>%
  as_tibble() %>%
  mutate(
    progress_window = as.character(progress_window)
  )

readr::write_csv(
  background_contrasts_by_sex,
  file.path(group_contrast_dir, "03_duration_sliding_window_background_group_by_window_medication_sex_contrasts.csv")
)

background_intervals_by_sex <- background_contrasts_by_sex %>%
  left_join(
    background_group_stats_by_sex %>%
      distinct(Sex, medication, progress_window, progress_window_index, progress_start, progress_end),
    by = c("Sex", "medication", "progress_window")
  ) %>%
  build_significance_intervals(group_vars = c("Sex", "medication"))

readr::write_csv(
  background_intervals_by_sex,
  file.path(group_contrast_dir, "04_duration_sliding_window_background_group_by_window_medication_sex_intervals.csv")
)
