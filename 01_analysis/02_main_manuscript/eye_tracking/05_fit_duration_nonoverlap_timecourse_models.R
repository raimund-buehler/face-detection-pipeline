# Manuscript section: Exploratory duration-based longitudinal analysis
# Analysis family: non-overlapping time-bin sensitivity check
# Primary input dataset(s): 00_data/derived/preprocessing/duration_sliding_window_data.csv
# Primary output(s): non-overlapping-bin omnibus tests, group contrasts, and attenuation contrasts

library(tidyverse)
library(glmmTMB)
library(emmeans)
library(car)
library(here)

beta_squeeze <- function(x) {
  n <- sum(!is.na(x))
  ((x * (n - 1)) + 0.5) / n
}

output_dir <- here(
  "02_outputs", "model_outputs", "main_manuscript", "eye_tracking", "timecourse_duration_nonoverlap"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- read_csv(
  here("00_data", "derived", "preprocessing", "duration_sliding_window_data.csv"),
  show_col_types = FALSE
) %>%
  filter(!is.na(window_fix_duration_prop), !is.na(medication), Sex %in% c("f", "m")) %>%
  # Select the 10%-wide windows that start at 0%, 10%, ..., 90%.
  # These are non-overlapping except for exact boundary timestamps.
  filter(abs((progress_start * 10) - round(progress_start * 10)) < 1e-8) %>%
  mutate(
    ID = Sub_ID,
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
    Group = factor(Group, levels = c("ASD", "CTRL")),
    session = factor(session),
    progress_bin = factor(
      sprintf("%02d", as.integer(round(progress_start * 10)) + 1),
      levels = sprintf("%02d", 1:10),
      ordered = TRUE
    ),
    progress_midpoint = progress_start + 0.05,
    window_fix_duration_prop_beta = beta_squeeze(window_fix_duration_prop)
  ) %>%
  filter(fix %in% c("Eyes", "Background"))

write_csv(
  df %>%
    distinct(progress_bin, progress_start, progress_midpoint, progress_end) %>%
    arrange(progress_bin),
  file.path(output_dir, "00_nonoverlap_bins.csv")
)

fit_block <- function(fix_label, panel_label) {
  df_block <- df %>%
    filter(
      fix == fix_label,
      if (panel_label == "Female") {
        Sex == "f"
      } else if (panel_label == "Male") {
        Sex == "m"
      } else {
        Sex %in% c("f", "m")
      }
    )

  model_formula <- if (panel_label == "Overall") {
    window_fix_duration_prop_beta ~ Group * progress_bin * medication + Sex + session + (1 | ID)
  } else {
    window_fix_duration_prop_beta ~ Group * progress_bin * medication + session + (1 | ID)
  }

  model <- glmmTMB(
    model_formula,
    family = beta_family(),
    data = df_block,
    control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
  )

  anova_table <- car::Anova(model, type = 3) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    as_tibble() %>%
    mutate(fix = fix_label, panel = panel_label)

  group_contrasts <- contrast(
    emmeans(model, ~ Group * medication),
    "revpairwise",
    by = "medication"
  )

  simple_contrasts <- summary(group_contrasts, infer = c(TRUE, TRUE), adjust = "none") %>%
    as_tibble() %>%
    mutate(fix = fix_label, panel = panel_label)

  attenuation_tests <- pairs(group_contrasts, by = NULL) %>%
    summary(infer = c(TRUE, TRUE), adjust = "none") %>%
    as_tibble() %>%
    mutate(fix = fix_label, panel = panel_label)

  list(
    anova = anova_table,
    simple_contrasts = simple_contrasts,
    attenuation_tests = attenuation_tests
  )
}

results <- list()
for (fix_label in c("Eyes", "Background")) {
  for (panel_label in c("Overall", "Female", "Male")) {
    results[[paste(fix_label, panel_label, sep = "_")]] <- fit_block(fix_label, panel_label)
  }
}

anova_results <- map_dfr(results, "anova")
simple_contrasts <- map_dfr(results, "simple_contrasts")
attenuation_tests <- map_dfr(results, "attenuation_tests")

write_csv(anova_results, file.path(output_dir, "01_nonoverlap_anova.csv"))
write_csv(simple_contrasts, file.path(output_dir, "02_nonoverlap_group_contrasts_by_medication.csv"))
write_csv(attenuation_tests, file.path(output_dir, "03_nonoverlap_attenuation_tests.csv"))

summary_lines <- c(
  "Non-overlapping time-bin sensitivity models",
  "",
  "Input:",
  "duration_sliding_window_data.csv, restricted to 10%-wide windows starting at 0%, 10%, ..., 90%",
  "",
  "Model:",
  "Overall: window_fix_duration_prop_beta ~ Group * progress_bin * medication + Sex + session + (1 | ID)",
  "Sex-stratified: window_fix_duration_prop_beta ~ Group * progress_bin * medication + session + (1 | ID)",
  "",
  "Purpose:",
  "Sensitivity check for medication attenuation tests without using heavily overlapping moving windows."
)

write_lines(summary_lines, file.path(output_dir, "00_model_notes.txt"))
