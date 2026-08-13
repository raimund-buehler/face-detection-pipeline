# Manuscript section: Exploratory sliding-window duration timecourse
# Analysis family: Eyes / Background timecourse models aligned to the preferred duration workflow
# Primary input dataset(s): 00_data/derived/preprocessing/duration_sliding_window_data.csv
# Primary output(s): saved model objects, omnibus tables, and reusable post-hoc results

library(tidyverse)
library(glmmTMB)
library(car)
library(emmeans)
library(here)

beta_squeeze <- function(x) {
  n <- sum(!is.na(x))
  ((x * (n - 1)) + 0.5) / n
}

write_anova_csv <- function(model, path) {
  car::Anova(model, type = 3) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    as_tibble() %>%
    write_csv(path)
}

save_model_bundle <- function(model, output_stub) {
  saveRDS(model, paste0(output_stub, ".rds"))
  capture.output(summary(model), file = paste0(output_stub, "_summary.txt"))
}

output_dir <- here(
  "02_outputs", "model_outputs", "main_manuscript", "eye_tracking", "timecourse_duration"
)
model_dir <- file.path(output_dir, "model_objects")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

df <- read_csv(
  here("00_data", "derived", "preprocessing", "duration_sliding_window_data.csv"),
  show_col_types = FALSE
) %>%
  filter(!is.na(window_fix_duration_prop), !is.na(medication), Sex %in% c("f", "m")) %>%
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
    Sex = factor(Sex),
    session = factor(session),
    progress_window = factor(
      progress_window,
      levels = unique(progress_window[order(progress_window_index)]),
      ordered = TRUE
    ),
    window_fix_duration_prop_beta = beta_squeeze(window_fix_duration_prop)
  ) %>%
  filter(fix %in% c("Eyes", "Background"))

fit_top_model <- function(fix_label) {
  df_fix <- filter(df, fix == fix_label)

  model <- glmmTMB(
    window_fix_duration_prop_beta ~ Group * progress_window + Group * Sex + session + (1 | ID),
    family = beta_family(),
    data = df_fix
  )

  stub <- file.path(model_dir, paste0("top_", str_to_lower(fix_label), "_group_progress_group_sex"))
  save_model_bundle(model, stub)
  write_anova_csv(model, file.path(output_dir, paste0("01_top_", str_to_lower(fix_label), "_anova.csv")))

  overall_contrasts <- summary(
    emmeans(model, pairwise ~ Group | progress_window, type = "response")$contrasts,
    infer = c(TRUE, TRUE)
  ) %>%
    as_tibble()

  female_contrasts <- summary(
    emmeans(model, pairwise ~ Group | progress_window * Sex, type = "response")$contrasts,
    infer = c(TRUE, TRUE)
  ) %>%
    as_tibble() %>%
    filter(Sex == "f")

  write_csv(
    overall_contrasts,
    file.path(output_dir, paste0("02_top_", str_to_lower(fix_label), "_group_contrasts_overall.csv"))
  )
  write_csv(
    female_contrasts,
    file.path(output_dir, paste0("03_top_", str_to_lower(fix_label), "_group_contrasts_female.csv"))
  )
}

fit_medication_model <- function(fix_label) {
  df_fix <- filter(df, fix == fix_label)

  model <- glmmTMB(
    window_fix_duration_prop_beta ~ Group * progress_window * medication + Group * Sex + session + (1 | ID),
    family = beta_family(),
    data = df_fix
  )

  stub <- file.path(model_dir, paste0("medication_", str_to_lower(fix_label), "_group_progress_medication_group_sex"))
  save_model_bundle(model, stub)
  write_anova_csv(model, file.path(output_dir, paste0("04_medication_", str_to_lower(fix_label), "_anova.csv")))

  overall_group_contrasts <- summary(
    emmeans(model, pairwise ~ Group | medication, type = "response")$contrasts,
    infer = c(TRUE, TRUE)
  ) %>%
    as_tibble()

  female_group_contrasts <- summary(
    emmeans(model, pairwise ~ Group | medication * Sex, type = "response")$contrasts,
    infer = c(TRUE, TRUE)
  ) %>%
    as_tibble() %>%
    filter(Sex == "f")

  female_emm <- emmeans(model, ~ Group * medication, by = "Sex", at = list(Sex = "f"))
  female_emm <- female_emm[female_emm@grid$Sex == "f", ]

  overall_attenuation <- pairs(
    contrast(emmeans(model, ~ Group * medication), "revpairwise", by = "medication"),
    by = NULL
  ) %>%
    summary(infer = c(TRUE, TRUE), adjust = "none") %>%
    as_tibble()

  female_attenuation <- pairs(
    contrast(female_emm, "revpairwise", by = "medication"),
    by = NULL
  ) %>%
    summary(infer = c(TRUE, TRUE), adjust = "none") %>%
    as_tibble()

  write_csv(
    overall_group_contrasts,
    file.path(output_dir, paste0("05_medication_", str_to_lower(fix_label), "_group_contrasts_overall.csv"))
  )
  write_csv(
    female_group_contrasts,
    file.path(output_dir, paste0("06_medication_", str_to_lower(fix_label), "_group_contrasts_female.csv"))
  )
  write_csv(
    overall_attenuation,
    file.path(output_dir, paste0("07_medication_", str_to_lower(fix_label), "_attenuation_tests_overall.csv"))
  )
  write_csv(
    female_attenuation,
    file.path(output_dir, paste0("08_medication_", str_to_lower(fix_label), "_attenuation_tests_female.csv"))
  )
}

fit_top_model("Eyes")
fit_top_model("Background")
fit_medication_model("Eyes")
fit_medication_model("Background")

notes <- c(
  "Exploratory sliding-window duration timecourse models",
  "",
  "Top models:",
  "window_fix_duration_prop_beta ~ Group * progress_window + Group * Sex + session + (1 | ID)",
  "",
  "Medication models:",
  "window_fix_duration_prop_beta ~ Group * progress_window * medication + Group * Sex + session + (1 | ID)",
  "",
  "Note:",
  "Because the medication models include Group * Sex but not Group * medication * Sex, female medication attenuation tests are algebraically identical to the overall attenuation tests."
)

write_lines(notes, file.path(output_dir, "00_model_notes.txt"))
