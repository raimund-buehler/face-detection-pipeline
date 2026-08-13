# Manuscript section: Exploratory cortisol analyses
# Analysis family: raw cortisol and fixation-duration association
# Primary input dataset(s): 00_data/derived/analysis/df_cortisol_merged.csv; 00_data/derived/analysis/df_cortisol_min_max.csv; 00_data/derived/preprocessing/duration_data.csv
# Primary output(s): current Figure 5 cortisol responder figure and full-sample exploratory responder figure using fixation-duration proportion

library(tidyverse)
library(glmmTMB)
library(car)
library(emmeans)
library(ggh4x)
library(patchwork)
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

output_dir <- here("02_outputs", "figures", "main", "current_story")
supp_output_dir <- here("02_outputs", "figures", "supplementary", "cortisol")
model_output_dir <- here(
  "02_outputs", "model_outputs", "main_manuscript", "cortisol",
  "duration_exploratory"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_output_dir, recursive = TRUE, showWarnings = FALSE)

output_png <- file.path(output_dir, "Figure_5_cortisol_duration.png")
output_t2t3_responder_full_png <- file.path(supp_output_dir, "Figure_5_cortisol_duration_t2t3_responder_full_sample_exploratory.png")
output_t2t3_responder_full_svg <- file.path(supp_output_dir, "figure5_cortisol_duration_t2t3_responder_full_sample_exploratory.svg")

medication_values <- c(
  "BOTH" = "#FDC010",
  "NAL" = "#C7B655",
  "OXT" = "#93AD7C",
  "PLA" = "#0B7A6B"
)
group_values <- c("ASD" = "#D55E00", "CTRL" = "#0072B2")
medication_levels_raw <- c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA")
medication_levels_plot <- c("PLA", "OXT", "NAL", "BOTH")
auci_response_threshold <- 1
t2t3_response_threshold <- 0.5

raw_cortisol <- read_csv(
  here("00_data", "derived", "analysis", "df_cortisol_merged.csv"),
  show_col_types = FALSE
) %>%
  filter(fix == "Eyes", !is.na(cortisol_1), !is.na(medication)) %>%
  mutate(
    medication = factor(
      medication,
      levels = medication_levels_raw,
      labels = c("BOTH", "NAL", "OXT", "PLA")
    ),
    medication = factor(medication, levels = medication_levels_plot),
    timepoint = factor(timepoint, levels = sort(unique(timepoint))),
    ID = factor(ID),
    Group = factor(Group, levels = c("ASD", "CTRL"))
  )

raw_summary <- raw_cortisol %>%
  group_by(medication, timepoint) %>%
  summarise(
    mean_cortisol = mean(cortisol_1, na.rm = TRUE),
    sem = sd(cortisol_1, na.rm = TRUE) / sqrt(sum(!is.na(cortisol_1))),
    n = sum(!is.na(cortisol_1)),
    .groups = "drop"
  )

participant_metadata <- read_csv(
  here("00_data", "derived", "analysis", "df_analysis_pub.csv"),
  show_col_types = FALSE
) %>%
  transmute(ID, Sex) %>%
  distinct() %>%
  filter(Sex %in% c("f", "m"))

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

duration_background <- read_csv(
  here("00_data", "derived", "preprocessing", "duration_data.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    ID = Sub_ID,
    session_norm = normalize_session(session),
    fix,
    percentage_fix_duration
  ) %>%
  filter(fix == "fix_on_background")

cortisol_duration <- read_csv(
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
    NAL_status = factor(
      if_else(medication %in% c("BOTH", "NAL"), "NAL-containing", "Non-NAL"),
      levels = c("Non-NAL", "NAL-containing")
    ),
    Group = factor(Group, levels = c("ASD", "CTRL")),
    ID = factor(ID),
    session = factor(session),
    T2T3 = MinMax,
    T2T3_z = as.numeric(scale(T2T3)),
    AUCi_z = as.numeric(scale(AUCi))
  ) %>%
  left_join(participant_metadata, by = "ID") %>%
  left_join(duration_eyes, by = c("ID", "session_norm")) %>%
  filter(
    !is.na(percentage_fix_duration),
    !is.na(T2T3_z),
    !is.na(AUCi_z),
    !is.na(Group),
    !is.na(medication),
    !is.na(Sex)
  ) %>%
  mutate(
    Sex = factor(Sex, levels = c("f", "m"), labels = c("Female", "Male")),
    percentage_fix_duration_beta = beta_squeeze(percentage_fix_duration),
    AUCi_response = factor(
      case_when(
        AUCi < -auci_response_threshold ~ "Decrease",
        AUCi > auci_response_threshold ~ "Increase",
        TRUE ~ "Non-response"
      ),
      levels = c("Decrease", "Non-response", "Increase")
    ),
    T2T3_response = factor(
      case_when(
        T2T3 < -t2t3_response_threshold ~ "Decrease",
        T2T3 > t2t3_response_threshold ~ "Increase",
        TRUE ~ "Non-response"
      ),
      levels = c("Decrease", "Non-response", "Increase")
    )
  )

cortisol_background_duration <- read_csv(
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
    T2T3 = MinMax,
    T2T3_z = as.numeric(scale(T2T3)),
    AUCi_z = as.numeric(scale(AUCi))
  ) %>%
  left_join(participant_metadata, by = "ID") %>%
  left_join(duration_background, by = c("ID", "session_norm")) %>%
  filter(
    !is.na(percentage_fix_duration),
    !is.na(T2T3_z),
    !is.na(AUCi_z),
    !is.na(Group),
    !is.na(medication),
    !is.na(Sex)
  ) %>%
  mutate(
    Sex = factor(Sex, levels = c("f", "m"), labels = c("Female", "Male")),
    percentage_fix_duration_beta = beta_squeeze(percentage_fix_duration)
  )

write_csv(
  cortisol_duration %>%
    summarise(
      n_observations = n(),
      n_participants = n_distinct(ID),
      t2t3_min = min(T2T3, na.rm = TRUE),
      t2t3_max = max(T2T3, na.rm = TRUE),
      eye_duration_min = min(percentage_fix_duration, na.rm = TRUE),
      eye_duration_max = max(percentage_fix_duration, na.rm = TRUE)
    ),
  file.path(model_output_dir, "00_duration_cortisol_sample_summary.csv")
)

t2t3_summary_model <- lmerTest::lmer(T2T3 ~ medication + Group + (1 | ID), data = cortisol_duration)
t2t3_duration_model <- glmmTMB(
  percentage_fix_duration_beta ~ Group * T2T3_z + medication + session + (1 | ID),
  family = beta_family(),
  data = cortisol_duration,
  control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
)
auci_duration_model <- glmmTMB(
  percentage_fix_duration_beta ~ Group * AUCi_z + medication + session + (1 | ID),
  family = beta_family(),
  data = cortisol_duration,
  control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
)
auci_trim_limits <- quantile(cortisol_duration$AUCi, probs = c(0.05, 0.95), na.rm = TRUE)
cortisol_duration_trim90 <- cortisol_duration %>%
  filter(AUCi >= auci_trim_limits[[1]], AUCi <= auci_trim_limits[[2]]) %>%
  mutate(
    AUCi_trim_z = as.numeric(scale(AUCi)),
    auci_bin = ntile(AUCi, 5)
  )
auci_duration_trim90_model <- glmmTMB(
  percentage_fix_duration_beta ~ Group * AUCi_trim_z + medication + session + (1 | ID),
  family = beta_family(),
  data = cortisol_duration_trim90,
  control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
)
auci_responder_model <- glmmTMB(
  percentage_fix_duration_beta ~ Group * AUCi_response * NAL_status + session + (1 | ID),
  family = beta_family(),
  data = cortisol_duration,
  control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
)
t2t3_responder_model <- glmmTMB(
  percentage_fix_duration_beta ~ Group * T2T3_response * NAL_status + session + (1 | ID),
  family = beta_family(),
  data = cortisol_duration,
  control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
)
t2t3_responder_full_model <- glmmTMB(
  percentage_fix_duration_beta ~ Group * T2T3_response + medication + session + (1 | ID),
  family = beta_family(),
  data = cortisol_duration,
  control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
)
auci_background_duration_model <- glmmTMB(
  percentage_fix_duration_beta ~ Group * AUCi_z + medication + session + (1 | ID),
  family = beta_family(),
  data = cortisol_background_duration,
  control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
)

write_model_anova(
  t2t3_summary_model,
  model_output_dir,
  "01_t2t3_reactivity_model_anova.csv",
  "t2t3_reactivity_medication_group"
)
write_csv(
  tidy_anova(t2t3_duration_model, "t2t3_eye_duration_model"),
  file.path(model_output_dir, "02_t2t3_eye_duration_model_anova.csv")
)
write_csv(
  tidy_anova(auci_duration_model, "auci_eye_duration_model"),
  file.path(model_output_dir, "03_auci_eye_duration_model_anova.csv")
)
write_csv(
  tidy_anova(auci_background_duration_model, "auci_background_duration_model"),
  file.path(model_output_dir, "06_auci_background_duration_model_anova.csv")
)
write_csv(
  tidy_anova(auci_duration_trim90_model, "auci_eye_duration_trim90_model"),
  file.path(model_output_dir, "08_auci_eye_duration_trim90_model_anova.csv")
)
write_csv(
  tidy_anova(auci_responder_model, "auci_eye_duration_responder_model"),
  file.path(model_output_dir, "10_auci_eye_duration_responder_model_anova.csv")
)
write_csv(
  tidy_anova(t2t3_responder_model, "t2t3_eye_duration_responder_model"),
  file.path(model_output_dir, "15_t2t3_eye_duration_responder_model_anova.csv")
)
write_csv(
  tidy_anova(t2t3_responder_full_model, "t2t3_eye_duration_responder_full_model"),
  file.path(model_output_dir, "22_t2t3_eye_duration_responder_full_model_anova.csv")
)
write_csv(
  broom::tidy(emtrends(t2t3_duration_model, ~ Group, var = "T2T3_z")) %>%
    mutate(
      model = "t2t3_eye_duration_model",
      p_display = format_p_value(p.value),
      sig = significance_stars(p.value)
    ) %>%
    relocate(model),
  file.path(model_output_dir, "04_t2t3_eye_duration_group_slopes.csv")
)
write_csv(
  broom::tidy(emtrends(auci_duration_model, ~ Group, var = "AUCi_z")) %>%
    mutate(
      model = "auci_eye_duration_model",
      p_display = format_p_value(p.value),
      sig = significance_stars(p.value)
    ) %>%
    relocate(model),
  file.path(model_output_dir, "05_auci_eye_duration_group_slopes.csv")
)
write_csv(
  broom::tidy(emtrends(auci_background_duration_model, ~ Group, var = "AUCi_z")) %>%
    mutate(
      model = "auci_background_duration_model",
      p_display = format_p_value(p.value),
      sig = significance_stars(p.value)
    ) %>%
    relocate(model),
  file.path(model_output_dir, "07_auci_background_duration_group_slopes.csv")
)
write_csv(
  broom::tidy(emtrends(auci_duration_trim90_model, ~ Group, var = "AUCi_trim_z")) %>%
    mutate(
      model = "auci_eye_duration_trim90_model",
      p_display = format_p_value(p.value),
      sig = significance_stars(p.value)
    ) %>%
    relocate(model),
  file.path(model_output_dir, "09_auci_eye_duration_trim90_group_slopes.csv")
)
write_csv(
  cortisol_duration %>%
    count(NAL_status, AUCi_response, Group, name = "n_sessions") %>%
    arrange(NAL_status, AUCi_response, Group),
  file.path(model_output_dir, "11_auci_response_counts_by_nal_group.csv")
)
write_csv(
  cortisol_duration %>%
    group_by(NAL_status, AUCi_response, Group) %>%
    summarise(
      n_sessions = n(),
      n_participants = n_distinct(ID),
      mean_eye_duration = mean(percentage_fix_duration, na.rm = TRUE),
      sd_eye_duration = sd(percentage_fix_duration, na.rm = TRUE),
      sem_eye_duration = sd_eye_duration / sqrt(n_sessions),
      mean_AUCi = mean(AUCi, na.rm = TRUE),
      .groups = "drop"
    ),
  file.path(model_output_dir, "12_auci_response_eye_duration_descriptives.csv")
)
write_csv(
  broom::tidy(emmeans(auci_responder_model, ~ Group | NAL_status * AUCi_response, type = "response")) %>%
    mutate(
      model = "auci_eye_duration_responder_model",
      p_display = NA_character_,
      sig = NA_character_
    ) %>%
    relocate(model),
  file.path(model_output_dir, "13_auci_response_group_emmeans.csv")
)
write_csv(
  broom::tidy(contrast(
    emmeans(auci_responder_model, ~ Group | NAL_status * AUCi_response, type = "response"),
    method = "revpairwise"
  )) %>%
    mutate(
      model = "auci_eye_duration_responder_model",
      p_display = format_p_value(p.value),
      sig = significance_stars(p.value)
    ) %>%
    relocate(model),
  file.path(model_output_dir, "14_auci_response_group_contrasts.csv")
)
write_csv(
  cortisol_duration %>%
    count(NAL_status, T2T3_response, Group, name = "n_sessions") %>%
    arrange(NAL_status, T2T3_response, Group),
  file.path(model_output_dir, "16_t2t3_response_counts_by_nal_group.csv")
)
write_csv(
  cortisol_duration %>%
    count(Sex, NAL_status, T2T3_response, Group, name = "n_sessions") %>%
    arrange(Sex, NAL_status, T2T3_response, Group),
  file.path(model_output_dir, "20_t2t3_response_counts_by_sex_nal_group.csv")
)
write_csv(
  cortisol_duration %>%
    count(T2T3_response, Group, name = "n_sessions") %>%
    arrange(T2T3_response, Group),
  file.path(model_output_dir, "23_t2t3_response_counts_full_sample_group.csv")
)
write_csv(
  cortisol_duration %>%
    group_by(NAL_status, T2T3_response, Group) %>%
    summarise(
      n_sessions = n(),
      n_participants = n_distinct(ID),
      mean_eye_duration = mean(percentage_fix_duration, na.rm = TRUE),
      sd_eye_duration = sd(percentage_fix_duration, na.rm = TRUE),
      sem_eye_duration = sd_eye_duration / sqrt(n_sessions),
      mean_T2T3 = mean(T2T3, na.rm = TRUE),
      .groups = "drop"
    ),
  file.path(model_output_dir, "17_t2t3_response_eye_duration_descriptives.csv")
)
write_csv(
  cortisol_duration %>%
    group_by(Sex, NAL_status, T2T3_response, Group) %>%
    summarise(
      n_sessions = n(),
      n_participants = n_distinct(ID),
      mean_eye_duration = mean(percentage_fix_duration, na.rm = TRUE),
      sd_eye_duration = sd(percentage_fix_duration, na.rm = TRUE),
      sem_eye_duration = sd_eye_duration / sqrt(n_sessions),
      mean_T2T3 = mean(T2T3, na.rm = TRUE),
      .groups = "drop"
    ),
  file.path(model_output_dir, "21_t2t3_response_eye_duration_descriptives_by_sex.csv")
)
write_csv(
  cortisol_duration %>%
    group_by(T2T3_response, Group) %>%
    summarise(
      n_sessions = n(),
      n_participants = n_distinct(ID),
      mean_eye_duration = mean(percentage_fix_duration, na.rm = TRUE),
      sd_eye_duration = sd(percentage_fix_duration, na.rm = TRUE),
      sem_eye_duration = sd_eye_duration / sqrt(n_sessions),
      mean_T2T3 = mean(T2T3, na.rm = TRUE),
      .groups = "drop"
    ),
  file.path(model_output_dir, "24_t2t3_response_eye_duration_descriptives_full_sample.csv")
)
write_csv(
  broom::tidy(contrast(
    emmeans(t2t3_responder_full_model, ~ Group | T2T3_response, type = "response"),
    method = "revpairwise"
  )) %>%
    mutate(
      model = "t2t3_eye_duration_responder_full_model",
      p_display = format_p_value(p.value),
      sig = significance_stars(p.value)
    ) %>%
    relocate(model),
  file.path(model_output_dir, "25_t2t3_response_group_contrasts_full_sample.csv")
)
write_csv(
  broom::tidy(emmeans(t2t3_responder_model, ~ Group | NAL_status * T2T3_response, type = "response")) %>%
    mutate(
      model = "t2t3_eye_duration_responder_model",
      p_display = NA_character_,
      sig = NA_character_
    ) %>%
    relocate(model),
  file.path(model_output_dir, "18_t2t3_response_group_emmeans.csv")
)
write_csv(
  broom::tidy(contrast(
    emmeans(t2t3_responder_model, ~ Group | NAL_status * T2T3_response, type = "response"),
    method = "revpairwise"
  )) %>%
    mutate(
      model = "t2t3_eye_duration_responder_model",
      p_display = format_p_value(p.value),
      sig = significance_stars(p.value)
    ) %>%
    relocate(model),
  file.path(model_output_dir, "19_t2t3_response_group_contrasts.csv")
)

t2t3_points <- cortisol_duration %>%
  group_by(ID, Group, medication) %>%
  summarise(
    T2T3 = mean(T2T3, na.rm = TRUE),
    .groups = "drop"
  )

auci_points <- cortisol_duration %>%
  group_by(ID, Group, medication) %>%
  summarise(
    AUCi = mean(AUCi, na.rm = TRUE),
    .groups = "drop"
  )

pred_grid <- expand_grid(
  T2T3_z = seq(
    min(cortisol_duration$T2T3_z, na.rm = TRUE),
    max(cortisol_duration$T2T3_z, na.rm = TRUE),
    length.out = 100
  ),
  Group = levels(cortisol_duration$Group),
  medication = "PLA",
  session = levels(cortisol_duration$session)[[1]],
  ID = cortisol_duration$ID[[1]]
) %>%
  mutate(
    medication = factor(medication, levels = medication_levels_plot),
    Group = factor(Group, levels = levels(cortisol_duration$Group)),
    session = factor(session, levels = levels(cortisol_duration$session)),
    ID = factor(ID, levels = levels(cortisol_duration$ID)),
    T2T3 = T2T3_z * attr(scale(cortisol_duration$T2T3), "scaled:scale") +
      attr(scale(cortisol_duration$T2T3), "scaled:center")
  )

pred <- predict(
  t2t3_duration_model,
  newdata = pred_grid,
  type = "response",
  se.fit = TRUE,
  re.form = NA
)
pred_df <- pred_grid %>%
  mutate(
    fit = pred$fit,
    se = pred$se.fit,
    lower = pmax(0, fit - 1.96 * se),
    upper = pmin(1, fit + 1.96 * se)
  )

auci_pred_grid <- expand_grid(
  AUCi_z = seq(
    min(cortisol_duration$AUCi_z, na.rm = TRUE),
    max(cortisol_duration$AUCi_z, na.rm = TRUE),
    length.out = 100
  ),
  Group = levels(cortisol_duration$Group),
  medication = "PLA",
  session = levels(cortisol_duration$session)[[1]],
  ID = cortisol_duration$ID[[1]]
) %>%
  mutate(
    medication = factor(medication, levels = medication_levels_plot),
    Group = factor(Group, levels = levels(cortisol_duration$Group)),
    session = factor(session, levels = levels(cortisol_duration$session)),
    ID = factor(ID, levels = levels(cortisol_duration$ID)),
    AUCi = AUCi_z * attr(scale(cortisol_duration$AUCi), "scaled:scale") +
      attr(scale(cortisol_duration$AUCi), "scaled:center")
  )

auci_pred <- predict(
  auci_duration_model,
  newdata = auci_pred_grid,
  type = "response",
  se.fit = TRUE,
  re.form = NA
)
auci_pred_df <- auci_pred_grid %>%
  mutate(
    fit = auci_pred$fit,
    se = auci_pred$se.fit,
    lower = pmax(0, fit - 1.96 * se),
    upper = pmin(1, fit + 1.96 * se)
  )

auci_trim_pred_grid <- expand_grid(
  AUCi_trim_z = seq(
    min(cortisol_duration_trim90$AUCi_trim_z, na.rm = TRUE),
    max(cortisol_duration_trim90$AUCi_trim_z, na.rm = TRUE),
    length.out = 100
  ),
  Group = levels(cortisol_duration_trim90$Group),
  medication = "PLA",
  session = levels(cortisol_duration_trim90$session)[[1]],
  ID = cortisol_duration_trim90$ID[[1]]
) %>%
  mutate(
    medication = factor(medication, levels = medication_levels_plot),
    Group = factor(Group, levels = levels(cortisol_duration_trim90$Group)),
    session = factor(session, levels = levels(cortisol_duration_trim90$session)),
    ID = factor(ID, levels = levels(cortisol_duration_trim90$ID)),
    AUCi = AUCi_trim_z * attr(scale(cortisol_duration_trim90$AUCi), "scaled:scale") +
      attr(scale(cortisol_duration_trim90$AUCi), "scaled:center")
  )

auci_trim_pred <- predict(
  auci_duration_trim90_model,
  newdata = auci_trim_pred_grid,
  type = "response",
  se.fit = TRUE,
  re.form = NA
)
auci_trim_pred_df <- auci_trim_pred_grid %>%
  mutate(
    fit = auci_trim_pred$fit,
    se = auci_trim_pred$se.fit,
    lower = pmax(0, fit - 1.96 * se),
    upper = pmin(1, fit + 1.96 * se)
  )

auci_binned_summary <- cortisol_duration_trim90 %>%
  group_by(Group, auci_bin) %>%
  summarise(
    AUCi = mean(AUCi, na.rm = TRUE),
    mean_duration = mean(percentage_fix_duration, na.rm = TRUE),
    sem_duration = sd(percentage_fix_duration, na.rm = TRUE) / sqrt(sum(!is.na(percentage_fix_duration))),
    n = sum(!is.na(percentage_fix_duration)),
    .groups = "drop"
  )

auci_response_summary <- cortisol_duration %>%
  group_by(NAL_status, AUCi_response, Group) %>%
  summarise(
    mean_duration = mean(percentage_fix_duration, na.rm = TRUE),
    sem_duration = sd(percentage_fix_duration, na.rm = TRUE) / sqrt(sum(!is.na(percentage_fix_duration))),
    n_sessions = sum(!is.na(percentage_fix_duration)),
    .groups = "drop"
  )

t2t3_response_summary <- cortisol_duration %>%
  group_by(NAL_status, T2T3_response, Group) %>%
  summarise(
    mean_duration = mean(percentage_fix_duration, na.rm = TRUE),
    sem_duration = sd(percentage_fix_duration, na.rm = TRUE) / sqrt(sum(!is.na(percentage_fix_duration))),
    n_sessions = sum(!is.na(percentage_fix_duration)),
    .groups = "drop"
  )

t2t3_response_sex_summary <- cortisol_duration %>%
  group_by(Sex, NAL_status, T2T3_response, Group) %>%
  summarise(
    mean_duration = mean(percentage_fix_duration, na.rm = TRUE),
    sem_duration = sd(percentage_fix_duration, na.rm = TRUE) / sqrt(sum(!is.na(percentage_fix_duration))),
    n_sessions = sum(!is.na(percentage_fix_duration)),
    .groups = "drop"
  )

t2t3_response_full_summary <- cortisol_duration %>%
  group_by(T2T3_response, Group) %>%
  summarise(
    mean_duration = mean(percentage_fix_duration, na.rm = TRUE),
    sem_duration = sd(percentage_fix_duration, na.rm = TRUE) / sqrt(sum(!is.na(percentage_fix_duration))),
    n_sessions = sum(!is.na(percentage_fix_duration)),
    .groups = "drop"
  )

auci_background_pred_grid <- expand_grid(
  AUCi_z = seq(
    min(cortisol_background_duration$AUCi_z, na.rm = TRUE),
    max(cortisol_background_duration$AUCi_z, na.rm = TRUE),
    length.out = 100
  ),
  Group = levels(cortisol_background_duration$Group),
  medication = "PLA",
  session = levels(cortisol_background_duration$session)[[1]],
  ID = cortisol_background_duration$ID[[1]]
) %>%
  mutate(
    medication = factor(medication, levels = medication_levels_plot),
    Group = factor(Group, levels = levels(cortisol_background_duration$Group)),
    session = factor(session, levels = levels(cortisol_background_duration$session)),
    ID = factor(ID, levels = levels(cortisol_background_duration$ID)),
    AUCi = AUCi_z * attr(scale(cortisol_background_duration$AUCi), "scaled:scale") +
      attr(scale(cortisol_background_duration$AUCi), "scaled:center")
  )

auci_background_pred <- predict(
  auci_background_duration_model,
  newdata = auci_background_pred_grid,
  type = "response",
  se.fit = TRUE,
  re.form = NA
)
auci_background_pred_df <- auci_background_pred_grid %>%
  mutate(
    fit = auci_background_pred$fit,
    se = auci_background_pred$se.fit,
    lower = pmax(0, fit - 1.96 * se),
    upper = pmin(1, fit + 1.96 * se)
  )

panel_a <- ggplot(raw_summary, aes(x = timepoint, y = mean_cortisol, color = medication, group = medication)) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.6) +
  geom_errorbar(aes(ymin = mean_cortisol - sem, ymax = mean_cortisol + sem), width = 0.08, linewidth = 0.35) +
  geom_point(
    data = tibble(
      medication = factor(medication_levels_plot, levels = medication_levels_plot),
      timepoint = factor("1", levels = levels(raw_summary$timepoint)),
      mean_cortisol = 0
    ),
    aes(x = timepoint, y = mean_cortisol, fill = medication),
    inherit.aes = FALSE,
    shape = 22,
    size = 0,
    alpha = 0,
    show.legend = TRUE
  ) +
  scale_color_manual(values = medication_values, breaks = medication_levels_plot, guide = "none") +
  scale_fill_manual(values = medication_values, breaks = medication_levels_plot, name = "Medication") +
  labs(
    title = "Raw Cortisol Across Timepoints",
    x = "Timepoint",
    y = "Cortisol (nmol/l)",
    fill = "Medication"
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    plot.title = element_text(size = 11, face = "bold", margin = margin(b = 2)),
    plot.margin = margin(2, 2, 2, 2)
  )

panel_b <- ggplot(cortisol_duration, aes(x = medication, y = T2T3)) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_violin(aes(fill = medication), alpha = 0.35, linewidth = 0.25, width = 0.75, color = NA, show.legend = FALSE) +
  geom_point(
    data = t2t3_points,
    aes(x = medication, y = T2T3, color = medication),
    inherit.aes = FALSE,
    position = position_jitter(width = 0.12, height = 0),
    alpha = 0.45,
    size = 0.8
  ) +
  stat_summary(fun = mean, geom = "point", color = "black", size = 1.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.14, color = "black", linewidth = 0.35) +
  scale_fill_manual(values = medication_values, breaks = medication_levels_plot, guide = "none") +
  scale_color_manual(values = medication_values, breaks = medication_levels_plot, guide = "none") +
  labs(
    title = "T2-T3 Cortisol Reactivity",
    x = "Medication",
    y = "T2-T3 (nmol/l)"
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    plot.title = element_text(size = 10, face = "bold", margin = margin(b = 2)),
    plot.margin = margin(2, 5, 2, 2)
  )

panel_c <- ggplot() +
  geom_point(
    data = cortisol_duration,
    aes(x = T2T3, y = percentage_fix_duration, color = Group),
    alpha = 0.35,
    size = 1.0,
    position = position_jitter(width = 0.03, height = 0)
  ) +
  geom_ribbon(
    data = pred_df,
    aes(x = T2T3, ymin = lower, ymax = upper, fill = Group),
    alpha = 0.16,
    color = NA
  ) +
  geom_line(
    data = pred_df,
    aes(x = T2T3, y = fit, color = Group),
    linewidth = 0.8
  ) +
  scale_color_manual(values = group_values) +
  scale_fill_manual(values = group_values) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Cortisol Reactivity and Eye Fixation",
    x = "T2-T3 Cortisol Reactivity",
    y = "% Eye Fixation Duration",
    color = "Group",
    fill = "Group"
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    plot.title = element_text(size = 11, face = "bold", margin = margin(b = 2)),
    plot.margin = margin(2, 2, 2, 2)
  )

panel_b_auci <- ggplot(cortisol_duration, aes(x = medication, y = AUCi, fill = medication)) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  geom_violin(alpha = 0.35, linewidth = 0.25, width = 0.75, color = NA) +
  geom_point(
    data = auci_points,
    aes(x = medication, y = AUCi, color = medication),
    inherit.aes = FALSE,
    position = position_jitter(width = 0.12, height = 0),
    alpha = 0.45,
    size = 0.8
  ) +
  stat_summary(fun = mean, geom = "point", color = "black", size = 1.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.14, color = "black", linewidth = 0.35) +
  scale_fill_manual(values = medication_values, breaks = medication_levels_plot) +
  scale_color_manual(values = medication_values, breaks = medication_levels_plot, guide = "none") +
  labs(
    title = "AUCi Cortisol Reactivity",
    x = "Medication",
    y = "AUCi"
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    plot.title = element_text(size = 11, face = "bold", margin = margin(b = 2)),
    plot.margin = margin(2, 2, 2, 2)
  )

panel_c_auci <- ggplot() +
  geom_point(
    data = cortisol_duration,
    aes(x = AUCi, y = percentage_fix_duration, color = Group),
    alpha = 0.35,
    size = 1.0,
    position = position_jitter(width = 0.03, height = 0)
  ) +
  geom_ribbon(
    data = auci_pred_df,
    aes(x = AUCi, ymin = lower, ymax = upper, fill = Group),
    alpha = 0.16,
    color = NA
  ) +
  geom_line(
    data = auci_pred_df,
    aes(x = AUCi, y = fit, color = Group),
    linewidth = 0.8
  ) +
  scale_color_manual(values = group_values) +
  scale_fill_manual(values = group_values) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "AUCi and Eye Fixation",
    x = "AUCi Cortisol Reactivity",
    y = "% Eye Fixation Duration",
    color = "Group",
    fill = "Group"
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    plot.title = element_text(size = 11, face = "bold", margin = margin(b = 2)),
    plot.margin = margin(2, 2, 2, 2)
  )

panel_c_auci_background <- ggplot() +
  geom_point(
    data = cortisol_background_duration,
    aes(x = AUCi, y = percentage_fix_duration, color = Group),
    alpha = 0.35,
    size = 1.0,
    position = position_jitter(width = 0.03, height = 0)
  ) +
  geom_ribbon(
    data = auci_background_pred_df,
    aes(x = AUCi, ymin = lower, ymax = upper, fill = Group),
    alpha = 0.16,
    color = NA
  ) +
  geom_line(
    data = auci_background_pred_df,
    aes(x = AUCi, y = fit, color = Group),
    linewidth = 0.8
  ) +
  scale_color_manual(values = group_values) +
  scale_fill_manual(values = group_values) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "AUCi and Background Fixation",
    x = "AUCi Cortisol Reactivity",
    y = "% Background Fixation Duration",
    color = "Group",
    fill = "Group"
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    plot.title = element_text(size = 11, face = "bold", margin = margin(b = 2)),
    plot.margin = margin(2, 2, 2, 2)
  )

panel_c_auci_robust <- ggplot() +
  geom_point(
    data = cortisol_duration,
    aes(x = AUCi, y = percentage_fix_duration, color = Group),
    alpha = 0.12,
    size = 0.8,
    position = position_jitter(width = 0.03, height = 0)
  ) +
  annotate(
    "rect",
    xmin = auci_trim_limits[[1]],
    xmax = auci_trim_limits[[2]],
    ymin = -Inf,
    ymax = Inf,
    fill = "grey50",
    alpha = 0.055
  ) +
  geom_ribbon(
    data = auci_trim_pred_df,
    aes(x = AUCi, ymin = lower, ymax = upper, fill = Group),
    alpha = 0.14,
    color = NA
  ) +
  geom_line(
    data = auci_trim_pred_df,
    aes(x = AUCi, y = fit, color = Group),
    linewidth = 0.8
  ) +
  geom_errorbar(
    data = auci_binned_summary,
    aes(
      x = AUCi,
      ymin = pmax(0, mean_duration - sem_duration),
      ymax = pmin(1, mean_duration + sem_duration),
      color = Group
    ),
    width = 0,
    linewidth = 0.35
  ) +
  geom_point(
    data = auci_binned_summary,
    aes(x = AUCi, y = mean_duration, fill = Group),
    shape = 21,
    color = "black",
    stroke = 0.25,
    size = 2.0
  ) +
  scale_color_manual(values = group_values) +
  scale_fill_manual(values = group_values) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "AUCi and Eye Fixation: Robust View",
    x = "AUCi Cortisol Reactivity",
    y = "% Eye Fixation Duration",
    color = "Group",
    fill = "Group"
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    plot.title = element_text(size = 11, face = "bold", margin = margin(b = 2)),
    plot.margin = margin(2, 2, 2, 2)
  )

panel_c_auci_responder <- ggplot(
  cortisol_duration,
  aes(x = AUCi_response, y = percentage_fix_duration, color = Group, fill = Group)
) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.10, jitter.height = 0, dodge.width = 0.48),
    alpha = 0.18,
    size = 0.75,
    stroke = 0
  ) +
  geom_errorbar(
    data = auci_response_summary,
    aes(
      y = mean_duration,
      ymin = pmax(0, mean_duration - sem_duration),
      ymax = pmin(1, mean_duration + sem_duration)
    ),
    position = position_dodge(width = 0.48),
    width = 0.12,
    linewidth = 0.35,
    color = "black"
  ) +
  geom_point(
    data = auci_response_summary,
    aes(y = mean_duration),
    position = position_dodge(width = 0.48),
    shape = 21,
    color = "black",
    stroke = 0.25,
    size = 2.1
  ) +
  facet_wrap(~ NAL_status, nrow = 1) +
  scale_x_discrete(labels = c("Decrease", "No response", "Increase")) +
  scale_color_manual(values = group_values) +
  scale_fill_manual(values = group_values) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Eye Fixation by AUCi Response",
    x = paste0("AUCi response (threshold +/-", auci_response_threshold, ")"),
    y = "% Eye Fixation Duration",
    color = "Group",
    fill = "Group"
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    plot.title = element_text(size = 11, face = "bold", margin = margin(b = 2)),
    plot.margin = margin(2, 2, 2, 2)
  )

panel_c_t2t3_responder <- ggplot(
  cortisol_duration,
  aes(x = T2T3_response, y = percentage_fix_duration, color = Group, fill = Group)
) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.10, jitter.height = 0, dodge.width = 0.48),
    alpha = 0.18,
    size = 0.75,
    stroke = 0
  ) +
  geom_errorbar(
    data = t2t3_response_summary,
    aes(
      y = mean_duration,
      ymin = pmax(0, mean_duration - sem_duration),
      ymax = pmin(1, mean_duration + sem_duration)
    ),
    position = position_dodge(width = 0.48),
    width = 0.12,
    linewidth = 0.35,
    color = "black"
  ) +
  geom_point(
    data = t2t3_response_summary,
    aes(y = mean_duration),
    position = position_dodge(width = 0.48),
    shape = 21,
    color = "black",
    stroke = 0.25,
    size = 2.1
  ) +
  facet_wrap(~ NAL_status, nrow = 1) +
  scale_x_discrete(labels = c("Decrease", "No response", "Increase")) +
  scale_color_manual(values = group_values) +
  scale_fill_manual(values = group_values) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Eye Fixation by T2-T3 Response",
    x = paste0("T2-T3 response (threshold +/-", t2t3_response_threshold, " nmol/l)"),
    y = "% Eye Fixation Duration",
    color = "Group",
    fill = "Group"
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    plot.title = element_text(size = 11, face = "bold", margin = margin(b = 2)),
    plot.margin = margin(2, 2, 2, 2)
  )

panel_t2t3_responder_sex <- ggplot(
  cortisol_duration,
  aes(x = T2T3_response, y = percentage_fix_duration, color = Group, fill = Group)
) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.10, jitter.height = 0, dodge.width = 0.48),
    alpha = 0.14,
    size = 0.65,
    stroke = 0
  ) +
  geom_errorbar(
    data = t2t3_response_sex_summary,
    aes(
      y = mean_duration,
      ymin = pmax(0, mean_duration - sem_duration),
      ymax = pmin(1, mean_duration + sem_duration)
    ),
    position = position_dodge(width = 0.48),
    width = 0.12,
    linewidth = 0.32,
    color = "black"
  ) +
  geom_point(
    data = t2t3_response_sex_summary,
    aes(y = mean_duration),
    position = position_dodge(width = 0.48),
    shape = 21,
    color = "black",
    stroke = 0.25,
    size = 1.8
  ) +
  geom_text(
    data = t2t3_response_sex_summary,
    aes(y = 0.98, label = n_sessions),
    position = position_dodge(width = 0.48),
    color = "grey30",
    size = 1.8,
    show.legend = FALSE
  ) +
  facet_grid(Sex ~ NAL_status) +
  scale_x_discrete(labels = c("Dec.", "No resp.", "Inc.")) +
  scale_color_manual(values = group_values) +
  scale_fill_manual(values = group_values) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Eye Fixation by T2-T3 Response, Sex, and NAL Status",
    x = paste0("T2-T3 response (+/-", t2t3_response_threshold, " nmol/l)"),
    y = "% Eye Fixation Duration",
    color = "Group",
    fill = "Group"
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    panel.spacing = grid::unit(0.55, "lines"),
    plot.title = element_text(size = 11, face = "bold", margin = margin(b = 2)),
    plot.margin = margin(2, 2, 2, 2)
  )

panel_c_t2t3_responder_full <- ggplot(
  cortisol_duration,
  aes(x = T2T3_response, y = percentage_fix_duration, color = Group, fill = Group)
) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.10, jitter.height = 0, dodge.width = 0.48),
    alpha = 0.18,
    size = 0.75,
    stroke = 0
  ) +
  geom_errorbar(
    data = t2t3_response_full_summary,
    aes(
      y = mean_duration,
      ymin = pmax(0, mean_duration - sem_duration),
      ymax = pmin(1, mean_duration + sem_duration)
    ),
    position = position_dodge(width = 0.48),
    width = 0.12,
    linewidth = 0.35,
    color = "black"
  ) +
  geom_point(
    data = t2t3_response_full_summary,
    aes(y = mean_duration),
    position = position_dodge(width = 0.48),
    shape = 21,
    color = "black",
    stroke = 0.25,
    size = 2.1
  ) +
  scale_x_discrete(labels = c("Decrease", "No response", "Increase")) +
  scale_color_manual(values = group_values) +
  scale_fill_manual(values = group_values) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Eye Fixation by T2-T3 Response",
    x = paste0("T2-T3 response (+/-", t2t3_response_threshold, " nmol/l)"),
    y = "% Eye Fixation Duration",
    color = "Group",
    fill = "Group"
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
    plot.title = element_text(size = 11, face = "bold", margin = margin(b = 2)),
    plot.margin = margin(2, 2, 2, 2)
  )

combined_plot <- (panel_a | panel_b | panel_c) +
  plot_layout(widths = c(1.05, 0.85, 1.15), guides = "collect") +
  plot_annotation(
    title = "Cortisol Reactivity and Eye-Directed Fixation Duration (Exploratory)",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 2)),
      plot.tag = element_text(size = 10, face = "plain"),
      plot.margin = margin(1, 1, 1, 1)
    )
  ) &
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.key = element_rect(fill = NA, color = NA)
  ) &
  guides(
    fill = guide_legend(
      override.aes = list(shape = 22, size = 3.5, alpha = 1, color = "black", linetype = 0)
    )
  )

combined_plot_auci <- (panel_a | panel_b_auci | panel_c_auci) +
  plot_layout(widths = c(1.05, 0.85, 1.15), guides = "collect") +
  plot_annotation(
    title = "AUCi Cortisol Reactivity and Eye-Directed Fixation Duration (Exploratory)",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 2)),
      plot.tag = element_text(size = 10, face = "plain"),
      plot.margin = margin(1, 1, 1, 1)
    )
  ) &
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.key = element_rect(fill = NA, color = NA)
  ) &
  guides(
    fill = guide_legend(
      override.aes = list(shape = 22, size = 3.5, alpha = 1, color = "black", linetype = 0)
    )
  )

combined_plot_auci_background <- (panel_a | panel_b_auci | panel_c_auci_background) +
  plot_layout(widths = c(1.05, 0.85, 1.15), guides = "collect") +
  plot_annotation(
    title = "AUCi Cortisol Reactivity and Background Fixation Duration (Exploratory)",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 2)),
      plot.tag = element_text(size = 10, face = "plain"),
      plot.margin = margin(1, 1, 1, 1)
    )
  ) &
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.key = element_rect(fill = NA, color = NA)
  ) &
  guides(
    fill = guide_legend(
      override.aes = list(shape = 22, size = 3.5, alpha = 1, color = "black", linetype = 0)
    )
  )

combined_plot_auci_robust <- (panel_a | panel_b_auci | panel_c_auci_robust) +
  plot_layout(widths = c(1.05, 0.85, 1.15), guides = "collect") +
  plot_annotation(
    title = "AUCi Cortisol Reactivity and Eye-Directed Fixation Duration: Robust View (Exploratory)",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 2)),
      plot.tag = element_text(size = 10, face = "plain"),
      plot.margin = margin(1, 1, 1, 1)
    )
  ) &
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.key = element_rect(fill = NA, color = NA)
  ) &
  guides(
    fill = guide_legend(
      override.aes = list(shape = 22, size = 3.5, alpha = 1, color = "black", linetype = 0)
    )
  )

combined_plot_auci_responder <- (panel_a | panel_b_auci | panel_c_auci_responder) +
  plot_layout(widths = c(0.90, 0.72, 1.58), guides = "collect") +
  plot_annotation(
    title = "AUCi Cortisol Response Groups and Eye-Directed Fixation Duration (Exploratory)",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 2)),
      plot.tag = element_text(size = 10, face = "plain"),
      plot.margin = margin(1, 1, 1, 1)
    )
  ) &
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.key = element_rect(fill = NA, color = NA)
  ) &
  guides(
    fill = guide_legend(
      override.aes = list(shape = 22, size = 3.5, alpha = 1, color = "black", linetype = 0)
    )
  )

combined_plot_t2t3_responder <- (panel_a | panel_b | panel_c_t2t3_responder) +
  plot_layout(widths = c(0.90, 0.72, 1.58), guides = "collect") +
  plot_annotation(
    title = "T2-T3 Cortisol Response Groups and Eye-Directed Fixation Duration (Exploratory)",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 2)),
      plot.tag = element_text(size = 10, face = "plain"),
      plot.margin = margin(1, 1, 1, 1)
    )
  ) &
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.key = element_rect(fill = NA, color = NA)
  ) &
  guides(
    fill = guide_legend(
      override.aes = list(shape = 22, size = 3.5, alpha = 1, color = "black", linetype = 0)
    )
  )

combined_plot_t2t3_responder_sex <- (panel_a | panel_b) / panel_t2t3_responder_sex +
  plot_layout(heights = c(0.78, 1.22), guides = "collect") +
  plot_annotation(
    title = "T2-T3 Cortisol Response Groups and Eye-Directed Fixation Duration by Sex (Exploratory)",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 2)),
      plot.tag = element_text(size = 10, face = "plain"),
      plot.margin = margin(1, 1, 1, 1)
    )
  ) &
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.key = element_rect(fill = NA, color = NA)
  ) &
  guides(
    fill = guide_legend(
      override.aes = list(shape = 22, size = 3.5, alpha = 1, color = "black", linetype = 0)
    )
  )

combined_plot_t2t3_responder_full <- (panel_a | panel_b | panel_c_t2t3_responder_full) +
  plot_layout(widths = c(0.98, 0.82, 1.30), guides = "collect") +
  plot_annotation(
    title = "T2-T3 Cortisol Response Groups and Eye-Directed Fixation Duration (Exploratory)",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0, margin = margin(b = 2)),
      plot.tag = element_text(size = 10, face = "plain"),
      plot.margin = margin(1, 1, 1, 1)
    )
  ) &
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.key = element_rect(fill = NA, color = NA)
  ) &
  guides(
    fill = guide_legend(
      override.aes = list(shape = 22, size = 3.5, alpha = 1, color = "black", linetype = 0)
    )
  )

ggsave(output_png, combined_plot_t2t3_responder_full, units = "mm", width = 200, height = 85, dpi = 300)
ggsave(output_t2t3_responder_full_svg, combined_plot_t2t3_responder_full, units = "mm", width = 200, height = 85)
ggsave(output_t2t3_responder_full_png, combined_plot_t2t3_responder_full, units = "mm", width = 200, height = 85, dpi = 300)
