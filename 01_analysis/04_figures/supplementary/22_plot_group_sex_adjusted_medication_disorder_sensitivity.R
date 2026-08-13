# Manuscript section: Supplementary sensitivity analyses
# Analysis family: adjusted Group x Sex gaze effect after ASD medication/comorbidity adjustment
# Primary input dataset(s): anonymized sensitivity datasets
# Primary output(s): adjusted Group x Sex profile plot for fixation-duration model

library(tidyverse)
library(glmmTMB)
library(emmeans)
library(here)

normalize_session <- function(x) {
  x %>%
    str_trim() %>%
    str_to_lower() %>%
    str_replace_all(" ", "_")
}

format_p <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "p < .001",
    TRUE ~ paste0("p = ", sub("^0", "", sprintf("%.3f", p)))
  )
}

read_anon_csv <- function(filename) {
  read_csv(
    here("00_data", "derived", "sensitivity", "medication_disorder_anonymized", filename),
    show_col_types = FALSE
  ) %>%
    rename_with(str_trim)
}

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

participant_groups <- session_metadata %>%
  distinct(anon_id, Group, Sex) %>%
  left_join(medication_type_groups, by = "anon_id") %>%
  left_join(disorder_type_groups, by = "anon_id") %>%
  mutate(
    medication_type_group = replace_na(medication_type_group, "No medication record"),
    disorder_type_group = replace_na(disorder_type_group, "No disorder record")
  )

df <- read_anon_csv("duration_data_anonymized.csv") %>%
  mutate(session_norm = normalize_session(session)) %>%
  left_join(session_metadata, by = c("anon_id", "session_norm")) %>%
  mutate(session = coalesce(session.y, session.x)) %>%
  left_join(
    participant_groups %>% select(anon_id, medication_type_group, disorder_type_group),
    by = "anon_id"
  ) %>%
  filter(!is.na(Group), Sex %in% c("f", "m"), !is.na(percentage_fix_duration)) %>%
  mutate(
    fix = factor(
      fix,
      levels = c("fix_on_eyes", "fix_on_mouth", "fix_on_face", "fix_on_background"),
      labels = c("Eyes", "Mouth", "Face", "Background")
    ),
    Group = factor(Group, levels = c("ASD", "CTRL")),
    Sex = factor(Sex, levels = c("f", "m"), labels = c("Female", "Male")),
    medication = factor(
      medication,
      levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA"),
      labels = c("BOTH", "NAL", "OXT", "PLA")
    ),
    session = factor(session),
    asd_med_ssri = as.integer(Group == "ASD" & medication_type_group == "SSRI/SNRI only"),
    asd_med_other = as.integer(Group == "ASD" & medication_type_group == "Other/multiple medication"),
    asd_dis_adhd = as.integer(Group == "ASD" & disorder_type_group == "ADHD only"),
    asd_dis_dep = as.integer(Group == "ASD" & disorder_type_group == "Depression only"),
    asd_dis_other = as.integer(Group == "ASD" & disorder_type_group == "Other/multiple disorder")
  )

adjusted_model <- glmmTMB(
  percentage_fix_duration ~ Group * fix * medication + Group * fix * Sex + session +
    (asd_med_ssri + asd_med_other + asd_dis_adhd + asd_dis_dep + asd_dis_other) * fix,
  family = beta_family(),
  data = df,
  control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
)

emm <- emmeans(
  adjusted_model,
  ~ Group | fix * Sex,
  type = "response",
  at = list(
    asd_med_ssri = 0,
    asd_med_other = 0,
    asd_dis_adhd = 0,
    asd_dis_dep = 0,
    asd_dis_other = 0
  )
)

plot_df <- as_tibble(emm) %>%
  filter(fix %in% c("Eyes", "Background")) %>%
  mutate(
    fix = factor(fix, levels = c("Eyes", "Background")),
    Sex = factor(Sex, levels = c("Female", "Male"))
  )

contrast_df <- broom::tidy(emmeans(
  adjusted_model,
  pairwise ~ Group | fix * Sex,
  type = "response",
  at = list(
    asd_med_ssri = 0,
    asd_med_other = 0,
    asd_dis_adhd = 0,
    asd_dis_dep = 0,
    asd_dis_other = 0
  )
)$contrasts) %>%
  filter(fix %in% c("Eyes", "Background")) %>%
  mutate(
    fix = factor(fix, levels = c("Eyes", "Background")),
    Sex = factor(Sex, levels = c("Female", "Male")),
    label = format_p(p.value)
  )

label_df <- plot_df %>%
  group_by(fix, Sex) %>%
  summarise(y = max(asymp.UCL, na.rm = TRUE) + 0.045, .groups = "drop") %>%
  left_join(contrast_df %>% select(fix, Sex, label), by = c("fix", "Sex")) %>%
  mutate(x = 1.5)

output_model_dir <- here(
  "02_outputs", "model_outputs", "main_manuscript", "eye_tracking",
  "group_aoi_medication_duration", "sensitivity", "medication_disorder"
)
dir.create(output_model_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(plot_df, file.path(output_model_dir, "20_adjusted_group_sex_plot_data.csv"))
write_csv(contrast_df, file.path(output_model_dir, "21_adjusted_group_sex_contrasts.csv"))

output_dir <- here("02_outputs", "figures", "supplementary", "medication_disorder_sensitivity")
dir.create(file.path(output_dir, "png"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "svg"), recursive = TRUE, showWarnings = FALSE)

group_colors <- c("ASD" = "#D55E00", "CTRL" = "#0072B2")

adjusted_plot <- ggplot(plot_df, aes(x = Group, y = response, color = Group)) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.12, linewidth = 0.55) +
  geom_line(aes(group = Sex), color = "grey50", linewidth = 0.75) +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.4
  ) +
  facet_grid(fix ~ Sex) +
  scale_color_manual(values = group_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(
    title = "Adjusted Group x Sex gaze pattern",
    x = NULL,
    y = "Adjusted fixation-duration proportion",
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0, size = 13),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10.5),
    panel.spacing = unit(6, "mm")
  )

ggsave(
  file.path(output_dir, "png", "adjusted_group_sex_gaze_sensitivity.png"),
  adjusted_plot,
  width = 7,
  height = 5,
  dpi = 300
)
ggsave(
  file.path(output_dir, "svg", "adjusted_group_sex_gaze_sensitivity.svg"),
  adjusted_plot,
  width = 7,
  height = 5
)

print(adjusted_plot)
