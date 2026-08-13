# Manuscript section: Supplementary sensitivity analyses
# Analysis family: ASD medication/comorbidity AOI profile visualization
# Primary input dataset(s): anonymized sensitivity datasets
# Primary output(s): supplementary sensitivity profile plot

library(tidyverse)
library(here)
library(patchwork)

normalize_session <- function(x) {
  x %>%
    str_trim() %>%
    str_to_lower() %>%
    str_replace_all(" ", "_")
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
  filter(Group == "ASD", Sex %in% c("f", "m"), !is.na(percentage_fix_duration)) %>%
  mutate(
    fix = factor(
      fix,
      levels = c("fix_on_eyes", "fix_on_mouth", "fix_on_face", "fix_on_background"),
      labels = c("Eyes", "Mouth", "Face", "Background")
    ),
    medication_type_group = factor(
      medication_type_group,
      levels = c("No medication record", "SSRI/SNRI only", "Other/multiple medication")
    ),
    disorder_type_group = factor(
      disorder_type_group,
      levels = c("No disorder record", "ADHD only", "Depression only", "Other/multiple disorder")
    )
  )

plot_df <- bind_rows(
  df %>%
    mutate(predictor = "Prescribed medication type", group = as.character(medication_type_group)),
  df %>%
    mutate(predictor = "Comorbidity type", group = as.character(disorder_type_group))
) %>%
  group_by(predictor, group, fix) %>%
  summarise(
    mean_duration = mean(percentage_fix_duration, na.rm = TRUE),
    sem = sd(percentage_fix_duration, na.rm = TRUE) / sqrt(n()),
    n_observations = n(),
    n_participants = n_distinct(anon_id),
    .groups = "drop"
  ) %>%
  mutate(
    lower = mean_duration - sem,
    upper = mean_duration + sem
  )

group_labels <- plot_df %>%
  distinct(predictor, group, n_participants) %>%
  mutate(group_label = paste0(group, " (N = ", n_participants, ")"))

plot_df <- plot_df %>%
  left_join(group_labels, by = c("predictor", "group", "n_participants"))

output_dir <- here("02_outputs", "figures", "supplementary", "medication_disorder_sensitivity")
dir.create(file.path(output_dir, "png"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "svg"), recursive = TRUE, showWarnings = FALSE)

write_csv(
  plot_df,
  here(
    "02_outputs", "model_outputs", "main_manuscript", "eye_tracking",
    "group_aoi_medication_duration", "sensitivity", "medication_disorder",
    "18_type_aoi_profile_plot_data.csv"
  )
)

make_profile_plot <- function(data, title, color_values) {
  ggplot(
    data,
    aes(x = fix, y = mean_duration, group = group_label, color = group_label)
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.12, linewidth = 0.55) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
    scale_color_manual(values = color_values) +
    labs(title = title, x = NULL, y = "Fixation-duration proportion", color = NULL) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 13),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = 9),
      legend.key.width = unit(0.55, "lines"),
      legend.spacing.x = unit(0.12, "cm"),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
}

medication_labels <- group_labels %>%
  filter(predictor == "Prescribed medication type") %>%
  arrange(match(group, c("No medication record", "SSRI/SNRI only", "Other/multiple medication")))

disorder_labels <- group_labels %>%
  filter(predictor == "Comorbidity type") %>%
  arrange(match(group, c("No disorder record", "ADHD only", "Depression only", "Other/multiple disorder")))

medication_colors <- setNames(
  c("#2F5D62", "#B35C44", "#6C6AA8")[seq_len(nrow(medication_labels))],
  medication_labels$group_label
)
disorder_colors <- setNames(
  c("#2F5D62", "#D08B2E", "#8F3F71", "#6C6AA8")[seq_len(nrow(disorder_labels))],
  disorder_labels$group_label
)

medication_plot <- make_profile_plot(
  plot_df %>% filter(predictor == "Prescribed medication type"),
  "Gaze by medication type",
  medication_colors
)

disorder_plot <- make_profile_plot(
  plot_df %>% filter(predictor == "Comorbidity type"),
  "Gaze by comorbidity type",
  disorder_colors
)

profile_plot <- medication_plot / disorder_plot + plot_layout(heights = c(1, 1))

ggsave(
  file.path(output_dir, "png", "medication_disorder_aoi_profiles.png"),
  profile_plot,
  width = 8,
  height = 7,
  dpi = 300
)
ggsave(
  file.path(output_dir, "svg", "medication_disorder_aoi_profiles.svg"),
  profile_plot,
  width = 8,
  height = 7
)

print(profile_plot)
