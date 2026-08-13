# Manuscript section: Supplementary sensitivity analyses
# Analysis family: ASD medication/comorbidity type and cortisol reactivity
# Primary input dataset(s): df_cortisol_min_max.csv; anonymized medication/disorder files; private ID lookup
# Primary output(s): supplementary cortisol sensitivity profile plot

library(tidyverse)
library(here)
library(patchwork)

read_anon_csv <- function(filename) {
  read_csv(
    here("00_data", "derived", "sensitivity", "medication_disorder_anonymized", filename),
    show_col_types = FALSE
  ) %>%
    rename_with(str_trim)
}

lookup_path <- here("00_data", "private", "id_lookup.csv")
if (!file.exists(lookup_path)) {
  stop(
    "Private lookup table not found at 00_data/private/id_lookup.csv.",
    call. = FALSE
  )
}

id_lookup <- read_csv(lookup_path, show_col_types = FALSE) %>%
  transmute(
    anon_id = as.integer(anon_id),
    ID = as.character(original_id)
  )

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

df <- read_csv(
  here("00_data", "derived", "analysis", "df_cortisol_min_max.csv"),
  show_col_types = FALSE
) %>%
  rename(T2T3 = MinMax) %>%
  mutate(ID = as.character(ID)) %>%
  left_join(id_lookup, by = "ID") %>%
  filter(Group == "ASD", !is.na(anon_id), !is.na(T2T3), !is.na(medication)) %>%
  left_join(medication_type_groups, by = "anon_id") %>%
  left_join(disorder_type_groups, by = "anon_id") %>%
  mutate(
    medication_type_group = replace_na(medication_type_group, "No medication record"),
    disorder_type_group = replace_na(disorder_type_group, "No disorder record"),
    medication = factor(
      medication,
      levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA"),
      labels = c("BOTH", "NAL", "OXT", "PLA")
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
  group_by(predictor, group, medication) %>%
  summarise(
    mean_t2t3 = mean(T2T3, na.rm = TRUE),
    sem = sd(T2T3, na.rm = TRUE) / sqrt(n()),
    n_observations = n(),
    n_participants = n_distinct(anon_id),
    .groups = "drop"
  ) %>%
  mutate(
    lower = mean_t2t3 - sem,
    upper = mean_t2t3 + sem
  )

group_labels <- plot_df %>%
  group_by(predictor, group) %>%
  summarise(n_participants_total = max(n_participants, na.rm = TRUE), .groups = "drop") %>%
  mutate(group_label = paste0(group, " (N = ", n_participants_total, ")"))

plot_df <- plot_df %>%
  left_join(group_labels, by = c("predictor", "group"))

output_model_dir <- here(
  "02_outputs", "model_outputs", "main_manuscript", "cortisol",
  "sensitivity", "medication_disorder"
)
dir.create(output_model_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(plot_df, file.path(output_model_dir, "06_cortisol_type_plot_data.csv"))

output_dir <- here("02_outputs", "figures", "supplementary", "medication_disorder_sensitivity")
dir.create(file.path(output_dir, "png"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "svg"), recursive = TRUE, showWarnings = FALSE)

make_cortisol_plot <- function(data, title, color_values) {
  ggplot(
    data,
    aes(x = medication, y = mean_t2t3, group = group_label, color = group_label)
  ) +
    geom_hline(yintercept = 0, color = "grey65", linewidth = 0.45) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.12, linewidth = 0.55) +
    scale_color_manual(values = color_values) +
    labs(title = title, x = NULL, y = "T2-T3 cortisol change", color = NULL) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 13),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = 9),
      legend.key.width = unit(0.55, "lines"),
      legend.spacing.x = unit(0.12, "cm")
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

medication_plot <- make_cortisol_plot(
  plot_df %>% filter(predictor == "Prescribed medication type"),
  "Cortisol by medication type",
  medication_colors
)

disorder_plot <- make_cortisol_plot(
  plot_df %>% filter(predictor == "Comorbidity type"),
  "Cortisol by comorbidity type",
  disorder_colors
)

cortisol_plot <- medication_plot / disorder_plot + plot_layout(heights = c(1, 1))

ggsave(
  file.path(output_dir, "png", "cortisol_medication_disorder_sensitivity.png"),
  cortisol_plot,
  width = 8,
  height = 7,
  dpi = 300
)
ggsave(
  file.path(output_dir, "svg", "cortisol_medication_disorder_sensitivity.svg"),
  cortisol_plot,
  width = 8,
  height = 7
)

print(cortisol_plot)
