# Manuscript section: Main eye-tracking figure (fixation duration)
# Analysis family: figure generation
# Primary input dataset(s): 00_data/derived/preprocessing/duration_data.csv; 00_data/derived/analysis/df_analysis_pub_prep.csv
# Primary output(s): duration-based Figure 2 using final no-random-effects post-hocs

library(tidyverse)
library(ggdist)
library(ggh4x)
library(patchwork)
library(here)

output_svg <- here("02_outputs", "figures", "main", "eye_tracking")
output_png <- here("02_outputs", "figures", "main", "eye_tracking")
dir.create(output_svg, recursive = TRUE, showWarnings = FALSE)
dir.create(output_png, recursive = TRUE, showWarnings = FALSE)

normalize_session <- function(x) {
  x %>%
    str_trim() %>%
    str_to_lower() %>%
    str_replace_all(" ", "_")
}

format_sig <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

group_offsets <- c("BOTH" = -0.30, "NAL" = -0.10, "OXT" = 0.10, "PLA" = 0.30)
medication_values <- c("BOTH" = "#FDC010", "NAL" = "#C7B655", "OXT" = "#93AD7C", "PLA" = "#0B7A6B")
group_values <- c("ASD" = "#D55E00", "CTRL" = "#0072B2")

session_metadata <- read_csv(
  here("00_data", "derived", "analysis", "df_analysis_pub_prep.csv"),
  show_col_types = FALSE
) %>%
  distinct(ID, session, Sex, Group, medication) %>%
  mutate(session_norm = normalize_session(session))

df <- read_csv(here("00_data", "derived", "preprocessing", "duration_data.csv"), show_col_types = FALSE) %>%
  mutate(ID = Sub_ID, session_norm = normalize_session(session)) %>%
  left_join(session_metadata, by = c("ID", "session_norm")) %>%
  mutate(session = coalesce(session.y, session.x)) %>%
  filter(!is.na(Group), !is.na(medication)) %>%
  mutate(
    fix = factor(
      fix,
      levels = c("fix_on_eyes", "fix_on_mouth", "fix_on_face", "fix_on_background"),
      labels = c("Eyes", "Mouth", "Face", "Background"),
      ordered = TRUE
    ),
    medication = factor(
      medication,
      levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA"),
      labels = c("BOTH", "NAL", "OXT", "PLA")
    ),
    Group = factor(Group, levels = c("ASD", "CTRL"))
  ) %>%
  select(-session_norm, -session.x, -session.y)

df_binary <- df %>%
  filter(Sex %in% c("f", "m")) %>%
  droplevels()

df_plot <- bind_rows(
  df_binary %>% mutate(Panel = "Overall"),
  df_binary %>% filter(Sex == "f") %>% mutate(Panel = "Female"),
  df_binary %>% filter(Sex == "m") %>% mutate(Panel = "Male")
) %>%
  mutate(
    Panel = factor(Panel, levels = c("Overall", "Female", "Male")),
    fix = factor(as.character(fix), levels = c("Eyes", "Background"), ordered = TRUE)
  ) %>%
  filter(!is.na(fix)) %>%
  droplevels()

group_only_summary <- df_plot %>%
  group_by(Panel, Group, fix) %>%
  summarise(
    mean_prop = mean(percentage_fix_duration, na.rm = TRUE),
    sem = sd(percentage_fix_duration, na.rm = TRUE) / sqrt(n()),
    lower = mean_prop - sem,
    upper = mean_prop + sem,
    .groups = "drop"
  )

med_split_summary <- df_plot %>%
  group_by(Panel, Group, fix, medication) %>%
  summarise(
    mean_prop = mean(percentage_fix_duration, na.rm = TRUE),
    sem = sd(percentage_fix_duration, na.rm = TRUE) / sqrt(n()),
    lower = mean_prop - sem,
    upper = mean_prop + sem,
    .groups = "drop"
  )

model_dir <- here("02_outputs", "model_outputs", "main_manuscript", "eye_tracking", "group_aoi_medication_duration")

overall_contrasts <- read_csv(
  file.path(model_dir, "03_group_by_fix.csv"),
  show_col_types = FALSE
) %>%
  transmute(Panel = "Overall", fix, sig = format_sig(p.value)) %>%
  filter(fix %in% c("Eyes", "Background"))

sex_contrasts <- read_csv(
  file.path(model_dir, "04_group_by_fix_sex.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    Panel = recode(Sex, "f" = "Female", "m" = "Male"),
    fix,
    sig = format_sig(p.value)
  ) %>%
  filter(fix %in% c("Eyes", "Background"))

group_only_sig <- bind_rows(overall_contrasts, sex_contrasts) %>%
  filter(sig != "") %>%
  left_join(
    group_only_summary %>%
      group_by(Panel, fix) %>%
      summarise(y = max(upper, na.rm = TRUE) + 0.055, .groups = "drop"),
    by = c("Panel", "fix")
  ) %>%
  mutate(
    Panel = factor(Panel, levels = c("Overall", "Female", "Male")),
    fix = factor(fix, levels = c("Eyes", "Background"), ordered = TRUE),
    x = 1.5,
    xstart = 1,
    xend = 2,
    ystem = y - 0.024
  )

med_overall_contrasts <- read_csv(
  file.path(model_dir, "05_group_by_fix_medication.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    Panel = "Overall",
    fix,
    medication = factor(
      medication,
      levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA"),
      labels = c("BOTH", "NAL", "OXT", "PLA")
    ),
    sig = format_sig(p.value)
  )

med_sex_contrasts <- read_csv(
  file.path(model_dir, "06_group_by_fix_medication_sex.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    Panel = recode(Sex, "f" = "Female", "m" = "Male"),
    fix,
    medication = factor(
      medication,
      levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA"),
      labels = c("BOTH", "NAL", "OXT", "PLA")
    ),
    sig = format_sig(p.value)
  )

med_split_sig <- bind_rows(med_overall_contrasts, med_sex_contrasts) %>%
  filter(fix %in% c("Eyes", "Background"), sig != "") %>%
  left_join(
    med_split_summary %>%
      group_by(Panel, fix, medication) %>%
      summarise(y = max(upper, na.rm = TRUE), .groups = "drop"),
    by = c("Panel", "fix", "medication")
  ) %>%
  group_by(Panel, fix) %>%
  arrange(medication, .by_group = TRUE) %>%
  mutate(
    Panel = factor(Panel, levels = c("Overall", "Female", "Male")),
    fix = factor(fix, levels = c("Eyes", "Background"), ordered = TRUE),
    y = max(y, na.rm = TRUE) + 0.065 + (row_number() - 1) * 0.08,
    ystem = y - 0.026,
    x = 1.5 + unname(group_offsets[as.character(medication)]),
    xstart = 1 + unname(group_offsets[as.character(medication)]),
    xend = 2 + unname(group_offsets[as.character(medication)]),
    color = unname(medication_values[as.character(medication)])
  ) %>%
  ungroup()

base_duration_theme <- theme_classic(base_size = 8) +
  theme(
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(color = "black", size = 9),
    plot.title = element_text(size = 10, face = "bold", hjust = 0),
    panel.spacing.y = unit(6, "mm"),
    panel.spacing.x = unit(5, "mm"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 6.5),
    legend.key.height = unit(2.3, "mm"),
    legend.key.width = unit(3.2, "mm"),
    legend.spacing.x = unit(1.0, "mm"),
    legend.key = element_blank()
  )

group_only_plot <- ggplot(
  group_only_summary,
  aes(x = Group, y = mean_prop, color = Group, fill = Group)
) +
  stat_halfeye(
    data = df_plot,
    aes(x = Group, y = percentage_fix_duration, fill = Group),
    adjust = 0.55,
    width = 0.55,
    .width = 0,
    point_color = NA,
    alpha = 0.38,
    slab_color = NA,
    normalize = "groups",
    scale = 0.40,
    justification = -0.28
  ) +
  geom_point(size = 2.1) +
  geom_segment(
    data = group_only_sig,
    aes(x = xstart, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.35
  ) +
  geom_segment(
    data = group_only_sig,
    aes(x = xstart, xend = xstart, y = ystem, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.35
  ) +
  geom_segment(
    data = group_only_sig,
    aes(x = xend, xend = xend, y = ystem, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.35
  ) +
  geom_text(
    data = group_only_sig,
    aes(x = x, y = y + 0.01, label = sig),
    inherit.aes = FALSE,
    size = 4.3
  ) +
  facet_grid2(Panel ~ fix, axes = "all", remove_labels = "all", switch = "y") +
  labs(
    title = "Overall Group Differences",
    x = "Group",
    y = "% Fixation Duration"
  ) +
  scale_color_manual(values = group_values) +
  scale_fill_manual(values = group_values) +
  base_duration_theme +
  theme(
    legend.position = "none",
    strip.placement = "outside"
  )

med_split_plot <- ggplot(
  med_split_summary,
  aes(x = Group, y = mean_prop, color = medication, fill = medication)
) +
  stat_halfeye(
    data = df_plot,
    aes(x = Group, y = percentage_fix_duration, fill = medication),
    adjust = 0.55,
    width = 0.55,
    .width = 0,
    point_color = NA,
    alpha = 0.45,
    slab_color = NA,
    normalize = "groups",
    scale = 0.40,
    position = position_dodge(width = 0.8),
    justification = -0.28
  ) +
  geom_point(
    position = position_dodge(width = 0.8),
    size = 1.9
  ) +
  geom_segment(
    data = med_split_sig,
    aes(x = xstart, xend = xend, y = y, yend = y, color = medication),
    inherit.aes = FALSE,
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  geom_segment(
    data = med_split_sig,
    aes(x = xstart, xend = xstart, y = ystem, yend = y, color = medication),
    inherit.aes = FALSE,
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  geom_segment(
    data = med_split_sig,
    aes(x = xend, xend = xend, y = ystem, yend = y, color = medication),
    inherit.aes = FALSE,
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  geom_text(
    data = med_split_sig,
    aes(x = x, y = y + 0.01, label = sig, color = medication),
    inherit.aes = FALSE,
    size = 4.0,
    show.legend = FALSE
  ) +
  facet_grid2(Panel ~ fix, axes = "all", remove_labels = "all", switch = "y") +
  labs(
    title = "Medication-Specific Differences",
    x = "Group",
    y = "% Fixation Duration"
  ) +
  scale_color_manual(values = medication_values) +
  scale_fill_manual(values = medication_values) +
  scale_x_discrete(expand = expansion(mult = c(0.18, 0.18))) +
  guides(
    color = "none",
    fill = guide_legend(
      override.aes = list(alpha = 1, shape = 15, size = 2.4, colour = NA, linewidth = 0),
      nrow = 1
    )
  ) +
  base_duration_theme +
  theme(strip.placement = "outside")

figure_duration <- group_only_plot / med_split_plot +
  plot_layout(heights = c(0.95, 1.05)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 12, face = "bold"),
      plot.tag.position = c(0.01, 0.995)
    )
  )

ggsave(file.path(output_svg, "Figure_2_fixation_duration.svg"), figure_duration, units = "mm", width = 180, height = 215)
ggsave(file.path(output_png, "Figure_2_fixation_duration.png"), figure_duration, units = "mm", width = 180, height = 215, dpi = 300)
