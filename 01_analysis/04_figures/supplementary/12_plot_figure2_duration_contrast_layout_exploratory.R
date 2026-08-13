# Manuscript section: Exploratory Figure 2 redesign (fixation duration)
# Analysis family: compact distribution panels plus main-model contrast panel

library(tidyverse)
library(ggdist)
library(ggh4x)
library(patchwork)
library(here)

output_dir <- here("02_outputs", "figures", "supplementary", "fixation_duration")
output_png <- file.path(output_dir, "Figure_2_fixation_duration_contrast_layout_exploratory.png")
output_svg <- file.path(output_dir, "figure2_fixation_duration_contrast_layout_exploratory.svg")
current_story_dir <- here("02_outputs", "figures", "main", "current_story")
current_story_png <- file.path(current_story_dir, "Figure_2_fixation_duration_contrast_layout.png")
current_story_svg <- file.path(current_story_dir, "Figure_2_fixation_duration_contrast_layout.svg")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(current_story_dir, recursive = TRUE, showWarnings = FALSE)

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

panel_levels <- c("Full Sample", "Female", "Male")
fix_levels <- c("Eyes", "Background")
med_levels <- c("BOTH", "NAL", "OXT", "PLA")
group_values <- c("ASD" = "#D55E00", "CTRL" = "#0072B2")
medication_values <- c("BOTH" = "#E69F00", "NAL" = "#CCB974", "OXT" = "#56B4A9", "PLA" = "#0B7A6B")
panel_values <- c("Full Sample" = "#303030", "Female" = "#B44E3A", "Male" = "#3B6C9E")

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
  filter(!is.na(Group), !is.na(medication), Sex %in% c("f", "m")) %>%
  mutate(
    fix = factor(
      fix,
      levels = c("fix_on_eyes", "fix_on_mouth", "fix_on_face", "fix_on_background"),
      labels = c("Eyes", "Mouth", "Face", "Background")
    ),
    medication = factor(
      medication,
      levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA"),
      labels = med_levels
    ),
    Group = factor(Group, levels = c("ASD", "CTRL"))
  ) %>%
  filter(fix %in% fix_levels)

df_plot <- bind_rows(
  df %>% mutate(Panel = "Full Sample"),
  df %>% filter(Sex == "f") %>% mutate(Panel = "Female"),
  df %>% filter(Sex == "m") %>% mutate(Panel = "Male")
) %>%
  mutate(
    Panel = factor(Panel, levels = panel_levels),
    fix = factor(as.character(fix), levels = fix_levels),
    medication = factor(medication, levels = med_levels)
  )

overall_positions <- tibble(
  Group = factor(c("ASD", "CTRL"), levels = c("ASD", "CTRL")),
  x = c(0.82, 1.18),
  side = c("left", "right")
)

df_overall <- df_plot %>%
  left_join(overall_positions, by = "Group")

df_overall_subject <- df_overall %>%
  group_by(Panel, fix, Group, x, ID) %>%
  summarise(
    percentage_fix_duration = mean(percentage_fix_duration, na.rm = TRUE),
    .groups = "drop"
  )

overall_summary <- df_overall_subject %>%
  group_by(Panel, fix, Group, x) %>%
  summarise(
    mean_prop = mean(percentage_fix_duration, na.rm = TRUE),
    sem = sd(percentage_fix_duration, na.rm = TRUE) / sqrt(sum(!is.na(percentage_fix_duration))),
    lower = pmax(0, mean_prop - 1.96 * sem),
    upper = pmin(1, mean_prop + 1.96 * sem),
    .groups = "drop"
  ) %>%
  mutate(fix = factor(as.character(fix), levels = fix_levels))

mean_box_summary <- df_overall_subject %>%
  group_by(Panel, fix, Group, x) %>%
  group_modify(~ {
    values <- .x$percentage_fix_duration
    values <- values[is.finite(values)]
    quartiles <- quantile(values, probs = c(0.25, 0.75), names = FALSE)
    iqr <- quartiles[[2]] - quartiles[[1]]

    tibble(
      mean_prop = mean(values),
      q1 = quartiles[[1]],
      q3 = quartiles[[2]],
      lower_whisker = min(values[values >= quartiles[[1]] - 1.5 * iqr]),
      upper_whisker = max(values[values <= quartiles[[2]] + 1.5 * iqr])
    )
  }) %>%
  ungroup() %>%
  mutate(fix = factor(as.character(fix), levels = fix_levels))

med_positions <- expand_grid(
  medication = factor(med_levels, levels = med_levels),
  Group = factor(c("ASD", "CTRL"), levels = c("ASD", "CTRL"))
) %>%
  mutate(
    med_x = as.numeric(medication),
    x = med_x + if_else(Group == "ASD", -0.13, 0.13),
    side = if_else(Group == "ASD", "left", "right")
  )

df_med <- df_plot %>%
  left_join(med_positions, by = c("medication", "Group"))

med_summary <- df_med %>%
  group_by(Panel, fix, medication, Group, x) %>%
  summarise(
    mean_prop = mean(percentage_fix_duration, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(fix = factor(as.character(fix), levels = fix_levels))

df_med_full_subject <- df_med %>%
  filter(Panel == "Full Sample") %>%
  group_by(ID, fix, medication, Group, x) %>%
  summarise(
    percentage_fix_duration = mean(percentage_fix_duration, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    fix = factor(as.character(fix), levels = fix_levels),
    medication = factor(medication, levels = med_levels)
  )

med_full_box_summary <- df_med_full_subject %>%
  group_by(fix, medication, Group, x) %>%
  group_modify(~ {
    values <- .x$percentage_fix_duration
    values <- values[is.finite(values)]
    quartiles <- quantile(values, probs = c(0.25, 0.75), names = FALSE)
    iqr <- quartiles[[2]] - quartiles[[1]]

    tibble(
      mean_prop = mean(values),
      q1 = quartiles[[1]],
      q3 = quartiles[[2]],
      lower_whisker = min(values[values >= quartiles[[1]] - 1.5 * iqr]),
      upper_whisker = max(values[values <= quartiles[[2]] + 1.5 * iqr])
    )
  }) %>%
  ungroup()

model_dir <- here("02_outputs", "model_outputs", "main_manuscript", "eye_tracking", "group_aoi_medication_duration")

overall_contrasts <- read_csv(file.path(model_dir, "03_group_by_fix.csv"), show_col_types = FALSE) %>%
  transmute(
    Panel = "Full Sample",
    fix,
    medication = "Overall",
    contrast,
    estimate,
    std.error,
    lower = estimate - 1.96 * std.error,
    upper = estimate + 1.96 * std.error,
    p.value
  )

sex_contrasts <- read_csv(file.path(model_dir, "04_group_by_fix_sex.csv"), show_col_types = FALSE) %>%
  transmute(
    Panel = recode(Sex, "f" = "Female", "m" = "Male"),
    fix,
    medication = "Overall",
    contrast,
    estimate,
    std.error,
    lower = estimate - 1.96 * std.error,
    upper = estimate + 1.96 * std.error,
    p.value
  )

med_overall_contrasts <- read_csv(file.path(model_dir, "05_group_by_fix_medication.csv"), show_col_types = FALSE) %>%
  transmute(
    Panel = "Full Sample",
    fix,
    medication = factor(
      medication,
      levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA"),
      labels = med_levels
    ) %>% as.character(),
    contrast,
    estimate,
    std.error,
    lower = estimate - 1.96 * std.error,
    upper = estimate + 1.96 * std.error,
    p.value
  )

med_sex_contrasts <- read_csv(file.path(model_dir, "06_group_by_fix_medication_sex.csv"), show_col_types = FALSE) %>%
  transmute(
    Panel = recode(Sex, "f" = "Female", "m" = "Male"),
    fix,
    medication = factor(
      medication,
      levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA"),
      labels = med_levels
    ) %>% as.character(),
    contrast,
    estimate,
    std.error,
    lower = estimate - 1.96 * std.error,
    upper = estimate + 1.96 * std.error,
    p.value
  )

contrast_df <- bind_rows(
  overall_contrasts,
  sex_contrasts,
  med_overall_contrasts,
  med_sex_contrasts
) %>%
  filter(fix %in% fix_levels) %>%
  mutate(
    Panel = factor(Panel, levels = panel_levels),
    fix = factor(fix, levels = fix_levels),
    medication = factor(medication, levels = c("Overall", med_levels)),
    medication_x = as.numeric(medication),
    significant = p.value < 0.05,
    sig = format_sig(p.value)
  )

panel_a_sig <- contrast_df %>%
  filter(medication == "Overall", significant) %>%
  left_join(
    overall_summary %>%
      group_by(Panel, fix) %>%
      summarise(y = max(upper, na.rm = TRUE) + 0.10, .groups = "drop"),
    by = c("Panel", "fix")
  ) %>%
  mutate(x = 1)

panel_dodge <- position_dodge(width = 0.55)

base_theme <- theme_classic(base_size = 8) +
  theme(
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(color = "black", size = 8.5),
    plot.title = element_text(size = 10, face = "bold", hjust = 0),
    panel.spacing.y = unit(5, "mm"),
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

mean_bar_width <- 0.075
mean_box_half_width <- 0.030
med_mean_box_half_width <- 0.055

panel_a <- ggplot() +
  stat_halfeye(
    data = df_overall_subject %>% filter(Group == "ASD"),
    aes(x = x, y = percentage_fix_duration, fill = Group),
    side = "left",
    adjust = 0.55,
    width = 0.62,
    .width = 0,
    point_color = NA,
    alpha = 0.46,
    slab_color = NA,
    normalize = "panels",
    scale = 0.62,
    show.legend = FALSE
  ) +
  stat_halfeye(
    data = df_overall_subject %>% filter(Group == "CTRL"),
    aes(x = x, y = percentage_fix_duration, fill = Group),
    side = "right",
    adjust = 0.55,
    width = 0.62,
    .width = 0,
    point_color = NA,
    alpha = 0.46,
    slab_color = NA,
    normalize = "panels",
    scale = 0.62,
    show.legend = FALSE
  ) +
  geom_errorbar(
    data = mean_box_summary,
    aes(
      x = x,
      ymin = lower_whisker,
      ymax = upper_whisker,
      color = Group
    ),
    width = 0.025,
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  geom_rect(
    data = mean_box_summary,
    aes(
      xmin = x - mean_box_half_width,
      xmax = x + mean_box_half_width,
      ymin = q1,
      ymax = q3,
      color = Group
    ),
    fill = "white",
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  geom_segment(
    data = mean_box_summary,
    aes(
      x = x - mean_box_half_width,
      xend = x + mean_box_half_width,
      y = mean_prop,
      yend = mean_prop,
      color = Group
    ),
    linewidth = 0.8,
    lineend = "round",
    show.legend = FALSE
  ) +
  geom_text(
    data = panel_a_sig,
    aes(x = x, y = y, label = sig),
    inherit.aes = FALSE,
    size = 3.8,
    fontface = "bold"
  ) +
  facet_grid2(Panel ~ fix, axes = "x", remove_labels = "x", switch = "y") +
  scale_x_continuous(
    breaks = c(0.82, 1.18),
    labels = c("ASD", "CTRL"),
    expand = expansion(mult = c(0.25, 0.25))
  ) +
  scale_color_manual(values = group_values) +
  scale_fill_manual(values = group_values) +
  coord_cartesian(xlim = c(0.42, 1.58), ylim = c(0, 1.1)) +
  labs(title = "Overall Group Differences", x = NULL, y = "% Fixation Duration") +
  base_theme +
  theme(
    legend.position = "none",
    strip.placement = "outside",
    axis.text.x = element_text(size = 7, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.25)
  ) +
  guides(
    color = "none",
    fill = "none"
  )

panel_c <- ggplot(
  contrast_df,
  aes(x = medication_x, y = estimate, color = Panel, shape = significant)
) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.25, linetype = "22") +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = panel_dodge,
    width = 0.14,
    linewidth = 0.4
  ) +
  geom_point(
    position = panel_dodge,
    size = 1.8,
    stroke = 0.7,
    fill = "white"
  ) +
  facet_grid(. ~ fix) +
  scale_x_continuous(
    breaks = seq_along(levels(contrast_df$medication)),
    labels = levels(contrast_df$medication),
    expand = expansion(mult = c(0.07, 0.07))
  ) +
  scale_color_manual(values = panel_values) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), guide = "none") +
  labs(
    title = "Group Contrasts (Exploratory)",
    x = NULL,
    y = "ASD - CTRL Contrast"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "top"
  ) +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 2.2)))

panel_b <- ggplot() +
  stat_halfeye(
    data = df_med %>% filter(Group == "ASD"),
    aes(x = x, y = percentage_fix_duration, fill = medication),
    side = "left",
    adjust = 0.55,
    width = 0.42,
    .width = 0,
    point_color = NA,
    alpha = 0.42,
    slab_color = NA,
    normalize = "panels",
    scale = 0.52
  ) +
  stat_halfeye(
    data = df_med %>% filter(Group == "CTRL"),
    aes(x = x, y = percentage_fix_duration, fill = medication),
    side = "right",
    adjust = 0.55,
    width = 0.42,
    .width = 0,
    point_color = NA,
    alpha = 0.42,
    slab_color = NA,
    normalize = "panels",
    scale = 0.52
  ) +
  geom_segment(
    data = med_summary,
    aes(
      x = x - mean_bar_width,
      xend = x + mean_bar_width,
      y = mean_prop,
      yend = mean_prop,
      color = medication
    ),
    linewidth = 0.85,
    lineend = "round"
  ) +
  facet_grid2(Panel ~ fix, axes = "x", remove_labels = "x", switch = "y") +
  scale_x_continuous(
    breaks = rep(seq_along(med_levels), each = 2) + rep(c(-0.13, 0.13), times = length(med_levels)),
    labels = rep(c("ASD", "CTRL"), times = length(med_levels)),
    minor_breaks = seq_along(med_levels),
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  scale_color_manual(values = medication_values) +
  scale_fill_manual(values = medication_values) +
  coord_cartesian(ylim = c(0, 1.1)) +
  labs(title = "Medication-Specific Differences (Exploratory)", x = NULL, y = "% Fixation Duration") +
  base_theme +
  theme(
    strip.placement = "outside",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_line(linewidth = 0.25),
    legend.position = "top"
  ) +
  guides(
    color = "none",
    fill = guide_legend(
      override.aes = list(alpha = 1, shape = 15, size = 2.4, colour = NA, linewidth = 0),
      nrow = 1
    )
  )

panel_c_descriptive <- ggplot() +
  geom_errorbar(
    data = med_full_box_summary,
    aes(
      x = x,
      ymin = lower_whisker,
      ymax = upper_whisker,
      color = Group
    ),
    width = 0.035,
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  geom_rect(
    data = med_full_box_summary,
    aes(
      xmin = x - med_mean_box_half_width,
      xmax = x + med_mean_box_half_width,
      ymin = q1,
      ymax = q3,
      color = Group
    ),
    fill = "white",
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  geom_segment(
    data = med_full_box_summary,
    aes(
      x = x - med_mean_box_half_width,
      xend = x + med_mean_box_half_width,
      y = mean_prop,
      yend = mean_prop,
      color = Group
    ),
    linewidth = 0.8,
    lineend = "round",
    show.legend = TRUE
  ) +
  facet_grid(. ~ fix) +
  scale_x_continuous(
    breaks = seq_along(med_levels),
    labels = med_levels,
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  scale_color_manual(values = group_values) +
  scale_fill_manual(values = group_values) +
  coord_cartesian(ylim = c(0, 1.1)) +
  labs(
    title = "Medication-Stratified Descriptives (Full Sample)",
    x = "Medication Condition",
    y = "% Fixation Duration",
    tag = "C"
  ) +
  base_theme +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.spacing = unit(0.5, "mm"),
    plot.tag = element_text(size = 11, face = "plain"),
    plot.tag.position = c(0.005, 0.995),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.25)
  ) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1.2), nrow = 1),
    fill = "none",
    size = "none"
  )

figure3_env <- new.env(parent = globalenv())
sys.source(
  here("01_analysis", "04_figures", "supplementary", "11_plot_nonoverlap_duration_timecourse_exploratory.R"),
  envir = figure3_env
)

figure3_panel <- figure3_env$top_plot +
  labs(title = "Timecourse Analysis", y = NULL, tag = "B") +
  theme(
    strip.text = element_text(color = "black", size = 8.5),
    strip.text.y.left = element_blank(),
    plot.title = element_text(size = 10, face = "bold", hjust = 0, margin = margin(b = 0, unit = "mm")),
    plot.tag = element_text(size = 11, face = "plain"),
    plot.tag.position = c(0.005, 0.995),
    panel.spacing.y = unit(5, "mm"),
    panel.spacing.x = unit(5, "mm"),
    legend.text = element_text(size = 6.5),
    legend.key.height = unit(1.8, "mm"),
    legend.key.width = unit(3.2, "mm"),
    legend.spacing.x = unit(1.0, "mm"),
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.spacing = unit(0.5, "mm"),
    strip.placement = "outside",
    legend.position = "top"
  ) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1.1), nrow = 1),
    fill = "none"
  )

panel_a <- panel_a +
  labs(tag = "A", subtitle = " ") +
  theme(
    plot.subtitle = element_text(size = 6.5, lineheight = 0.8, margin = margin(t = 0.25, b = 0.25, unit = "mm")),
    plot.tag = element_text(size = 11, face = "plain"),
    plot.tag.position = c(0.005, 0.995)
  )

figure_body <- (panel_a | figure3_panel) +
  plot_layout(widths = c(0.44, 0.56))

figure_exploratory <- figure_body / panel_c_descriptive +
  plot_layout(heights = c(2.15, 0.8))

ggsave(output_svg, figure_exploratory, units = "mm", width = 183, height = 170)
ggsave(output_png, figure_exploratory, units = "mm", width = 183, height = 170, dpi = 300)
ggsave(current_story_svg, figure_exploratory, units = "mm", width = 183, height = 170)
ggsave(current_story_png, figure_exploratory, units = "mm", width = 183, height = 170, dpi = 300)
