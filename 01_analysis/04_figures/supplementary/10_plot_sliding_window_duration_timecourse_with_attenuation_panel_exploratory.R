# Manuscript section: Exploratory duration-based longitudinal analysis
# Analysis family: reduced combined timecourse figure

library(tidyverse)
library(ggh4x)
library(patchwork)
library(glmmTMB)
library(emmeans)
library(here)

input_path <- here("00_data", "derived", "preprocessing", "duration_sliding_window_data.csv")
output_png <- here(
  "02_outputs", "figures", "supplementary", "timecourse", "png",
  "Supp_Figure_duration_sliding_window_with_attenuation_panel_exploratory.png"
)
output_svg <- here(
  "02_outputs", "figures", "supplementary", "timecourse", "svg",
  "duration_sliding_window_with_attenuation_panel_exploratory.svg"
)
panel_levels <- c("Full Sample", "Female")
fix_levels <- c("Eyes", "Background")
group_values <- c("ASD" = "#D55E00", "CTRL" = "#0072B2")
interval_colors <- c("ASD_gt_CTRL" = group_values[["ASD"]], "ASD_lt_CTRL" = group_values[["CTRL"]])
plot_x_min <- 0.05
plot_x_max <- 0.95
plot_x_breaks <- seq(plot_x_min, plot_x_max, length.out = 5)
plot_x_labels <- c("0%", "25%", "50%", "75%", "100%")

beta_squeeze <- function(x) {
  n <- sum(!is.na(x))
  ((x * (n - 1)) + 0.5) / n
}

build_significance_intervals <- function(contrast_df, group_vars) {
  effect_col <- intersect(c("estimate", "odds.ratio", "ratio"), names(contrast_df))
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
      .groups = "drop"
    ) %>%
    filter(is.finite(interval_start), is.finite(interval_end))
}

df <- read_csv(input_path, show_col_types = FALSE) %>%
  filter(!is.na(window_fix_duration_prop), !is.na(medication), Sex %in% c("f", "m")) %>%
  mutate(
    ID = factor(Sub_ID),
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
    subject_session = interaction(ID, session, drop = TRUE),
    progress_window = factor(
      progress_window,
      levels = unique(progress_window[order(progress_window_index)]),
      ordered = TRUE
    ),
    window_fix_duration_prop_beta = beta_squeeze(window_fix_duration_prop)
  ) %>%
  filter(fix %in% fix_levels)

prepare_top_stats <- function(data, panel_label) {
  data %>%
    group_by(fix, Group, progress_window_index, progress_midpoint, progress_start, progress_end) %>%
    summarise(
      mean_window_fix_duration_prop = mean(window_fix_duration_prop, na.rm = TRUE),
      sem = sd(window_fix_duration_prop, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(Panel = factor(panel_label, levels = panel_levels))
}

prepare_placebo_stats <- function(data, panel_label) {
  data %>%
    filter(medication == "PLA") %>%
    group_by(fix, Group, progress_window_index, progress_midpoint, progress_start, progress_end) %>%
    summarise(
      mean_window_fix_duration_prop = mean(window_fix_duration_prop, na.rm = TRUE),
      sem = sd(window_fix_duration_prop, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(Panel = factor(panel_label, levels = panel_levels))
}

fit_top_model <- function(data, fix_label, panel_label) {
  df_fix <- data %>% filter(fix == fix_label)
  model_formula <- if (panel_label == "Full Sample") {
    window_fix_duration_prop_beta ~ Group * progress_window + Sex +
      (1 | ID) + (1 | subject_session)
  } else {
    window_fix_duration_prop_beta ~ Group * progress_window +
      (1 | ID) + (1 | subject_session)
  }
  model <- glmmTMB(
    model_formula,
    family = beta_family(),
    data = df_fix,
    control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
  )
  list(model = model, data = df_fix)
}

fit_medication_model <- function(data, fix_label, panel_label) {
  df_fix <- if (panel_label == "Female") {
    data %>% filter(fix == fix_label, Sex == "f")
  } else {
    data %>% filter(fix == fix_label)
  }
  model_formula <- if (panel_label == "Female") {
    window_fix_duration_prop_beta ~ Group * progress_window * medication + session +
      (1 | ID) + (1 | subject_session)
  } else {
    window_fix_duration_prop_beta ~ Group * progress_window * medication + Sex + session +
      (1 | ID) + (1 | subject_session)
  }
  model <- glmmTMB(
    model_formula,
    family = beta_family(),
    data = df_fix,
    control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
  )
  list(model = model, data = df_fix)
}

build_top_panel_intervals <- function(model, df_fix, fix_label, panel_label) {
  emmeans(model, pairwise ~ Group | progress_window, type = "response")$contrasts %>%
    summary() %>%
    as_tibble() %>%
    mutate(progress_window = as.character(progress_window)) %>%
    left_join(
      df_fix %>% distinct(progress_window, progress_window_index, progress_start, progress_end),
      by = "progress_window"
    ) %>%
    filter(!is.na(progress_start), !is.na(progress_end), !is.na(progress_window_index)) %>%
    build_significance_intervals(group_vars = character()) %>%
    mutate(
      Panel = factor(panel_label, levels = panel_levels),
      fix = factor(fix_label, levels = fix_levels)
    )
}

build_placebo_intervals <- function(model, df_fix, fix_label, panel_label) {
  emmeans(model, pairwise ~ Group | progress_window * medication, type = "response")$contrasts %>%
    summary() %>%
    as_tibble() %>%
    mutate(progress_window = as.character(progress_window)) %>%
    filter(medication == "PLA") %>%
    left_join(
      df_fix %>% distinct(progress_window, progress_window_index, progress_start, progress_end),
      by = "progress_window"
    ) %>%
    filter(!is.na(progress_start), !is.na(progress_end), !is.na(progress_window_index)) %>%
    build_significance_intervals(group_vars = character()) %>%
    mutate(
      Panel = factor(panel_label, levels = panel_levels),
      fix = factor(fix_label, levels = fix_levels)
    )
}

overall_top <- prepare_top_stats(df, "Full Sample")
female_top <- prepare_top_stats(filter(df, Sex == "f"), "Female")
top_stats <- bind_rows(overall_top, female_top) %>% mutate(fix = factor(fix, levels = fix_levels))

overall_placebo <- prepare_placebo_stats(df, "Full Sample")
female_placebo <- prepare_placebo_stats(filter(df, Sex == "f"), "Female")
placebo_stats <- bind_rows(overall_placebo, female_placebo) %>% mutate(fix = factor(fix, levels = fix_levels))

top_eyes_overall <- fit_top_model(df, "Eyes", "Full Sample")
top_eyes_female <- fit_top_model(filter(df, Sex == "f"), "Eyes", "Female")
top_bg_overall <- fit_top_model(df, "Background", "Full Sample")
top_bg_female <- fit_top_model(filter(df, Sex == "f"), "Background", "Female")

top_intervals <- bind_rows(
  build_top_panel_intervals(top_eyes_overall$model, top_eyes_overall$data, "Eyes", "Full Sample"),
  build_top_panel_intervals(top_eyes_female$model, top_eyes_female$data, "Eyes", "Female"),
  build_top_panel_intervals(top_bg_overall$model, top_bg_overall$data, "Background", "Full Sample"),
  build_top_panel_intervals(top_bg_female$model, top_bg_female$data, "Background", "Female")
)

top_interval_positions <- top_stats %>%
  group_by(Panel, fix) %>%
  summarise(interval_y = max(mean_window_fix_duration_prop + sem, na.rm = TRUE) + 0.03, .groups = "drop")

top_intervals <- top_intervals %>%
  mutate(
    interval_start = pmax(interval_start, plot_x_min),
    interval_end = pmin(interval_end, plot_x_max)
  ) %>%
  filter(interval_start < interval_end) %>%
  left_join(top_interval_positions, by = c("Panel", "fix"))

pla_eyes_overall <- fit_medication_model(df, "Eyes", "Full Sample")
pla_eyes_female <- fit_medication_model(df, "Eyes", "Female")
pla_bg_overall <- fit_medication_model(df, "Background", "Full Sample")
pla_bg_female <- fit_medication_model(df, "Background", "Female")

placebo_intervals <- bind_rows(
  build_placebo_intervals(pla_eyes_overall$model, pla_eyes_overall$data, "Eyes", "Full Sample"),
  build_placebo_intervals(pla_eyes_female$model, pla_eyes_female$data, "Eyes", "Female"),
  build_placebo_intervals(pla_bg_overall$model, pla_bg_overall$data, "Background", "Full Sample"),
  build_placebo_intervals(pla_bg_female$model, pla_bg_female$data, "Background", "Female")
)

placebo_interval_positions <- placebo_stats %>%
  group_by(Panel, fix) %>%
  summarise(interval_y = max(mean_window_fix_duration_prop + sem, na.rm = TRUE) + 0.03, .groups = "drop")

placebo_intervals <- placebo_intervals %>%
  mutate(
    interval_start = pmax(interval_start, plot_x_min),
    interval_end = pmin(interval_end, plot_x_max)
  ) %>%
  filter(interval_start < interval_end) %>%
  left_join(placebo_interval_positions, by = c("Panel", "fix"))

base_timecourse_theme <- theme_classic(base_size = 8) +
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

make_base_plot <- function(plot_df, title_text, show_legend = TRUE) {
  ggplot(
    plot_df,
    aes(x = progress_midpoint, y = mean_window_fix_duration_prop, group = Group, color = Group)
  ) +
    geom_line(linewidth = 0.4) +
    geom_ribbon(
      aes(ymin = mean_window_fix_duration_prop - sem, ymax = mean_window_fix_duration_prop + sem, fill = Group),
      alpha = 0.45,
      linewidth = 0.15
    ) +
    scale_x_continuous(
      limits = c(plot_x_min, plot_x_max),
      breaks = plot_x_breaks,
      labels = plot_x_labels,
      expand = expansion(mult = c(0, 0))
    ) +
    labs(x = "Normalized Interview Progress", y = "% Fixation Duration", title = title_text) +
    scale_color_manual(values = group_values) +
    scale_fill_manual(values = group_values) +
    base_timecourse_theme +
    theme(legend.position = if (show_legend) "top" else "none") +
    facet_grid2(Panel ~ fix, axes = "all", remove_labels = "all", switch = "y") +
    guides(
      color = guide_legend(override.aes = list(fill = group_values, color = group_values, alpha = 1, shape = 15, size = 3.2)),
      fill = "none"
    )
}

top_plot <- make_base_plot(top_stats, "Overall Group Differences", show_legend = TRUE) +
  geom_segment(
    data = top_intervals %>% filter(direction == "ASD_gt_CTRL"),
    aes(x = interval_start, xend = interval_end, y = interval_y, yend = interval_y),
    inherit.aes = FALSE,
    color = interval_colors[["ASD_gt_CTRL"]],
    linewidth = 2.1,
    lineend = "round"
  ) +
  geom_segment(
    data = top_intervals %>% filter(direction == "ASD_lt_CTRL"),
    aes(x = interval_start, xend = interval_end, y = interval_y, yend = interval_y),
    inherit.aes = FALSE,
    color = interval_colors[["ASD_lt_CTRL"]],
    linewidth = 2.1,
    lineend = "round"
  )

placebo_plot <- make_base_plot(placebo_stats, "Group Differences in the Placebo Condition", show_legend = FALSE) +
  geom_segment(
    data = placebo_intervals %>% filter(direction == "ASD_gt_CTRL"),
    aes(x = interval_start, xend = interval_end, y = interval_y, yend = interval_y),
    inherit.aes = FALSE,
    color = interval_colors[["ASD_gt_CTRL"]],
    linewidth = 2.1,
    lineend = "round"
  ) +
  geom_segment(
    data = placebo_intervals %>% filter(direction == "ASD_lt_CTRL"),
    aes(x = interval_start, xend = interval_end, y = interval_y, yend = interval_y),
    inherit.aes = FALSE,
    color = interval_colors[["ASD_lt_CTRL"]],
    linewidth = 2.1,
    lineend = "round"
  )

combined_plot <- top_plot / placebo_plot +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Timecourse Analysis (Exploratory)",
    tag_levels = "A",
    theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
  ) &
  theme(
    plot.tag = element_text(size = 10, face = "plain"),
    plot.tag.position = c(0, 1)
  )

ggsave(output_svg, combined_plot, units = "mm", width = 183, height = 155)
ggsave(output_png, combined_plot, units = "mm", width = 183, height = 155)
