# Manuscript section: Exploratory duration-based longitudinal analysis
# Analysis family: non-overlapping time-bin visualization with matched inference strips

library(tidyverse)
library(ggh4x)
library(patchwork)
library(glmmTMB)
library(emmeans)
library(here)

input_path <- here("00_data", "derived", "preprocessing", "duration_sliding_window_data.csv")
output_png <- here(
  "02_outputs", "figures", "supplementary", "timecourse", "png",
  "Supp_Figure_duration_nonoverlap_timecourse_exploratory.png"
)
output_svg <- here(
  "02_outputs", "figures", "supplementary", "timecourse", "svg",
  "duration_nonoverlap_timecourse_exploratory.svg"
)
current_story_dir <- here("02_outputs", "figures", "main", "current_story")
current_story_png <- file.path(current_story_dir, "Figure_3_timecourse_nonoverlap.png")
current_story_svg <- file.path(current_story_dir, "Figure_3_timecourse_nonoverlap.svg")
current_story_inset_png <- file.path(current_story_dir, "Figure_3_timecourse_nonoverlap_inset.png")
contrast_output <- here(
  "02_outputs", "model_outputs", "main_manuscript", "eye_tracking",
  "timecourse_duration_nonoverlap_random_structure", "05_nonoverlap_exploratory_window_contrasts.csv"
)

dir.create(dirname(output_png), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_svg), recursive = TRUE, showWarnings = FALSE)
dir.create(current_story_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(contrast_output), recursive = TRUE, showWarnings = FALSE)

panel_levels <- c("Full Sample", "Female", "Male")
fix_levels <- c("Eyes", "Background")
group_values <- c("ASD" = "#D55E00", "CTRL" = "#0072B2")
interval_colors <- c("ASD_gt_CTRL" = group_values[["ASD"]], "ASD_lt_CTRL" = group_values[["CTRL"]])

beta_squeeze <- function(x) {
  n <- sum(!is.na(x))
  ((x * (n - 1)) + 0.5) / n
}

build_significance_intervals <- function(contrast_df, group_vars) {
  effect_col <- intersect(c("estimate", "odds.ratio", "ratio"), names(contrast_df))
  effect_col <- effect_col[[1]]

  sig_df <- contrast_df %>%
    filter(!is.na(p.value_fdr), p.value_fdr < 0.05)

  if (nrow(sig_df) == 0) {
    return(tibble(
      direction = character(),
      run_id = integer(),
      interval_start = numeric(),
      interval_end = numeric()
    ))
  }

  sig_df %>%
    arrange(across(all_of(c(group_vars, "progress_bin_index")))) %>%
    mutate(
      direction = case_when(
        str_detect(contrast, "^CTRL\\s*/\\s*ASD|^CTRL\\s*-\\s*ASD") & .data[[effect_col]] > 1 ~ "ASD_lt_CTRL",
        str_detect(contrast, "^CTRL\\s*/\\s*ASD|^CTRL\\s*-\\s*ASD") & .data[[effect_col]] < 1 ~ "ASD_gt_CTRL",
        str_detect(contrast, "^ASD\\s*/\\s*CTRL|^ASD\\s*-\\s*CTRL") & .data[[effect_col]] > 1 ~ "ASD_gt_CTRL",
        str_detect(contrast, "^ASD\\s*/\\s*CTRL|^ASD\\s*-\\s*CTRL") & .data[[effect_col]] < 1 ~ "ASD_lt_CTRL",
        TRUE ~ "no_difference"
      )
    ) %>%
    filter(direction != "no_difference") %>%
    group_by(across(all_of(group_vars)), direction) %>%
    mutate(run_id = cumsum(c(TRUE, diff(progress_bin_index) != 1))) %>%
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
  filter(abs((progress_start * 10) - round(progress_start * 10)) < 1e-8) %>%
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
    progress_bin = factor(
      sprintf("%02d", as.integer(round(progress_start * 10)) + 1),
      levels = sprintf("%02d", 1:10),
      ordered = TRUE
    ),
    progress_bin_index = as.integer(progress_bin),
    progress_midpoint = progress_start + 0.05,
    subject_session = interaction(ID, session, drop = TRUE),
    window_fix_duration_prop_beta = beta_squeeze(window_fix_duration_prop)
  ) %>%
  filter(fix %in% fix_levels)

prepare_timecourse_stats <- function(data, panel_label, placebo_only = FALSE) {
  plot_data <- if (placebo_only) {
    data %>% filter(medication == "PLA")
  } else {
    data
  }

  plot_data %>%
    group_by(fix, Group, progress_bin, progress_bin_index, progress_midpoint, progress_start, progress_end) %>%
    summarise(
      mean_window_fix_duration_prop = mean(window_fix_duration_prop, na.rm = TRUE),
      sem = sd(window_fix_duration_prop, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      Panel = factor(panel_label, levels = panel_levels),
      fix = factor(fix, levels = fix_levels)
    )
}

fit_preferred_model <- function(data, fix_label) {
  df_fix <- data %>% filter(fix == fix_label)

  if (fix_label == "Eyes") {
    model_formula <- window_fix_duration_prop_beta ~
      Group * progress_bin * medication + Group * Sex * progress_bin + session +
      (1 | ID) + (1 | subject_session)
  } else {
    model_formula <- window_fix_duration_prop_beta ~
      Group * progress_bin * medication + Group * Sex * progress_bin + session +
      (1 | ID) + diag(0 + medication | ID)
  }

  glmmTMB(
    model_formula,
    family = beta_family(),
    data = df_fix,
    control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
  )
}

build_window_contrasts <- function(model, source_data, fix_label, panel_label, figure_panel, at_values = list()) {
  emmeans(
    model,
    ~ Group | progress_bin,
    at = at_values,
    type = "response"
  ) %>%
    contrast(method = "revpairwise") %>%
    summary() %>%
    as_tibble() %>%
    mutate(progress_bin = as.character(progress_bin)) %>%
    left_join(
      source_data %>%
        distinct(progress_bin, progress_bin_index, progress_start, progress_end),
      by = "progress_bin"
    ) %>%
    filter(!is.na(progress_start), !is.na(progress_end), !is.na(progress_bin_index)) %>%
    mutate(
      Figure_Panel = figure_panel,
      Panel = factor(panel_label, levels = panel_levels),
      fix = factor(fix_label, levels = fix_levels)
    )
}

top_stats <- bind_rows(
  prepare_timecourse_stats(df, "Full Sample", placebo_only = FALSE),
  prepare_timecourse_stats(filter(df, Sex == "f"), "Female", placebo_only = FALSE),
  prepare_timecourse_stats(filter(df, Sex == "m"), "Male", placebo_only = FALSE)
)

overall_models <- list(
  eyes = fit_preferred_model(df, "Eyes"),
  background = fit_preferred_model(df, "Background")
)

window_contrasts <- bind_rows(
  build_window_contrasts(overall_models$eyes, df %>% filter(fix == "Eyes"), "Eyes", "Full Sample", "A"),
  build_window_contrasts(overall_models$eyes, df %>% filter(fix == "Eyes"), "Eyes", "Female", "A", at_values = list(Sex = "f")),
  build_window_contrasts(overall_models$eyes, df %>% filter(fix == "Eyes"), "Eyes", "Male", "A", at_values = list(Sex = "m")),
  build_window_contrasts(overall_models$background, df %>% filter(fix == "Background"), "Background", "Full Sample", "A"),
  build_window_contrasts(overall_models$background, df %>% filter(fix == "Background"), "Background", "Female", "A", at_values = list(Sex = "f")),
  build_window_contrasts(overall_models$background, df %>% filter(fix == "Background"), "Background", "Male", "A", at_values = list(Sex = "m"))
) %>%
  group_by(Figure_Panel, Panel, fix) %>%
  mutate(
    p.value_fdr = p.adjust(p.value, method = "BH"),
    significant_fdr = !is.na(p.value_fdr) & p.value_fdr < 0.05
  ) %>%
  ungroup()

write_csv(window_contrasts, contrast_output)

intervals <- window_contrasts %>%
  group_by(Figure_Panel, Panel, fix) %>%
  group_modify(~ build_significance_intervals(.x, group_vars = character())) %>%
  ungroup()

interval_positions <- bind_rows(
  top_stats %>% mutate(Figure_Panel = "A")
) %>%
  group_by(Figure_Panel, Panel, fix) %>%
  summarise(interval_y = max(mean_window_fix_duration_prop + sem, na.rm = TRUE) + 0.03, .groups = "drop")

intervals <- intervals %>%
  mutate(
    interval_start = pmax(interval_start, 0.05),
    interval_end = pmin(interval_end, 0.95)
  ) %>%
  filter(interval_start < interval_end) %>%
  left_join(interval_positions, by = c("Figure_Panel", "Panel", "fix"))

base_theme <- theme_classic(base_size = 8) +
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

make_timecourse_plot <- function(plot_df, figure_panel, title_text, show_legend = TRUE) {
  panel_intervals <- intervals %>% filter(Figure_Panel == figure_panel)

  ggplot(
    plot_df,
    aes(x = progress_midpoint, y = mean_window_fix_duration_prop, group = Group, color = Group)
  ) +
    geom_ribbon(
      aes(ymin = mean_window_fix_duration_prop - sem, ymax = mean_window_fix_duration_prop + sem, fill = Group),
      alpha = 0.35,
      linewidth = 0.15,
      show.legend = FALSE
    ) +
    geom_line(linewidth = 0.45) +
    geom_segment(
      data = panel_intervals %>% filter(direction == "ASD_gt_CTRL"),
      aes(x = interval_start, xend = interval_end, y = interval_y, yend = interval_y),
      inherit.aes = FALSE,
      color = interval_colors[["ASD_gt_CTRL"]],
      linewidth = 1.6,
      lineend = "round"
    ) +
    geom_segment(
      data = panel_intervals %>% filter(direction == "ASD_lt_CTRL"),
      aes(x = interval_start, xend = interval_end, y = interval_y, yend = interval_y),
      inherit.aes = FALSE,
      color = interval_colors[["ASD_lt_CTRL"]],
      linewidth = 1.6,
      lineend = "round"
    ) +
    scale_x_continuous(
      limits = c(0.05, 0.95),
      breaks = seq(0.05, 0.95, length.out = 5),
      labels = c("0%", "25%", "50%", "75%", "100%"),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_color_manual(values = group_values) +
    scale_fill_manual(values = group_values) +
    labs(
      x = "Normalized Interview Progress",
      y = "% Fixation Duration",
      title = title_text
    ) +
    base_theme +
    theme(
      legend.position = if (show_legend) "top" else "none",
      strip.placement = "outside"
    ) +
    facet_grid2(Panel ~ fix, axes = "x", remove_labels = "x", switch = "y") +
    guides(
      color = guide_legend(override.aes = list(linewidth = 1.1)),
      fill = "none"
    )
}

top_plot <- make_timecourse_plot(
  top_stats,
  figure_panel = "A",
  title_text = "Overall Group Differences",
  show_legend = TRUE
)

combined_plot <- top_plot +
  plot_annotation(
    title = "Timecourse Analysis",
    tag_levels = "A",
    theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
  ) &
  theme(
    plot.tag = element_text(size = 10, face = "plain"),
    plot.tag.position = c(0, 1)
  )

inset_plot <- top_plot +
  labs(title = "Timecourse Analysis") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0))

ggsave(output_svg, combined_plot, units = "mm", width = 183, height = 150)
ggsave(output_png, combined_plot, units = "mm", width = 183, height = 150, dpi = 300)
ggsave(current_story_svg, combined_plot, units = "mm", width = 183, height = 150)
ggsave(current_story_png, combined_plot, units = "mm", width = 183, height = 150, dpi = 300)
ggsave(current_story_inset_png, inset_plot, units = "mm", width = 183, height = 150, dpi = 300)
