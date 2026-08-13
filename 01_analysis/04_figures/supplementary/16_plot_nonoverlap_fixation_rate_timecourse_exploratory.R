library(tidyverse)
library(ggh4x)
library(patchwork)
library(lme4)
library(lmerTest)
library(emmeans)
library(here)

input_path <- here("00_data", "derived", "preprocessing", "fixation_rate_sliding_window_data.csv")
output_png <- here(
  "02_outputs", "figures", "supplementary", "fixation_rate_analogs_current",
  "Figure_S_fixrate_3_current_layout.png"
)
dir.create(dirname(output_png), recursive = TRUE, showWarnings = FALSE)

panel_levels <- c("Overall", "Female")
fix_levels <- c("Eyes", "Background")
group_values <- c("ASD" = "#D55E00", "CTRL" = "#0072B2")
interval_colors <- c("ASD_gt_CTRL" = group_values[["ASD"]], "ASD_lt_CTRL" = group_values[["CTRL"]])

build_significance_intervals <- function(contrast_df) {
  sig_df <- contrast_df %>% filter(!is.na(p.value), p.value < 0.05)
  if (nrow(sig_df) == 0) {
    return(tibble(direction = character(), run_id = integer(), interval_start = numeric(), interval_end = numeric()))
  }

  sig_df %>%
    arrange(progress_bin_index) %>%
    mutate(
      direction = case_when(
        str_detect(contrast, "^CTRL\\s*/\\s*ASD|^CTRL\\s*-\\s*ASD") & estimate > 0 ~ "ASD_lt_CTRL",
        str_detect(contrast, "^CTRL\\s*/\\s*ASD|^CTRL\\s*-\\s*ASD") & estimate < 0 ~ "ASD_gt_CTRL",
        str_detect(contrast, "^ASD\\s*/\\s*CTRL|^ASD\\s*-\\s*CTRL") & estimate > 0 ~ "ASD_gt_CTRL",
        str_detect(contrast, "^ASD\\s*/\\s*CTRL|^ASD\\s*-\\s*CTRL") & estimate < 0 ~ "ASD_lt_CTRL",
        TRUE ~ "no_difference"
      )
    ) %>%
    filter(direction != "no_difference") %>%
    group_by(direction) %>%
    mutate(run_id = cumsum(c(TRUE, diff(progress_bin_index) != 1))) %>%
    group_by(direction, run_id) %>%
    summarise(
      interval_start = min(progress_start, na.rm = TRUE),
      interval_end = max(progress_end, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(is.finite(interval_start), is.finite(interval_end))
}

df <- read_csv(input_path, show_col_types = FALSE) %>%
  filter(!is.na(window_fix_rate), is.finite(window_fix_rate), !is.na(medication), Sex %in% c("f", "m")) %>%
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
    window_fix_rate_sqrt = sqrt(window_fix_rate)
  ) %>%
  filter(fix %in% fix_levels)

prepare_timecourse_stats <- function(data, panel_label, placebo_only = FALSE) {
  plot_data <- if (placebo_only) data %>% filter(medication == "PLA") else data
  plot_data %>%
    group_by(fix, Group, progress_bin, progress_bin_index, progress_midpoint, progress_start, progress_end) %>%
    summarise(
      sem = sd(window_fix_rate_sqrt, na.rm = TRUE) / sqrt(n()),
      mean_window_fix_rate = mean(window_fix_rate_sqrt, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Panel = factor(panel_label, levels = panel_levels),
      fix = factor(fix, levels = fix_levels)
    )
}

fit_preferred_model <- function(data, fix_label) {
  df_fix <- data %>% filter(fix == fix_label)
  lmer(
    window_fix_rate_sqrt ~ Group * progress_bin * medication + Group * Sex * progress_bin + session + (1 | ID) + (1 | subject_session),
    data = df_fix
  )
}

build_window_contrasts <- function(model, source_data, fix_label, panel_label, figure_panel, at_values = list()) {
  emmeans(model, ~ Group | progress_bin, at = at_values) %>%
    contrast(method = "revpairwise") %>%
    summary() %>%
    as_tibble() %>%
    mutate(progress_bin = as.character(progress_bin)) %>%
    left_join(
      source_data %>% distinct(progress_bin, progress_bin_index, progress_start, progress_end),
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
  prepare_timecourse_stats(df, "Overall", placebo_only = FALSE),
  prepare_timecourse_stats(filter(df, Sex == "f"), "Female", placebo_only = FALSE)
)

placebo_stats <- bind_rows(
  prepare_timecourse_stats(df, "Overall", placebo_only = TRUE),
  prepare_timecourse_stats(filter(df, Sex == "f"), "Female", placebo_only = TRUE)
)

overall_models <- list(
  eyes = fit_preferred_model(df, "Eyes"),
  background = fit_preferred_model(df, "Background")
)

window_contrasts <- bind_rows(
  build_window_contrasts(overall_models$eyes, df %>% filter(fix == "Eyes"), "Eyes", "Overall", "A"),
  build_window_contrasts(overall_models$eyes, df %>% filter(fix == "Eyes"), "Eyes", "Female", "A", at_values = list(Sex = "f")),
  build_window_contrasts(overall_models$background, df %>% filter(fix == "Background"), "Background", "Overall", "A"),
  build_window_contrasts(overall_models$background, df %>% filter(fix == "Background"), "Background", "Female", "A", at_values = list(Sex = "f")),
  build_window_contrasts(overall_models$eyes, df %>% filter(fix == "Eyes", medication == "PLA"), "Eyes", "Overall", "B", at_values = list(medication = "PLA")),
  build_window_contrasts(overall_models$eyes, df %>% filter(fix == "Eyes", Sex == "f", medication == "PLA"), "Eyes", "Female", "B", at_values = list(medication = "PLA", Sex = "f")),
  build_window_contrasts(overall_models$background, df %>% filter(fix == "Background", medication == "PLA"), "Background", "Overall", "B", at_values = list(medication = "PLA")),
  build_window_contrasts(overall_models$background, df %>% filter(fix == "Background", Sex == "f", medication == "PLA"), "Background", "Female", "B", at_values = list(medication = "PLA", Sex = "f"))
)

intervals <- window_contrasts %>%
  group_by(Figure_Panel, Panel, fix) %>%
  group_modify(~ build_significance_intervals(.x)) %>%
  ungroup()

interval_positions <- bind_rows(
  top_stats %>% mutate(Figure_Panel = "A"),
  placebo_stats %>% mutate(Figure_Panel = "B")
) %>%
  group_by(Figure_Panel, Panel, fix) %>%
  summarise(interval_y = max(mean_window_fix_rate + sem, na.rm = TRUE) + 0.08, .groups = "drop")

intervals <- intervals %>%
  mutate(interval_start = pmax(interval_start, 0.05), interval_end = pmin(interval_end, 0.95)) %>%
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
    legend.title = element_blank(),
    legend.text = element_text(size = 6.5),
    legend.key = element_blank()
  )

make_timecourse_plot <- function(plot_df, figure_panel, title_text, show_legend = TRUE) {
  panel_intervals <- intervals %>% filter(Figure_Panel == figure_panel)

  ggplot(plot_df, aes(x = progress_midpoint, y = mean_window_fix_rate, group = Group, color = Group)) +
    geom_ribbon(aes(ymin = mean_window_fix_rate - sem, ymax = mean_window_fix_rate + sem, fill = Group), alpha = 0.35, linewidth = 0.15) +
    geom_line(linewidth = 0.45) +
    geom_point(size = 1.1) +
    geom_segment(
      data = panel_intervals %>% filter(direction == "ASD_gt_CTRL"),
      aes(x = interval_start, xend = interval_end, y = interval_y, yend = interval_y),
      inherit.aes = FALSE, color = interval_colors[["ASD_gt_CTRL"]], linewidth = 2.1, lineend = "round"
    ) +
    geom_segment(
      data = panel_intervals %>% filter(direction == "ASD_lt_CTRL"),
      aes(x = interval_start, xend = interval_end, y = interval_y, yend = interval_y),
      inherit.aes = FALSE, color = interval_colors[["ASD_lt_CTRL"]], linewidth = 2.1, lineend = "round"
    ) +
    facet_grid2(Panel ~ fix, axes = "all", remove_labels = "x", switch = "y") +
    scale_x_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0%", "25%", "50%", "75%", "100%"),
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_color_manual(values = group_values) +
    scale_fill_manual(values = group_values) +
    labs(title = title_text, x = "Normalized Interview Progress", y = "sqrt(Fixations/s)") +
    base_theme +
    theme(
      legend.position = if (show_legend) "top" else "none",
      axis.text.x = element_text(size = 8)
    ) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2.4, shape = 15)), fill = "none")
}

top_plot <- make_timecourse_plot(top_stats, "A", "Overall Group Differences", show_legend = TRUE)
placebo_plot <- make_timecourse_plot(placebo_stats, "B", "Group Differences in the Placebo Condition (Exploratory)", show_legend = FALSE)

figure <- top_plot / placebo_plot +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Timecourse Analysis",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.tag = element_text(size = 11, face = "plain")
    )
  )

ggsave(output_png, figure, units = "mm", width = 183, height = 220, dpi = 300)
