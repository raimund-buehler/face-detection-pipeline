# Manuscript section: Exploratory dimensional trait figure
# Analysis family: AQ/SIAS dimensional visualization
# Primary input dataset(s): 00_data/derived/preprocessing/duration_data.csv; 00_data/derived/analysis/df_analysis_pub.csv
# Primary output(s): exploratory two-panel AQ simple-slope and AQ-SIAS correlation figure

library(tidyverse)
library(glmmTMB)
library(emmeans)
library(ggh4x)
library(patchwork)
library(here)

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

output_dir <- here("02_outputs", "figures", "supplementary", "social_traits")
current_story_dir <- here("02_outputs", "figures", "main", "current_story")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(current_story_dir, recursive = TRUE, showWarnings = FALSE)

output_png <- file.path(output_dir, "Figure_4_dimensional_aq_exploratory.png")
output_svg <- file.path(output_dir, "figure4_dimensional_aq_exploratory.svg")
output_panel_a_png <- file.path(output_dir, "Figure_4_dimensional_aq_panel_A_exploratory.png")
output_panel_a_svg <- file.path(output_dir, "figure4_dimensional_aq_panel_A_exploratory.svg")
output_panel_b_png <- file.path(output_dir, "Figure_4_dimensional_aq_panel_B_exploratory.png")
output_panel_b_svg <- file.path(output_dir, "figure4_dimensional_aq_panel_B_exploratory.svg")
output_with_hist_png <- file.path(output_dir, "Figure_4_dimensional_aq_with_histogram_exploratory.png")
output_with_hist_svg <- file.path(output_dir, "figure4_dimensional_aq_with_histogram_exploratory.svg")
output_with_density_png <- file.path(output_dir, "Figure_4_dimensional_aq_with_panel_b_density_exploratory.png")
output_with_density_svg <- file.path(output_dir, "figure4_dimensional_aq_with_panel_b_density_exploratory.svg")
current_story_png <- file.path(current_story_dir, "Figure_4_dimensional_aq.png")
current_story_svg <- file.path(current_story_dir, "Figure_4_dimensional_aq.svg")
current_story_combined_png <- file.path(current_story_dir, "Figure_4_dimensional_aq_cortisol.png")
current_story_combined_svg <- file.path(current_story_dir, "Figure_4_dimensional_aq_cortisol.svg")

group_values <- c("ASD" = "#D55E00", "CTRL" = "#0072B2")
use_aoi_colors <- FALSE
aoi_values_color <- c("Eyes" = "#A23B72", "Background" = "#5E6B32")
aoi_values_gray <- c("Eyes" = "#4A4A4A", "Background" = "#4A4A4A")
aoi_values <- if (isTRUE(use_aoi_colors)) aoi_values_color else aoi_values_gray
panel_levels <- c("Full Sample", "Female")
fix_levels <- c("Eyes", "Background")

session_metadata <- read_csv(
  here("00_data", "derived", "analysis", "df_analysis_pub.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    ID,
    session,
    session_norm = normalize_session(session),
    Sex,
    Group,
    medication,
    AQ = `AQ-K Score`,
    SIAS = `SIAS Score`
  ) %>%
  distinct()

df <- read_csv(
  here("00_data", "derived", "preprocessing", "duration_data.csv"),
  show_col_types = FALSE
) %>%
  mutate(ID = Sub_ID, session_norm = normalize_session(session)) %>%
  left_join(session_metadata, by = c("ID", "session_norm")) %>%
  mutate(session = coalesce(session.y, session.x)) %>%
  filter(!is.na(Group), Sex %in% c("f", "m"), !is.na(AQ)) %>%
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
    Sex = factor(Sex, levels = c("f", "m")),
    session = factor(session),
    AQ_z = as.numeric(scale(AQ)),
    percentage_fix_duration_beta = beta_squeeze(percentage_fix_duration)
  ) %>%
  select(ID, session, Sex, Group, medication, fix, percentage_fix_duration, percentage_fix_duration_beta, AQ, SIAS, AQ_z)

aq_center <- attr(scale(df$AQ), "scaled:center")
aq_scale <- attr(scale(df$AQ), "scaled:scale")

aq_model <- glmmTMB(
  percentage_fix_duration_beta ~ AQ_z * fix * medication + AQ_z * fix * Sex + session,
  family = beta_family(),
  data = df,
  control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
)

aq_raw_range <- seq(
  floor(min(df$AQ, na.rm = TRUE)),
  ceiling(max(df$AQ, na.rm = TRUE)),
  length.out = 80
)
aq_range <- (aq_raw_range - aq_center) / aq_scale

pred_full <- emmeans(
  aq_model,
  ~ AQ_z | fix,
  at = list(AQ_z = aq_range, fix = fix_levels),
  type = "response"
) %>%
  summary(infer = c(TRUE, TRUE)) %>%
  as_tibble() %>%
  mutate(Panel = "Full Sample")

pred_female <- emmeans(
  aq_model,
  ~ AQ_z | fix,
  at = list(AQ_z = aq_range, fix = fix_levels, Sex = "f"),
  type = "response"
) %>%
  summary(infer = c(TRUE, TRUE)) %>%
  as_tibble() %>%
  mutate(Panel = "Female")

pred_df <- bind_rows(pred_full, pred_female) %>%
  filter(fix %in% fix_levels) %>%
  mutate(
    AQ = AQ_z * aq_scale + aq_center,
    Panel = factor(Panel, levels = panel_levels),
    fix = factor(fix, levels = fix_levels)
  )

points_full <- df %>%
  filter(fix %in% fix_levels) %>%
  group_by(ID, fix, AQ_z) %>%
  summarise(
    AQ = mean(AQ, na.rm = TRUE),
    percentage_fix_duration = mean(percentage_fix_duration, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Panel = "Full Sample")

points_female <- df %>%
  filter(fix %in% fix_levels, Sex == "f") %>%
  group_by(ID, fix, AQ_z) %>%
  summarise(
    AQ = mean(AQ, na.rm = TRUE),
    percentage_fix_duration = mean(percentage_fix_duration, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Panel = "Female")

points_df <- bind_rows(points_full, points_female) %>%
  mutate(
    Panel = factor(Panel, levels = panel_levels),
    fix = factor(fix, levels = fix_levels)
  )

slope_labels <- tibble(
  Panel = factor(c("Full Sample", "Full Sample", "Female", "Female"), levels = panel_levels),
  fix = factor(c("Eyes", "Background", "Eyes", "Background"), levels = fix_levels),
  label = c(
    "z = -1.70, p = .089",
    "z = 3.15, p = .002",
    "z = -2.47, p = .013",
    "z = 3.60, p < .001"
  ),
  x = min(aq_raw_range) + 0.08 * diff(range(aq_raw_range)),
  y = c(0.94, 0.94, 0.94, 0.94)
)

panel_a <- ggplot() +
  geom_point(
    data = points_df,
    aes(x = AQ, y = percentage_fix_duration, color = fix),
    alpha = 0.22,
    size = 0.9,
    position = position_jitter(width = 0.035, height = 0)
  ) +
  geom_ribbon(
    data = pred_df,
    aes(x = AQ, ymin = asymp.LCL, ymax = asymp.UCL, fill = fix),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(
    data = pred_df,
    aes(x = AQ, y = response, color = fix),
    linewidth = 0.85
  ) +
  geom_text(
    data = slope_labels,
    aes(x = x, y = y, label = label, color = fix),
    hjust = 0,
    size = 2.45,
    show.legend = FALSE
  ) +
  facet_grid2(Panel ~ fix, axes = "all", remove_labels = "x", switch = "y") +
  scale_color_manual(values = aoi_values, name = "AOI") +
  scale_fill_manual(values = aoi_values, name = "AOI") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = range(aq_raw_range), ylim = c(0, 1), clip = "off") +
  labs(
    title = "AQ Simple Slopes",
    x = "AQ",
    y = "% Fixation Duration"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title.position = "plot",
    aspect.ratio = 1,
    legend.position = "none",
    legend.justification = "left",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 9, face = "plain", margin = margin(1, 1, 1, 1)),
    plot.title = element_text(size = 11, face = "bold", hjust = 0, margin = margin(b = 2, l = 8)),
    panel.spacing = grid::unit(0.25, "lines"),
    plot.margin = margin(1, 1, 1, 1)
  )

df_cor <- session_metadata %>%
  distinct(ID, .keep_all = TRUE) %>%
  filter(!is.na(AQ), !is.na(SIAS)) %>%
  mutate(
    Sex_symbol = recode(Sex, "f" = "\u2640", "m" = "\u2642"),
    Sex_label = factor(recode(Sex, "f" = "Female", "m" = "Male"), levels = c("Female", "Male"))
  )

aq_axis_max <- ceiling(max(df_cor$AQ, na.rm = TRUE) / 5) * 5
sias_axis_max <- ceiling(max(df_cor$SIAS, na.rm = TRUE) / 10) * 10

density_height <- sias_axis_max * 0.025
density_lower_limit <- -sias_axis_max * 0.16
density_lanes <- tibble(
  Group = factor(c("CTRL", "CTRL", "ASD", "ASD"), levels = names(group_values)),
  Sex_label = factor(c("Female", "Male", "Female", "Male"), levels = c("Female", "Male")),
  density_base = -sias_axis_max * c(0.140, 0.105, 0.070, 0.035),
  density_label = c("CTRL F", "CTRL M", "ASD F", "ASD M")
)

aq_density_df <- df_cor %>%
  group_by(Group, Sex_label) %>%
  group_modify(~ {
    if (sum(!is.na(.x$AQ)) < 2) {
      return(tibble(AQ = numeric(), density_scaled = numeric()))
    }
    density_est <- density(.x$AQ, from = 0, to = aq_axis_max, n = 160, na.rm = TRUE)
    tibble(
      AQ = density_est$x,
      density_scaled = density_est$y / max(density_est$y, na.rm = TRUE)
    )
  }) %>%
  ungroup() %>%
  left_join(density_lanes, by = c("Group", "Sex_label")) %>%
  mutate(density_y = density_base + density_scaled * density_height)

aq_density_means <- df_cor %>%
  group_by(Group, Sex_label) %>%
  summarise(mean_aq = mean(AQ, na.rm = TRUE), .groups = "drop") %>%
  left_join(density_lanes, by = c("Group", "Sex_label"))

aq_sex_contrasts <- read_csv(
  here(
    "02_outputs", "model_outputs", "main_manuscript", "eye_tracking",
    "group_aoi_medication_duration", "sensitivity", "aq_adjusted_group_sex",
    "04_aq_sex_contrasts_by_group.csv"
  ),
  show_col_types = FALSE
)

asd_sex_contrast <- aq_sex_contrasts %>%
  filter(Group == "ASD", contrast == "f - m") %>%
  slice(1)

asd_density_bracket <- density_lanes %>%
  filter(Group == "ASD", Sex_label %in% c("Female", "Male")) %>%
  left_join(aq_density_means, by = c("Group", "Sex_label", "density_base", "density_label")) %>%
  summarise(
    x_low = min(mean_aq),
    x_high = max(mean_aq),
    y = max(density_base + density_height * 1.72),
    y_tip = max(density_base + density_height * 1.25),
    label = if_else(nrow(asd_sex_contrast) > 0, asd_sex_contrast$sig[[1]], "*"),
    .groups = "drop"
  )

cor_test <- cor.test(df_cor$AQ, df_cor$SIAS, method = "pearson")
cor_label <- paste0(
  "r = ", sprintf("%.2f", unname(cor_test$estimate)),
  "\n95% CI [", sprintf("%.2f", cor_test$conf.int[[1]]), ", ",
  sprintf("%.2f", cor_test$conf.int[[2]]), "]\np < .001"
)

panel_b <- ggplot(df_cor, aes(x = AQ, y = SIAS)) +
  geom_smooth(method = "lm", color = "#2F6F62", fill = "#2F6F62", alpha = 0.18, linewidth = 0.75) +
  geom_point(aes(color = Group), alpha = 0, size = 2.2, show.legend = TRUE) +
  geom_text(aes(color = Group, label = Sex_symbol), alpha = 0.78, size = 3.1, show.legend = FALSE) +
  annotate(
    "text",
    x = min(df_cor$AQ, na.rm = TRUE) + 0.05 * diff(range(df_cor$AQ, na.rm = TRUE)),
    y = max(df_cor$SIAS, na.rm = TRUE) - 0.05 * diff(range(df_cor$SIAS, na.rm = TRUE)),
    label = cor_label,
    hjust = 0,
    vjust = 1,
    size = 3.0,
    color = "#303030"
  ) +
  scale_color_manual(values = group_values, name = "Group") +
  scale_x_continuous(limits = c(0, aq_axis_max), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, sias_axis_max), expand = expansion(mult = c(0, 0))) +
  coord_cartesian(xlim = c(0, aq_axis_max), ylim = c(0, sias_axis_max), expand = FALSE, clip = "off") +
  labs(
    title = "AQ and Social Anxiety",
    x = "AQ",
    y = "SIAS"
  ) +
  theme_classic(base_size = 9) +
  theme(
    aspect.ratio = 1,
    legend.position = c(0.5, 1.12),
    legend.justification = c(0.5, 1),
    legend.direction = "horizontal",
    legend.background = element_rect(fill = scales::alpha("white", 0.82), color = NA),
    legend.margin = margin(0, 2, 0, 2),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    axis.ticks.length = grid::unit(1.7, "pt"),
    plot.title = element_text(size = 11, face = "bold", hjust = 0, margin = margin(b = 2, l = 8)),
    plot.margin = margin(1, 1, 1, 1)
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, shape = 15, size = 3.2)))

panel_b_density <- panel_b +
  geom_ribbon(
    data = aq_density_df,
    aes(x = AQ, ymin = density_base, ymax = density_y, fill = Group, group = interaction(Group, Sex_label)),
    inherit.aes = FALSE,
    alpha = 0.14,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line(
    data = aq_density_df,
    aes(x = AQ, y = density_y, color = Group, linetype = Sex_label, group = interaction(Group, Sex_label)),
    inherit.aes = FALSE,
    linewidth = 0.35,
    alpha = 0.8,
    show.legend = FALSE
  ) +
  geom_segment(
    data = aq_density_means,
    aes(
      x = mean_aq,
      xend = mean_aq,
      y = density_base,
      yend = density_base + density_height * 1.08,
      color = Group
    ),
    inherit.aes = FALSE,
    linewidth = 0.55,
    show.legend = FALSE
  ) +
  geom_text(
    data = density_lanes,
    aes(x = aq_axis_max * 0.875, y = density_base + density_height * 0.45, label = density_label, color = Group),
    inherit.aes = FALSE,
    hjust = 1,
    size = 2.05,
    show.legend = FALSE
  ) +
  geom_segment(
    data = asd_density_bracket,
    aes(x = x_low, xend = x_high, y = y, yend = y),
    inherit.aes = FALSE,
    color = group_values[["ASD"]],
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  geom_segment(
    data = asd_density_bracket,
    aes(x = x_low, xend = x_low, y = y, yend = y_tip),
    inherit.aes = FALSE,
    color = group_values[["ASD"]],
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  geom_segment(
    data = asd_density_bracket,
    aes(x = x_high, xend = x_high, y = y, yend = y_tip),
    inherit.aes = FALSE,
    color = group_values[["ASD"]],
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  geom_text(
    data = asd_density_bracket,
    aes(x = (x_low + x_high) / 2, y = y + density_height * 0.25, label = label),
    inherit.aes = FALSE,
    color = group_values[["ASD"]],
    hjust = 0.5,
    vjust = 0,
    size = 2.8,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = group_values) +
  scale_linetype_manual(values = c("Female" = "solid", "Male" = "22")) +
  scale_y_continuous(
    breaks = seq(0, sias_axis_max, by = 20),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(
    xlim = c(0, aq_axis_max),
    ylim = c(density_lower_limit, sias_axis_max),
    expand = FALSE,
    clip = "off"
  ) +
  guides(fill = "none", linetype = "none")

aq_hist_means <- df_cor %>%
  group_by(Group, Sex_label) %>%
  summarise(mean_aq = mean(AQ, na.rm = TRUE), .groups = "drop")

panel_c <- ggplot(df_cor, aes(x = AQ, fill = Group, color = Group)) +
  geom_histogram(binwidth = 3, boundary = 0, alpha = 0.28, linewidth = 0.35) +
  geom_vline(
    data = aq_hist_means,
    aes(xintercept = mean_aq, color = Group),
    linewidth = 0.8
  ) +
  facet_grid(Group ~ Sex_label) +
  scale_fill_manual(values = group_values) +
  scale_color_manual(values = group_values) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  coord_cartesian(xlim = c(0, aq_axis_max), expand = FALSE) +
  labs(
    title = "AQ Distribution by Group and Sex",
    x = "AQ",
    y = "Count"
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    strip.background = element_blank(),
    strip.text = element_text(size = 8, face = "plain", margin = margin(1, 1, 1, 1)),
    panel.spacing = grid::unit(0.35, "lines"),
    plot.title = element_text(size = 11, face = "bold", hjust = 0, margin = margin(b = 2, l = 8)),
    plot.margin = margin(1, 1, 1, 1)
  )

combined_plot <- panel_a + panel_b +
  plot_layout(widths = c(1.08, 1.0)) +
  plot_annotation(
    title = "Dimensional Trait Analysis (Exploratory)",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(b = 1)),
      plot.tag = element_text(size = 10, face = "plain"),
      plot.margin = margin(0, 0, 0, 0)
    )
  )

combined_plot_with_hist <- (panel_a + panel_b + plot_layout(widths = c(1.08, 1.0))) / panel_c +
  plot_layout(heights = c(1.0, 0.38)) +
  plot_annotation(
    title = "Dimensional Trait Analysis (Exploratory)",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(b = 1)),
      plot.tag = element_text(size = 10, face = "plain"),
      plot.margin = margin(0, 0, 0, 0)
    )
  )

combined_plot_with_density <- panel_a + panel_b_density +
  plot_layout(widths = c(1.08, 1.0)) +
  plot_annotation(
    title = "Dimensional Trait Analysis (Exploratory)",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(b = 1)),
      plot.tag = element_text(size = 10, face = "plain"),
      plot.margin = margin(0, 0, 0, 0)
    )
  )

cortisol_env <- new.env(parent = globalenv())
sys.source(
  here(
    "01_analysis", "04_figures", "main", "cortisol",
    "05_plot_raw_cortisol_group_timecourse.R"
  ),
  envir = cortisol_env
)

cortisol_panel <- cortisol_env$figure_5 +
  scale_x_discrete(expand = expansion(add = 0.35)) +
  scale_y_continuous(
    limits = c(1, 6),
    breaks = 1:6,
    expand = expansion(mult = c(0, 0))
  ) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 11, face = "bold", hjust = 0, margin = margin(b = 2, l = 8)),
    legend.position = "right",
    legend.justification = "top",
    legend.direction = "vertical",
    legend.box.spacing = grid::unit(0, "pt"),
    legend.margin = margin(5, 0, 0, 4),
    legend.text = element_text(size = 7),
    strip.text = element_text(size = 8),
    panel.spacing.x = grid::unit(7, "pt"),
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    axis.title.x = element_text(size = 8, margin = margin(t = 2)),
    axis.title.y = element_text(size = 8, margin = margin(r = 2)),
    axis.text = element_text(size = 7),
    plot.margin = margin(0, 8, 0, 18)
  )

cortisol_panel_free <- patchwork::free(
  cortisol_panel,
  type = "space",
  side = "l"
)

panel_b_density_tall <- panel_b_density +
  theme(
    aspect.ratio = NULL,
    plot.title.position = "plot",
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.margin = margin(0, 0, 1, 0),
    legend.box.spacing = grid::unit(0, "pt"),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 11, face = "bold", hjust = 0, margin = margin(b = 1, l = 8))
  )

current_story_combined <- wrap_plots(
  A = panel_a,
  B = cortisol_panel_free,
  C = panel_b_density_tall,
  design = c(
    area(t = 1, l = 1, b = 1, r = 1),
    area(t = 2, l = 1, b = 2, r = 1),
    area(t = 1, l = 2, b = 2, r = 2)
  ),
  widths = c(0.80, 1.20),
  heights = c(1.0, 0.48)
) +
  plot_annotation(
    title = "Dimensional Trait and Cortisol Analysis (Exploratory)",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(b = 1)),
      plot.tag = element_text(size = 10, face = "plain"),
      plot.margin = margin(0, 0, 0, 0)
    )
  )

ggsave(output_svg, combined_plot, units = "mm", width = 168, height = 112)
ggsave(output_png, combined_plot, units = "mm", width = 168, height = 112, dpi = 300)
ggsave(current_story_svg, current_story_combined, units = "mm", width = 178, height = 118)
ggsave(current_story_png, current_story_combined, units = "mm", width = 178, height = 118, dpi = 300)
ggsave(current_story_combined_svg, current_story_combined, units = "mm", width = 178, height = 118)
ggsave(current_story_combined_png, current_story_combined, units = "mm", width = 178, height = 118, dpi = 300)
ggsave(output_with_density_svg, combined_plot_with_density, units = "mm", width = 168, height = 112)
ggsave(output_with_density_png, combined_plot_with_density, units = "mm", width = 168, height = 112, dpi = 300)
ggsave(output_with_hist_svg, combined_plot_with_hist, units = "mm", width = 168, height = 145)
ggsave(output_with_hist_png, combined_plot_with_hist, units = "mm", width = 168, height = 145, dpi = 300)
ggsave(output_panel_a_svg, panel_a, units = "mm", width = 105, height = 125)
ggsave(output_panel_a_png, panel_a, units = "mm", width = 105, height = 125, dpi = 300)
ggsave(output_panel_b_svg, panel_b, units = "mm", width = 95, height = 95)
ggsave(output_panel_b_png, panel_b, units = "mm", width = 95, height = 95, dpi = 300)
