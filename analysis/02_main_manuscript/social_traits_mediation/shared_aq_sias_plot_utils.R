# Manuscript section: Shared helper for dimensional and mediation analyses
# Analysis family: shared helper
# Original source path: scripts/publication/analysis/aq_sias_plot_utils.R
# Primary input dataset(s): manuscript analysis tables and mediation panel assets
# Primary output(s): in-memory AQ/SIAS analysis table and composite panel objects
# Known TODOs: review dependence on pre-rendered mediation image assets; formalize canonical input table provenance
# Scientific logic note: copied from source without changing scientific logic

library(here)
library(tidyverse)
library(ggh4x)
library(ggsignif)
library(patchwork)
library(magick)

if (!exists("apply_custom_settings")) {
  source(here("01_analysis", "04_figures", "shared", "custom_plot_settings.R"))
}

prepare_aq_sias_data <- function(df = NULL) {
  if (is.null(df)) {
    df <- read_csv(here("00_data", "derived", "analysis", "df_analysis_pub.csv"))
  }

  df <- df %>%
    mutate(
      fix = factor(
        fix,
        levels = c("fix_on_eyes", "fix_on_mouth", "fix_on_face", "fix_on_background"),
        labels = c("Eyes", "Mouth", "Face", "Background"),
        ordered = TRUE
      )
    ) %>%
    dplyr::select(
      ID, session, Sex, Group, `AQ-K Score`, `SIAS Score`, `IRI Score`, Soz_Scales,
      medication, fix, rate_all, accuracy_degrees
    ) %>%
    rename(
      AQ = `AQ-K Score`,
      SIAS = `SIAS Score`,
      IRI = `IRI Score`
    ) %>%
    distinct(ID, session, fix, .keep_all = TRUE) %>%
    ungroup()

  df$sqrt_rate_all <- sqrt(df$rate_all)

  df
}

create_fixation_plot <- function(df,
                                 score_var,
                                 palette,
                                 base_size,
                                 axis_text_size,
                                 axis_title_size,
                                 title_size,
                                 legend_text_size,
                                 strip_text_size,
                                 panel_spacing,
                                 annotation_size) {
  score_sym <- rlang::sym(score_var)

  df_summary <- df %>%
    group_by(ID, fix) %>%
    reframe(mean_rate_all = mean(sqrt_rate_all), !!score_sym := .data[[score_var]]) %>%
    distinct(ID, fix, .keep_all = TRUE)

  ggplot(
    df_summary,
    aes(x = !!score_sym, y = mean_rate_all, color = fix, fill = fix)
  ) +
    geom_point(
      position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.6),
      alpha = 0.6,
      size = 1
    ) +
    geom_smooth(method = lm, se = TRUE, alpha = 0.15, linewidth = 0.35) +
    labs(
      y = "Fixations / Second",
      x = score_var
    ) +
    ggh4x::facet_wrap2(~fix, scales = "fixed", axes = "all", remove_labels = "all") +
    apply_custom_settings(
      values = palette,
      base_size = base_size,
      axis_text_size = axis_text_size,
      axis_title_size = axis_title_size,
      title_size = title_size,
      legend_text_size = legend_text_size,
      strip_text_size = strip_text_size,
      panel_spacing = panel_spacing
    )
}

build_aq_sias_panels <- function(df = NULL,
                                 palette = c("#FDC010", "#C7B655", "#93AD7C", "#0B7A6B"),
                                 base_size = 8,
                                 axis_text_size = 7,
                                 axis_title_size = 8,
                                 title_size = 9,
                                 legend_text_size = 7,
                                 strip_text_size = 7,
                                 panel_spacing = grid::unit(0.5, "lines"),
                                 annotation_size = 3.2) {
  df_prepared <- prepare_aq_sias_data(df)
  df_cor <- df_prepared %>% distinct(ID, .keep_all = TRUE) %>% filter(!is.na(AQ))

  plot_aq_k <- create_fixation_plot(
    df_prepared,
    "AQ",
    palette,
    base_size,
    axis_text_size,
    axis_title_size,
    title_size,
    legend_text_size,
    strip_text_size,
    panel_spacing,
    annotation_size
  ) +
    labs(title = "Autism Spectrum Quotient") +
    theme(legend.position = "none") +
    geom_text(
      data = subset(df_prepared, fix == "Eyes"),
      aes(x = 18, y = 1.7, label = "r = -0.15*"),
      color = palette[1],
      size = annotation_size
    ) +
    geom_text(
      data = subset(df_prepared, fix == "Background"),
      aes(x = 18, y = 1.7, label = "r = 0.19**"),
      color = palette[4],
      size = annotation_size
    )

  plot_sias <- create_fixation_plot(
    df_prepared,
    "SIAS",
    palette,
    base_size,
    axis_text_size,
    axis_title_size,
    title_size,
    legend_text_size,
    strip_text_size,
    panel_spacing,
    annotation_size
  ) +
    labs(title = "Social Interaction Anxiety") +
    theme(legend.position = "none") +
    geom_text(
      data = subset(df_prepared, fix == "Background"),
      aes(x = 38, y = 1.65, label = "r = 0.22***"),
      color = palette[4],
      size = annotation_size
    )

  corr_sias_aq <- ggplot(df_cor, aes(x = AQ, y = SIAS)) +
    geom_smooth(method = "lm", color = palette[3], fill = palette[3], alpha = 0.2, linewidth = 0.35) +
    geom_point(position = position_jitter(width = 0.2), color = palette[3], alpha = 0.6, size = 1) +
    labs(title = "Social Anxiety and AQ", x = "Autism Spectrum Quotient", y = "SIAS") +
    apply_custom_settings(
      values = palette[3:4],
      base_size = base_size,
      axis_text_size = axis_text_size,
      axis_title_size = axis_title_size,
      title_size = title_size,
      legend_text_size = legend_text_size,
      strip_text_size = strip_text_size,
      panel_spacing = panel_spacing
    ) +
    theme(legend.position = "none") +
    geom_text(
      aes(x = 19, y = 15, label = "r = 0.85***"),
      color = palette[3],
      size = annotation_size
    )

  SIAS_by_group <- ggplot(df_cor, aes(x = Group, fill = Group, color = Group, y = SIAS)) +
    geom_violin(alpha = 0.5, linewidth = 0.35) +
    geom_point(
      position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0.1),
      alpha = 0.6,
      size = 1
    ) +
    stat_summary(fun = mean, geom = "point", size = 1, color = "black") +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.05, color = "black", linewidth = 0.35) +
    labs(title = "Social Anxiety by Group", y = "SIAS", x = NULL) +
    apply_custom_settings(
      values = palette[c(4, 1)],
      base_size = base_size,
      axis_text_size = axis_text_size,
      axis_title_size = axis_title_size,
      title_size = title_size,
      legend_text_size = legend_text_size,
      strip_text_size = strip_text_size,
      panel_spacing = panel_spacing
    ) +
    theme(legend.position = "none") +
    geom_signif(
      comparisons = list(c("ASD", "CTRL")),
      annotations = "***",
      map_signif_level = TRUE,
      y_position = 64,
      tip_length = 0.03,
      color = "black",
      size = 0.35
    ) +
    geom_text(
      aes(x = 1.2, y = 15, label = "d = 2.07"),
      color = "black",
      size = annotation_size
    )

  list(
    plots = list(
      plot_aq_k = plot_aq_k,
      plot_sias = plot_sias,
      corr_sias_aq = corr_sias_aq,
      SIAS_by_group = SIAS_by_group
    ),
    data = list(
      df = df_prepared,
      df_cor = df_cor
    )
  )
}

load_mediation_grobs <- function(scale_width = "1500", trim = TRUE) {
  png_plot_ASD <- image_read(here("02_outputs", "figures", "assets", "Mediation_plot_ASD.png"))
  png_plot_AQ <- image_read(here("02_outputs", "figures", "assets", "Mediation_plot_AQ.png"))

  if (isTRUE(trim)) {
    png_plot_ASD <- image_trim(png_plot_ASD)
    png_plot_AQ <- image_trim(png_plot_AQ)
  }

  if (!is.null(scale_width)) {
    png_plot_ASD <- image_scale(png_plot_ASD, scale_width)
    png_plot_AQ <- image_scale(png_plot_AQ, scale_width)
  }

  list(
    asd = wrap_elements(grid::rasterGrob(as.raster(png_plot_ASD), interpolate = TRUE)),
    aq = wrap_elements(grid::rasterGrob(as.raster(png_plot_AQ), interpolate = TRUE))
  )
}

assemble_compact_aq_sias_mediation <- function(panels,
                                                mediation_grobs = load_mediation_grobs(),
                                                base_size = 8) {
  layout <- "
AABBCC
DDEEFF
"

  wrap_plots(
    A = panels$SIAS_by_group,
    B = panels$plot_aq_k,
    C = panels$plot_sias,
    D = panels$corr_sias_aq,
    E = mediation_grobs$asd,
    F = mediation_grobs$aq,
    design = layout,
    widths = c(0.8, 0.5, 0.7, 0.7, 0.7, 0.7),
    heights = c(1, 0.75)
  ) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = base_size + 2))
}
