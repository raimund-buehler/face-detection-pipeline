# Manuscript section: Cortisol analyses
# Analysis family: figure generation
# Original source path: legacy/cortisol_reports/markdown_source_bundle/4_cortisol_maxmin.Rmd
# Primary input dataset(s): data/derived/analysis/df_cortisol_min_max.csv
# Primary output(s): in-memory T2T3 summary and predictor panels for manuscript figure assembly
# Known TODOs: no standalone export path was defined in the original markdown workflow
# Scientific logic note: plotting code lives directly in this topical script so plot edits can be made here

library(tidyverse)
library(emmeans)
library(ggh4x)
library(here)
library(lmerTest)
source(here("01_analysis", "04_figures", "shared", "custom_plot_settings.R"))
source(here("01_analysis", "04_figures", "shared", "custom_theme_settings_publication.R"))

build_interaction_plot <- function(model, panel_df, cortisol_predictor) {
  emm <- emtrends(model, ~ medication * Group, var = cortisol_predictor)
  emm_df <- as.data.frame(emm)
  colnames(emm_df)[colnames(emm_df) == paste0(cortisol_predictor, ".trend")] <- "trend"

  emm_df <- emm_df %>%
    mutate(
      lower.CL = trend - SE * qt(0.975, df),
      upper.CL = trend + SE * qt(0.975, df)
    )

  ggplot(emm_df, aes_string(x = "medication", y = "trend", group = "Group", color = "Group")) +
    geom_line(size = 0.25) +
    geom_ribbon(
      aes_string(ymin = "lower.CL", ymax = "upper.CL", fill = "Group"),
      alpha = 0.4,
      size = 0.25,
      color = NA
    ) +
    geom_point(
      data = panel_df,
      aes_string(x = cortisol_predictor, y = "sqrt_rate_all", color = "Group"),
      inherit.aes = FALSE,
      position = position_jitter(width = 0.1, height = 0),
      alpha = 0.4,
      size = 1
    ) +
    theme_minimal() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      x = "Medication Condition",
      y = paste("Predicted Effect of", cortisol_predictor, "on Eye Gaze"),
      title = "Predicted Slopes of Cortisol Effect by Group and Medication"
    ) +
    facet_wrap2(~Group, scales = "fixed", axes = "all", remove_labels = "all") +
    apply_custom_settings(values = c("#FDC010", "#0B7A6B"))
}

df <- read_csv(here("00_data", "derived", "analysis", "df_cortisol_min_max.csv")) %>%
  rename(T2T3 = MinMax, T2T3_scaled = MinMax_scaled)

df_long <- df %>%
  select(Group, medication, T2T3, AUCi) %>%
  pivot_longer(
    cols = c(T2T3, AUCi),
    names_to = "reactivity_type",
    values_to = "reactivity_value"
  ) %>%
  mutate(
    medication_simple = factor(
      medication,
      levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA"),
      labels = c("BOTH", "NAL", "OXT", "PLA")
    )
  )

strip_reactivity <- function(labels) {
  sub("^,\\s*", "", labels)
}

reactivity_barplot <- ggplot(df_long, aes(x = medication_simple, y = reactivity_value, fill = Group)) +
  geom_bar(
    stat = "summary",
    fun = "mean",
    position = position_dodge(width = 0.7),
    width = 0.7,
    alpha = 0.6
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.25,
    position = position_dodge(width = 0.7),
    linewidth = 0.2
  ) +
  scale_x_discrete(limits = c("PLA", "OXT", "NAL", "BOTH")) +
  geom_hline(yintercept = 0, color = "grey30") +
  facet_wrap2(
    ~ Group * reactivity_type,
    axes = "margins",
    remove_labels = "all",
    labeller = labeller(
      Group = function(x) rep("", length(x)),
      reactivity_type = function(x) strip_reactivity(x)
    ),
    strip = strip_themed(
      background_x = element_blank(),
      background_y = element_blank()
    )
  ) +
  labs(title = "Cortisol Reactivity", x = "Medication", y = "Reactivity Value") +
  publication_palette(aesthetics = "fill") +
  coord_flip()

model_t2t3 <- lmer(sqrt_rate_all ~ medication * Group * T2T3 + (1 | ID), data = df)
interaction_t2t3_plot <- build_interaction_plot(model_t2t3, df, "T2T3")

reactivity_barplot
interaction_t2t3_plot
