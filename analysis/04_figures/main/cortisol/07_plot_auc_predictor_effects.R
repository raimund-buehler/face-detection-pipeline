# Manuscript section: Cortisol analyses
# Analysis family: figure generation
# Original source path: legacy/cortisol_reports/markdown_source_bundle/3_cortisol_auc_pred.Rmd
# Primary input dataset(s): data/derived/analysis/df_cortisol_merged_with_auc.csv
# Primary output(s): in-memory AUC predictor panels for manuscript figure assembly
# Known TODOs: no standalone export path was defined in the original markdown workflow
# Scientific logic note: plotting code lives directly in this topical script so plot edits can be made here

library(tidyverse)
library(emmeans)
library(ggh4x)
library(here)
library(lmerTest)
source(here("01_analysis", "04_figures", "shared", "custom_plot_settings.R"))

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

df_session <- read_csv(here("00_data", "derived", "analysis", "df_cortisol_merged_with_auc.csv"))

model_aucg <- lmer(sqrt_rate_all ~ medication * Group * AUCg + (1 | ID), data = df_session)
interaction_aucg_plot <- build_interaction_plot(model_aucg, df_session, "AUCg") +
  geom_text(
    aes(label = "r = 0.24**"),
    x = 3.5,
    y = 2.8,
    color = "#0B7A6B",
    size = 4,
    data = df_session %>% filter(Group == "CTRL", medication == "PLA/PLA"),
    inherit.aes = FALSE
  )

model_auci <- lmer(sqrt_rate_all ~ medication * Group * AUCi + (1 | ID), data = df_session)
interaction_auci_plot <- build_interaction_plot(model_auci, df_session, "AUCi")

interaction_aucg_plot
interaction_auci_plot
