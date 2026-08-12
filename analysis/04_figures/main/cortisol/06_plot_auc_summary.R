# Manuscript section: Cortisol analyses
# Analysis family: figure generation
# Original source path: legacy/cortisol_reports/markdown_source_bundle/2_cortisol_auc.Rmd
# Primary input dataset(s): data/derived/analysis/df_cortisol_merged_with_auc.csv
# Primary output(s): in-memory AUC summary panel for manuscript figure assembly
# Known TODOs: no standalone export path was defined in the original markdown workflow
# Scientific logic note: plotting code lives directly in this topical script so plot edits can be made here

library(tidyverse)
library(ggh4x)
library(ggsignif)
library(here)
source(here("01_analysis", "04_figures", "shared", "custom_plot_settings.R"))

df <- read_csv(here("00_data", "derived", "analysis", "df_cortisol_merged_with_auc.csv"))

df_long <- df %>%
  select(Group, medication, AUCg, AUCi) %>%
  pivot_longer(cols = c(AUCg, AUCi), names_to = "AUC_type", values_to = "AUC_value")

auc_summary_plot <- ggplot(df_long, aes(x = medication, y = AUC_value, fill = AUC_type)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.9), width = 0.7, alpha = 0.6) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1, position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 0, color = "grey30") +
  facet_wrap2(~Group * AUC_type, axes = "all", remove_labels = "all") +
  labs(title = "AUCg and AUCi by Medication", x = "Medication", y = "AUC Value") +
  apply_custom_settings(values = c("#C7B655", "#93AD7C")) +
  coord_flip() +
  geom_signif(
    comparisons = list(c("NAL/OXT", "PLA/PLA"), c("NAL/PLA", "PLA/OXT"), c("NAL/PLA", "PLA/PLA")),
    annotations = c("*", "*", "***"),
    y_position = c(11, 9, 7),
    tip_length = 0.01,
    vjust = 0.5,
    data = df_long %>% filter(AUC_type == "AUCg", Group == "CTRL")
  ) +
  geom_signif(
    comparisons = list(c("NAL/OXT", "PLA/OXT"), c("NAL/PLA", "PLA/OXT")),
    annotations = c("*", "*"),
    y_position = c(0, 2),
    tip_length = 0.01,
    vjust = 0.5,
    data = df_long %>% filter(AUC_type == "AUCi", Group == "ASD")
  )

auc_summary_plot
