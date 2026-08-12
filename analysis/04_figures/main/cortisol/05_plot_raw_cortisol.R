# Manuscript section: Cortisol analyses
# Analysis family: figure generation
# Original source path: legacy/cortisol_reports/markdown_source_bundle/1_cortisol_raw.Rmd
# Primary input dataset(s): data/derived/analysis/df_cortisol_merged.csv
# Primary output(s): in-memory raw cortisol panel for manuscript figure assembly
# Known TODOs: no standalone export path was defined in the original markdown workflow
# Scientific logic note: plotting code lives directly in this topical script so plot edits can be made here

library(tidyverse)
library(ggh4x)
library(here)
source(here("01_analysis", "04_figures", "shared", "custom_plot_settings.R"))

df <- read_csv(here("00_data", "derived", "analysis", "df_cortisol_merged.csv")) %>%
  filter(fix == "Eyes") %>%
  mutate(timepoint = as.factor(timepoint))

data_summary <- df %>%
  group_by(Group, timepoint, medication) %>%
  summarise(
    mean_cortisol_1 = mean(cortisol_1, na.rm = TRUE),
    sem = sd(cortisol_1, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(timepoint = as.factor(timepoint))

data_summary$medication_simple <- factor(
  data_summary$medication,
  levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA"),
  labels = c("BOTH", "NAL", "OXT", "PLA")
)

raw_plot <- ggplot(data_summary, aes(x = timepoint, y = mean_cortisol_1, color = Group, group = Group)) +
  geom_line(linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = mean_cortisol_1 - sem, ymax = mean_cortisol_1 + sem),
    width = 0.1,
    linewidth = 0.3
  ) +
  facet_wrap2(~medication_simple, scales = "fixed", axes = "all", remove_labels = "all") +
  labs(
    title = "Raw Cortisol",
    x = "Timepoint",
    y = "Cortisol Level (nmol/l)"
  ) +
  apply_custom_settings(values = c("#FDC010", "#0B7A6B"))

raw_plot
