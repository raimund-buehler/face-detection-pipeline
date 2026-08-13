# Manuscript section: Main eye-tracking model
# Analysis family: group x AOI x medication
# Original source path: scripts/publication/analysis/overall_model_pub.R
# Primary input dataset(s): data/derived/analysis/df_analysis_pub_prep.csv
# Primary output(s): model summaries and emmeans contrasts for the main eye-tracking result
# Known TODOs: df_analysis_pub_prep.csv provenance is unresolved; figure generation moved to analysis/04_figures/main/eye_tracking/01_plot_group_aoi_medication.R
# Scientific logic note: scientific logic and model formulas are unchanged from source; plotting was separated for structural cleanup

library(tidyverse)
library(lmerTest)
library(emmeans)
library(effectsize)
library(performance)
library(here)
source(here("01_analysis", "shared", "model_output_utils.R"))

df <- read_csv(here("00_data", "derived", "analysis", "df_analysis_pub_prep.csv"))

df$fix <- factor(
  df$fix,
  levels = c("Eyes", "Mouth", "Face", "Background"),
  labels = c("Eyes", "Mouth", "Face", "Background"),
  ordered = TRUE
)

df$medication_simple <- factor(
  df$medication,
  levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA"),
  labels = c("BOTH", "NAL", "OXT", "PLA")
)

model <- lmer(
  sqrt_rate_all ~ Group * fix * medication + Sex + session + (1 | ID),
  data = df
)

output_dir <- here("02_outputs", "model_outputs", "main_manuscript", "eye_tracking", "group_aoi_medication")

summary(model)
anova(model)
write_model_anova(model, output_dir, "01_main_model_anova.csv", "main_group_fix_medication")

emmeans(model, pairwise ~ Sex)
emmeans(model, pairwise ~ fix)

df %>%
  group_by(fix) %>%
  summarise(
    Mean = mean(sqrt_rate_all, na.rm = TRUE),
    SD = sd(sqrt_rate_all, na.rm = TRUE),
    N = n()
  )

emmeans(model, pairwise ~ Group, by = c("fix"))

# Eyes
t_to_d(-2.036, 285)

# Background
t_to_d(2.289, 285)

emmeans(model, pairwise ~ Group, by = c("fix", "medication"))
emmeans(model, pairwise ~ medication, by = c("fix"))
emmeans(model, pairwise ~ medication, by = c("fix", "Group"))

df %>%
  filter(fix == "Eyes") %>%
  distinct(ID, session, .keep_all = TRUE) %>%
  select(ID, session, Group, medication, rate_all, accuracy_degrees) %>%
  arrange(-rate_all)

df %>%
  filter(fix == "Eyes") %>%
  distinct(ID, session, .keep_all = TRUE) %>%
  select(ID, session, Group, medication, rate_all, accuracy_degrees) %>%
  filter(accuracy_degrees < 1.5) %>%
  arrange(rate_all) %>%
  print(n = 20)

df_eyes <- df %>%
  filter(fix == "Eyes")

model_eyes <- lmer(
  sqrt_rate_all ~ Group * medication + Sex + session + (1 | ID),
  data = df_eyes
)

summary(model_eyes)
anova(model_eyes)
write_model_anova(model_eyes, output_dir, "02_eyes_only_model_anova.csv", "eyes_only_group_medication")
