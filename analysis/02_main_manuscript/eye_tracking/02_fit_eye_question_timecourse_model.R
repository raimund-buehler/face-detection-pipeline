# Manuscript section: Timecourse and question-type analysis
# Analysis family: fixations on eyes over interview questions / open vs closed questions
# Original source path: scripts/publication/analysis/within_questions.R
# Primary input dataset(s): data/derived/analysis/df_analysis_pub.csv
# Primary output(s): model summaries, emmeans contrasts, and plot-ready summary tables for the timecourse result
# Known TODOs: confirm whether df_analysis_pub.csv or a prep variant is the intended canonical input; figure generation moved to analysis/04_figures/main/eye_tracking/02_plot_eye_question_timecourse.R
# Scientific logic note: scientific logic and model formulas are unchanged from source; plotting was separated for structural cleanup

library(tidyverse)
library(lmerTest)
library(emmeans)
library(ggh4x)
library(here)
source(here("01_analysis", "shared", "model_output_utils.R"))

df <- read_csv(here("00_data", "derived", "analysis", "df_analysis_pub.csv"))

df <- df %>% filter(!is.na(rate_q))
df <- df %>% filter(!is.na(medication))

df <- df %>%
  select(ID, session, Sex, Group, medication, label, fix, rate_q, accuracy_degrees) %>%
  distinct(ID, session, fix, label, .keep_all = TRUE) %>%
  ungroup()

df$sqrt_rate_q <- sqrt(df$rate_q)

df$fix <- factor(
  df$fix,
  levels = c("fix_on_eyes", "fix_on_mouth", "fix_on_face", "fix_on_background"),
  labels = c("Eyes", "Mouth", "Face", "Background")
)

ordered_labels <- c("Intro", paste("Question", 1:13))

df$label <- factor(
  df$label,
  levels = ordered_labels,
  labels = c("Intro", paste(1:13)),
  ordered = TRUE
)

df <- df %>% filter(fix == "Eyes", label != "Intro")

df <- df %>%
  mutate(q_type = ifelse((df$label %in% paste(1:6)), "closed", "open"))

q_type_model <- lmer(
  data = df,
  sqrt_rate_q ~ Group * q_type * medication + Sex + session + (1 | ID/q_type)
)

output_dir <- here("02_outputs", "model_outputs", "main_manuscript", "eye_tracking", "timecourse_questions")

summary(q_type_model)
anova(q_type_model)
write_model_anova(q_type_model, output_dir, "01_question_type_model_anova.csv", "question_type_model")

emmeans(q_type_model, pairwise ~ Group)
emmeans(q_type_model, pairwise ~ Group | medication)
emmeans(q_type_model, pairwise ~ Group | q_type)
emmeans(q_type_model, pairwise ~ Group | q_type + medication)
emmeans(q_type_model, pairwise ~ q_type)
emmeans(q_type_model, pairwise ~ q_type | Group)
emmeans(q_type_model, pairwise ~ q_type | medication)

emm_options(pbkrtest.limit = 3046)
emmeans(q_type_model, pairwise ~ q_type | medication * Group)

label_model <- lmer(data = df, sqrt_rate_q ~ Group * label * medication + (1 | ID))

summary(label_model)
anova(label_model)
write_model_anova(label_model, output_dir, "02_label_timecourse_model_anova.csv", "label_timecourse_model")

emmeans(label_model, pairwise ~ Group)
emmeans(label_model, pairwise ~ Group | medication)
emmeans(label_model, pairwise ~ Group, by = c("label", "medication"))

timecourse_group_stats <- df %>%
  group_by(Group, medication, label) %>%
  summarise(
    mean_sqrt_rate_q = mean(sqrt_rate_q, na.rm = TRUE),
    sem = sd(sqrt_rate_q, na.rm = TRUE) / sqrt(n()),
    lower_ci = mean_sqrt_rate_q - qt(0.975, df = n() - 1) * sem,
    upper_ci = mean_sqrt_rate_q + qt(0.975, df = n() - 1) * sem,
    .groups = "drop"
  )

timecourse_emm_results <- emmeans(label_model, pairwise ~ Group | label * medication)
timecourse_contrasts <- summary(timecourse_emm_results$contrasts)

timecourse_contrasts$stars <- ifelse(
  timecourse_contrasts$`p.value` < 0.01,
  "**",
  ifelse(timecourse_contrasts$`p.value` < 0.05, "*", "")
)

timecourse_contrasts$label_medication <- paste(
  timecourse_contrasts$label,
  timecourse_contrasts$medication,
  sep = "_"
)

timecourse_group_stats$label_medication <- paste(
  timecourse_group_stats$label,
  timecourse_group_stats$medication,
  sep = "_"
)

timecourse_group_stats <- merge(
  timecourse_group_stats,
  timecourse_contrasts[, c("label_medication", "stars")],
  by = "label_medication",
  all.x = TRUE
)
