# Manuscript section: Supplementary / control analyses
# Analysis family: menstrual cycle analysis
# Original source path: scripts/publication/analysis/menst_analysis.r
# Primary input dataset(s): data/derived/analysis/df_analysis_pub.csv
# Primary output(s): outputs/figures/supplementary/menst_phases.svg; outputs/figures/supplementary/Supp_Figure_1.png
# Known TODOs: phase provenance depends on the external menstrual-prep step
# Scientific logic note: plotting logic was separated from the original mixed script

library(tidyverse)
library(here)
source(here("01_analysis", "04_figures", "shared", "custom_plot_settings.R"))

normalize_session <- function(x) {
  x %>%
    str_trim() %>%
    str_to_lower() %>%
    str_replace_all(" ", "_")
}

session_metadata <- read_csv(
  here("00_data", "derived", "analysis", "df_analysis_pub.csv"),
  show_col_types = FALSE
) %>%
  distinct(ID, session, phase, Sex, Group, medication) %>%
  mutate(session_norm = normalize_session(session))

df <- read_csv(
  here("00_data", "derived", "preprocessing", "duration_data.csv"),
  show_col_types = FALSE
) %>%
  mutate(ID = Sub_ID, session_norm = normalize_session(session)) %>%
  left_join(session_metadata, by = c("ID", "session_norm")) %>%
  mutate(session = coalesce(session.y, session.x)) %>%
  filter(!is.na(phase), Sex == "f") %>%
  mutate(
    fix = factor(
      fix,
      levels = c("fix_on_eyes", "fix_on_mouth", "fix_on_face", "fix_on_background"),
      labels = c("Eyes", "Mouth", "Face", "Background"),
      ordered = TRUE
    ),
    phase = factor(phase, levels = c("Follicular", "Luteal")),
    Group = factor(Group),
    medication = factor(
      medication,
      levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA")
    )
  ) %>%
  select(ID, session, Sex, Group, medication, phase, fix, percentage_fix_duration)

df_menst <- df %>%
  ungroup() %>%
  filter(!is.na(fix), fix %in% c("Eyes", "Background"))

summary_df <- df_menst %>%
  group_by(Group, fix, phase) %>%
  summarise(
    mean_var = mean(percentage_fix_duration, na.rm = TRUE),
    se = sd(percentage_fix_duration, na.rm = TRUE) / sqrt(n()),
    n_observations = n(),
    n_participants = n_distinct(ID),
    .groups = "drop"
  ) %>%
  mutate(
    lower = mean_var - se,
    upper = mean_var + se
  )

phase_values <- c("Follicular" = "#2F5D62", "Luteal" = "#B35C44")

base_phase_theme <- theme_classic(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(color = "black", size = 10.5, face = "bold"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.45),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9),
    legend.key.height = unit(3.0, "mm"),
    legend.key.width = unit(4.0, "mm"),
    legend.spacing.x = unit(1.4, "mm"),
    legend.key = element_blank(),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    panel.spacing.x = unit(5, "mm"),
    panel.spacing.y = unit(6, "mm")
  )

menst_phase_plot <- ggplot(
  summary_df,
  aes(x = Group, y = mean_var, color = phase, group = phase)
) +
  geom_line(linewidth = 0.9, position = position_dodge(width = 0.25)) +
  geom_point(size = 2.7, position = position_dodge(width = 0.25)) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.12,
    linewidth = 0.55,
    position = position_dodge(width = 0.25)
  ) +
  facet_wrap(~fix, scales = "fixed", nrow = 1) +
  labs(
    title = "Gaze by menstrual phase",
    y = "Fixation-duration proportion",
    x = NULL,
    color = "Menstrual phase"
  ) +
  scale_color_manual(values = phase_values) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.65), expand = expansion(mult = c(0.02, 0.06))) +
  base_phase_theme

write_csv(
  summary_df,
  here("02_outputs", "model_outputs", "supplementary", "menstrual_cycle", "08_phase_plot_data.csv")
)

participant_counts <- df_menst %>%
  distinct(ID, Group, phase) %>%
  count(Group, phase, name = "n_participants")

write_csv(
  participant_counts,
  here("02_outputs", "model_outputs", "supplementary", "menstrual_cycle", "09_phase_participant_counts.csv")
)

menst_phase_plot

ggsave(
  here("02_outputs", "figures", "supplementary", "menstrual_cycle", "svg", "menst_phases.svg"),
  menst_phase_plot,
  units = "mm",
  width = 180,
  height = 95
)
ggsave(
  here("02_outputs", "figures", "supplementary", "menstrual_cycle", "png", "Supp_Figure_1.png"),
  menst_phase_plot,
  width = 7.1,
  height = 3.75,
  dpi = 600
)
