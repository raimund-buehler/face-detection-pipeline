# Supplementary figure: T2-T3 cortisol reactivity and gaze robustness
# Shows the initial Group x T2-T3 association and its sensitivity to trimming.

library(tidyverse)
library(here)
library(scales)
library(svglite)

source(here("01_analysis", "shared", "model_output_utils.R"))

normalize_session <- function(x) {
  x %>%
    str_trim() %>%
    str_to_lower() %>%
    str_replace_all(" ", "_")
}

format_p_plot <- function(p) {
  ifelse(
    is.na(p), "NA",
    ifelse(p < .001, "< .001", sub("^0", "", sprintf("%.3f", p)))
  )
}

fig_dir <- here("02_outputs", "figures", "supplementary", "cortisol")
robustness_dir <- here(
  "02_outputs", "model_outputs", "main_manuscript", "cortisol",
  "duration_exploratory", "robustness"
)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(robustness_dir, recursive = TRUE, showWarnings = FALSE)

medication_levels_raw <- c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA")
medication_levels_plot <- c("PLA", "OXT", "NAL", "BOTH")

duration_eyes <- read_csv(
  here("00_data", "derived", "preprocessing", "duration_data.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    ID = Sub_ID,
    session_norm = normalize_session(session),
    fix,
    percentage_fix_duration
  ) %>%
  filter(fix == "fix_on_eyes")

df <- read_csv(
  here("00_data", "derived", "analysis", "df_cortisol_min_max.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    session_norm = normalize_session(session),
    medication = factor(
      medication,
      levels = medication_levels_raw,
      labels = c("BOTH", "NAL", "OXT", "PLA")
    ),
    medication = factor(medication, levels = medication_levels_plot),
    Group = factor(Group, levels = c("ASD", "CTRL")),
    ID = factor(ID),
    T2T3 = MinMax
  ) %>%
  left_join(duration_eyes, by = c("ID", "session_norm")) %>%
  filter(
    !is.na(percentage_fix_duration),
    !is.na(T2T3),
    !is.na(Group),
    !is.na(medication)
  )

t2t3_q <- quantile(df$T2T3, probs = c(.025, .05, .95, .975), na.rm = TRUE)

plot_df <- bind_rows(
  df %>%
    mutate(
      trim_model = "t2t3_raw_full",
      trim_panel = "Untrimmed"
    ),
  df %>%
    filter(T2T3 >= t2t3_q[["2.5%"]], T2T3 <= t2t3_q[["97.5%"]]) %>%
    mutate(
      trim_model = "t2t3_trim_2_5_97_5",
      trim_panel = "Trimmed outer 2.5%"
    ),
  df %>%
    filter(T2T3 >= t2t3_q[["5%"]], T2T3 <= t2t3_q[["95%"]]) %>%
    mutate(
      trim_model = "t2t3_trim_5_95",
      trim_panel = "Trimmed outer 5%"
    )
) %>%
  mutate(
    trim_panel = factor(
      trim_panel,
      levels = c("Untrimmed", "Trimmed outer 2.5%", "Trimmed outer 5%")
    )
  )

stats_df <- read_csv(
  file.path(robustness_dir, "06_interaction_robustness_summary.csv"),
  show_col_types = FALSE
) %>%
  filter(model %in% unique(plot_df$trim_model)) %>%
  transmute(
    trim_model = model,
    label = paste0("p = ", format_p_plot(interaction_p))
  )

annotation_df <- plot_df %>%
  distinct(trim_model, trim_panel) %>%
  left_join(stats_df, by = "trim_model") %>%
  mutate(
    x = min(df$T2T3, na.rm = TRUE) +
      .04 * diff(range(df$T2T3, na.rm = TRUE)),
    y = .62
  )

write_csv(
  plot_df %>%
    select(ID, Group, session, medication, T2T3, percentage_fix_duration, trim_model, trim_panel),
  file.path(robustness_dir, "07_t2t3_trim_gaze_plot_data.csv")
)

group_cols <- c("ASD" = "#D55E00", "CTRL" = "#0072B2")

p <- ggplot(plot_df, aes(x = T2T3, y = percentage_fix_duration, color = Group)) +
  geom_point(alpha = .42, size = 1.6, stroke = 0) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = .95,
    alpha = .16,
    show.legend = FALSE
  ) +
  geom_label(
    data = annotation_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.6,
    lineheight = .98,
    color = "grey15",
    fill = "white",
    linewidth = 0,
    label.padding = unit(.16, "lines")
  ) +
  facet_wrap(~ trim_panel, nrow = 1) +
  scale_color_manual(values = group_cols) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(.01, .03))
  ) +
  scale_x_continuous(expand = expansion(mult = c(.03, .03))) +
  coord_cartesian(ylim = c(0, .66)) +
  labs(
    title = "T2-T3 cortisol reactivity and eye gaze: trimming sensitivity",
    x = "T2-T3 cortisol change",
    y = "Eye fixation-duration proportion",
    color = "Group"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "grey15"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.width = unit(.8, "lines"),
    legend.key.height = unit(.8, "lines"),
    panel.spacing.x = unit(1.1, "lines"),
    plot.margin = margin(8, 10, 8, 10)
  ) +
  guides(
    color = guide_legend(
      override.aes = list(alpha = 1, size = 3)
    )
  )

ggsave(
  file.path(fig_dir, "Figure_5_t2t3_eye_gaze_trim_robustness.png"),
  p,
  width = 10.5,
  height = 3.8,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(fig_dir, "figure5_t2t3_eye_gaze_trim_robustness.svg"),
  p,
  width = 10.5,
  height = 3.8,
  bg = "white"
)
