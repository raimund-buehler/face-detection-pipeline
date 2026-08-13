library(tidyverse)
library(ggdist)
library(ggh4x)
library(patchwork)
library(glmmTMB)
library(emmeans)
library(here)

output_dir <- here("02_outputs", "figures", "supplementary", "fixation_duration")
output_png <- file.path(output_dir, "Figure_2_fixation_duration_exclude_worst10_accuracy_sensitivity.png")
output_panelA_png <- file.path(output_dir, "Figure_2A_fixation_duration_exclude_worst10_accuracy_sensitivity.png")
output_csv <- here(
  "02_outputs", "model_outputs", "supplementary", "quality_checks",
  "10_figure2_duration_exclude_worst10_accuracy_contrasts.csv"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)

normalize_session <- function(x) {
  x %>% str_trim() %>% str_to_lower() %>% str_replace_all(" ", "_")
}

beta_squeeze <- function(x) {
  n <- sum(!is.na(x))
  ((x * (n - 1)) + 0.5) / n
}

format_sig <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

panel_levels <- c("Full Sample", "Female", "Male")
fix_levels <- c("Eyes", "Background")
med_levels <- c("BOTH", "NAL", "OXT", "PLA")
group_values <- c("ASD" = "#D55E00", "CTRL" = "#0072B2")
medication_values <- c("BOTH" = "#E69F00", "NAL" = "#CCB974", "OXT" = "#56B4A9", "PLA" = "#0B7A6B")
panel_values <- c("Full Sample" = "#303030", "Female" = "#B44E3A", "Male" = "#3B6C9E")

session_metadata <- read_csv(
  here("00_data", "derived", "analysis", "df_analysis_pub.csv"),
  show_col_types = FALSE
) %>%
  distinct(ID, session, Sex, Group, medication, accuracy_degrees) %>%
  mutate(session_norm = normalize_session(session))

worst10 <- session_metadata %>%
  filter(!is.na(accuracy_degrees)) %>%
  arrange(desc(accuracy_degrees)) %>%
  slice(1:10) %>%
  transmute(ID, session_norm)

df <- read_csv(here("00_data", "derived", "preprocessing", "duration_data.csv"), show_col_types = FALSE) %>%
  mutate(ID = Sub_ID, session_norm = normalize_session(session)) %>%
  left_join(session_metadata, by = c("ID", "session_norm")) %>%
  anti_join(worst10, by = c("ID", "session_norm")) %>%
  mutate(session = coalesce(session.y, session.x)) %>%
  filter(!is.na(Group), !is.na(medication), Sex %in% c("f", "m")) %>%
  mutate(
    fix = factor(
      fix,
      levels = c("fix_on_eyes", "fix_on_mouth", "fix_on_face", "fix_on_background"),
      labels = c("Eyes", "Mouth", "Face", "Background")
    ),
    medication = factor(
      medication,
      levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA"),
      labels = med_levels
    ),
    Group = factor(Group, levels = c("ASD", "CTRL")),
    session = factor(session),
    Sex = factor(Sex, levels = c("f", "m")),
    percentage_fix_duration_beta = beta_squeeze(percentage_fix_duration)
  ) %>%
  filter(fix %in% fix_levels)

model <- glmmTMB(
  percentage_fix_duration_beta ~ Group * fix * medication + Group * fix * Sex + session,
  family = beta_family(),
  data = df
)

get_contrasts <- function(by_vars, panel_fun, med_fun = NULL) {
  emm <- emmeans(model, as.formula(paste("~ Group |", paste(by_vars, collapse = " * "))))
  ctab <- summary(contrast(emm, method = "revpairwise")) %>% as_tibble()
  ctab %>%
    mutate(
      Panel = panel_fun(.),
      medication = if (is.null(med_fun)) "Overall" else med_fun(.)
    ) %>%
    transmute(
      Panel,
      fix,
      medication,
      contrast,
      estimate,
      std.error = SE,
      lower = estimate - 1.96 * SE,
      upper = estimate + 1.96 * SE,
      p.value
    )
}

overall_contrasts <- get_contrasts(
  by_vars = c("fix"),
  panel_fun = function(x) rep("Full Sample", nrow(x))
)

sex_contrasts <- get_contrasts(
  by_vars = c("fix", "Sex"),
  panel_fun = function(x) recode(as.character(x$Sex), "f" = "Female", "m" = "Male")
)

med_overall_contrasts <- get_contrasts(
  by_vars = c("fix", "medication"),
  panel_fun = function(x) rep("Full Sample", nrow(x)),
  med_fun = function(x) as.character(x$medication)
)

med_sex_contrasts <- get_contrasts(
  by_vars = c("fix", "medication", "Sex"),
  panel_fun = function(x) recode(as.character(x$Sex), "f" = "Female", "m" = "Male"),
  med_fun = function(x) as.character(x$medication)
)

contrast_df <- bind_rows(
  overall_contrasts,
  sex_contrasts,
  med_overall_contrasts,
  med_sex_contrasts
) %>%
  filter(fix %in% fix_levels) %>%
  mutate(
    Panel = factor(Panel, levels = panel_levels),
    fix = factor(fix, levels = fix_levels),
    medication = factor(medication, levels = c("Overall", med_levels)),
    medication_x = as.numeric(medication),
    significant = p.value < 0.05,
    sig = format_sig(p.value)
  )

write_csv(contrast_df, output_csv)

df_plot <- bind_rows(
  df %>% mutate(Panel = "Full Sample"),
  df %>% filter(Sex == "f") %>% mutate(Panel = "Female"),
  df %>% filter(Sex == "m") %>% mutate(Panel = "Male")
) %>%
  mutate(
    Panel = factor(Panel, levels = panel_levels),
    fix = factor(as.character(fix), levels = fix_levels),
    medication = factor(medication, levels = med_levels)
  )

summary_plot <- df_plot %>%
  group_by(Panel, fix, Group) %>%
  summarise(
    mean_prop = mean(percentage_fix_duration, na.rm = TRUE),
    sem = sd(percentage_fix_duration, na.rm = TRUE) / sqrt(n()),
    n_observations = n(),
    n_participants = n_distinct(ID),
    .groups = "drop"
  ) %>%
  mutate(
    lower = mean_prop - sem,
    upper = mean_prop + sem
  )

write_csv(
  summary_plot,
  here(
    "02_outputs", "model_outputs", "supplementary", "quality_checks",
    "11_figure2_duration_exclude_worst10_accuracy_plot_data.csv"
  )
)

sig_plot <- contrast_df %>%
  filter(medication == "Overall", significant) %>%
  left_join(
    summary_plot %>%
      group_by(Panel, fix) %>%
      summarise(y = max(upper, na.rm = TRUE) + 0.055, .groups = "drop"),
    by = c("Panel", "fix")
  ) %>%
  mutate(x = 1.5)

base_theme <- theme_classic(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(color = "black", size = 10.5, face = "bold"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0),
    panel.spacing.y = unit(6, "mm"),
    panel.spacing.x = unit(5, "mm"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.key = element_blank(),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10)
  )

figure_sens <- ggplot(
  summary_plot,
  aes(x = Group, y = mean_prop, color = Group, group = 1)
) +
  geom_line(color = "grey55", linewidth = 0.75) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.12, linewidth = 0.55) +
  geom_text(
    data = sig_plot,
    aes(x = x, y = y, label = sig),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  facet_grid(Panel ~ fix, switch = "y") +
  scale_color_manual(values = group_values) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.65), expand = expansion(mult = c(0.02, 0.06))) +
  labs(
    title = "Group gaze pattern after excluding 10 worst-accuracy sessions",
    x = NULL,
    y = "Fixation-duration proportion"
  ) +
  base_theme +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90),
    plot.title = element_text(size = 12, face = "bold", hjust = 0)
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 0.9, size = 2.8)))

ggsave(output_png, figure_sens, units = "mm", width = 170, height = 135, dpi = 600)
ggsave(output_panelA_png, figure_sens, units = "mm", width = 170, height = 135, dpi = 600)
