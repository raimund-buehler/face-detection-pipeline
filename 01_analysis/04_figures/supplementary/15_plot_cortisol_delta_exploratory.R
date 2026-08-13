library(tidyverse)
library(here)

plot_dir <- here("02_outputs", "figures", "supplementary", "cortisol", "delta_exploratory")
data_dir <- here("02_outputs", "model_outputs", "main_manuscript", "cortisol", "delta_exploratory")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

group_values <- c("ASD" = "#D55E00", "CTRL" = "#0072B2")

make_delta_plot <- function(data, title, subtitle, x_label, y_label, out_file) {
  p <- ggplot(data, aes(delta_t2t3, delta_gaze, color = Group, shape = Sex)) +
    geom_hline(yintercept = 0, linetype = 3, color = "grey70") +
    geom_vline(xintercept = 0, linetype = 3, color = "grey70") +
    geom_point(size = 2.4, alpha = 0.9) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
    facet_wrap(~ fix, scales = "free_y") +
    scale_color_manual(values = group_values) +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = y_label,
      color = NULL,
      shape = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold")
    )

  ggsave(file.path(plot_dir, out_file), p, width = 9, height = 4.8, dpi = 300)
}

prepare_delta_df <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    filter(Sex %in% c("f", "m")) %>%
    mutate(
      Group = factor(Group, levels = c("ASD", "CTRL")),
      Sex = factor(Sex, levels = c("f", "m"), labels = c("Female", "Male")),
      fix = factor(fix, levels = c("Eyes", "Background"))
    )
}

nal_df <- prepare_delta_df(file.path(data_dir, "34_double_delta_nal_pla_with_covariates_dataset.csv"))
both_df <- prepare_delta_df(file.path(data_dir, "35_double_delta_both_pla_with_covariates_dataset.csv"))
both_clean_df <- prepare_delta_df(file.path(data_dir, "48_both_pla_with_covariates_dataset_extremes_removed.csv"))

make_delta_plot(
  data = nal_df,
  title = "NAL - PLA Double-Delta: Cortisol Change vs Gaze Change",
  subtitle = "Points are participants; lines are group-wise linear fits. Shapes denote sex.",
  x_label = "Delta T2-T3 cortisol (NAL - PLA)",
  y_label = "Delta fixation-duration proportion (NAL - PLA)",
  out_file = "Figure_5_nal_pla_double_delta_exploratory.png"
)

make_delta_plot(
  data = both_df,
  title = "BOTH - PLA Double-Delta: Cortisol Change vs Gaze Change",
  subtitle = "Points are participants; lines are group-wise linear fits. Shapes denote sex.",
  x_label = "Delta T2-T3 cortisol (BOTH - PLA)",
  y_label = "Delta fixation-duration proportion (BOTH - PLA)",
  out_file = "Figure_5_both_pla_double_delta_exploratory.png"
)

make_delta_plot(
  data = both_clean_df,
  title = "BOTH - PLA Double-Delta After Removing Influential Cases",
  subtitle = "Participants flagged by Cook's distance in either AOI were excluded. Shapes denote sex.",
  x_label = "Delta T2-T3 cortisol (BOTH - PLA)",
  y_label = "Delta fixation-duration proportion (BOTH - PLA)",
  out_file = "Figure_5_both_pla_double_delta_extremes_removed.png"
)
