library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggdist)
library(ggh4x)
library(patchwork)
library(here)

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

output_dir <- here("02_outputs", "figures", "supplementary", "fixation_rate_analogs_current")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_png <- file.path(output_dir, "Figure_S_fixrate_2_current_layout.png")

df <- read_csv(here("00_data", "derived", "analysis", "df_analysis_pub.csv"), show_col_types = FALSE) %>%
  filter(Sex %in% c("f", "m")) %>%
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
    ID = factor(ID),
    session = factor(session),
    Sex = factor(Sex, levels = c("f", "m"))
  ) %>%
  filter(fix %in% fix_levels, !is.na(sqrt(rate_all))) %>%
  mutate(sqrt_rate_all = sqrt(rate_all))

model <- lmer(sqrt_rate_all ~ Group * fix * medication * Sex + session + (1 | ID), data = df)

overall_contrasts <- emmeans(model, ~ Group | fix, type = "response") %>%
  contrast(method = "revpairwise") %>%
  summary(infer = TRUE) %>%
  as_tibble() %>%
  mutate(Panel = "Full Sample", medication = "Overall")

sex_contrasts <- emmeans(model, ~ Group | fix * Sex, type = "response") %>%
  contrast(method = "revpairwise") %>%
  summary(infer = TRUE) %>%
  as_tibble() %>%
  mutate(
    Panel = recode(as.character(Sex), "f" = "Female", "m" = "Male"),
    medication = "Overall"
  )

med_overall_contrasts <- emmeans(model, ~ Group | fix * medication, type = "response") %>%
  contrast(method = "revpairwise") %>%
  summary(infer = TRUE) %>%
  as_tibble() %>%
  mutate(Panel = "Full Sample")

med_sex_contrasts <- emmeans(model, ~ Group | fix * medication * Sex, type = "response") %>%
  contrast(method = "revpairwise") %>%
  summary(infer = TRUE) %>%
  as_tibble() %>%
  mutate(Panel = recode(as.character(Sex), "f" = "Female", "m" = "Male"))

contrast_df <- bind_rows(
  overall_contrasts %>% select(Panel, fix, medication, contrast, estimate, SE, p.value),
  sex_contrasts %>% select(Panel, fix, medication, contrast, estimate, SE, p.value),
  med_overall_contrasts %>% select(Panel, fix, medication, contrast, estimate, SE, p.value),
  med_sex_contrasts %>% select(Panel, fix, medication, contrast, estimate, SE, p.value)
) %>%
  mutate(
    Panel = factor(Panel, levels = panel_levels),
    fix = factor(as.character(fix), levels = fix_levels),
    medication = factor(as.character(medication), levels = c("Overall", med_levels)),
    lower = estimate - 1.96 * SE,
    upper = estimate + 1.96 * SE,
    medication_x = as.numeric(medication),
    significant = p.value < 0.05,
    sig = format_sig(p.value)
  )

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

overall_positions <- tibble(
  Group = factor(c("ASD", "CTRL"), levels = c("ASD", "CTRL")),
  x = c(0.88, 1.12)
)
df_overall <- df_plot %>% left_join(overall_positions, by = "Group")
overall_summary <- df_overall %>%
  group_by(Panel, fix, Group, x) %>%
  summarise(mean_rate = mean(sqrt_rate_all, na.rm = TRUE), .groups = "drop")

med_positions <- expand_grid(
  medication = factor(med_levels, levels = med_levels),
  Group = factor(c("ASD", "CTRL"), levels = c("ASD", "CTRL"))
) %>%
  mutate(
    med_x = as.numeric(medication),
    x = med_x + if_else(Group == "ASD", -0.13, 0.13)
  )
df_med <- df_plot %>% left_join(med_positions, by = c("medication", "Group"))
med_summary <- df_med %>%
  group_by(Panel, fix, medication, Group, x) %>%
  summarise(mean_rate = mean(sqrt_rate_all, na.rm = TRUE), .groups = "drop")

panel_a_sig <- contrast_df %>%
  filter(medication == "Overall", significant) %>%
  left_join(
    overall_summary %>% group_by(Panel, fix) %>% summarise(y = max(mean_rate, na.rm = TRUE) + 0.22, .groups = "drop"),
    by = c("Panel", "fix")
  ) %>%
  mutate(x = 1)

panel_b_sig <- contrast_df %>%
  filter(significant) %>%
  mutate(
    panel_offset = recode(as.character(Panel), "Full Sample" = -0.183, "Female" = 0, "Male" = 0.183),
    x_star = medication_x + panel_offset - 0.025,
    y = if_else(estimate >= 0, upper + 0.22, lower - 0.22)
  )

panel_c_sig <- contrast_df %>%
  filter(medication != "Overall", significant) %>%
  left_join(
    med_summary %>% group_by(Panel, fix, medication) %>% summarise(y = max(mean_rate, na.rm = TRUE) + 0.22, .groups = "drop"),
    by = c("Panel", "fix", "medication")
  ) %>%
  mutate(x = as.numeric(factor(as.character(medication), levels = med_levels)) - 0.01)

base_theme <- theme_classic(base_size = 8) +
  theme(
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(color = "black", size = 8.5),
    plot.title = element_text(size = 10, face = "bold", hjust = 0),
    panel.spacing.y = unit(5, "mm"),
    panel.spacing.x = unit(5, "mm"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 6.5),
    legend.key = element_blank()
  )

mean_bar_width <- 0.075

panel_a <- ggplot() +
  stat_halfeye(
    data = df_overall %>% filter(Group == "ASD"),
    aes(x = x, y = sqrt_rate_all, fill = Group),
    side = "left", adjust = 1.35, width = 0.62, .width = 0,
    point_color = NA, alpha = 0.46, slab_color = NA, normalize = "panels", scale = 0.62
  ) +
  stat_halfeye(
    data = df_overall %>% filter(Group == "CTRL"),
    aes(x = x, y = sqrt_rate_all, fill = Group),
    side = "right", adjust = 1.35, width = 0.62, .width = 0,
    point_color = NA, alpha = 0.46, slab_color = NA, normalize = "panels", scale = 0.62
  ) +
  geom_segment(
    data = overall_summary,
    aes(x = x - mean_bar_width, xend = x + mean_bar_width, y = mean_rate, yend = mean_rate, color = Group),
    linewidth = 0.9, lineend = "round"
  ) +
  geom_text(data = panel_a_sig, aes(x = x, y = y, label = sig), inherit.aes = FALSE, size = 3.8, fontface = "bold") +
  facet_grid2(Panel ~ fix, axes = "x", remove_labels = "x", switch = "y") +
  scale_x_continuous(breaks = c(0.88, 1.12), labels = c("ASD", "CTRL"), expand = expansion(mult = c(0.25, 0.25))) +
  scale_color_manual(values = group_values) +
  scale_fill_manual(values = group_values) +
  coord_cartesian(xlim = c(0.42, 1.58)) +
  labs(title = "Overall Group Differences", x = NULL, y = "sqrt(Fixations/s)") +
  base_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_line(linewidth = 0.25)) +
  guides(color = guide_legend(override.aes = list(fill = group_values, color = group_values, alpha = 1, shape = 15, size = 2.4)), fill = "none")

panel_dodge <- position_dodge(width = 0.55)
panel_b <- ggplot(contrast_df, aes(x = medication_x, y = estimate, color = Panel, shape = significant)) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.25, linetype = "22") +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = panel_dodge, width = 0.14, linewidth = 0.4) +
  geom_point(position = panel_dodge, size = 1.8, stroke = 0.7, fill = "white") +
  geom_text(data = panel_b_sig, aes(x = x_star, y = y, label = sig, color = Panel), inherit.aes = FALSE, size = 3.0, fontface = "bold", show.legend = FALSE) +
  facet_grid(. ~ fix) +
  scale_x_continuous(breaks = seq_along(levels(contrast_df$medication)), labels = levels(contrast_df$medication), expand = expansion(mult = c(0.07, 0.07))) +
  scale_color_manual(values = panel_values) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), guide = "none") +
  labs(title = "Group Contrasts (Exploratory)", x = NULL, y = "ASD - CTRL Contrast") +
  base_theme +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  guides(color = guide_legend(override.aes = list(shape = 16, size = 2.2)))

panel_c <- ggplot() +
  stat_halfeye(
    data = df_med %>% filter(Group == "ASD"),
    aes(x = x, y = sqrt_rate_all, fill = medication),
    side = "left", adjust = 1.35, width = 0.42, .width = 0,
    point_color = NA, alpha = 0.42, slab_color = NA, normalize = "panels", scale = 0.52
  ) +
  stat_halfeye(
    data = df_med %>% filter(Group == "CTRL"),
    aes(x = x, y = sqrt_rate_all, fill = medication),
    side = "right", adjust = 1.35, width = 0.42, .width = 0,
    point_color = NA, alpha = 0.42, slab_color = NA, normalize = "panels", scale = 0.52
  ) +
  geom_segment(
    data = med_summary,
    aes(x = x - mean_bar_width, xend = x + mean_bar_width, y = mean_rate, yend = mean_rate, color = medication),
    linewidth = 0.85, lineend = "round"
  ) +
  geom_text(data = panel_c_sig, aes(x = x, y = y, label = sig, color = medication), inherit.aes = FALSE, size = 3.5, fontface = "bold", show.legend = FALSE) +
  facet_grid2(Panel ~ fix, axes = "x", remove_labels = "x", switch = "y") +
  scale_x_continuous(
    breaks = rep(seq_along(med_levels), each = 2) + rep(c(-0.13, 0.13), times = length(med_levels)),
    labels = rep(c("ASD", "CTRL"), times = length(med_levels)),
    minor_breaks = seq_along(med_levels),
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  scale_color_manual(values = medication_values) +
  scale_fill_manual(values = medication_values) +
  labs(title = "Medication-Specific Differences (Exploratory)", x = NULL, y = "sqrt(Fixations/s)") +
  base_theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.ticks.x = element_line(linewidth = 0.25)) +
  guides(color = "none", fill = guide_legend(override.aes = list(alpha = 1, shape = 15, size = 2.4, colour = NA, linewidth = 0), nrow = 1))

figure <- (panel_a | panel_b) / panel_c +
  plot_layout(heights = c(0.83, 1.17), widths = c(0.40, 0.60), guides = "keep") +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 11, face = "plain"), plot.tag.position = c(0.005, 0.995))
  )

ggsave(output_png, figure, units = "mm", width = 183, height = 220, dpi = 300)
