library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggh4x)
library(patchwork)
library(here)

output_dir <- here("02_outputs", "figures", "supplementary", "fixation_rate_analogs_current")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_png <- file.path(output_dir, "Figure_S_fixrate_4_current_layout.png")

group_values <- c("ASD" = "#D55E00", "CTRL" = "#0072B2")
aoi_values <- c("Eyes" = "#4A4A4A", "Background" = "#4A4A4A")
panel_levels <- c("Full Sample", "Female")
fix_levels <- c("Eyes", "Background")

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
      labels = c("BOTH", "NAL", "OXT", "PLA")
    ),
    Sex = factor(Sex, levels = c("f", "m")),
    Group = factor(Group, levels = c("ASD", "CTRL")),
    ID = factor(ID),
    session = factor(session),
    AQ = .data[["AQ-K Score"]],
    SIAS = .data[["SIAS Score"]],
    AQ_z = as.numeric(scale(.data[["AQ-K Score"]]))
  ) %>%
  filter(!is.na(AQ), !is.na(SIAS), fix %in% fix_levels, !is.na(rate_all)) %>%
  mutate(sqrt_rate_all = sqrt(rate_all))

aq_center <- attr(scale(df$AQ), "scaled:center")
aq_scale <- attr(scale(df$AQ), "scaled:scale")

aq_model <- lmer(sqrt_rate_all ~ AQ_z * fix * medication + AQ_z * fix * Sex + session + (1 | ID), data = df)

aq_raw_range <- seq(floor(min(df$AQ, na.rm = TRUE)), ceiling(max(df$AQ, na.rm = TRUE)), length.out = 80)
aq_range <- (aq_raw_range - aq_center) / aq_scale

pred_full <- emmeans(aq_model, ~ AQ_z | fix, at = list(AQ_z = aq_range, fix = fix_levels)) %>%
  summary(infer = c(TRUE, TRUE)) %>%
  as_tibble() %>%
  mutate(Panel = "Full Sample")

pred_female <- emmeans(aq_model, ~ AQ_z | fix, at = list(AQ_z = aq_range, fix = fix_levels, Sex = "f")) %>%
  summary(infer = c(TRUE, TRUE)) %>%
  as_tibble() %>%
  mutate(Panel = "Female")

pred_df <- bind_rows(pred_full, pred_female) %>%
  filter(fix %in% fix_levels)

if ("lower.CL" %in% names(pred_df)) {
  pred_df <- pred_df %>% mutate(lower = lower.CL)
} else {
  pred_df <- pred_df %>% mutate(lower = asymp.LCL)
}
if ("upper.CL" %in% names(pred_df)) {
  pred_df <- pred_df %>% mutate(upper = upper.CL)
} else {
  pred_df <- pred_df %>% mutate(upper = asymp.UCL)
}
if ("emmean" %in% names(pred_df)) {
  pred_df <- pred_df %>% mutate(fit = emmean)
} else {
  pred_df <- pred_df %>% mutate(fit = response)
}

pred_df <- pred_df %>%
  mutate(
    AQ = AQ_z * aq_scale + aq_center,
    Panel = factor(Panel, levels = panel_levels),
    fix = factor(fix, levels = fix_levels)
  )

points_df <- bind_rows(
  df %>% mutate(Panel = "Full Sample"),
  df %>% filter(Sex == "f") %>% mutate(Panel = "Female")
) %>%
  filter(fix %in% fix_levels) %>%
  group_by(Panel, ID, fix, AQ) %>%
  summarise(sqrt_rate_all = mean(sqrt_rate_all, na.rm = TRUE), .groups = "drop") %>%
  mutate(Panel = factor(Panel, levels = panel_levels), fix = factor(fix, levels = fix_levels))

slope_stats <- bind_rows(
  emtrends(aq_model, ~ fix, var = "AQ_z") %>% summary(infer = TRUE) %>% as_tibble() %>% mutate(Panel = "Full Sample"),
  emtrends(aq_model, ~ fix, var = "AQ_z", at = list(Sex = "f")) %>% summary(infer = TRUE) %>% as_tibble() %>% mutate(Panel = "Female")
) %>%
  mutate(
    label = case_when(
      p.value < 0.001 ~ sprintf("z = %.2f, p < .001", z.ratio),
      TRUE ~ sprintf("z = %.2f, p = %s", z.ratio, sub("^0", "", sprintf("%.3f", p.value)))
    ),
    Panel = factor(Panel, levels = panel_levels),
    fix = factor(fix, levels = fix_levels),
    x = min(aq_raw_range) + 0.08 * diff(range(aq_raw_range)),
    y = max(points_df$sqrt_rate_all, na.rm = TRUE) * 0.93
  )

panel_a <- ggplot() +
  geom_point(data = points_df, aes(x = AQ, y = sqrt_rate_all, color = fix), alpha = 0.22, size = 0.9, position = position_jitter(width = 0.035, height = 0)) +
  geom_ribbon(data = pred_df, aes(x = AQ, ymin = lower, ymax = upper, fill = fix), alpha = 0.18, color = NA) +
  geom_line(data = pred_df, aes(x = AQ, y = fit, color = fix), linewidth = 0.85) +
  geom_text(data = slope_stats, aes(x = x, y = y, label = label, color = fix), hjust = 0, size = 2.45, show.legend = FALSE) +
  facet_grid2(Panel ~ fix, axes = "all", remove_labels = "x", switch = "y") +
  scale_color_manual(values = aoi_values) +
  scale_fill_manual(values = aoi_values) +
  coord_cartesian(xlim = range(aq_raw_range), clip = "off") +
  labs(title = "AQ Simple Slopes", x = "AQ", y = "sqrt(Fixations/s)") +
  theme_classic(base_size = 9) +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 9),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    panel.spacing = unit(0.25, "lines")
  )

df_cor <- df %>%
  distinct(ID, Group, AQ, SIAS, Sex) %>%
  mutate(Sex_symbol = recode(Sex, "f" = "\u2640", "m" = "\u2642"), Sex_label = factor(recode(Sex, "f" = "Female", "m" = "Male"), levels = c("Female", "Male")))

corr_test <- cor.test(df_cor$AQ, df_cor$SIAS)

panel_b <- ggplot(df_cor, aes(x = AQ, y = SIAS)) +
  geom_smooth(method = "lm", se = TRUE, color = "#2E7D6E", fill = "#2E7D6E", alpha = 0.20, linewidth = 1.25) +
  geom_point(aes(color = Group), size = 0, alpha = 0) +
  geom_text(aes(label = Sex_symbol, color = Group), size = 16/.pt, alpha = 0.85, show.legend = FALSE) +
  annotate(
    "text",
    x = min(df_cor$AQ, na.rm = TRUE) + 0.03 * diff(range(df_cor$AQ, na.rm = TRUE)),
    y = max(df_cor$SIAS, na.rm = TRUE) - 0.04 * diff(range(df_cor$SIAS, na.rm = TRUE)),
    hjust = 0, vjust = 1,
    label = sprintf("r = %.2f\n95%% CI [%.2f, %.2f]\np < .001", corr_test$estimate, corr_test$conf.int[1], corr_test$conf.int[2]),
    size = 5.4/.pt, color = "#333333"
  ) +
  scale_color_manual(values = group_values) +
  labs(title = "AQ and Social Anxiety", x = "AQ", y = "SIAS") +
  theme_classic(base_size = 9) +
  theme(
    aspect.ratio = 1,
    legend.position = "top",
    legend.title = element_blank(),
    strip.background = element_blank(),
    plot.title = element_text(size = 11, face = "bold", hjust = 0)
  ) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 7, alpha = 1)))

figure <- panel_a + panel_b +
  plot_layout(widths = c(1.05, 0.95), guides = "collect") +
  plot_annotation(
    title = "Dimensional Trait Analysis (Exploratory)",
    tag_levels = "A",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0), plot.tag = element_text(size = 16, face = "plain"))
  )

ggsave(output_png, figure, units = "mm", width = 183, height = 120, dpi = 300)
