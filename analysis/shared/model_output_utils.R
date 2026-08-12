significance_stars <- function(p_values) {
  case_when(
    is.na(p_values) ~ "",
    p_values < 0.001 ~ "***",
    p_values < 0.01 ~ "**",
    p_values < 0.05 ~ "*",
    p_values < 0.1 ~ ".",
    TRUE ~ ""
  )
}

format_p_value <- function(p_values) {
  case_when(
    is.na(p_values) ~ NA_character_,
    p_values < 0.001 ~ "<0.001",
    TRUE ~ formatC(p_values, format = "f", digits = 3)
  )
}

prepare_anova_table <- function(model, model_name) {
  anova_df <- anova(model) %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    as_tibble()
  prepare_anova_table_from_df(anova_df, model_name)
}

prepare_anova_table_from_df <- function(anova_df, model_name) {
  p_col <- intersect(c("Pr(>F)", "Pr(>Chisq)", "Pr(>Chi)", "p.value"), names(anova_df))
  if (length(p_col) == 0) {
    anova_df <- anova_df %>% mutate(p_value = NA_real_)
  } else {
    anova_df <- anova_df %>% rename(p_value = all_of(p_col[[1]]))
  }

  anova_df %>%
    mutate(
      model = model_name,
      p_display = format_p_value(p_value),
      sig = significance_stars(p_value),
      significant_0_05 = !is.na(p_value) & p_value < 0.05,
      highlight = case_when(
        significant_0_05 ~ "significant",
        !is.na(p_value) & p_value < 0.1 ~ "trend",
        TRUE ~ ""
      ),
      p_display = if_else(sig == "", p_display, paste0(p_display, " ", sig))
    ) %>%
    relocate(model, term, p_value, p_display, sig, significant_0_05, highlight)
}

write_model_anova <- function(model, output_dir, filename, model_name) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  table <- prepare_anova_table(model, model_name)
  readr::write_csv(table, file.path(output_dir, filename))
  table
}

write_anova_table_from_df <- function(anova_df, output_dir, filename, model_name) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  table <- prepare_anova_table_from_df(anova_df, model_name)
  readr::write_csv(table, file.path(output_dir, filename))
  table
}

write_lm_coefficients <- function(model, output_dir, filename, model_name) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  coef_table <- summary(model)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    as_tibble()

  p_col <- intersect(c("Pr(>|t|)", "Pr(>|z|)", "p.value"), names(coef_table))
  if (length(p_col) == 0) {
    coef_table <- coef_table %>% mutate(p_value = NA_real_)
  } else {
    coef_table <- coef_table %>% rename(p_value = all_of(p_col[[1]]))
  }

  coef_table <- coef_table %>%
    mutate(
      model = model_name,
      p_display = format_p_value(p_value),
      sig = significance_stars(p_value),
      significant_0_05 = !is.na(p_value) & p_value < 0.05,
      highlight = case_when(
        significant_0_05 ~ "significant",
        !is.na(p_value) & p_value < 0.1 ~ "trend",
        TRUE ~ ""
      ),
      p_display = if_else(sig == "", p_display, paste0(p_display, " ", sig))
    ) %>%
    relocate(model, term, p_value, p_display, sig, significant_0_05, highlight)

  readr::write_csv(coef_table, file.path(output_dir, filename))
  coef_table
}

write_htest_result <- function(test_result, output_dir, filename, test_name) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  estimate_value <- if (!is.null(test_result$estimate)) unname(test_result$estimate[[1]]) else NA_real_
  statistic_value <- if (!is.null(test_result$statistic)) unname(test_result$statistic[[1]]) else NA_real_
  parameter_value <- if (!is.null(test_result$parameter)) unname(test_result$parameter[[1]]) else NA_real_
  conf_low <- if (!is.null(test_result$conf.int)) test_result$conf.int[1] else NA_real_
  conf_high <- if (!is.null(test_result$conf.int)) test_result$conf.int[2] else NA_real_
  p_value <- test_result$p.value
  method_name <- test_result$method
  p_display <- format_p_value(p_value)
  sig <- significance_stars(p_value)
  significant_0_05 <- !is.na(p_value) & p_value < 0.05
  highlight <- case_when(
    significant_0_05 ~ "significant",
    !is.na(p_value) & p_value < 0.1 ~ "trend",
    TRUE ~ ""
  )
  p_display_label <- if_else(sig == "", p_display, paste0(p_display, " ", sig))

  result <- tibble(
    test = test_name,
    statistic = statistic_value,
    parameter = parameter_value,
    p_value = p_value,
    p_display = p_display_label,
    estimate = estimate_value,
    conf_low = conf_low,
    conf_high = conf_high,
    method = method_name,
    sig = sig,
    significant_0_05 = significant_0_05,
    highlight = highlight
  )

  readr::write_csv(result, file.path(output_dir, filename))
  result
}

write_text_output <- function(lines, output_dir, filename) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, file.path(output_dir, filename))
}
