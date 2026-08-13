# Manuscript section: Shared figure/report formatting helper
# Analysis family: shared helper
# Original source path: scripts/publication/format_markdown/format_emmeans.R
# Primary input dataset(s): emmeans outputs from manuscript analyses
# Primary output(s): formatted marginal means text/tables for markdown reporting
# Known TODOs: verify which markdown reports still actively depend on this helper
# Scientific logic note: copied from source without changing scientific logic

library(broom)
library(scales)
library(effectsize)

library(effectsize)

format_emmeans <- function(emmeans_result, contrast_lvl, at = NULL, statistic = "full", effect_size = NULL) {
  if (!"contrasts" %in% names(emmeans_result)) {
    stop("The provided object does not contain contrasts. Please use an emmeans object.")
  }
  
  tidy_emmeans <- tidy(emmeans_result$contrasts)
  
  if (!is.null(at)) {
    for (factor_name in names(at)) {
      factor_value <- at[[factor_name]]
      if (factor_name %in% colnames(tidy_emmeans)) {
        tidy_emmeans <- tidy_emmeans %>% filter(!!sym(factor_name) == factor_value)
      } else {
        stop(paste("Factor", factor_name, "not found in the data. Available factors:", paste(colnames(tidy_emmeans), collapse = ", ")))
      }
    }
  }
  
  result <- tidy_emmeans %>% filter(contrast == contrast_lvl)
  if (nrow(result) != 1) {
    stop("Please refine your filters. Available contrasts:", paste(unique(tidy_emmeans$contrast), collapse = ", "))
  }
  
  t_ratio <- result$statistic
  df <- result$df
  p_value <- result[[if ("p.value" %in% colnames(result)) "p.value" else "adj.p.value"]]
  formatted_p <- scales::pvalue(p_value)
  
  # Check if the formatted p-value contains "<" or ">", and adjust the output accordingly
  if (grepl("<|>", formatted_p)) {
    p_string <- paste0("_p_ ", formatted_p)  # No '=' if there is < or >
  } else {
    p_string <- paste0("_p_ = ", formatted_p)  # Include '=' if exact p-value
  }
  
  
  # Ensure effect_size has a default value if NULL
  effect_size <- ifelse(is.null(effect_size), "", effect_size)
  
  # Calculate effect size only if specified
  effect_size_string <- switch(effect_size,
                               "d" = paste0(", _d_ = ", round(t_to_d(t_ratio, df = df)$d, 2)),
                               "r" = paste0(", _r_ = ", round(t_to_r(t_ratio, df = df)$r, 2)),
                               "")
  
  # Helper function for formatting output
  format_output <- function(stat, t, df, p, es) {
    switch(stat,
           "full" = paste0("_t_(", round(df, 0), ") = ", round(t, 2), ", ", p_string, es),
           "p" = p_string,
           "t" = paste0("_t_(", round(df, 0), ") = ", round(t, 2), es),
           stop("Invalid statistic type specified. Use 'full', 'p', or 't'."))
  }
  
  # Return formatted output
  format_output(statistic, t_ratio, df, p_string, effect_size_string)
}




