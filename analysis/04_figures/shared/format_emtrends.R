# Manuscript section: Shared figure/report formatting helper
# Analysis family: shared helper
# Original source path: scripts/publication/format_markdown/format_emtrends.R
# Primary input dataset(s): emtrends outputs from manuscript analyses
# Primary output(s): formatted trend summaries for markdown reporting
# Known TODOs: verify which markdown reports still actively depend on this helper
# Scientific logic note: copied from source without changing scientific logic

library(broom)
library(scales)

format_emtrends <- function(emtrends_result, variable, at = NULL, statistic = "full", effect_size = NULL) {
  tidy_emtrends <- tidy(emtrends_result)
  
  trend_column <- paste0(variable, ".trend")
  if (!(trend_column %in% colnames(tidy_emtrends))) {
    stop(paste("Trend column", trend_column, "not found in the data. Available columns:", paste(colnames(tidy_emtrends), collapse = ", ")))
  }
  
  if (!is.null(at)) {
    for (factor_name in names(at)) {
      factor_value <- at[[factor_name]]
      if (factor_name %in% colnames(tidy_emtrends)) {
        tidy_emtrends <- tidy_emtrends %>% filter(!!sym(factor_name) == factor_value)
      } else {
        stop(paste("Factor", factor_name, "not found in the data. Available factors:", paste(colnames(tidy_emtrends), collapse = ", ")))
      }
    }
  }
  
  if (nrow(tidy_emtrends) != 1) {
    stop("Please refine your filters.")
  }
  
  t_ratio <- tidy_emtrends$statistic
  df <- tidy_emtrends$df
  p_value <- tidy_emtrends[[if ("p.value" %in% colnames(tidy_emtrends)) "p.value" else "adj.p.value"]]
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
  
  # Use the same helper function for formatting output
  format_output <- function(stat, t, df, p, es) {
    switch(stat,
           "full" = paste0("_t_(", round(df, 0), ") = ", round(t, 2), ", ", p_string, es),
           "p" = p_string,
           "t" = paste0("_t_(", round(df, 0), ") = ", round(t, 2), es),
           stop("Invalid statistic type specified. Use 'full', 'p', or 't'."))
  }
  
  # Return formatted output
  format_output(statistic, t_ratio, df, formatted_p, effect_size_string)
}

