# Manuscript section: Shared figure/report formatting helper
# Analysis family: shared helper
# Original source path: scripts/publication/format_markdown/format_anova.R
# Primary input dataset(s): ANOVA tables produced by manuscript analyses
# Primary output(s): formatted ANOVA text/tables for markdown reporting
# Known TODOs: verify which markdown reports still actively depend on this helper
# Scientific logic note: copied from source without changing scientific logic

library(scales)

format_anova <- function(anova_result, term, statistic = "full") {
  # Extract the relevant row for the specified term
  df1 <- anova_result[term, "NumDF"]
  df2 <- anova_result[term, "DenDF"]
  F_value <- anova_result[term, "F value"]
  p_value <- anova_result[term, "Pr(>F)"]
  
  # Use the scales package to format p-values
  formatted_p <- scales::pvalue(p_value)
  
  # Check if the formatted p-value contains "<" or ">", and adjust the output accordingly
  if (grepl("<|>", formatted_p)) {
    p_string <- paste0("_p_ ", formatted_p)  # No '=' if there is < or >
  } else {
    p_string <- paste0("_p_ = ", formatted_p)  # Include '=' if exact p-value
  }
  
  # Build the result based on the desired statistic
  if (statistic == "full") {
    # Italicize the F and return the full string
    formatted_string <- paste0("_F_(", round(df1, 0), ", ", round(df2, 0), 
                               ") = ", round(F_value, 2), ", ", p_string)
  } else if (statistic == "p") {
    # Return only the formatted p-value
    formatted_string <- p_string
  } else if (statistic == "F") {
    # Return only the F-statistic and degrees of freedom
    formatted_string <- paste0("_F_(", round(df1, 0), ", ", round(df2, 0), 
                               ") = ", round(F_value, 2))
  } else {
    stop("Invalid statistic type specified. Use 'full', 'p', or 'F'.")
  }
  
  return(formatted_string)
}
