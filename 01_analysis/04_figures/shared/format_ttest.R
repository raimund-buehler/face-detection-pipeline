# Manuscript section: Shared figure/report formatting helper
# Analysis family: shared helper
# Original source path: scripts/publication/format_markdown/format_ttest.R
# Primary input dataset(s): t-test outputs from manuscript analyses
# Primary output(s): formatted t-test text/tables for markdown reporting
# Known TODOs: verify whether this helper is still used in the cleaned manuscript workflow
# Scientific logic note: copied from source without changing scientific logic

# Function to format t-test objects
format_ttest <- function(t_test_result, statistic = "full") {
  # Extract t value, degrees of freedom, and p-value
  t_value <- t_test_result$statistic
  df <- t_test_result$parameter
  p_value <- t_test_result$p.value
  
  # Use the scales package to format p-values
  library(scales)
  formatted_p <- scales::pvalue(p_value)
  
  # Check if the formatted p-value contains "<" or ">", and adjust the output accordingly
  if (grepl("<|>", formatted_p)) {
    p_string <- paste0("_p_ ", formatted_p)  # No '=' if there is < or >
  } else {
    p_string <- paste0("_p_ = ", formatted_p)  # Include '=' if exact p-value
  }
  
  # Build the result based on the desired statistic
  if (statistic == "full") {
    # Italicize t and return the full string
    formatted_string <- paste0("_t_(", round(df, 2), ") = ", round(t_value, 2), ", ", p_string)
  } else if (statistic == "p") {
    # Return only the formatted p-value
    formatted_string <- p_string
  } else if (statistic == "t") {
    # Return only the t-statistic and degrees of freedom
    formatted_string <- paste0("_t_(", round(df, 2), ") = ", round(t_value, 2))
  } else {
    stop("Invalid statistic type specified. Use 'full', 'p', or 't'.")
  }
  
  return(formatted_string)
}
