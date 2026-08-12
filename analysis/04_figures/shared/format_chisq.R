# Manuscript section: Shared figure/report formatting helper
# Analysis family: shared helper
# Original source path: scripts/publication/format_markdown/format_chisq.R
# Primary input dataset(s): chi-square test outputs
# Primary output(s): formatted chi-square text/tables for markdown reporting
# Known TODOs: verify whether this helper is still used in the cleaned manuscript workflow
# Scientific logic note: copied from source without changing scientific logic

library(scales)

format_chisq <- function(chisq_result, statistic = "full") {
  # Extract chi-square value, degrees of freedom, and p-value
  X_squared <- chisq_result$statistic
  df <- chisq_result$parameter
  p_value <- chisq_result$p.value
  
  # Use the scales package to format p-values
  formatted_p <- scales::pvalue(p_value)
  
  # Check if the formatted p-value contains "<" or ">", and adjust the output accordingly
  if (grepl("<|>", formatted_p)) {
    p_string <- paste0("_p_ ", formatted_p)  # No '=' if there is < or >
  } else {
    p_string <- paste0("_p_ = ", formatted_p)  # Include '=' if exact p-value
  }
  
  # Use Greek letter χ and superscript ² for chi-squared
  chi_symbol <- "\u03c7"  # Greek chi
  squared_symbol <- "\u00b2"  # Superscript 2
  
  # Build the result based on the desired statistic
  if (statistic == "full") {
    # Format the string with proper symbols
    formatted_string <- paste0(chi_symbol, squared_symbol, "(", round(df, 0), ") = ", round(X_squared, 2), ", ", p_string)
  } else if (statistic == "p") {
    # Return only the formatted p-value
    formatted_string <- p_string
  } else if (statistic == "X2") {
    # Return only the chi-square statistic and degrees of freedom
    formatted_string <- paste0(chi_symbol, "^", squared_symbol, "(", round(df, 0), ") = ", round(X_squared, 2))
  } else {
    stop("Invalid statistic type specified. Use 'full', 'p', or 'X2'.")
  }
  
  return(formatted_string)
}
