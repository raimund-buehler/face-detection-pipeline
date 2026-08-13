# Manuscript section: Results organization
# Analysis family: model output summary
# Original source path: repository structural cleanup
# Primary input dataset(s): outputs/model_outputs/**/*.csv
# Primary output(s): outputs/model_outputs/model_results_overview.csv and outputs/model_outputs/model_results_overview.md
# Known TODOs: extend if additional model-output object types are added beyond CSV tables
# Scientific logic note: this script only aggregates saved model outputs; it does not refit or change any analyses

library(tidyverse)
library(here)

model_output_root <- here("02_outputs", "model_outputs")
output_files <- list.files(
  model_output_root,
  pattern = "[.]csv$",
  recursive = TRUE,
  full.names = TRUE
)

output_files <- output_files[!basename(output_files) %in% c("model_results_overview.csv")]

read_model_output <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  relative_path <- sub(paste0("^", model_output_root, "/?"), "", path)

  if (!"model" %in% names(df)) {
    df <- df %>% mutate(model = NA_character_)
  }

  if (!"term" %in% names(df)) {
    if ("test" %in% names(df)) {
      df <- df %>% rename(term = test)
    } else {
      df <- df %>% mutate(term = NA_character_)
    }
  }

  if (!"p_value" %in% names(df)) {
    df <- df %>% mutate(p_value = NA_real_)
  }

  if (!"p_display" %in% names(df)) {
    df <- df %>% mutate(p_display = NA_character_)
  }
  df <- df %>% mutate(p_display = as.character(p_display))

  if (!"sig" %in% names(df)) {
    df <- df %>% mutate(sig = NA_character_)
  }

  if (!"highlight" %in% names(df)) {
    df <- df %>% mutate(highlight = NA_character_)
  }

  if (!"significant_0_05" %in% names(df)) {
    df <- df %>% mutate(significant_0_05 = FALSE)
  }

  df %>%
    mutate(file = relative_path) %>%
    select(file, model, term, p_value, p_display, sig, highlight, everything())
}

overview <- purrr::map_dfr(output_files, read_model_output) %>%
  arrange(file, model, desc(significant_0_05 %||% FALSE), p_value)

readr::write_csv(
  overview,
  here("02_outputs", "model_outputs", "model_results_overview.csv")
)

overview_md <- overview %>%
  mutate(
    model = replace_na(model, ""),
    term = replace_na(term, ""),
    p_display = replace_na(p_display, ""),
    highlight = replace_na(highlight, ""),
    row_label = if_else(highlight == "significant", paste0("**", term, "**"), term)
  ) %>%
  select(file, model, row_label, p_display, highlight)

md_lines <- c(
  "# Model Results Overview",
  "",
  "| File | Model | Term | p | Highlight |",
  "| --- | --- | --- | --- | --- |",
  purrr::pmap_chr(
    overview_md,
    \(file, model, row_label, p_display, highlight) {
      paste0(
        "| ", file, " | ", model, " | ", row_label, " | ", p_display, " | ", highlight, " |"
      )
    }
  )
)

writeLines(md_lines, here("02_outputs", "model_outputs", "model_results_overview.md"))
