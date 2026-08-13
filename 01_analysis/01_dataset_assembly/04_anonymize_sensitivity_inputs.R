# Manuscript section: Sensitivity-analysis data preparation
# Analysis family: private ID anonymization helper
# Primary input dataset(s): configured eye-gaze, medication, and disorder CSV files
# Primary output(s): anonymized CSV files using a shared numeric participant key
#
# This script is intentionally path-driven. Keep the lookup table outside the
# repository, or in a gitignored private directory, because it contains the
# original ID to anonymous ID crosswalk.

library(tidyverse)
library(here)

private_config_paths <- c(
  here("00_data", "private", "anonymization_config.R"),
  here("00_data", "private", "anonymization_config_template.R")
)
private_config_path <- private_config_paths[file.exists(private_config_paths)][1]
if (!is.na(private_config_path)) {
  source(private_config_path)
}

parse_cli_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) {
      next
    }
    key <- sub("^--", "", sub("=.*$", "", arg))
    value <- sub("^[^=]*=", "", arg)
    out[[key]] <- value
  }
  out
}

split_paths <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(character())
  }

  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(character())
  }

  if (length(x) > 1) {
    return(str_trim(x) %>% discard(~ !nzchar(.x)))
  }

  if (!nzchar(x)) {
    return(character())
  }

  str_split(x, ",", simplify = FALSE)[[1]] %>%
    str_trim() %>%
    discard(~ !nzchar(.x))
}

is_absolute_or_home_path <- function(path) {
  startsWith(path, .Platform$file.sep) || startsWith(path, "~")
}

resolve_path <- function(path) {
  if (is_absolute_or_home_path(path)) {
    return(path.expand(path))
  }

  here(path)
}

coalesce_config <- function(cli_value, env_name, config_name, default = "") {
  if (!is.null(cli_value) && !is.na(cli_value) && nzchar(cli_value)) {
    return(cli_value)
  }

  env_value <- Sys.getenv(env_name, unset = "")
  if (nzchar(env_value)) {
    return(env_value)
  }

  if (exists(config_name, inherits = TRUE)) {
    config_value <- get(config_name, inherits = TRUE)
    if (!is.function(config_value) && !is.null(config_value) && length(config_value) > 0) {
      return(config_value)
    }
  }

  default
}

read_csv_auto <- function(path) {
  first_line <- readLines(path, n = 1, warn = FALSE)
  delimiter <- if (length(first_line) > 0 && str_count(first_line, ";") > str_count(first_line, ",")) {
    ";"
  } else {
    ","
  }

  data <- read_delim(path, delim = delimiter, show_col_types = FALSE, trim_ws = TRUE)
  names(data) <- str_trim(names(data))
  data
}

first_existing_id_column <- function(data, candidates = c("ID", "Id", "id", "Sub_ID", "participant_id", "participant", "subject_id")) {
  id_col <- intersect(candidates, names(data))
  if (length(id_col) == 0) {
    stop(
      "Could not find an ID column. Expected one of: ",
      paste(candidates, collapse = ", "),
      call. = FALSE
    )
  }
  id_col[[1]]
}

read_ids <- function(path, id_candidates) {
  data <- read_csv_auto(path)
  id_col <- first_existing_id_column(data, id_candidates)
  data %>%
    transmute(original_id = as.character(.data[[id_col]])) %>%
    filter(!is.na(original_id), nzchar(original_id)) %>%
    distinct()
}

create_id_lookup <- function(paths, lookup_path, id_candidates) {
  lookup <- paths %>%
    map_dfr(read_ids, id_candidates = id_candidates) %>%
    distinct(original_id) %>%
    arrange(original_id) %>%
    mutate(anon_id = row_number()) %>%
    select(anon_id, original_id)

  dir.create(dirname(lookup_path), recursive = TRUE, showWarnings = FALSE)
  write_csv(lookup, lookup_path)
  lookup
}

load_or_create_lookup <- function(paths, lookup_path, id_candidates, create_if_missing) {
  if (file.exists(lookup_path)) {
    return(read_csv_auto(lookup_path) %>%
      mutate(
        anon_id = as.integer(anon_id),
        original_id = as.character(original_id)
      ))
  }

  if (!isTRUE(create_if_missing)) {
    stop(
      "Lookup table does not exist and create_lookup is FALSE: ",
      lookup_path,
      call. = FALSE
    )
  }

  create_id_lookup(paths, lookup_path, id_candidates)
}

anonymize_table <- function(input_path, output_path, lookup, id_candidates) {
  data <- read_csv_auto(input_path)
  id_col <- first_existing_id_column(data, id_candidates)

  anonymized <- data %>%
    mutate(original_id = as.character(.data[[id_col]])) %>%
    left_join(lookup, by = "original_id") %>%
    relocate(anon_id, .before = 1)

  missing_lookup <- anonymized %>%
    filter(is.na(anon_id), !is.na(original_id), nzchar(original_id)) %>%
    distinct(original_id) %>%
    nrow()

  if (missing_lookup > 0) {
    stop(
      "Found ", missing_lookup, " IDs in ", basename(input_path),
      " that are absent from the lookup table. Recreate or update the lookup before anonymizing.",
      call. = FALSE
    )
  }

  anonymized <- anonymized %>%
    select(-all_of(id_col), -original_id)

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  write_csv(anonymized, output_path)

  tibble(
    input = input_path,
    output = output_path,
    rows = nrow(anonymized),
    unique_anon_ids = n_distinct(anonymized$anon_id, na.rm = TRUE)
  )
}

make_output_path <- function(input_path, output_dir) {
  file.path(output_dir, paste0(tools::file_path_sans_ext(basename(input_path)), "_anonymized.csv"))
}

args <- parse_cli_args()

lookup_path <- coalesce_config(args$lookup, "SENS_LOOKUP_PATH", "lookup_path")
if (!nzchar(lookup_path)) {
  stop(
    "Provide a private lookup path in 00_data/private/anonymization_config.R, ",
    "or with --lookup=/private/path/id_lookup.csv, ",
    "or with the SENS_LOOKUP_PATH environment variable.",
    call. = FALSE
  )
}
lookup_path <- resolve_path(lookup_path)

default_eye_paths <- c(
  here("00_data", "derived", "preprocessing", "duration_data.csv"),
  here("00_data", "derived", "analysis", "df_analysis_pub_prep.csv")
) %>%
  keep(file.exists)

eye_paths <- split_paths(coalesce_config(args$eye, "SENS_EYE_PATHS", "eye_paths", default = ""))
eye_paths <- if (length(eye_paths) == 0) default_eye_paths else map_chr(eye_paths, resolve_path)

medication_paths <- split_paths(coalesce_config(args$meds, "SENS_MEDICATION_PATHS", "medication_paths", default = "")) %>%
  map_chr(resolve_path)
disorder_paths <- split_paths(coalesce_config(args$disorders, "SENS_DISORDER_PATHS", "disorder_paths", default = "")) %>%
  map_chr(resolve_path)

input_paths <- c(eye_paths, medication_paths, disorder_paths)
if (length(input_paths) == 0) {
  stop("No input files configured.", call. = FALSE)
}
missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs) > 0) {
  stop("Input file(s) not found:\n", paste(missing_inputs, collapse = "\n"), call. = FALSE)
}

output_dir <- resolve_path(coalesce_config(
  args$`out-dir`,
  "SENS_ANON_OUTPUT_DIR",
  "output_dir",
  default = file.path("00_data", "derived", "sensitivity", "medication_disorder_anonymized")
))

create_lookup_if_missing <- as.logical(coalesce_config(
  args$`create-lookup`,
  "SENS_CREATE_LOOKUP",
  "create_lookup",
  default = "TRUE"
))
id_candidates <- c("ID", "Id", "id", "Sub_ID", "participant_id", "participant", "subject_id")

lookup <- load_or_create_lookup(
  paths = input_paths,
  lookup_path = lookup_path,
  id_candidates = id_candidates,
  create_if_missing = create_lookup_if_missing
)

summary <- input_paths %>%
  map_dfr(~ anonymize_table(
    input_path = .x,
    output_path = make_output_path(.x, output_dir),
    lookup = lookup,
    id_candidates = id_candidates
  ))

write_csv(summary, file.path(output_dir, "00_anonymization_summary.csv"))

message("Anonymization complete.")
message("Lookup path was used but not copied into the output directory.")
message("Anonymized files written to: ", output_dir)
message("Rows written by file:")
walk2(summary$output, summary$rows, ~ message("  ", basename(.x), ": ", .y))
