# Machine-specific paths live in config/paths.R (not committed).
source(here::here("config", "paths.R"))
library(tidyverse)
library(here)
library(stringr)
library(stringi)

base_directory <- recordings_root
window_width <- 0.10
window_step <- 0.025
window_half_width <- window_width / 2
window_midpoints <- seq(window_half_width, 1 - window_half_width, by = window_step)
window_labels <- sprintf("%02d", seq_along(window_midpoints))

window_bounds <- tibble(
  progress_window = window_labels,
  progress_midpoint = window_midpoints,
  progress_start = pmax(window_midpoints - window_half_width, 0),
  progress_end = pmin(window_midpoints + window_half_width, 1),
  progress_window_index = seq_along(window_midpoints)
)

analysis_sessions <- read_csv(here("00_data", "derived", "analysis", "df_analysis_pub.csv"), show_col_types = FALSE) %>%
  distinct(ID, session, Sex, Group, medication, accuracy_degrees) %>%
  mutate(
    ID = stri_trans_general(ID, "Latin-ASCII"),
    session = str_replace_all(session, "_", " "),
    session = str_to_title(session)
  )

desired_order <- paste("Question", 1:13)
label_map <- c(
  "1" = "Question 1", "2" = "Question 2", "3" = "Question 3", "4" = "Question 4",
  "5" = "Question 5", "6" = "Question 6", "7" = "Question 7", "8" = "Question 8",
  "9" = "Question 9", "0" = "Question 10", "q" = "Question 11", "w" = "Question 12", "e" = "Question 13"
)

extract_frame_timestamp <- function(timestamp_string, which = c("first", "last")) {
  position <- match.arg(which)
  if (is.na(timestamp_string) || !nzchar(timestamp_string)) return(NA_real_)

  values <- str_extract_all(timestamp_string, "-?\\d+\\.?\\d*")[[1]] %>% as.double()
  if (length(values) == 0 || all(is.na(values))) return(NA_real_)
  non_missing_idx <- which(!is.na(values))

  if (position == "first") values[non_missing_idx[1]] else values[tail(non_missing_idx, 1)]
}

read_result_with_labels <- function(result_path, marker_path, sub_id, session) {
  message("Fix-rate sliding window: ", sub_id, " / ", session)

  result <- read_csv(result_path, show_col_types = FALSE) %>%
    mutate(id_fix = as.integer(id_fix))

  marker <- read_csv(marker_path, show_col_types = FALSE, col_types = cols(label = col_character())) %>%
    select(index, label) %>%
    rename(frame_id = index) %>%
    mutate(frame_id = as.double(frame_id)) %>%
    mutate(label = recode(label, !!!label_map, .default = label)) %>%
    filter(label %in% desired_order) %>%
    arrange(frame_id) %>%
    mutate(label = factor(label, levels = desired_order, ordered = TRUE))

  result %>%
    left_join(marker, by = "frame_id") %>%
    fill(label, .direction = "down") %>%
    mutate(
      label = fct_explicit_na(label, "Intro"),
      label = factor(label, levels = c("Intro", desired_order), ordered = TRUE),
      fix_on_face = if_else((fix_on_face == TRUE & fix_on_mouth == TRUE) |
        (fix_on_face == TRUE & fix_on_eyes == TRUE) |
        (fix_on_face == FALSE), FALSE, TRUE),
      fix_on_background = if_else((fix_on_face == FALSE & fix_on_eyes == FALSE & fix_on_mouth == FALSE), TRUE, FALSE)
    )
}

compute_sliding_window_fix_rate <- function(sub_id, session, result_marked) {
  fixations_long <- result_marked %>%
    filter(label != "Intro") %>%
    distinct(id_fix, .keep_all = TRUE) %>%
    mutate(pupil_timestamp = map_dbl(gaze_tmstmps_frame, extract_frame_timestamp, which = "first")) %>%
    arrange(pupil_timestamp, frame_id) %>%
    pivot_longer(
      cols = c(fix_on_face, fix_on_eyes, fix_on_mouth, fix_on_background),
      names_to = "fix",
      values_to = "value"
    ) %>%
    filter(value == TRUE) %>%
    mutate(pupil_timestamp = as.double(pupil_timestamp)) %>%
    filter(!is.na(pupil_timestamp))

  if (nrow(fixations_long) == 0) return(NULL)

  interview_start <- min(fixations_long$pupil_timestamp, na.rm = TRUE)
  interview_end <- max(fixations_long$pupil_timestamp, na.rm = TRUE)
  interview_duration_sec <- interview_end - interview_start

  if (!is.finite(interview_duration_sec) || interview_duration_sec <= 0) return(NULL)

  fixations_long <- fixations_long %>%
    mutate(
      progress = (pupil_timestamp - interview_start) / interview_duration_sec,
      progress = pmin(pmax(progress, 0), 1)
    )

  window_data <- purrr::map_dfr(seq_len(nrow(window_bounds)), function(i) {
    current_window <- window_bounds[i, ]
    window_fixations <- fixations_long %>%
      filter(progress >= current_window$progress_start, progress <= current_window$progress_end)

    window_duration_sec <- (current_window$progress_end - current_window$progress_start) * interview_duration_sec
    if (!is.finite(window_duration_sec) || window_duration_sec <= 0) return(NULL)

    if (nrow(window_fixations) == 0) {
      return(expand_grid(
        progress_window = current_window$progress_window,
        fix = c("fix_on_eyes", "fix_on_mouth", "fix_on_face", "fix_on_background")
      ) %>%
        mutate(
          window_fix_count = 0L,
          window_total_fix_count = 0L,
          window_duration_sec = window_duration_sec,
          progress_midpoint = current_window$progress_midpoint,
          progress_start = current_window$progress_start,
          progress_end = current_window$progress_end,
          progress_window_index = current_window$progress_window_index
        ))
    }

    window_total_fix_count <- nrow(window_fixations)

    window_fixations %>%
      count(fix, name = "window_fix_count") %>%
      right_join(tibble(fix = c("fix_on_eyes", "fix_on_mouth", "fix_on_face", "fix_on_background")), by = "fix") %>%
      mutate(
        progress_window = current_window$progress_window,
        progress_midpoint = current_window$progress_midpoint,
        progress_start = current_window$progress_start,
        progress_end = current_window$progress_end,
        progress_window_index = current_window$progress_window_index,
        window_fix_count = replace_na(window_fix_count, 0L),
        window_total_fix_count = window_total_fix_count,
        window_duration_sec = window_duration_sec
      )
  })

  window_data %>%
    mutate(
      window_fix_rate = if_else(window_duration_sec > 0, window_fix_count / window_duration_sec, NA_real_),
      Sub_ID = sub_id,
      session = str_replace_all(session, " ", "_") %>% str_to_lower()
    ) %>%
    select(
      Sub_ID, session, progress_window, progress_window_index, progress_midpoint,
      progress_start, progress_end, fix, window_fix_count, window_total_fix_count,
      window_duration_sec, window_fix_rate
    )
}

subject_directories <- list.dirs(base_directory, full.names = TRUE, recursive = FALSE)
all_data <- list()

for (subject_dir in subject_directories) {
  subject_id <- basename(subject_dir)
  if (!subject_id %in% analysis_sessions$ID) next

  session_dirs <- list.dirs(subject_dir, full.names = TRUE, recursive = FALSE)
  for (session_dir in session_dirs) {
    session_label <- basename(session_dir)
    if (!any(analysis_sessions$ID == subject_id & analysis_sessions$session == session_label)) next

    result_path <- file.path(session_dir, "000/exports/000/result.csv")
    marker_path <- file.path(session_dir, "000/exports/000/annotations.csv")
    if (!(file.exists(result_path) && file.exists(marker_path))) next

    session_data <- tryCatch(
      {
        marked <- read_result_with_labels(result_path, marker_path, subject_id, session_label)
        compute_sliding_window_fix_rate(subject_id, session_label, marked)
      },
      error = function(e) {
        message("Skipping session: ", session_dir)
        message("Reason: ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(session_data)) {
      all_data[[paste(subject_id, session_label)]] <- session_data
    }
  }
}

fix_rate_sliding_window_data <- bind_rows(all_data) %>%
  left_join(
    analysis_sessions %>% mutate(session = str_replace_all(session, " ", "_") %>% str_to_lower()),
    by = c("Sub_ID" = "ID", "session")
  ) %>%
  relocate(Sex, Group, medication, accuracy_degrees, .after = session)

write_csv(
  fix_rate_sliding_window_data,
  here("00_data", "derived", "preprocessing", "fixation_rate_sliding_window_data.csv")
)
