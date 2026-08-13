# Manuscript section: Preprocessing utility for AOI and marker transformation
# Analysis family: preprocessing helper
# Original source path: scripts/transform_data.R
# Primary input dataset(s): raw pupil positions; fixation result export; annotation markers
# Primary output(s): transformed fixation table with AOI labels and timestamps
# Known TODOs: review question-label recoding assumptions; document marker irregularity handling
# Scientific logic note: copied from source without changing scientific logic

library(tidyverse)

transform_data <- function(time, result, marker) {
  time <-
    time %>%
    filter(method == "pye3d 0.3.0 real-time") %>%
    select(frame_id = world_index, eye_id, pupil_timestamp, diameter_3d, confidence) %>%
    mutate(
      pupil_timestamp = round(pupil_timestamp - first(pupil_timestamp), 2),
      eye_id = as.factor(eye_id)
    ) %>%
    # Timestamps for distinct frames only
    group_by(frame_id, eye_id) %>%
    summarize(
      pupil_timestamp = mean(pupil_timestamp),
      diameter_3d = mean(diameter_3d, na.rm = TRUE),
      confidence = first(confidence)
    )

  # Reset marker

  # marker <- sub_data[[1]]["marker"]
  # marker <- marker[[1]]

  levels <- c(
    "Question 1" = "1",
    "Question 2" = "2",
    "Question 3" = "3",
    "Question 4" = "4",
    "Question 5" = "5",
    "Question 6" = "6",
    "Question 7" = "7",
    "Question 8" = "8",
    "Question 9" = "9",
    "Question 10" = "0",
    "Question 11" = "q",
    "Question 12" = "w",
    "Question 13" = "e",
    "Question 10" = "10",
    "Question 11" = "11",
    "Question 12" = "12",
    "Question 13" = "13",
    "Question 1" = "Question 1",
    "Question 2" = "Question 2",
    "Question 3" = "Question 3",
    "Question 4" = "Question 4",
    "Question 5" = "Question 5",
    "Question 6" = "Question 6",
    "Question 7" = "Question 7",
    "Question 8" = "Question 8",
    "Question 9" = "Question 9",
    "Question 10" = "Question 10",
    "Question 11" = "Question 11",
    "Question 12" = "Question 12",
    "Question 13" = "Question 13"
  )

  levels <- fct_recode(marker$label, !!!levels)

  marker$label <- levels

  marker$frame_id <- as.double(marker$frame_id)

  desired_order <- c(paste("Question", 1:13))

  marker <- marker %>%
    filter(label %in% desired_order)

  marker <- marker %>%
    mutate(label = factor(label, levels = desired_order, ordered = TRUE))

  desired_factor <- factor(desired_order, levels = desired_order, ordered = T)

  filter_label <- function(label, marker) {
    irregular_indices <- which(label != desired_factor)
    if (length(label) > 13) {
      marker <- marker %>%
        filter(!(row_number() == first(irregular_indices)))

      label <- marker$label

      return(filter_label(label, marker))
    }
    return(marker)
  }

  irregular_indices <- which(marker$label != desired_factor)

  if (length(irregular_indices) > 0 & length(marker$label) >= 13) {
    # Check if length of label is too long
    marker <- filter_label(marker$label, marker)
    irregular_indices <- which(marker$label != desired_factor)
    # Check if the values of "label" in the marker dataframe are not as expected
    if (!identical(marker$label, desired_factor) & length(marker$label) == 13) {
      # If the values are not as expected, replace them in the marker dataframe
      marker$label[first(irregular_indices):last(irregular_indices)] <-
        desired_order[first(irregular_indices):last(irregular_indices)]
    }
  }

  # new_labels <- fct_match(marker$label, !!!levels)
  #
  # labels_old <- c(marker %>% distinct(label) %>% pull())
  # labels_new <- factor(c(paste("Question", 1:12)), levels = c(paste("Question", 1:12)), ordered = T)
  #
  # replacements <- setNames(labels_old, labels_new)
  #
  # marker$label <- recode(marker$label, !!!replacements)

  result_marked <-
    # Join question markers and result file
    left_join(result, marker, by = "frame_id") %>%
    fill(label, .direction = "down") %>%
    mutate(
      label = fct_explicit_na(label, "Intro"),
      label = factor(label, levels = c("Intro", paste("Question", 1:13)), ordered = T)
    ) %>%
    # convert columns
    mutate(
      id_fix = as.integer(id_fix),
      conf_fix = as.double(conf_fix)
    ) %>%
    # making AOI's independent (fix on eye/mouth doesn't count as fix on face)
    mutate(fix_on_face = if_else((fix_on_face == T & fix_on_mouth == T) |
      (fix_on_face == T & fix_on_eyes == T) |
      (fix_on_face == F), F, T)) %>%
    # select needed columns
    select(label, frame_id, id_fix, duration_fix, conf_fix, fix_on_face, fix_on_eyes, fix_on_mouth) %>%
    # Fill with NA's where no fix is present
    mutate_at(c("fix_on_face", "fix_on_eyes", "fix_on_mouth"), ~ ifelse(is.na(id_fix), NA, .)) %>%
    # create fix_on_background column
    mutate(
      "fix_on_background" =
        ifelse((fix_on_face == F & fix_on_eyes == F & fix_on_mouth == F), T, F)
    ) %>%
    relocate(fix_on_background, .after = "fix_on_mouth")

  # Join result file and timestamps/diameter
  result_merged <- left_join(result_marked, time, by = "frame_id")

  result_fix <-
    result_merged %>%
    # filter(!id_fix == "False") %>%
    group_by(id_fix) %>%
    mutate(
      fix_on_eyes = mean(fix_on_eyes),
      fix_on_mouth = mean(fix_on_mouth),
      fix_on_face = mean(fix_on_face),
      fix_on_background = mean(fix_on_background)
    ) %>%
    distinct(id_fix, frame_id, .keep_all = T) %>%
    mutate(
      fix_on_eyes = if_else(fix_on_eyes >= 0.5,
        TRUE,
        FALSE
      ),
      fix_on_mouth = if_else(fix_on_mouth >= 0.5,
        TRUE,
        FALSE
      ),
      fix_on_face = if_else(fix_on_face >= 0.5 & !(fix_on_eyes >= 0.5) & !(fix_on_mouth >= 0.5),
        TRUE,
        FALSE
      ),
      fix_on_background = if_else(fix_on_background >= 0.5,
        TRUE,
        FALSE
      )
    )


  return(result_fix)
}
