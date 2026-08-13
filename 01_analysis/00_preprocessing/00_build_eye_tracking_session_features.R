# Machine-specific paths live in config/paths.R (not committed).
source(here::here("config", "paths.R"))
# Manuscript section: Preprocessing and session-level eye-tracking feature assembly
# Analysis family: preprocessing
# Original source path: scripts/timecourse_all_preprocess.R
# Primary input dataset(s): raw eye-tracking export folders under result_csv.nosync
# Primary output(s): data/derived/preprocessing/full_data.rds; data/derived/preprocessing/duration_data.csv; legacy df_analysis/full_data.csv path reference
# Known TODOs: replace absolute raw-data paths; document current canonical output path(s); review legacy df_analysis references
# Scientific logic note: copied from source without changing scientific logic

library(tidyverse)
library(here)

source(here("01_analysis", "00_preprocessing", "shared_read_data.R"))
source(here("01_analysis", "00_preprocessing", "shared_transform_data.R"))
source(here("01_analysis", "00_preprocessing", "shared_plot_sub.R"))
source(here("01_analysis", "00_preprocessing", "shared_analyze_sub.R"))

# Directory where the subject directories are located
base_directory <- recordings_root


# List all subject directories within the base directory
subject_directories <- list.dirs(base_directory, full.names = TRUE, recursive = FALSE)

# Initialize a list to store the data for all subjects
all_subjects_data <- list()

# Loop through each subject directory
for (subject_dir in subject_directories) {
  subject_id <- basename(subject_dir)
  session_dirs <- list.dirs(subject_dir, full.names = TRUE, recursive = FALSE)

  for (session_dir in session_dirs) {
    session_label <- basename(session_dir)

    if (any(excluded_sessions$Sub_ID == subject_id & excluded_sessions$session == session_label)) {
      message("Skipping excluded session: ", subject_id, " / ", session_label)
      next
    }

    # Check if the subject directory contains the required files

    time_path <- file.path(session_dir, "000/exports/000/pupil_positions.csv")
    result_path <- file.path(session_dir, "000/exports/000/result.csv")
    marker_path <- file.path(session_dir, "000/exports/000/annotations.csv")

    # Check if the files exist for the current subject
    if (file.exists(time_path) &&
      file.exists(result_path) &&
      file.exists(marker_path)) {
      processed_session <- tryCatch(
        {
          sub_data <- read_data(time_path, result_path, marker_path)
          result_fix <- transform_data(sub_data$time, sub_data$result, sub_data$marker)
          analyze_sub(sub_data$Sub_ID, sub_data$session, result_fix)
        },
        error = function(e) {
          message("Skipping session due to preprocessing error: ", session_dir)
          message("Reason: ", conditionMessage(e))
          NULL
        }
      )

      if (!is.null(processed_session)) {
        all_subjects_data[[subject_id %>% paste(session_label)]] <- processed_session
      }
    }
  }
}

# Combine data for all subjects into a single data frame
combined_data <- do.call(rbind, all_subjects_data)

# Delete Rownames
rownames(combined_data) <- NULL

folder_path <- recordings_root

# check if all subs from folders are in df
folder_names <- list.files(path = folder_path, full.names = FALSE, all.files = FALSE)
sub_ids <- combined_data %>%
  distinct(Sub_ID) %>%
  pull()
setdiff(folder_names, sub_ids)
# (any folders without matching data are listed by the setdiff above)

# Subs for which not all markers are present
not_all_markers <-
  combined_data %>%
  group_by(Sub_ID, session) %>%
  summarise(row_count = n()) %>%
  ungroup() %>%
  group_by(row_count) %>%
  summarise(count = n())

# 10 missing completely
missings <- combined_data %>% filter(is.na(rate_all))

missings <- combined_data %>%
  filter(is.na(rate_all)) %>%
  distinct(Sub_ID, session)

duration_data <- combined_data %>%
  distinct(Sub_ID, session, fix, .keep_all = T) %>%
  select(Sub_ID, session, fix, percentage_fix_duration, total_fix_duration, rate_all)

duration_data %>%
  group_by(Sub_ID, session) %>%
  summarise(sum = sum(percentage_fix_duration, na.rm = T))

write.csv(duration_data, here("00_data", "derived", "preprocessing", "duration_data.csv"), fileEncoding = "UTF-8")

write_rds(combined_data, here("00_data", "derived", "preprocessing", "full_data.rds"))
write.csv(combined_data, here("00_data", "derived", "preprocessing", "full_data.csv"), fileEncoding = "UTF-8")

# High Fixation on eyes Subjects
combined_data %>%
  filter(rate_all > 3 & fix == "fix_on_eyes") %>%
  distinct(Sub_ID, session, .keep_all = T) %>%
  select(Sub_ID, session, rate_all)
# stbä01, s2

# Low Fixation on eyes Subjects
combined_data %>%
  filter(rate_all < 0.1 & fix == "fix_on_eyes") %>%
  distinct(Sub_ID, session, .keep_all = T) %>%
  select(Sub_ID, session, rate_all) %>%
  arrange(rate_all)
# (inspect the low-fixation sessions listed above)

# COMPARE WITH REGULARLY READ IN DATA

df_old <- read_csv(here("00_data", "derived", "preprocessing", "df_all.csv"))

df_old <- df_old %>%
  select(Sub_ID, session, AOI, norm_n) %>%
  rename("fix" = AOI)

df_old$fix <- recode(df_old$fix, Eyes = "fix_on_eyes", Face = "fix_on_face", Mouth = "fix_on_mouth", Non_Face = "fix_on_background")

df_new <- combined_data %>% distinct(Sub_ID, session, fix, .keep_all = T)

merge_oldnew <-
  left_join(df_new, df_old, by = c("Sub_ID", "session", "fix"))

# Check differences
merge_oldnew <- merge_oldnew %>% mutate("diff" = rate_all - norm_n)

merge_oldnew %>% filter(diff == max(diff, na.rm = T))

# Differences present, but not major, cluster around zero
hist(merge_oldnew$diff)
