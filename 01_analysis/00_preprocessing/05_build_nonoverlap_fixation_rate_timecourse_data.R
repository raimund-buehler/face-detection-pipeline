library(tidyverse)
library(here)

output_path <- here("00_data", "derived", "preprocessing", "fixation_rate_nonoverlap_bins_data.csv")
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

df <- read_csv(here("00_data", "derived", "preprocessing", "full_data.csv"), show_col_types = FALSE) %>%
  transmute(
    Sub_ID,
    session,
    fix,
    label,
    rate_q,
    time_question
  ) %>%
  filter(!is.na(rate_q), is.finite(rate_q))

session_meta <- read_csv(here("00_data", "derived", "analysis", "df_analysis_pub.csv"), show_col_types = FALSE) %>%
  distinct(ID, session, Sex, Group, medication, accuracy_degrees) %>%
  mutate(session_norm = str_replace_all(str_to_lower(session), " ", "_"))

df <- df %>%
  mutate(
    session_norm = str_replace_all(str_to_lower(session), " ", "_"),
    q_num = parse_number(label)
  ) %>%
  filter(!is.na(q_num), q_num >= 1, q_num <= 13) %>%
  mutate(
    q_start = (q_num - 1) / 13,
    q_end = q_num / 13
  ) %>%
  left_join(
    session_meta %>% mutate(session_norm = str_replace_all(str_to_lower(session), " ", "_")),
    by = c("Sub_ID" = "ID", "session_norm"),
    suffix = c("", "_meta")
  ) %>%
  filter(!is.na(Group), !is.na(Sex), !is.na(medication))

bins <- tibble(
  progress_bin_index = 1:10,
  progress_start = seq(0, 0.9, by = 0.1),
  progress_end = seq(0.1, 1.0, by = 0.1),
  progress_midpoint = progress_start + 0.05,
  progress_bin = sprintf("%02d", progress_bin_index)
)

expanded <- df %>%
  crossing(bins) %>%
  mutate(
    overlap_start = pmax(q_start, progress_start),
    overlap_end = pmin(q_end, progress_end),
    overlap = pmax(0, overlap_end - overlap_start)
  ) %>%
  filter(overlap > 0)

rate_bins <- expanded %>%
  group_by(
    Sub_ID, session_norm, Sex, Group, medication, accuracy_degrees,
    progress_bin, progress_bin_index, progress_midpoint, progress_start, progress_end, fix
  ) %>%
  summarise(
    window_fix_rate = weighted.mean(rate_q, w = overlap, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(session = session_norm)

write_csv(rate_bins, output_path)
