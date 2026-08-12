# Manuscript section: Preprocessing utility for subject/session feature extraction
# Analysis family: preprocessing helper
# Original source path: scripts/analyze_sub.R
# Primary input dataset(s): transformed per-session eye-tracking fixation data
# Primary output(s): per-subject/session rate, slope, and duration summaries
# Known TODOs: trace downstream canonical consumers; keep synchronized with preprocessing wrapper if edited later
# Scientific logic note: copied from source without changing scientific logic

library(tidyverse)
library(broom)

analyze_sub <- function(Sub_ID, session, result_fix) {
  time_interval_size <- 1

  # RATES
  # (define result_a, equivalent to result_plot, but without integer for fix)
  result_a <- result_fix %>%
    pivot_longer(
      cols = fix_on_face:fix_on_background,
      names_to = "fix", values_to = "value"
    ) %>%
    mutate(fix = factor(fix, levels = unique(fix))) %>%
    group_by(fix) %>%
    mutate(rolling_average = rollmean(value, k = 25, fill = NA, align = "right")) %>%
    mutate(change_rate = c(rep(NA, 2), diff(rolling_average, lag = 2))) %>%
    mutate(
      time_interval =
        floor(pupil_timestamp / time_interval_size) * time_interval_size
    ) %>%
    group_by(time_interval, fix) %>%
    distinct(id_fix, .keep_all = T) %>%
    mutate(rate = sum(value, na.rm = T) / time_interval_size)

  # Step 1: Filter to only rows where the fixation is active (`value == TRUE`)
  filtered_result_a <- result_a %>%
    filter(value == TRUE)

  # Step 2: Get distinct fixations (`id_fix`) for each `fix`, ensuring each fixation is counted only once
  distinct_fixations <- filtered_result_a %>%
    distinct(id_fix, .keep_all = TRUE)

  # Step 3: Calculate the total fixation duration across all fixations
  total_fix_duration <- sum(as.double(distinct_fixations$duration_fix), na.rm = TRUE)

  question_duration_totals <- distinct_fixations %>%
    group_by(label) %>%
    summarise(
      total_fix_duration_q = sum(as.double(duration_fix), na.rm = TRUE),
      .groups = "drop"
    )

  question_fix_duration <- distinct_fixations %>%
    group_by(fix, label) %>%
    summarise(
      question_fix_duration = sum(as.double(duration_fix), na.rm = TRUE),
      .groups = "drop"
    )

  result_a <- result_a %>%
    filter(value == T) %>%
    group_by(fix) %>%
    mutate(
      percentage_fix_duration = sum(as.double(duration_fix), na.rm = T) / total_fix_duration,
      total_fix_duration = total_fix_duration
    )
  # distinct(fix, .keep_all = T) %>%
  # select(fix, percentage_fix_duration) %>%
  # print()

  # result_a %>% ungroup() %>% distinct(fix, percentage_fix_duration) %>% summarise(sum = sum(percentage_fix_duration, na.rm = T)) %>% print()

  # Calculate rates
  rates <-
    result_a %>%
    ungroup() %>%
    group_by(fix) %>%
    distinct(id_fix, .keep_all = T) %>%
    mutate(
      rate_all = sum(value, na.rm = T) / max(pupil_timestamp, na.rm = T)
    ) %>%
    group_by(fix, label) %>%
    mutate(time_question = max(pupil_timestamp, na.rm = T) - min(pupil_timestamp, na.rm = T)) %>%
    mutate(rate_q = sum(value, na.rm = T) / time_question) %>%
    ungroup() %>%
    distinct(label, fix, .keep_all = T) %>%
    left_join(question_fix_duration, by = c("fix", "label")) %>%
    left_join(question_duration_totals, by = "label") %>%
    mutate(
      duration_q_prop = if_else(total_fix_duration_q > 0, question_fix_duration / total_fix_duration_q, NA_real_),
      duration_q_per_sec = if_else(time_question > 0, question_fix_duration / time_question, NA_real_)
    ) %>%
    select(
      fix, label, time_question, rate_all, percentage_fix_duration, total_fix_duration, rate_q,
      question_fix_duration, total_fix_duration_q, duration_q_prop, duration_q_per_sec
    )

  # SLOPES
  # Entire timecourse
  slopes_all <-
    result_a %>%
    group_by(fix) %>%
    do(fit = lm(rate ~ time_interval, data = .)) %>%
    mutate(slope_all = tidy(fit)$estimate[2]) %>%
    ungroup()

  # By Question
  slopes_q <-
    result_a %>%
    group_by(fix, label) %>%
    do(fit = lm(rate ~ time_interval, data = .)) %>%
    mutate(slope_q = tidy(fit)$estimate[2]) %>%
    ungroup()

  # Slope for first five seconds after marker
  slopes_five <-
    result_a %>%
    group_by(label) %>%
    mutate(marker_time = min(time_interval, na.rm = T)) %>%
    filter(time_interval <= marker_time + 5) %>%
    # fill NA in time_interval
    fill(time_interval, .direction = "up") %>%
    group_by(fix, label) %>%
    # fit slope with lm
    do(fit = lm(rate ~ time_interval, data = .)) %>%
    mutate(slope_five = round(tidy(fit)$estimate[2], 6)) %>%
    ungroup()

  # Merge Slopes
  slopes_merged <-
    left_join(slopes_all, slopes_q, by = c("fix")) %>%
    select(label, fix, slope_all, slope_q)

  slopes_merged <-
    left_join(slopes_merged, slopes_five, by = c("fix", "label")) %>%
    select(-fit)

  df <-
    left_join(slopes_merged, rates, by = c("fix", "label")) %>%
    mutate(
      Sub_ID = Sub_ID,
      session = session
    ) %>%
    select(Sub_ID, session, everything())

  return(df)
}
