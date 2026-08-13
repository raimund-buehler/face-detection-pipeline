## Manuscript section: supplementary / control tables
## Analysis family: session interval summary table
## Original source path: legacy/data_archive/original_layout/data_questionnaire/Medication/date_analysis.r
## Primary input dataset(s): external DATA_DIR/part_man_sess.csv, data/raw/questionnaires/participants_1812 (1).xlsx
## Primary output(s): outputs/tables/medication_controls/session_interval_summary_overall.csv, outputs/tables/medication_controls/session_interval_summary_by_group.csv
## Known TODOs: session date source still depends on DATA_DIR
## Scientific logic unchanged from source except for canonical file-path cleanup and explicit table export

library(tidyverse)
library(readxl)
library(lubridate)
library(here)

test_dates <- read_csv(file.path(Sys.getenv("DATA_DIR"), "part_man_sess.csv")) %>%
    rename(
        ID = "Personal ID",
        date_1 = "1. Date",
        date_2 = "2. Date",
        date_3 = "3. Date",
        date_4 = "4. Date",
    ) %>%
    select(ID, date_1, date_2, date_3, date_4)

group_gender <- read_xlsx(here("00_data", "raw", "questionnaires", "participants_1812 (1).xlsx")) %>% rename(ID = "Personal ID")

# whhich ID is in test_dates but not in group_gender
setdiff(test_dates$ID, group_gender$ID)

dates <- left_join(group_gender, test_dates, by = "ID")

# Convert the date columns to date objects
dates <- dates %>%
    mutate(
        date_1 = ymd(date_1),
        date_2 = ymd(date_2),
        date_3 = ymd(date_3),
        date_4 = ymd(date_4)
    )

# Correct date
dates <- dates %>%
    mutate(
        # Step 2: Correct the date where year is 2320
        date_2 = if_else(year(date_2) == 2320,
            ymd(paste0("2023", substr(date_2, 5, 10))),
            date_2
        )
    )

intervals <- dates %>%
    rowwise() %>%
    mutate(
        interval_1_2 = as.numeric(difftime(date_2, date_1, units = "days")),
        interval_2_3 = as.numeric(difftime(date_3, date_2, units = "days")),
        interval_3_4 = as.numeric(difftime(date_4, date_3, units = "days"))
    ) %>%
    ungroup()

intervals %>% select(ID, interval_1_2:interval_3_4)

results <- intervals %>%
    summarise(
        mean_interval = mean(c(interval_1_2, interval_2_3, interval_3_4), na.rm = TRUE),
        std_interval = sd(c(interval_1_2, interval_2_3, interval_3_4), na.rm = TRUE)
    )

print(results)

results_grouped <- intervals %>%
    group_by(Group) %>%
    summarise(
        mean_interval = mean(c(interval_1_2, interval_2_3, interval_3_4), na.rm = TRUE),
        std_interval = sd(c(interval_1_2, interval_2_3, interval_3_4), na.rm = TRUE)
    )

print(results_grouped)

write_csv(results, here("02_outputs", "tables", "medication_controls", "session_interval_summary_overall.csv"))
write_csv(results_grouped, here("02_outputs", "tables", "medication_controls", "session_interval_summary_by_group.csv"))
