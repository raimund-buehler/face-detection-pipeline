# Manuscript section: Supplementary / control analyses
# Analysis family: menstrual cycle dataset preparation
# Original source path: data_questionnaire/Medication/menstrual_phase.r
# Primary input dataset(s): menstrual-cycle questionnaire exports and participant management files
# Primary output(s): menst_phase.csv and related menstrual-phase tables
# Known TODOs: several data-location assumptions remain external; clarify canonical output location relative to manuscript datasets
# Scientific logic note: copied from source without changing scientific logic

library(tidyverse)
library(readxl)
library(here)

df <- read_csv(file.path(Sys.getenv("DATA_DIR"), "sosci_data.csv"))

df <-
    df %>%
    rename(
        Gender = "Male or female",
        ID = "Codeeingabe: Code",
        cycle_day = "Datum Zyklus: (Tag) ... ",
        cycle_month = "Datum Zyklus: (Monat) ... 2022"
    )

menst <-
    df %>%
    filter(
        QUESTNNR == "PreScr_new",
    ) %>%
    rename(Group = "ASD diagnosis") %>%
    select(ID, Group, Gender, Pille_Verhütung, Zyklustage, cycle_day, cycle_month) %>%
    mutate(ID = as.character(tolower(ID)))

test_dates <- read_csv(file.path(Sys.getenv("DATA_DIR"), "part_man_sess.csv")) %>%
    rename(
        ID = "Personal ID",
        date_1 = "1. Date",
        date_2 = "2. Date",
        date_3 = "3. Date",
        date_4 = "4. Date",
    ) %>%
    select(ID, date_1, date_2, date_3, date_4)

menst <- left_join(test_dates, menst, by = "ID") %>% filter(Gender == "weiblich (Zyklusfragen..)")

# Convert the date_1 column to a date object
menst <- menst %>%
    mutate(
        date_1 = dmy(date_1),
        date_2 = dmy(date_2),
        date_3 = dmy(date_3),
        date_4 = dmy(date_4)
    )

# Extract the year and month from date_1
menst <- menst %>%
    mutate(
        year = year(date_1),
        month_1 = month(date_1),
        cycle_month = month(as.numeric(cycle_month))
    )

# Adjust the year based on the cycle_month
menst <- menst %>%
    mutate(
        adjusted_year = ifelse(cycle_month > month_1, year - 1, year)
    )

# Create a full date from cycle_day, cycle_month, and adjusted_year
menst <- menst %>%
    mutate(
        cycle_date = make_date(adjusted_year, cycle_month, cycle_day)
    ) %>%
    select(-c(cycle_day, cycle_month, year, month_1, adjusted_year))


menst <- menst %>% filter(!is.na(cycle_date))

# Remove duplicates that are exactly the same
menst <- menst %>% distinct()

# Keep Zyklustage 28 for dahö96 (more reasonable entry)
menst <- menst %>%
    arrange(ID, desc(cycle_date))

menst <- menst %>%
    group_by(ID) %>%
    mutate(
        Zyklustage = if_else(
            ID == "dahö96" & row_number() == 1, # Replace "your_conflict_id" with the actual ID
            lead(Zyklustage), # Use the Zyklustage from the less recent entry
            Zyklustage
        )
    ) %>%
    ungroup()

# Keep most recent date for other duplicates
menst <- menst %>%
    arrange(ID, desc(cycle_date)) %>%
    distinct(ID, .keep_all = TRUE)


# check duplicate IDs
duplicate_IDs <- menst %>%
    group_by(ID) %>%
    filter(n() > 1) %>%
    pull(ID) %>%
    unique()

duplicate_IDs

dups <- menst %>%
    filter(ID %in% duplicate_IDs)

# Convert Zyklustage to numeric
menst$Zyklustage_num <- as.numeric(str_extract(menst$Zyklustage, "\\d+"))

# Check which Zyklustage are above 40 (3 IDs)
menst %>%
    filter(Zyklustage_num >= 40) %>%
    select(ID, Zyklustage_num)

# Correct date
menst <- menst %>%
    mutate(
        # Step 2: Correct the date where year is 2320
        date_2 = if_else(year(date_2) == 2320,
            ymd(paste0("2023", substr(date_2, 5, 10))),
            date_2
        )
    )

# Calculate the cycle day for each testing date
menst <- menst %>%
    mutate(
        actual_diff_1 = as.numeric(difftime(date_1, cycle_date, units = "days")),
        effective_cycle_day_1 = (actual_diff_1 %% Zyklustage_num) + 1,
        actual_diff_2 = as.numeric(difftime(date_2, cycle_date, units = "days")),
        effective_cycle_day_2 = (actual_diff_2 %% Zyklustage_num) + 1,
        actual_diff_3 = as.numeric(difftime(date_3, cycle_date, units = "days")),
        effective_cycle_day_3 = (actual_diff_3 %% Zyklustage_num) + 1,
        actual_diff_4 = as.numeric(difftime(date_4, cycle_date, units = "days")),
        effective_cycle_day_4 = (actual_diff_4 %% Zyklustage_num) + 1
    )

menst_no_filter <- menst

# check which cycle_day_1 are > 100
menst %>%
    filter(actual_diff_1 > 120) %>%
    select(ID, actual_diff_1)

# 4 particpants
# filter them out
menst <- menst %>%
    filter(actual_diff_1 <= 120)

# Define a function to categorize the menstrual phase
get_menstrual_phase <- function(cycle_day) {
    if (is.na(cycle_day) || !is.numeric(cycle_day)) {
        return(NA) # Return NA if the input is invalid or missing
    }

    if (cycle_day >= 1 & cycle_day <= 5) {
        return("Menstrual")
    } else if (cycle_day >= 6 & cycle_day <= 13) {
        return("Follicular")
    } else if (cycle_day >= 14 & cycle_day <= 16) {
        return("Ovulatory")
    } else if (cycle_day >= 17) {
        return("Luteal")
    } else {
        return(NA) # Return NA if cycle_day is out of expected range
    }
}

# Apply this function to each cycle day to determine the phase
menst <- menst %>%
    mutate(
        phase_1 = sapply(effective_cycle_day_1, get_menstrual_phase),
        phase_2 = sapply(effective_cycle_day_2, get_menstrual_phase),
        phase_3 = sapply(effective_cycle_day_3, get_menstrual_phase),
        phase_4 = sapply(effective_cycle_day_4, get_menstrual_phase)
    )

write_csv(menst, file.path(Sys.getenv("DATA_DIR"), "menst_all.csv"))

menst %>%
    group_by(Pille_Verhütung) %>%
    count()

# Save the data
menst_phase <- menst %>% select(ID, phase_1, phase_2, phase_3, phase_4)

write_csv(menst_phase, here("00_data", "derived", "assembly", "menst_phase.csv"))

part_man <- read_xlsx(here("00_data", "raw", "questionnaires", "participants_1812 (1).xlsx")) %>% filter(Gender == "f")

part_man_ID <- part_man$"Personal ID"

part_man_ID

# which part_man_ID are not in menst

missing <- setdiff(part_man_ID, menst_no_filter$ID)

df_missing <- df %>%
    filter(ID %in% missing) %>%
    filter(QUESTNNR == "PreScr_new")
