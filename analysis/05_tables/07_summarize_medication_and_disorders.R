## Manuscript section: supplementary / control tables
## Analysis family: ASD medication and comorbidity summary tables
## Original source path: legacy/data_archive/original_layout/data_questionnaire/Medication/medication_and_disorders.r
## Primary input dataset(s): external DATA_DIR/sosci_data.csv, external DATA_DIR/df_meds.csv, external DATA_DIR/df_disorders.csv, data/derived/assembly/exp_stat.csv
## Primary output(s): outputs/tables/medication_controls/asd_medication_type_counts.csv, outputs/tables/medication_controls/asd_medication_count_by_person.csv, outputs/tables/medication_controls/asd_disorder_counts.csv
## Known TODOs: depends on manually cleaned external `df_meds.csv` and `df_disorders.csv`
## Scientific logic unchanged from source except for canonical file-path cleanup and explicit table export

library(tidyverse)
library(here)
library(janitor)
library(stringdist)

df <- read_csv(file.path(Sys.getenv("DATA_DIR"), "sosci_data.csv"))
IDs <- read_csv(here("00_data", "derived", "assembly", "exp_stat.csv"), locale = locale(encoding = "UTF-8"))

IDs <- IDs %>% pull(ID)

df <-
    df %>%
    filter(
        QUESTNNR == "PreScr_new",
    ) %>%
    rename(Group = "ASD diagnosis") %>%
    mutate(ID = as.character(tolower(ID)))

df <- df %>%
    filter(Group == "Ja") %>%
    filter(ID %in% IDs)

df <- df %>%
    rename(
        "Antipsychotika" = Meds,
        meds_n = "Meds Name: Number of mentions",
        meds_1 = "Meds Name: Mention 1",
        meds_2 = "Meds Name: Mention 2",
        meds_3 = "Meds Name: Mention 3",
        disorders = Disturbs,
        disorders_n = "Disorders 2: Number of mentions",
        disorders_1 = "Disorders 2: Mention 1",
        disorders_2 = "Disorders 2: Mention 2",
        disorders_3 = "Disorders 2: Mention 3",
    ) %>%
    select(ID, Medikamente, Antipsychotika, meds_n, meds_1, meds_2, meds_3, disorders, disorders_n, disorders_1, disorders_2, disorders_3)

df <- df %>% clean_names()
# MEDS

df_meds <- df %>%
    filter(medikamente == "Ja") %>%
    select(id, meds_n:meds_3)

df_meds <- df_meds %>%
    pivot_longer(cols = c(meds_1, meds_2, meds_3), names_to = "meds", values_to = "meds_name") %>%
    filter(!is.na(meds_n) & !is.na(meds_name))

# write_csv(df_meds, file.path(Sys.getenv("DATA_DIR"), "df_meds.csv"))

# Manually cleaned data

df_meds <- read_csv(file.path(Sys.getenv("DATA_DIR"), "df_meds.csv"))

df_meds %>%
    distinct(id) %>%
    nrow()
# 20 on meds

df_meds %>%
    group_by(type) %>%
    summarise(n = n()) %>%
    arrange(desc(n))

medication_type_counts <- df_meds %>%
    group_by(type) %>%
    summarise(n = n()) %>%
    arrange(desc(n))

medication_count_by_person <- df_meds %>%
    group_by(meds_n) %>%
    distinct(id) %>%
    summarise(n = n())

df_disorders <- df %>%
    filter(disorders == "Ja") %>%
    select(id, disorders_n:disorders_3)

df_disorders <- df_disorders %>%
    pivot_longer(cols = c(disorders_1, disorders_2, disorders_3), names_to = "disorders", values_to = "name") %>%
    filter(!is.na(disorders) & !is.na(name))

# write_csv(df_disorders, file.path(Sys.getenv("DATA_DIR"), "df_disorders.csv"))

# MANUAL CLEANING
df_disorders <- read_csv(file.path(Sys.getenv("DATA_DIR"), "df_disorders.csv"))

df_disorders %>%
    distinct(id) %>%
    nrow()
# 23 with disorders

df_disorders %>%
    group_by(name_clean) %>%
    summarise(n = n()) %>%
    arrange(desc(n))

disorder_counts <- df_disorders %>%
    group_by(name_clean) %>%
    summarise(n = n()) %>%
    arrange(desc(n))

write_csv(medication_type_counts, here("02_outputs", "tables", "medication_controls", "asd_medication_type_counts.csv"))
write_csv(medication_count_by_person, here("02_outputs", "tables", "medication_controls", "asd_medication_count_by_person.csv"))
write_csv(disorder_counts, here("02_outputs", "tables", "medication_controls", "asd_disorder_counts.csv"))
