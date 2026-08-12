# Manuscript section: Dataset assembly for manuscript core analysis table
# Analysis family: dataset assembly
# Original source path: scripts/merge_pm_df.R
# Primary input dataset(s): participant management files; accuracy summaries; questionnaire files; duration-derived analysis table
# Primary output(s): data/derived/assembly/df_full_merge.csv; data/derived/preprocessing/duration_data.csv
# Known TODOs: several path assumptions remain legacy; verify intended canonical write paths and RDS/CSV naming
# Scientific logic note: copied from source without changing scientific logic

library(tidyverse)
library(stringi)
library(here)
library(readxl)

part_man <- read_csv(file.path(Sys.getenv("DATA_DIR"), "Particpant Management - Overall Info.csv"))
session_info <- read_csv(file.path(Sys.getenv("DATA_DIR"), "Particpant Management - Session Info.csv"))

part_man <-
  part_man %>%
  select(`Personal ID`, `Male/Female`, Group) %>%
  rename(
    "Sub_ID" = `Personal ID`,
    "Sex" = `Male/Female`,
    "Group" = Group
  )

write_csv2(part_man, here("00_data", "derived", "intermediate", "part_man.csv"))

# transformed in excel

part_man <- read_csv(here("00_data", "raw", "participant_management", "part_man_trans.csv"))
participant_registry <- read_csv(here("00_data", "raw", "participant_management", "participants_1812.xlsx - Tabelle1.csv")) %>%
  select(`Personal ID`, Gender) %>%
  rename(Sub_ID = `Personal ID`, Sex_registry = Gender) %>%
  mutate(Sub_ID = stri_trans_general(Sub_ID, "Latin-ASCII")) %>%
  filter(!is.na(Sub_ID)) %>%
  distinct(Sub_ID, .keep_all = TRUE)
df <- read.csv(here("00_data", "derived", "preprocessing", "full_data.csv"))

# duration
df <- read_csv(here("00_data", "derived", "preprocessing", "duration_data.csv"))

# adjust Sub_ID for Umlauts
df$Sub_ID <- stri_trans_general(df$Sub_ID, "Latin-ASCII")
part_man$Sub_ID <- stri_trans_general(part_man$Sub_ID, "Latin-ASCII")
part_man <- part_man %>%
  left_join(participant_registry, by = "Sub_ID") %>%
  mutate(Sex = coalesce(Sex_registry, Sex)) %>%
  select(-Sex_registry)

df_group <- left_join(df, part_man, by = "Sub_ID")

df_group %>%
  filter(is.na(Group)) %>%
  distinct(Sub_ID)

session_info <-
  session_info %>%
  select(
    `Personal ID`, `1. Pill`, `1. Spray`, `2. Pill`, `2. Spray`,
    `3. Pill`, `3. Spray`, `4. Pill`, `4. Spray`
  )

write_csv2(session_info, here("00_data", "derived", "intermediate", "session_info.csv"))

# did transforming (unblind and wide to long) in excel

sessions_fixed <- read_csv(here("00_data", "raw", "participant_management", "sessions_fixed.csv"))
sessions_fixed$Sub_ID <- stri_trans_general(sessions_fixed$Sub_ID, "Latin-ASCII")

df_group_session <- left_join(df_group, sessions_fixed, by = c("Sub_ID", "session"))

df_group_session %>%
  filter(is.na(medication)) %>%
  distinct(Sub_ID)

# compare with actually tested participants

session_ids <- read_csv(here("00_data", "raw", "participant_management", "participants_1812.xlsx - Tabelle1.csv"))

tested <- session_ids %>%
  slice(0:79) %>%
  pull(`Personal ID`)
in_data <- df_group_session %>%
  distinct(Sub_ID) %>%
  pull(Sub_ID)

tested[which(!tested %in% in_data)]

write_csv(df_group_session, here("00_data", "derived", "intermediate", "full_data_session.csv"))
write_rds(df_group_session, here("00_data", "derived", "intermediate", "full_data_session.rds"))

# merge with accuracy

acc <- read_rds(here("00_data", "derived", "quality", "acc_df.rds"))
acc$Sub_ID <- stri_trans_general(acc$Sub_ID, "Latin-ASCII")

df_group_session <-
  df_group_session %>%
  mutate(session = str_extract(session, "(Session\\s[0-9]+)") %>%
    str_replace_all(" ", "_") %>%
    str_to_lower())

df_group_session_acc <- left_join(df_group_session, acc, by = c("Sub_ID", "session"))

# merge with questionnaires

quests <- read_xlsx(here("00_data", "raw", "questionnaires", "Fragebögen_ausgewertet.xlsx"))
panas <- read_xlsx(here("00_data", "raw", "questionnaires", "PANAS_ausgewertet_final_corrected.xlsx"))
quests$Sub_ID <- stri_trans_general(quests$Sub_ID, "Latin-ASCII")
panas$Sub_ID <- stri_trans_general(panas$Sub_ID, "Latin-ASCII")

quests %>%
  filter(is.na(`AQ-K Score`)) %>%
  distinct(Sub_ID)

anti_join(df_group_session_acc, quests, by = "Sub_ID") %>% distinct(Sub_ID)

df_quest_merge <- left_join(df_group_session_acc, quests, by = "Sub_ID")

df_quest_merge %>%
  filter(is.na(`AQ-K Score`)) %>%
  distinct(Sub_ID)

panas <-
  panas %>%
  mutate(session = str_extract(session, "(Session\\s[0-9]+)") %>%
    str_replace_all(" ", "_") %>%
    str_to_lower())

df_full_merge <- left_join(df_quest_merge, panas, by = c("Sub_ID", "session"))

anti_join(df_full_merge, panas, by = c("Sub_ID", "session")) %>% distinct(Sub_ID, session)

anti_join(df_full_merge, df_quest_merge, by = c("Sub_ID", "session")) %>% distinct(Sub_ID, session)

panas_id <- panas %>% distinct(Sub_ID, session)
df_quest_id <- df_quest_merge %>% distinct(Sub_ID, session)

setdiff(panas_id, df_quest_id)
setdiff(df_quest_id, panas_id)

# write data

write_csv(df_full_merge, here("00_data", "derived", "assembly", "df_full_merge.csv"))
write_rds(df_full_merge, here("00_data", "derived", "assembly", "df_full_merge.rds"))
