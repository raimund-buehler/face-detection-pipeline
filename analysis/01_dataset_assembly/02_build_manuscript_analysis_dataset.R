# Manuscript section: Dataset assembly for manuscript-ready analysis dataset
# Analysis family: dataset assembly
# Original source path: legacy/scripts_archive/original_scripts/publication/merging/merge_aq_menst.r
# Primary input dataset(s): data/derived/assembly/df_full_merge.csv; data/derived/assembly/menst_phase.csv; data/derived/assembly/exp_stat.csv; data/derived/questionnaire_scoring/AQ_final_data.csv
# Primary output(s): data/derived/analysis/df_analysis_pub.csv
# Known TODOs: provenance for downstream df_analysis_pub_prep.csv remains unresolved
# Scientific logic note: copied from source without changing scientific logic

library(tidyverse)
library(here)
# DATA PREP

df <- read_csv(here("00_data", "derived", "assembly", "df_full_merge.csv")) %>%
    rename(ID = Sub_ID) %>%
    select(-last_col())
menst <- read_csv(here("00_data", "derived", "assembly", "menst_phase.csv")) %>% rename_with(~ str_replace(., "phase", "session"))
exp_stat <- read_csv(here("00_data", "derived", "assembly", "exp_stat.csv"), locale = locale(encoding = "UTF-8"))
AQ_final_data <- read_csv(here("00_data", "derived", "questionnaire_scoring", "AQ_final_data.csv"))


exp_stat %>%
    group_by(Experimenter) %>%
    summarise(n = n())

# turn menst into long format
menst_long <- menst %>%
    pivot_longer(cols = -ID, names_to = "session", values_to = "phase") %>%
    filter(!is.na(phase))

# merge df and menst by ID and session
df_menst <- left_join(df, menst_long, by = c("ID", "session"))

df_exp_stat <- left_join(df_menst, exp_stat, by = "ID")

df_exp_stat <- df_exp_stat %>% filter(!(Status %in% c("excluded", "side effects")))

df_exp_stat %>%
    group_by(ID, session) %>%
    distinct(ID, session) %>%
    summarise(n = n())

# merge AQ data
df_exp_stat <- left_join(df_exp_stat, AQ_final_data, by = "ID")

write_csv(df_exp_stat, here("00_data", "derived", "analysis", "df_analysis_pub.csv"))
