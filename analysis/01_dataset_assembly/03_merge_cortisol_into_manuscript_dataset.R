# Manuscript section: Dataset assembly for cortisol manuscript analyses
# Analysis family: dataset assembly
# Original source path: legacy/scripts_archive/original_scripts/publication/merging/merge_cortisol.r
# Primary input dataset(s): data/derived/analysis/df_analysis_pub.csv; data/derived/analysis/df_analysis_pub_prep.csv; data/derived/assembly/cortisol_joined.csv
# Primary output(s): data/derived/analysis/df_cortisol_merged.csv
# Known TODOs: provenance for df_analysis_pub_prep.csv is unresolved; cortisol source harmonization remains pending
# Scientific logic note: copied from source without changing scientific logic

library(tidyverse)
library(lmerTest)
library(here)

df_analysis_pub <- read_csv(here("00_data", "derived", "analysis", "df_analysis_pub.csv"))
df_analysis_pub_prep <- read_csv(here("00_data", "derived", "analysis", "df_analysis_pub_prep.csv"))

df_cortisol <- read_csv(here("00_data", "derived", "assembly", "cortisol_joined.csv")) %>%
    rename(ID = Sub_ID) %>%
    select(-lab_notes, -cortisol_2, -medication, -cv, -Sex, -Group)

# check which IDs in the cortisol data are not in the df_analysis_pub data
df_cortisol %>%
    anti_join(df_analysis_pub, by = "ID") %>%
    distinct(ID)

# repair known ID typos (mapping supplied in config/paths.R)
df_cortisol <- df_cortisol %>% mutate(ID = ifelse(ID %in% names(id_corrections), id_corrections[ID], ID))

# make session lowercase and replace whitespace with underscore
df_cortisol <- df_cortisol %>%
    mutate(session = tolower(session)) %>%
    mutate(session = str_replace_all(session, " ", "_"))

df_analysis_pub_prep_cort <- left_join(df_analysis_pub_prep, df_cortisol, by = c("ID", "session"))

df_analysis_pub_prep %>%
    count(ID, session) # Original row counts per session

df_analysis_pub_prep_cort %>%
    count(ID, session)

# no duplicates
df_duplicates <- df_analysis_pub_prep_cort %>%
    group_by(across(-c(cortisol_1, timepoint))) %>%
    filter(n() > 3)

write_csv(df_analysis_pub_prep_cort, here("00_data", "derived", "analysis", "df_cortisol_merged.csv"))
