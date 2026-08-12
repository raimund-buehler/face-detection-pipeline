## Manuscript section: manuscript / supplementary tables
## Analysis family: PANAS summary table
## Original source path: legacy/data_archive/original_layout/data_questionnaire/PANAS.r
## Primary input dataset(s): data/raw/questionnaires/PANAS_ausgewertet_final_corrected.xlsx, data/derived/intermediate/merged.csv, data/derived/analysis/df_analysis_pub_prep.csv
## Primary output(s): outputs/tables/questionnaires/panas_summary_table.xlsx
## Known TODOs: inferential PANAS model code is preserved below the table build for provenance
## Scientific logic unchanged from source except for canonical file-path cleanup

library(tidyverse)
library(readxl)
library(writexl)
library(here)
library(lmerTest)
library(emmeans)

# PANAS

panas <- read_xlsx(here("00_data", "raw", "questionnaires", "PANAS_ausgewertet_final_corrected.xlsx"))

# rename column Sub_ID to ID in panas and fill NA's in Group with ASD
panas <- panas %>% rename(ID = Sub_ID, pos = `mean Positive Affect`, neg = `mean Negative Affect`)

# read data from overall_model analysis (without excludes)
overall_df <- read_csv(here("00_data", "derived", "analysis", "df_analysis_pub_prep.csv"))

#check excludes and remove
overall_df %>% distinct(ID, session) %>% summarise(n())
excludes <- anti_join(panas, overall_df, by = "ID") %>% distinct(ID, session)
panas <- panas %>% filter(!ID %in% excludes$ID)
panas %>% distinct(ID) %>% summarise(n())

merged <- read_csv2(here("00_data", "derived", "intermediate", "merged.csv"))

# join panas and merged data to get group information (select only group and ID)
panas <- left_join(panas %>% select(-5), merged %>% select(ID, Group), by = "ID") %>% mutate(Group = ifelse(is.na(Group), "ASD", Group))

panas <- panas %>%
  mutate(session = tolower(gsub(" ", "_", session)))

overall_merge <- overall_df %>% distinct(ID, session, medication)

# join medication
panas <- left_join(panas, overall_merge %>% select(ID, session, medication), by = c("ID", "session"))

panas <- panas %>% filter(!is.na(medication))

panas_t <- panas %>%
    pivot_longer(
        cols = starts_with("mean"),
        names_to = "Affect_Type",
        values_to = "Score"
    ) %>%
    mutate(Affect_Type = str_remove(Affect_Type, "^mean ")) %>%
    group_by(session, Affect_Type) %>%
    summarise(
        t_test = list(t.test(Score ~ Group, na.rm = TRUE)),
        .groups = "drop"
    ) %>%
    mutate(
        t_value = map_dbl(t_test, ~ .x$statistic),
        p_value = map_dbl(t_test, ~ .x$p.value)
    ) %>%
    select(-t_test) %>%
    arrange(Affect_Type, session)

panas <-
    panas %>%
    group_by(Group, session) %>%
    summarise(across(c("mean Positive Affect", "mean Negative Affect"), list(mean = mean, sd = sd), na.rm = TRUE), .groups = "drop")

# Pivot longer for Positive and Negative Affect
panas_long <- panas %>%
    pivot_longer(
        cols = starts_with("mean"), # Select columns that start with "mean"
        names_to = c("Affect_Type", ".value"), # Split names into Affect_Type and mean/sd
        names_pattern = "mean (.*)_(.*)" # Specify the pattern to split the names
    ) %>%
    # sort by Affect_Type
    arrange(Affect_Type, session)

# pivot wider (for group)
panas_wide <- panas_long %>%
    pivot_wider(names_from = Group, values_from = c(mean, sd)) %>%
    select(session, Affect_Type, contains("ASD"), contains("CTRL"), everything())

# merge panas_wide with panas_t
panas_result <- panas_wide %>%
    left_join(panas_t, by = c("session", "Affect_Type")) %>%
    arrange(Affect_Type, session)

panas_result

write_xlsx(panas_result, here("02_outputs", "tables", "questionnaires", "panas_summary_table.xlsx"))


## By Medication

# Fit a mixed model (random intercept for subjects)
# Positive Affect
model <- lmer(pos ~ medication * Group * session + (1 | ID), data = panas)

# Show model summary
summary(model)
anova(model)

emmeans(model, pairwise ~ Group)
emmeans(model, pairwise ~ Group | session)
emmeans(model, pairwise ~ Group | medication)

# Fit a mixed model (random intercept for subjects)
# Negative Affect
model <- lmer(neg ~ medication * Group * session + (1 | ID), data = panas)

# Show model summary
summary(model)
anova(model)

emmeans(model, pairwise ~ Group)
emmeans(model, pairwise ~ Group | session)
emmeans(model, pairwise ~ Group | medication)
