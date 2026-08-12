## Manuscript section: manuscript / supplementary tables
## Analysis family: Sympa summary table
## Original source path: legacy/data_archive/original_layout/data_questionnaire/Sympa.r
## Primary input dataset(s): legacy/data_archive/original_layout/data_questionnaire/Sympa_raw.xlsx, data/raw/questionnaires/participants_1812 (1).xlsx, data/derived/analysis/df_analysis_pub_prep.csv
## Primary output(s): data/derived/intermediate/Sympa_merged_clean.xlsx, outputs/tables/questionnaires/sympa_summary_table.xlsx
## Known TODOs: Sympa raw source still lives in the archived original layout
## Scientific logic unchanged from source except for canonical file-path cleanup

library(tidyverse)
library(readxl)
library(writexl)
library(lmerTest)
library(here)

sympa <- read_xlsx(here("legacy", "data_archive", "original_layout", "data_questionnaire", "Sympa_raw.xlsx")) %>%
    rename(ID = SF29_01, attr = SF73, good = SF74, nice = SF75, openminded = SF76, session = SF77) %>%
    select(ID, session, attr, good, nice, openminded) %>%
    mutate(ID = tolower(ID)) %>%
    arrange(ID, session)

# Apply the transformation to columns 3 to 6
sympa[3:6] <- lapply(sympa[3:6], function(col) {
    # Convert the factor to a character vector
    col <- as.character(col)

    # Extract the numeric part or assign NA for "[NA] Not answered"
    col <- ifelse(col == "[NA] Not answered",
        NA,
        as.numeric(gsub(".*\\[(\\d+)\\].*", "\\1", col))
    )

    return(col)
})

part <- read_xlsx(here("00_data", "raw", "questionnaires", "participants_1812 (1).xlsx")) %>%
    rename(ID = `Personal ID`) %>%
    arrange(ID) %>%
    select(ID, Group, Gender, Sessions, "Session status")

merged <- left_join(sympa, part, by = "ID")

write_xlsx(merged, here("00_data", "derived", "intermediate", "Sympa_merged_clean.xlsx"))

merged <- read_xlsx(here("00_data", "derived", "intermediate", "Sympa_merged_clean.xlsx")) %>% arrange(Group, ID, session) %>% mutate(ID = ifelse(ID %in% names(id_corrections), id_corrections[ID], ID))

merged %>% distinct(ID) %>% summarise(n())

overall_df <- read_csv(here("00_data", "derived", "analysis", "df_analysis_pub_prep.csv"))

#check excludes and remove
overall_df %>% distinct(ID) %>% summarise(n())
excludes <- anti_join(merged, overall_df, by = "ID") %>% distinct(ID)
merged <- merged %>% filter(!ID %in% excludes$ID)
merged %>% distinct(ID) %>% summarise(n())

means <- merged %>%
    group_by(Group) %>%
    summarise(across(attr:openminded, list(mean = mean, sd = sd), na.rm = TRUE), .groups = "drop")

# perform tests

# Create a list to store the models
lmm_models <- list()

# Loop through columns 3 to 6
for (i in 3:6) {
    # Define the formula for the LMM
    formula <- as.formula(paste0("merged[[", i, "]] ~ Group + (1 | ID)"))

    # Fit the LMM
    model <- lmer(formula, data = merged)

    # Store the model in the list
    lmm_models[[colnames(merged)[i]]] <- model
}

# Display the summary of each model
lapply(lmm_models, summary)

# Extract the t-values and p-values from the models
t_values <- lapply(lmm_models, function(model) {
    t_value <- summary(model)$coefficients[2, "t value"]
    p_value <- summary(model)$coefficients[2, "Pr(>|t|)"]
    return(c(t_value, p_value))
})

t_values <- as_tibble(t_values)

t <- t(t_values)

t <- as.data.frame(t)

colnames(t) <- c("t_value", "p_value")
t$Variable <- rownames(t)

t <- as_tibble(t) %>% mutate(t_value = as.numeric(t_value), p_value = as.numeric(p_value))


df <- merged %>% select(Group, attr:openminded)

means_sds <- df %>%
    group_by(Group) %>%
    summarise(across(everything(), list(mean = mean, sd = sd), na.rm = TRUE), .groups = "drop") %>%
    pivot_longer(-Group, names_to = c("Variable", ".value"), names_pattern = "(.*)_(.*)") %>%
    pivot_wider(names_from = Group, values_from = c(mean, sd))

result <- means_sds %>%
    left_join(t, by = "Variable") %>%
    arrange(Variable) %>%
    # reorder columns, ASD before CTRLs
    select(Variable, contains("ASD"), contains("CTRL"), everything())

result

write_xlsx(result, here("02_outputs", "tables", "questionnaires", "sympa_summary_table.xlsx"))

           
