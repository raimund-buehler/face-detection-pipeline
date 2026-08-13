# Machine-specific paths live in config/paths.R (not committed).
source(here::here("config", "paths.R"))
## Manuscript section: manuscript / supplementary tables
## Analysis family: questionnaire summary table
## Original source path: legacy/data_archive/original_layout/data_questionnaire/Numeric_quest.r
## Primary input dataset(s): data/raw/questionnaires/Fragebögen_ausgewertet.xlsx, data/raw/questionnaires/participants_1812 (1).xlsx, data/derived/questionnaire_scoring/AQ_final_data.csv, data/derived/analysis/df_analysis_pub_prep.csv
## Primary output(s): data/derived/intermediate/merged.csv, outputs/tables/questionnaires/questionnaire_summary_table.xlsx
## Known TODOs: legacy-style descriptive formatting is preserved; this remains a table-builder rather than a core assembly step
## Scientific logic unchanged from source except for canonical file-path cleanup

library(tidyverse)
library(purrr)
library(broom)
library(janitor)
library(readxl)
library(writexl)
library(here)
source(here("01_analysis", "04_figures", "shared", "format_ttest.R"))

quest <- read_xlsx(here("00_data", "raw", "questionnaires", "Fragebögen_ausgewertet.xlsx"))
part <- read_xlsx(here("00_data", "raw", "questionnaires", "participants_1812 (1).xlsx"))
aq <- read_csv(here("00_data", "derived", "questionnaire_scoring", "AQ_final_data.csv")) %>% mutate(ID = ifelse(ID %in% names(id_corrections), id_corrections[ID], ID))

# rename column Personal ID to ID in part and Sub_ID to ID in quest
part <- part %>%
    rename(ID = `Personal ID`) %>%
    mutate(ID = ifelse(ID %in% names(id_corrections), id_corrections[ID], ID))
quest <- quest %>%
    rename(ID = Sub_ID) %>%
    mutate(ID = ifelse(ID %in% names(id_corrections), id_corrections[ID], ID))

# Merge the two data frames
merged <- left_join(part, quest, by = "ID")
merged <- left_join(merged, aq, by = "ID") %>% arrange(ID)

# save the merged data frame as excel file
write_excel_csv2(merged, here("00_data", "derived", "intermediate", "merged.csv"))

# read data from overall_model analysis (without excludes)
overall_df <- read_csv(here("00_data", "derived", "analysis", "df_analysis_pub_prep.csv"))

#check excludes and remove
overall_df %>% distinct(ID) %>% summarise(n())
excludes <- anti_join(merged, overall_df, by = "ID") %>% distinct(ID)
merged <- merged %>% filter(!ID %in% excludes$ID)

# check if the columns are numeric
numeric_columns <- sapply(merged, is.numeric)

# select the numeric columns using dplyr
df <- merged %>%
    select(Group, which(numeric_columns), -Sessions) %>%
    rename(AQOld = `AQ-K Score`)

df <- df %>%
    clean_names()

# check nas
df %>% filter(aq_old != overall_aq)
df %>% filter(overall_aq != aq_old)

df %>%
    summarise(
        na_in_aq_old = sum(is.na(aq_old)),
        na_in_overall_aq = sum(is.na(overall_aq))
    )


# Perform t-tests
t_test_results <- df %>%
    reframe(across(-group, ~ {
        t_test <- t.test(. ~ group, na.rm = TRUE)
        tibble(
            t_value = t_test$statistic,
            p_value = t_test$p.value
        )
    }, .names = "{.col}")) %>%
    pivot_longer(cols = everything(), names_to = "Variable_Metric", values_to = "Value") %>%
    separate(Variable_Metric, into = c("Variable", "Metric"), sep = "-", fill = "right") %>%
    pivot_wider(names_from = Metric, values_from = Value, names_sep = ".", names_repair = "unique") %>%
    unnest(everything(), names_sep = ".")


# Calculate means and standard deviations
means_sds <- df %>%
    group_by(group) %>%
    summarise(across(everything(), list(mean = mean, sd = sd), na.rm = TRUE), .groups = "drop") %>%
    pivot_longer(-group, names_to = c("Variable", ".value"), names_pattern = "(.*)_(.*)") %>%
    pivot_wider(names_from = group, values_from = c(mean, sd))

# Combine t-test results with means and SDs
result <- means_sds %>%
    left_join(t_test_results, by = "Variable") %>%
    arrange(Variable) %>%
    # reorder columns, ASD before CTRLs
    select(Variable, contains("ASD"), contains("CTRL"), everything())

# Step 1: Replace any missing or empty column names with generic names
names(result) <- make.names(names(result), unique = TRUE)

# Step 2: Rename the first column explicitly to "variable"
names(result)[1] <- "variable"

#correct names and order

result <- result %>%
  filter(!variable == "overall_aq") %>% 
  mutate(
    variable = case_when(
      variable == "age" ~ "Age",
      variable == "aq_old" ~ "AQ - Autism Spectrum Quotient",
      variable == "comm_rec" ~ "AQ - Communication and Reciprocity",
      variable == "soz_spont" ~ "AQ - Social Interaction and Spontaneity",
      variable == "imag" ~ "AQ - Imagination",
      variable == "soz_scales" ~ "AQ - Social Scales",
      variable == "bdi_a_score" ~ "BDI-A - Beck Depression Inventory",
      variable == "iq" ~ "Verbal IQ",
      variable == "iri_score" ~ "IRI - Interpersonal Reactivity Index",
      variable == "iri_subscale_empathetic_concern" ~ "IRI - Empathetic Concern",
      variable == "iri_subscale_fantasy" ~ "IRI - Fantasy",
      variable == "iri_subscale_personal_distress" ~ "IRI - Personal Distress",
      variable == "iri_subscale_perspective_taking" ~ "IRI - Perspective Taking",
      variable == "sias_score" ~ "SIAS - Social Interaction Anxiety Scale",
      variable == "stai_t_score" ~ "STAI-t - State-Trait Anxiety Inventory (trait)",
      TRUE ~ variable  # Keep the original value if no match
    )
  )

desired_order <- c(
  "Age",
  "AQ - Autism Spectrum Quotient",
  "AQ - Communication and Reciprocity",
  "AQ - Social Interaction and Spontaneity",
  "AQ - Imagination",
  "AQ - Social Scales",
  "BDI-A - Beck Depression Inventory",
  "Verbal IQ",
  "IRI - Interpersonal Reactivity Index",
  "IRI - Empathetic Concern",
  "IRI - Fantasy",
  "IRI - Personal Distress",
  "IRI - Perspective Taking",
  "SIAS - Social Interaction Anxiety Scale",
  "STAI-t - State-Trait Anxiety Inventory (trait)"
)

result <- result %>%
  mutate(order = match(variable, desired_order)) %>%  # Add a helper column with the order value
  arrange(order) %>%  # Arrange rows based on the order column
  select(-order)  # Remove the helper column

# View the final result
print(result)

write_xlsx(result, here("02_outputs", "tables", "questionnaires", "questionnaire_summary_table.xlsx"))
