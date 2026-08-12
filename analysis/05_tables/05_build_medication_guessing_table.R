## Manuscript section: supplementary / control tables
## Analysis family: medication guessing chi-square table
## Original source path: legacy/data_archive/original_layout/data_questionnaire/Medication/chisq_med.r
## Primary input dataset(s): legacy/data_archive/original_layout/data_questionnaire/Medication/merged_chi_square.csv, legacy/data_archive/original_layout/data_questionnaire/Medication/sessions_fixed.csv
## Primary output(s): outputs/tables/medication_controls/medication_guessing_chisq_overall.csv, outputs/tables/medication_controls/medication_guessing_chisq_by_medication.csv
## Known TODOs: this script still depends on archived medication-questionnaire exports
## Scientific logic unchanged from source except for canonical file-path cleanup and explicit table export

library(tidyverse)
library(readxl)
library(here)

merged_chi_square <- read_csv(here("legacy", "data_archive", "original_layout", "data_questionnaire", "Medication", "merged_chi_square.csv")) %>%
    rename(ID = Code, session = Session) %>%
    mutate(ID = ifelse(ID %in% names(id_corrections), id_corrections[ID], ID))

sessions_fixed <- read_csv(here("legacy", "data_archive", "original_layout", "data_questionnaire", "Medication", "sessions_fixed.csv")) %>% rename(ID = Sub_ID)

sessions_fixed$session <- as.integer(gsub("Session ", "", sessions_fixed$session))

merged_chi_square <- merged_chi_square %>%
    left_join(sessions_fixed, by = c("ID", "session"))

# Aggregate the data across all sessions for Pille
pille_table <- table(merged_chi_square$`Guessed Pille`, merged_chi_square$`Received Pill`)

# Aggregate the data across all sessions for Nebulizer
nebulizer_table <- table(merged_chi_square$`Guessed Nebulizer`, merged_chi_square$`Received Spray`)

# Perform chi-square test on the aggregated data for Pille
pille_test <- chisq.test(pille_table)

# Perform chi-square test on the aggregated data for Nebulizer
nebulizer_test <- chisq.test(nebulizer_table)

# Output the results
print(pille_test)
print(nebulizer_test)

overall_results <- tibble(
    measure = c("pille", "nebulizer"),
    statistic = c(unname(pille_test$statistic), unname(nebulizer_test$statistic)),
    parameter = c(unname(pille_test$parameter), unname(nebulizer_test$parameter)),
    p_value = c(pille_test$p.value, nebulizer_test$p.value),
    method = c(pille_test$method, nebulizer_test$method)
)

## OVER SESSIONS
sessions <- unique(merged_chi_square$session)

# OVER MEDICATION
sessions <- unique(merged_chi_square$medication)

# Initialize lists to store results
pille_results <- list()
nebulizer_results <- list()

# Loop through each session
for (cond in sessions) {
    # Filter the data for the current session
    session_data <- subset(merged_chi_square, medication == cond)

    # Create contingency tables for Pille and Nebulizer
    pille_table <- table(session_data$`Guessed Pille`, session_data$`Received Pill`)
    nebulizer_table <- table(session_data$`Guessed Nebulizer`, session_data$`Received Spray`)

    # Perform chi-square tests
    pille_test <- chisq.test(pille_table)
    nebulizer_test <- chisq.test(nebulizer_table)

    # Store results in the lists
    pille_results[[paste("Session", cond)]] <- pille_test
    nebulizer_results[[paste("Session", cond)]] <- nebulizer_test
}

# Print results for each session
print("Pille Results by Session")
print(pille_results)

print("Nebulizer Results by Session")
print(nebulizer_results)

by_medication_results <- bind_rows(
    lapply(names(pille_results), function(cond) {
        tibble(
            medication = cond,
            measure = "pille",
            statistic = unname(pille_results[[cond]]$statistic),
            parameter = unname(pille_results[[cond]]$parameter),
            p_value = pille_results[[cond]]$p.value,
            method = pille_results[[cond]]$method
        )
    }),
    lapply(names(nebulizer_results), function(cond) {
        tibble(
            medication = cond,
            measure = "nebulizer",
            statistic = unname(nebulizer_results[[cond]]$statistic),
            parameter = unname(nebulizer_results[[cond]]$parameter),
            p_value = nebulizer_results[[cond]]$p.value,
            method = nebulizer_results[[cond]]$method
        )
    })
)

write_csv(overall_results, here("02_outputs", "tables", "medication_controls", "medication_guessing_chisq_overall.csv"))
write_csv(by_medication_results, here("02_outputs", "tables", "medication_controls", "medication_guessing_chisq_by_medication.csv"))
