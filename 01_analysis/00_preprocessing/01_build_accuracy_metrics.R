## Manuscript section: preprocessing
## Analysis family: eye-tracking quality / accuracy metrics
## Original source path: legacy/scripts_archive/original_scripts/precision:accuracy.R
## Primary input dataset(s): legacy/data_archive/original_layout/accuracy/data/*.json
## Primary output(s): data/derived/quality/acc_df.rds, data/derived/quality/acc_df.csv
## Known TODOs: raw accuracy JSON files still live in the archived original layout; script has only been path-cleaned
## Scientific logic unchanged from source except for output-path cleanup

library(tidyverse)
library(jsonlite)

filenames <- list.files(
  "legacy/data_archive/original_layout/accuracy/data",
  full.names = TRUE
)

acc_df <- lapply(filenames, fromJSON) %>% bind_rows()

acc_df <- acc_df %>%
  mutate(
    Sub_ID = str_extract(recording, "(?<=recordings_sorted\\/)[^\\/]+"),
    session = str_extract(recording, "(Session\\s[0-9]+)") %>%
      str_replace_all(" ", "_") %>%
      str_to_lower()
  ) %>% 
  unnest_wider(accuracy, names_sep = "_") %>%
  unnest_wider(precision, names_sep = "_") %>% 
  # Move 'subject_id' and 'session_id' to the front
  select(Sub_ID, session, accuracy_degrees, precision_degrees) %>% 
  arrange(Sub_ID, session)

write_rds(acc_df, "data/derived/quality/acc_df.rds")
write.csv(acc_df, "data/derived/quality/acc_df.csv", row.names = FALSE)

#visualize acc
hist(acc_df$accuracy_degrees)
mean(acc_df$accuracy_degrees)
sd(acc_df$accuracy_degrees)

#bad recordings
acc_df %>% filter(accuracy_degrees > 3)

#best recordings
acc_df %>% filter(accuracy_degrees < 1)

#visualize precision
hist(acc_df$precision_degrees)
mean(acc_df$precision_degrees)
sd(acc_df$precision_degrees)

#Filter outliers in accuracy

# Calculate the IQR for accuracy_degrees
Q1 <- quantile(acc_df$accuracy_degrees, 0.25)
Q3 <- quantile(acc_df$accuracy_degrees, 0.75)
IQR <- Q3 - Q1

# Define the upper bound for what you consider a normal range
upper_bound <- Q3 + 0.5 * IQR

# Filter the dataframe to exclude values above this upper bound
acc_df %>% filter(!accuracy_degrees <= upper_bound)
