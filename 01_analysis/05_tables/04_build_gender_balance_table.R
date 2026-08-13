## Manuscript section: manuscript / supplementary tables
## Analysis family: gender balance table
## Original source path: legacy/data_archive/original_layout/data_questionnaire/gender.r
## Primary input dataset(s): data/derived/analysis/df_analysis_pub_prep.csv
## Primary output(s): outputs/tables/questionnaires/gender_balance_counts.csv, outputs/tables/questionnaires/gender_balance_tests.csv
## Known TODOs: this remains a promoted historical check script rather than a manuscript-critical pipeline step
## Scientific logic unchanged from source except for canonical file-path cleanup and using the analyzed sample directly

library(tidyverse)
library(here)

df <- read_csv(here("00_data", "derived", "analysis", "df_analysis_pub_prep.csv")) %>%
  distinct(ID, Group, Sex) %>%
  filter(Sex %in% c("f", "m"))

gender_table <- df %>%
  count(Group, Sex, name = "n") %>%
  pivot_wider(names_from = Sex, values_from = n, values_fill = 0) %>%
  arrange(Group)

chi_square_matrix <- as.matrix(gender_table[, -1])
rownames(chi_square_matrix) <- gender_table$Group

chisq_corrected <- chisq.test(chi_square_matrix)
chisq_uncorrected <- chisq.test(chi_square_matrix, correct = FALSE)
fisher_result <- fisher.test(chi_square_matrix)

print(gender_table)
print(chisq_corrected)
print(chisq_uncorrected)
print(fisher_result)

write_csv(
  gender_table,
  here("02_outputs", "tables", "questionnaires", "gender_balance_counts.csv")
)

write_csv(
  tibble(
    test = c("pearson_chisq_yates", "pearson_chisq_uncorrected", "fisher_exact"),
    statistic = c(
      unname(chisq_corrected$statistic),
      unname(chisq_uncorrected$statistic),
      NA_real_
    ),
    parameter = c(
      unname(chisq_corrected$parameter),
      unname(chisq_uncorrected$parameter),
      NA_real_
    ),
    p_value = c(
      chisq_corrected$p.value,
      chisq_uncorrected$p.value,
      fisher_result$p.value
    ),
    estimate = c(
      NA_real_,
      NA_real_,
      unname(fisher_result$estimate)
    ),
    conf_low = c(
      NA_real_,
      NA_real_,
      fisher_result$conf.int[1]
    ),
    conf_high = c(
      NA_real_,
      NA_real_,
      fisher_result$conf.int[2]
    ),
    method = c(
      chisq_corrected$method,
      chisq_uncorrected$method,
      fisher_result$method
    )
  ),
  here("02_outputs", "tables", "questionnaires", "gender_balance_tests.csv")
)
