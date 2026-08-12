# Manuscript section: Cortisol analyses
# Analysis family: AUC dataset assembly
# Original source path: scripts/publication/analysis/cortisol_analysis/3_cortisol_AUC_mean.r
# Primary input dataset(s): data/derived/analysis/df_cortisol_merged.csv
# Primary output(s): data/derived/analysis/df_with_auc.csv; data/derived/analysis/df_cortisol_merged_with_auc.csv
# Known TODOs: confirm whether df_with_auc.csv remains needed once downstream consumers are fully standardized
# Scientific logic note: extracted from the existing AUC preparation logic without changing the calculation

library(tidyverse)
library(pracma)
library(here)

calculate_auc <- function(df_session) {
  timepoints <- df_session$timepoint
  cortisol_levels <- df_session$cortisol_1

  aucg <- trapz(timepoints, cortisol_levels)
  baseline <- cortisol_levels[1]
  auci <- trapz(timepoints, cortisol_levels - baseline)

  data.frame(AUCg = aucg, AUCi = auci)
}

df <- read_csv(here("00_data", "derived", "analysis", "df_cortisol_merged.csv")) %>%
  filter(fix == "Eyes")

auc_df <- df %>%
  group_by(ID, session) %>%
  do(calculate_auc(.))

df_with_auc <- left_join(df, auc_df, by = c("ID", "session"))

write_csv(df_with_auc, here("00_data", "derived", "analysis", "df_with_auc.csv"))

df_cortisol_merged_with_auc <- df_with_auc %>%
  group_by(ID, Group, session, medication) %>%
  summarise(
    rate_all = first(rate_all),
    sqrt_rate_all = first(sqrt_rate_all),
    AUCg = first(AUCg),
    AUCi = first(AUCi),
    .groups = "drop"
  )

df_cortisol_merged_with_auc$AUCg_scaled <- datawizard::standardise(df_cortisol_merged_with_auc$AUCg)
df_cortisol_merged_with_auc$AUCi_scaled <- datawizard::standardise(df_cortisol_merged_with_auc$AUCi)

write_csv(
  df_cortisol_merged_with_auc,
  here("00_data", "derived", "analysis", "df_cortisol_merged_with_auc.csv")
)
