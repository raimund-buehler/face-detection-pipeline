# Manuscript section: Cortisol analyses
# Analysis family: T2T3 / MinMax cortisol reactivity preparation
# Original source path: legacy/cortisol_reports/markdown_source_bundle/4_cortisol_minmax_prep.R
# Primary input dataset(s): data/derived/analysis/df_with_auc.csv
# Primary output(s): data/derived/analysis/df_cortisol_min_max.csv
# Known TODOs: confirm whether T2T3 and MinMax should be harmonized to a single variable name across the repo
# Scientific logic note: extracted from the manuscript preprocessing script without changing the calculation

library(here)
library(tidyverse)

calculate_minmax <- function(df_session) {
  if (any(df_session$timepoint == 2) & any(df_session$timepoint == 3)) {
    cortisol_tp2 <- df_session$cortisol_1[df_session$timepoint == 2]
    cortisol_tp3 <- df_session$cortisol_1[df_session$timepoint == 3]
    minmax <- cortisol_tp3 - cortisol_tp2

    return(data.frame(MinMax = minmax))
  }

  data.frame(MinMax = NA)
}

df <- read_csv(here("00_data", "derived", "analysis", "df_with_auc.csv")) %>%
  filter(fix == "Eyes")

minmax_df <- df %>%
  group_by(ID, session) %>%
  do(calculate_minmax(.))

df_with_minmax <- left_join(df, minmax_df, by = c("ID", "session"))

df_cortisol_min_max <- df_with_minmax %>%
  group_by(ID, Group, session, medication) %>%
  summarise(
    rate_all = first(rate_all),
    sqrt_rate_all = first(sqrt_rate_all),
    AUCg = first(AUCg),
    AUCi = first(AUCi),
    MinMax = first(MinMax),
    .groups = "drop"
  )

df_cortisol_min_max$AUCg_scaled <- datawizard::standardise(df_cortisol_min_max$AUCg)
df_cortisol_min_max$AUCi_scaled <- datawizard::standardise(df_cortisol_min_max$AUCi)
df_cortisol_min_max$MinMax_scaled <- datawizard::standardise(df_cortisol_min_max$MinMax)

write_csv(
  df_cortisol_min_max,
  here("00_data", "derived", "analysis", "df_cortisol_min_max.csv")
)
