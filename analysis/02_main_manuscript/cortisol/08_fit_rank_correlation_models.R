# Manuscript section: Cortisol analyses
# Analysis family: rank correlations
# Original source path: legacy/cortisol_reports/markdown_source_bundle/5_rank_correlation.rmd
# Primary input dataset(s): data/derived/analysis/df_cortisol_min_max.csv
# Primary output(s): Spearman and clustered Kendall correlation results
# Known TODOs: source markdown had this section disabled; keeping it active here for traceable code-only access
# Scientific logic note: extracted from the manuscript markdown workflow without changing the analysis steps

library(tidyverse)
library(coin)
library(here)
source(here("01_analysis", "shared", "model_output_utils.R"))

df <- read_csv(here("00_data", "derived", "analysis", "df_cortisol_min_max.csv")) %>%
  rename(T2T3 = MinMax, T2T3_scaled = MinMax_scaled)

auc_subject <- df %>%
  group_by(ID, Group) %>%
  summarise(
    AUCi = mean(AUCi, na.rm = TRUE),
    eye_rate = mean(rate_all, na.rm = TRUE),
    .groups = "drop"
  )

t2t3_subject <- df %>%
  group_by(ID, Group) %>%
  summarise(
    T2T3 = mean(T2T3, na.rm = TRUE),
    eye_rate = mean(rate_all, na.rm = TRUE),
    .groups = "drop"
  )

cor.test(
  auc_subject$AUCi[auc_subject$Group == "ASD"],
  auc_subject$eye_rate[auc_subject$Group == "ASD"],
  method = "spearman"
)
output_dir <- here("02_outputs", "model_outputs", "main_manuscript", "cortisol", "rank_correlations")

write_htest_result(
  cor.test(
    auc_subject$AUCi[auc_subject$Group == "ASD"],
    auc_subject$eye_rate[auc_subject$Group == "ASD"],
    method = "spearman"
  ),
  output_dir,
  "01_spearman_auci_eye_rate_asd.csv",
  "spearman_auci_eye_rate_asd"
)

write_htest_result(
  cor.test(
    auc_subject$AUCi[auc_subject$Group == "CTRL"],
    auc_subject$eye_rate[auc_subject$Group == "CTRL"],
    method = "spearman"
  ),
  output_dir,
  "02_spearman_auci_eye_rate_ctrl.csv",
  "spearman_auci_eye_rate_ctrl"
)

write_htest_result(
  cor.test(
    t2t3_subject$T2T3[t2t3_subject$Group == "ASD"],
    t2t3_subject$eye_rate[t2t3_subject$Group == "ASD"],
    method = "spearman"
  ),
  output_dir,
  "03_spearman_t2t3_eye_rate_asd.csv",
  "spearman_t2t3_eye_rate_asd"
)

write_htest_result(
  cor.test(
    t2t3_subject$T2T3[t2t3_subject$Group == "CTRL"],
    t2t3_subject$eye_rate[t2t3_subject$Group == "CTRL"],
    method = "spearman"
  ),
  output_dir,
  "04_spearman_t2t3_eye_rate_ctrl.csv",
  "spearman_t2t3_eye_rate_ctrl"
)

df_subject <- df %>%
  group_by(ID, Group) %>%
  summarise(
    T2T3 = mean(T2T3, na.rm = TRUE),
    AUCi = mean(AUCi, na.rm = TRUE),
    eye_rate = mean(rate_all, na.rm = TRUE),
    .groups = "drop"
  )

ctrl_subject <- df_subject %>% filter(Group == "CTRL")
asd_subject <- df_subject %>% filter(Group == "ASD")

write_htest_result(
  cor.test(ctrl_subject$AUCi, ctrl_subject$eye_rate, method = "kendall"),
  output_dir,
  "05_kendall_auci_eye_rate_ctrl.csv",
  "kendall_auci_eye_rate_ctrl"
)
write_htest_result(
  cor.test(asd_subject$AUCi, asd_subject$eye_rate, method = "kendall"),
  output_dir,
  "06_kendall_auci_eye_rate_asd.csv",
  "kendall_auci_eye_rate_asd"
)
write_htest_result(
  cor.test(ctrl_subject$T2T3, ctrl_subject$eye_rate, method = "kendall"),
  output_dir,
  "07_kendall_t2t3_eye_rate_ctrl.csv",
  "kendall_t2t3_eye_rate_ctrl"
)
write_htest_result(
  cor.test(asd_subject$T2T3, asd_subject$eye_rate, method = "kendall"),
  output_dir,
  "08_kendall_t2t3_eye_rate_asd.csv",
  "kendall_t2t3_eye_rate_asd"
)

ctrl_df <- df %>%
  filter(Group == "CTRL") %>%
  mutate(ID = as.factor(ID)) %>%
  group_by(ID) %>%
  filter(n() > 1) %>%
  ungroup()

kendall_ctrl_t2t3 <- independence_test(
  rate_all ~ T2T3 | ID,
  data = ctrl_df,
  teststat = "quad",
  distribution = approximate(B = 10000)
)

kendall_ctrl_auci <- independence_test(
  rate_all ~ AUCi | ID,
  data = ctrl_df,
  teststat = "quad",
  distribution = approximate(B = 10000)
)

asd_df <- df %>%
  filter(Group == "ASD") %>%
  mutate(ID = as.factor(ID)) %>%
  group_by(ID) %>%
  filter(n() > 1) %>%
  ungroup()

kendall_asd_t2t3 <- independence_test(
  rate_all ~ T2T3 | ID,
  data = asd_df,
  teststat = "quad",
  distribution = approximate(B = 10000)
)

kendall_asd_auci <- independence_test(
  rate_all ~ AUCi | ID,
  data = asd_df,
  teststat = "quad",
  distribution = approximate(B = 10000)
)

kendall_ctrl_t2t3
kendall_ctrl_auci
kendall_asd_t2t3
kendall_asd_auci
