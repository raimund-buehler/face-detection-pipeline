# Manuscript section: Supplementary / control analyses
# Analysis family: session coverage check
# Original source path: legacy/scripts_archive/original_scripts/publication/analysis/session_check.R
# Primary input dataset(s): data/derived/analysis/df_analysis_pub_prep.csv
# Primary output(s): session-level inspection output in console
# Known TODOs: limited scope script; decide later whether to keep as QC-only or formal supplementary check
# Scientific logic note: copied from source without changing scientific logic

library(ggridges)
library(ggdist)
library(rlang)
library(here)
library(ggh4x)
library(ggsignif)
library(effectsize)
library(performance)
source(here("01_analysis", "04_figures", "shared", "custom_plot_settings.R"))

# DATA PREP

df <- read_csv(here("00_data", "derived", "analysis", "df_analysis_pub_prep.csv"))

df %>% distinct(ID, session)
