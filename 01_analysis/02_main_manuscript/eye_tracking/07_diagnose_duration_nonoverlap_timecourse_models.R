# Manuscript section: Exploratory duration-based longitudinal analysis
# Analysis family: non-overlapping time-bin model diagnostics
# Primary input dataset(s): saved glmmTMB model objects from random-structure sensitivity check
# Primary output(s): DHARMa and performance diagnostics for selected timecourse models

library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(performance)
library(here)

beta_squeeze <- function(x) {
  n <- sum(!is.na(x))
  ((x * (n - 1)) + 0.5) / n
}

output_dir <- here(
  "02_outputs", "model_outputs", "main_manuscript", "eye_tracking",
  "timecourse_duration_nonoverlap_random_structure"
)
model_dir <- file.path(output_dir, "model_objects")
diag_dir <- file.path(output_dir, "diagnostics")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

model_index <- read_csv(
  file.path(model_dir, "00_saved_model_index.csv"),
  show_col_types = FALSE
)

df <- read_csv(
  here("00_data", "derived", "preprocessing", "duration_sliding_window_data.csv"),
  show_col_types = FALSE
) %>%
  filter(!is.na(window_fix_duration_prop), !is.na(medication), Sex %in% c("f", "m")) %>%
  filter(abs((progress_start * 10) - round(progress_start * 10)) < 1e-8) %>%
  mutate(
    ID = factor(Sub_ID),
    fix = factor(
      fix,
      levels = c("fix_on_eyes", "fix_on_mouth", "fix_on_face", "fix_on_background"),
      labels = c("Eyes", "Mouth", "Face", "Background")
    ),
    medication = factor(
      medication,
      levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA"),
      labels = c("BOTH", "NAL", "OXT", "PLA")
    ),
    Group = factor(Group, levels = c("ASD", "CTRL")),
    Sex = factor(Sex),
    session = factor(session),
    progress_bin = factor(
      sprintf("%02d", as.integer(round(progress_start * 10)) + 1),
      levels = sprintf("%02d", 1:10),
      ordered = TRUE
    ),
    subject_session = interaction(ID, session, drop = TRUE),
    window_fix_duration_prop_beta = beta_squeeze(window_fix_duration_prop)
  ) %>%
  filter(fix %in% c("Eyes", "Background"))

run_diagnostics <- function(model_row, n_sim = 500) {
  model_name <- paste(model_row$fix, tools::file_path_sans_ext(model_row$file), sep = "_")
  model_data <- df %>% filter(fix == model_row$fix)
  formula_text <- str_replace(
    model_row$formula,
    "^y ~",
    "window_fix_duration_prop_beta ~"
  )

  fit_model <- glmmTMB(
    as.formula(formula_text),
    family = beta_family(),
    data = model_data,
    control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
  )

  sim_res <- DHARMa::simulateResiduals(fit_model, n = n_sim, plot = FALSE)

  png(file.path(diag_dir, paste0(model_name, "_01_dharma_overview.png")), width = 1600, height = 1200, res = 150)
  plot(sim_res)
  dev.off()

  png(file.path(diag_dir, paste0(model_name, "_02_fitted_vs_observed.png")), width = 1600, height = 1200, res = 150)
  plot(
    fitted(fit_model),
    model_data$window_fix_duration_prop_beta,
    xlab = "Fitted values",
    ylab = "Observed fixation-duration proportion",
    pch = 16,
    col = scales::alpha("#0B7A6B", 0.25)
  )
  abline(0, 1, col = "red", lwd = 2)
  dev.off()

  overdispersion <- performance::check_overdispersion(fit_model)

  tibble(
    model = model_name,
    fix = model_row$fix,
    saved_file = model_row$file,
    formula = formula_text,
    n_observations = nrow(model_data),
    aic = AIC(fit_model),
    bic = BIC(fit_model),
    dharma_uniformity_p = DHARMa::testUniformity(sim_res)$p.value,
    dharma_dispersion_p = DHARMa::testDispersion(sim_res)$p.value,
    dharma_outlier_p = DHARMa::testOutliers(sim_res)$p.value,
    dharma_zero_inflation_p = DHARMa::testZeroInflation(sim_res)$p.value,
    check_singularity = as.logical(performance::check_singularity(fit_model)),
    overdispersion_ratio = overdispersion$dispersion_ratio,
    overdispersion_p = overdispersion$p_value
  )
}

diagnostic_summary <- map_dfr(
  seq_len(nrow(model_index)),
  ~ run_diagnostics(model_index[.x, ])
)

write_csv(diagnostic_summary, file.path(diag_dir, "00_timecourse_diagnostic_summary.csv"))

summary_lines <- c(
  "Non-overlapping timecourse model diagnostics",
  "",
  paste("DHARMa simulations per model:", 500),
  "",
  apply(
    diagnostic_summary,
    1,
    function(x) {
      paste0(
        x[["model"]],
        ": uniformity p = ", sprintf("%.3f", as.numeric(x[["dharma_uniformity_p"]])),
        ", dispersion p = ", sprintf("%.3f", as.numeric(x[["dharma_dispersion_p"]])),
        ", outlier p = ", sprintf("%.3f", as.numeric(x[["dharma_outlier_p"]])),
        ", zero-inflation p = ", sprintf("%.3f", as.numeric(x[["dharma_zero_inflation_p"]])),
        ", overdispersion ratio = ", sprintf("%.3f", as.numeric(x[["overdispersion_ratio"]])),
        ", overdispersion p = ", sprintf("%.3f", as.numeric(x[["overdispersion_p"]]))
      )
    }
  )
)

writeLines(summary_lines, file.path(diag_dir, "00_timecourse_diagnostic_summary.txt"))
