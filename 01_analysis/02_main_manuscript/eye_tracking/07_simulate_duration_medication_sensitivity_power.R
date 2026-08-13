## Manuscript section: Main eye-tracking model sensitivity / detectability
## Analysis family: simulation-based sensitivity analysis for medication null findings
## Primary input dataset(s): 00_data/derived/preprocessing/duration_data.csv; 00_data/derived/analysis/df_analysis_pub_prep.csv
## Primary output(s): benchmark timings, smoke-test sensitivity estimates, optional power curve

library(tidyverse)
library(glmmTMB)
library(car)
library(here)
source(here("01_analysis", "shared", "model_output_utils.R"))

args <- commandArgs(trailingOnly = TRUE)

get_arg_value <- function(flag, default = NULL) {
  prefix <- paste0("--", flag, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) {
    return(default)
  }
  sub(prefix, "", hit[[1]])
}

normalize_session <- function(x) {
  x %>%
    str_trim() %>%
    str_to_lower() %>%
    str_replace_all(" ", "_")
}

clamp_beta <- function(x, eps = 1e-4) {
  pmin(pmax(x, eps), 1 - eps)
}

simulate_beta_response <- function(mu, phi) {
  rbeta(length(mu), shape1 = mu * phi, shape2 = (1 - mu) * phi)
}

extract_anova_p <- function(model, term) {
  out <- tryCatch(car::Anova(model, type = 3), error = function(e) NULL)
  if (is.null(out) || !(term %in% rownames(out))) {
    return(NA_real_)
  }
  p_col <- intersect(c("Pr(>Chisq)", "Pr(>Chi)", "Pr(>F)", "p.value"), colnames(out))
  if (length(p_col) == 0) {
    return(NA_real_)
  }
  as.numeric(out[term, p_col[[1]]])
}

output_subdir <- get_arg_value("output_subdir", "power_simulation")
output_dir <- here(
  "02_outputs", "model_outputs", "main_manuscript", "eye_tracking",
  "group_aoi_medication_duration", "sensitivity", output_subdir
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

session_metadata <- read_csv(
  here("00_data", "derived", "analysis", "df_analysis_pub_prep.csv"),
  show_col_types = FALSE
) %>%
  distinct(ID, session, Sex, Group, medication) %>%
  mutate(session_norm = normalize_session(session))

df <- read_csv(
  here("00_data", "derived", "preprocessing", "duration_data.csv"),
  show_col_types = FALSE
) %>%
  mutate(ID = Sub_ID, session_norm = normalize_session(session)) %>%
  left_join(session_metadata, by = c("ID", "session_norm")) %>%
  mutate(session = coalesce(session.y, session.x)) %>%
  filter(!is.na(Group), Sex %in% c("f", "m")) %>%
  mutate(
    fix = factor(
      fix,
      levels = c("fix_on_eyes", "fix_on_mouth", "fix_on_face", "fix_on_background"),
      labels = c("Eyes", "Mouth", "Face", "Background"),
      ordered = TRUE
    ),
    medication = factor(
      medication,
      levels = c("NAL/OXT", "NAL/PLA", "PLA/OXT", "PLA/PLA"),
      labels = c("BOTH", "NAL", "OXT", "PLA")
    ),
    Group = factor(Group),
    Sex = factor(Sex),
    session = factor(session)
  ) %>%
  select(ID, session, Sex, Group, medication, fix, percentage_fix_duration)

formula_main <- percentage_fix_duration ~ Group * fix * medication + Group * fix * Sex + session
formula_null <- percentage_fix_duration ~ Group * fix + Group * medication + fix * medication + Group * fix * Sex + session
focal_term <- "Group:fix:medication"

fit_main <- glmmTMB(
  formula_main,
  family = beta_family(),
  data = df
)

fit_null <- glmmTMB(
  formula_null,
  family = beta_family(),
  data = df
)

mu_base <- fitted(fit_null, type = "response")
phi_base <- predict(fit_main, type = "disp")

simulate_scenario_data <- function(response_shift) {
  active_med <- df$medication %in% c("BOTH", "NAL", "OXT")
  asd_rows <- df$Group == "ASD"
  eyes_rows <- df$fix == "Eyes"
  background_rows <- df$fix == "Background"

  mu_sim <- mu_base
  mu_sim[asd_rows & active_med & eyes_rows] <- mu_sim[asd_rows & active_med & eyes_rows] + response_shift
  mu_sim[asd_rows & active_med & background_rows] <- mu_sim[asd_rows & active_med & background_rows] - response_shift
  mu_sim <- clamp_beta(mu_sim)

  df %>%
    mutate(
      percentage_fix_duration = simulate_beta_response(mu_sim, phi_base),
      scenario_shift = response_shift
    )
}

refit_and_extract <- function(sim_df) {
  fit <- tryCatch(
    glmmTMB(formula_main, family = beta_family(), data = sim_df),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(tibble(
      converged = FALSE,
      pdHess = FALSE,
      p_group_fix_medication = NA_real_,
      p_medication = NA_real_,
      error = conditionMessage(fit)
    ))
  }

  tibble(
    converged = fit$fit$convergence == 0,
    pdHess = isTRUE(fit$sdr$pdHess),
    p_group_fix_medication = extract_anova_p(fit, focal_term),
    p_medication = extract_anova_p(fit, "medication"),
    error = NA_character_
  )
}

benchmark_fit <- system.time(
  glmmTMB(formula_main, family = beta_family(), data = df)
)

benchmark_sim_refit <- system.time({
  sim_df_benchmark <- simulate_scenario_data(response_shift = 0.05)
  refit_and_extract(sim_df_benchmark)
})

benchmark_df <- tibble(
  step = c("fit_empirical_model", "simulate_and_refit_shift_0.05"),
  elapsed_sec = c(unname(benchmark_fit["elapsed"]), unname(benchmark_sim_refit["elapsed"])),
  user_sec = c(unname(benchmark_fit["user.self"]), unname(benchmark_sim_refit["user.self"])),
  system_sec = c(unname(benchmark_fit["sys.self"]), unname(benchmark_sim_refit["sys.self"]))
)

set.seed(20260507)
nsim <- as.integer(get_arg_value("nsim", "100"))
scenario_arg <- get_arg_value("shifts", "0,0.02,0.05,0.08,0.12,0.16,0.20,0.24")
scenario_grid <- tibble(
  scenario_shift = as.numeric(str_split(scenario_arg, ",", simplify = TRUE))
) %>%
  filter(!is.na(scenario_shift)) %>%
  arrange(scenario_shift)

simulation_results <- map_dfr(seq_len(nrow(scenario_grid)), function(i) {
  shift <- scenario_grid$scenario_shift[[i]]
  map_dfr(seq_len(nsim), function(rep_id) {
    sim_df <- simulate_scenario_data(shift)
    stats <- refit_and_extract(sim_df)
    stats %>%
      mutate(
        scenario_shift = shift,
        simulation = rep_id
      ) %>%
      relocate(scenario_shift, simulation)
  })
}) %>%
  mutate(
    scenario_label = case_when(
      scenario_shift == 0 ~ "Null baseline",
      TRUE ~ paste0("+/-", scales::percent(scenario_shift, accuracy = 1), " shift")
    )
  )

simulation_summary <- simulation_results %>%
  group_by(scenario_shift) %>%
  summarise(
    successful_fits = sum(converged & pdHess, na.rm = TRUE),
    total_runs = n(),
    detect_rate_group_fix_medication = mean(p_group_fix_medication < 0.05, na.rm = TRUE),
    detect_rate_medication = mean(p_medication < 0.05, na.rm = TRUE),
    median_p_group_fix_medication = median(p_group_fix_medication, na.rm = TRUE),
    median_p_medication = median(p_medication, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    scenario_label = case_when(
      scenario_shift == 0 ~ "Null baseline",
      TRUE ~ paste0("+/-", scales::percent(scenario_shift, accuracy = 1), " shift")
    )
  )

threshold_candidates <- simulation_summary %>%
  filter(scenario_shift > 0)

threshold_row <- threshold_candidates %>%
  filter(detect_rate_group_fix_medication >= 0.8) %>%
  slice_head(n = 1)

if (nrow(threshold_row) == 0) {
  threshold_summary <- tibble(
    target = "Group x AOI x medication",
    target_power = 0.8,
    estimated_shift = NA_real_,
    lower_grid_shift = max(threshold_candidates$scenario_shift, na.rm = TRUE),
    upper_grid_shift = NA_real_,
    note = "80% power not reached within simulated shift grid"
  )
} else {
  upper_shift <- threshold_row$scenario_shift[[1]]
  upper_power <- threshold_row$detect_rate_group_fix_medication[[1]]
  lower_row <- threshold_candidates %>%
    filter(scenario_shift < upper_shift) %>%
    slice_tail(n = 1)

  if (nrow(lower_row) == 0) {
    interpolated_shift <- upper_shift
    lower_shift <- NA_real_
  } else {
    lower_shift <- lower_row$scenario_shift[[1]]
    lower_power <- lower_row$detect_rate_group_fix_medication[[1]]
    interpolated_shift <- approx(
      x = c(lower_power, upper_power),
      y = c(lower_shift, upper_shift),
      xout = 0.8,
      ties = "ordered"
    )$y
  }

  threshold_summary <- tibble(
    target = "Group x AOI x medication",
    target_power = 0.8,
    estimated_shift = interpolated_shift,
    lower_grid_shift = lower_shift,
    upper_grid_shift = upper_shift,
    note = "Linear interpolation between nearest simulated grid points"
  )
}

power_plot <- ggplot(simulation_summary, aes(x = scenario_shift, y = detect_rate_group_fix_medication)) +
  geom_line(color = "#0072B2", linewidth = 0.9) +
  geom_point(color = "#0072B2", size = 2) +
  scale_x_continuous(
    breaks = scenario_grid$scenario_shift,
    labels = scales::label_percent(accuracy = 1)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::label_percent(accuracy = 1)
  ) +
  labs(
    title = "Detectability of Group x AOI x Medication Effects",
    subtitle = "Null baseline from reduced model; active-medication shifts applied in ASD for Eyes (+) and Background (-)",
    x = "Injected response-scale shift per active medication condition",
    y = "Detection rate for Group x AOI x Medication (p < .05)"
  ) +
  theme_classic(base_size = 10)

reviewer_table <- simulation_summary %>%
  transmute(
    `Scenario` = scenario_label,
    `Injected shift` = scales::percent(scenario_shift, accuracy = 1),
    `Successful fits` = paste0(successful_fits, "/", total_runs),
    `Detection rate: Group x AOI x medication` = scales::percent(detect_rate_group_fix_medication, accuracy = 1),
    `Median p: Group x AOI x medication` = sprintf("%.3f", median_p_group_fix_medication),
    `Detection rate: medication main effect` = scales::percent(detect_rate_medication, accuracy = 1),
    `Median p: medication main effect` = sprintf("%.3f", median_p_medication)
  )

write_csv(benchmark_df, file.path(output_dir, "01_benchmark_timings.csv"))
write_csv(simulation_results, file.path(output_dir, "02_simulation_results.csv"))
write_csv(simulation_summary, file.path(output_dir, "03_simulation_summary.csv"))
write_csv(reviewer_table, file.path(output_dir, "05_reviewer_summary_table.csv"))
write_csv(threshold_summary, file.path(output_dir, "06_estimated_shift_for_80_power.csv"))
ggsave(
  filename = file.path(output_dir, "04_power_curve.png"),
  plot = power_plot,
  width = 8,
  height = 5,
  dpi = 160
)

notes <- c(
  "Simulation-based sensitivity analysis for the main fixation-duration beta model",
  "",
  "Goal:",
  "Assess how often the omnibus Group x AOI x medication term is detected under plausible active-medication attenuation scenarios.",
  "",
  "Model:",
  "percentage_fix_duration ~ Group * fix * medication + Group * fix * Sex + session",
  "",
  "Simulation design:",
  paste0(
    "- Null baseline: fitted means from a reduced beta model without the Group x AOI x medication interaction."
  ),
  paste0(
    "- Reduced-model AIC = ", round(AIC(fit_null), 3),
    "; full-model AIC = ", round(AIC(fit_main), 3), "."
  ),
  "- For ASD rows under active medications (BOTH, NAL, OXT), the response mean was shifted upward for Eyes and downward for Background by the same absolute amount.",
  "- This creates a pharmacological attenuation pattern relative to placebo while leaving the rest of the design unchanged.",
  paste0("- Simulated scenarios: ", paste(scales::percent(scenario_grid$scenario_shift, accuracy = 1), collapse = ', '), "."),
  paste0("- Simulations per scenario: ", nsim, "."),
  "",
  "Interpretation:",
  "- This is a sensitivity / detectability analysis, not observed power based on the realized p-value.",
  "- The 0% scenario estimates the false-positive rate for the focal omnibus term under the reduced-model null.",
  "- The 80% threshold is estimated by interpolation across the simulated grid and should be treated as approximate."
)
write_text_output(notes, output_dir, "00_readme.txt")

print(benchmark_df)
print(simulation_summary)
print(threshold_summary)
