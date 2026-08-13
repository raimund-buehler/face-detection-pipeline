## Manuscript section: Main eye-tracking model sensitivity / detectability
## Analysis family: simulation-based sensitivity analysis for hypothesized medication x AOI effects
## Primary input dataset(s): 00_data/derived/preprocessing/duration_data.csv; 00_data/derived/analysis/df_analysis_pub_prep.csv
## Primary output(s): benchmark timings, targeted simulation estimates for AOI-specific medication effects

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

output_subdir <- get_arg_value("output_subdir", "power_simulation_medication_aoi_effect")
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
formula_null <- percentage_fix_duration ~ Group * fix * Sex + session
focal_term <- "fix:medication"
interaction_term <- "Group:fix:medication"

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
phi_base <- predict(fit_null, type = "disp")

simulate_scenario_data <- function(max_shift) {
  both_rows <- df$medication == "BOTH"
  oxt_rows <- df$medication == "OXT"
  nal_rows <- df$medication == "NAL"
  eyes_rows <- df$fix == "Eyes"
  background_rows <- df$fix == "Background"

  mu_sim <- mu_base

  # Hypothesized ordered medication effect:
  # BOTH has the largest effect; OXT and NAL each have half that magnitude.
  mu_sim[both_rows & eyes_rows] <- mu_sim[both_rows & eyes_rows] + max_shift
  mu_sim[both_rows & background_rows] <- mu_sim[both_rows & background_rows] - max_shift

  mu_sim[oxt_rows & eyes_rows] <- mu_sim[oxt_rows & eyes_rows] + (max_shift / 2)
  mu_sim[oxt_rows & background_rows] <- mu_sim[oxt_rows & background_rows] - (max_shift / 2)

  mu_sim[nal_rows & eyes_rows] <- mu_sim[nal_rows & eyes_rows] + (max_shift / 2)
  mu_sim[nal_rows & background_rows] <- mu_sim[nal_rows & background_rows] - (max_shift / 2)

  mu_sim <- clamp_beta(mu_sim)

  df %>%
    mutate(
      percentage_fix_duration = simulate_beta_response(mu_sim, phi_base),
      max_shift = max_shift
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
      p_fix_medication = NA_real_,
      p_group_fix_medication = NA_real_,
      error = conditionMessage(fit)
    ))
  }

  tibble(
    converged = fit$fit$convergence == 0,
    pdHess = isTRUE(fit$sdr$pdHess),
    p_fix_medication = extract_anova_p(fit, focal_term),
    p_group_fix_medication = extract_anova_p(fit, interaction_term),
    error = NA_character_
  )
}

benchmark_fit <- system.time(
  glmmTMB(formula_main, family = beta_family(), data = df)
)

benchmark_sim_refit <- system.time({
  sim_df_benchmark <- simulate_scenario_data(max_shift = 0.08)
  refit_and_extract(sim_df_benchmark)
})

benchmark_df <- tibble(
  step = c("fit_empirical_model", "simulate_and_refit_max_shift_0.08"),
  elapsed_sec = c(unname(benchmark_fit["elapsed"]), unname(benchmark_sim_refit["elapsed"])),
  user_sec = c(unname(benchmark_fit["user.self"]), unname(benchmark_sim_refit["user.self"])),
  system_sec = c(unname(benchmark_fit["sys.self"]), unname(benchmark_sim_refit["sys.self"]))
)

set.seed(20260507)
nsim <- as.integer(get_arg_value("nsim", "500"))
scenario_arg <- get_arg_value("shifts", "0,0.08,0.12,0.16")
scenario_grid <- tibble(
  max_shift = as.numeric(str_split(scenario_arg, ",", simplify = TRUE))
) %>%
  filter(!is.na(max_shift)) %>%
  arrange(max_shift)

simulation_results <- map_dfr(seq_len(nrow(scenario_grid)), function(i) {
  shift <- scenario_grid$max_shift[[i]]
  map_dfr(seq_len(nsim), function(rep_id) {
    sim_df <- simulate_scenario_data(shift)
    stats <- refit_and_extract(sim_df)
    stats %>%
      mutate(
        max_shift = shift,
        simulation = rep_id
      ) %>%
      relocate(max_shift, simulation)
  })
}) %>%
  mutate(
    scenario_label = case_when(
      max_shift == 0 ~ "Null baseline",
      TRUE ~ paste0("BOTH +/-", scales::percent(max_shift, accuracy = 1), ", OXT/NAL +/-", scales::percent(max_shift / 2, accuracy = 1))
    )
  )

simulation_summary <- simulation_results %>%
  group_by(max_shift) %>%
  summarise(
    successful_fits = sum(converged & pdHess, na.rm = TRUE),
    total_runs = n(),
    detect_rate_fix_medication = mean(p_fix_medication < 0.05, na.rm = TRUE),
    detect_rate_group_fix_medication = mean(p_group_fix_medication < 0.05, na.rm = TRUE),
    median_p_fix_medication = median(p_fix_medication, na.rm = TRUE),
    median_p_group_fix_medication = median(p_group_fix_medication, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    scenario_label = case_when(
      max_shift == 0 ~ "Null baseline",
      TRUE ~ paste0("BOTH +/-", scales::percent(max_shift, accuracy = 1), ", OXT/NAL +/-", scales::percent(max_shift / 2, accuracy = 1))
    )
  )

reviewer_table <- simulation_summary %>%
  transmute(
    `Scenario` = scenario_label,
    `Maximum injected shift (BOTH)` = scales::percent(max_shift, accuracy = 1),
    `Successful fits` = paste0(successful_fits, "/", total_runs),
    `Detection rate: AOI x medication` = scales::percent(detect_rate_fix_medication, accuracy = 1),
    `Median p: AOI x medication` = sprintf("%.3f", median_p_fix_medication),
    `Detection rate: Group x AOI x medication` = scales::percent(detect_rate_group_fix_medication, accuracy = 1),
    `Median p: Group x AOI x medication` = sprintf("%.3f", median_p_group_fix_medication)
  )

power_plot <- ggplot(simulation_summary, aes(x = max_shift, y = detect_rate_fix_medication)) +
  geom_line(color = "#0072B2", linewidth = 0.9) +
  geom_point(color = "#0072B2", size = 2) +
  scale_x_continuous(
    breaks = scenario_grid$max_shift,
    labels = scales::label_percent(accuracy = 1)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::label_percent(accuracy = 1)
  ) +
  labs(
    title = "Detectability of Hypothesized AOI-specific Medication Effects",
    subtitle = "Ordered pattern: BOTH > OXT/NAL > PLA, imposed equally across groups",
    x = "Maximum injected shift in BOTH condition",
    y = "Detection rate for AOI x Medication (p < .05)"
  ) +
  theme_classic(base_size = 10)

write_csv(benchmark_df, file.path(output_dir, "01_benchmark_timings.csv"))
write_csv(simulation_results, file.path(output_dir, "02_simulation_results.csv"))
write_csv(simulation_summary, file.path(output_dir, "03_simulation_summary.csv"))
write_csv(reviewer_table, file.path(output_dir, "04_reviewer_summary_table.csv"))
ggsave(
  filename = file.path(output_dir, "05_power_curve.png"),
  plot = power_plot,
  width = 8,
  height = 5,
  dpi = 160
)

notes <- c(
  "Simulation-based sensitivity analysis for hypothesized AOI-specific medication effects",
  "",
  "Goal:",
  "Assess detectability of the preregistered-style medication hypothesis in the primary beta model.",
  "",
  "Model:",
  "percentage_fix_duration ~ Group * fix * medication + Group * fix * Sex + session",
  "",
  "Simulation design:",
  "- Null baseline: fitted means and precision from a reduced beta model without medication terms.",
  paste0(
    "- Reduced-model AIC = ", round(AIC(fit_null), 3),
    "; full-model AIC = ", round(AIC(fit_main), 3), "."
  ),
  "- Medication effect was imposed equally across groups.",
  "- In BOTH, eye-directed fixation was increased and background-directed fixation decreased by the full shift amount.",
  "- In OXT and NAL, eye-directed fixation was increased and background-directed fixation decreased by half that amount.",
  "- Placebo was left unchanged.",
  paste0("- Simulated scenarios: ", paste(scales::percent(scenario_grid$max_shift, accuracy = 1), collapse = ', '), "."),
  paste0("- Simulations per scenario: ", nsim, "."),
  "",
  "Interpretation:",
  "- The focal omnibus target is AOI x medication, because the hypothesized medication effect was AOI-specific rather than a global medication main effect.",
  "- The Group x AOI x medication term is also tracked to confirm that the injected effect was not group-specific."
)
write_text_output(notes, output_dir, "00_readme.txt")

print(benchmark_df)
print(simulation_summary)
