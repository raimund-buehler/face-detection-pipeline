# Manuscript section: Exploratory duration-based longitudinal analysis
# Analysis family: non-overlapping time-bin random-structure sensitivity check
# Primary input dataset(s): 00_data/derived/preprocessing/duration_sliding_window_data.csv
# Primary output(s): model comparison and omnibus tests across random structures

library(tidyverse)
library(glmmTMB)
library(car)
library(here)

beta_squeeze <- function(x) {
  n <- sum(!is.na(x))
  ((x * (n - 1)) + 0.5) / n
}

output_dir <- here(
  "02_outputs", "model_outputs", "main_manuscript", "eye_tracking",
  "timecourse_duration_nonoverlap_random_structure"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- read_csv(
  here("00_data", "derived", "preprocessing", "duration_sliding_window_data.csv"),
  show_col_types = FALSE
) %>%
  filter(!is.na(window_fix_duration_prop), !is.na(medication), Sex %in% c("f", "m")) %>%
  # Select the 10%-wide windows that start at 0%, 10%, ..., 90%.
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

candidate_random <- tibble(
  random_structure = c(
    "id_intercept",
    "id_plus_subject_session",
    "id_plus_medication_diag",
    "id_plus_session_diag",
    "id_plus_subject_session_plus_medication_diag"
  ),
  random_rhs = c(
    "(1 | ID)",
    "(1 | ID) + (1 | subject_session)",
    "(1 | ID) + diag(0 + medication | ID)",
    "(1 | ID) + diag(0 + session | ID)",
    "(1 | ID) + (1 | subject_session) + diag(0 + medication | ID)"
  )
)

fit_candidate <- function(fix_label, random_structure, random_rhs) {
  df_block <- df %>% filter(fix == fix_label)
  formula_text <- paste(
    "window_fix_duration_prop_beta ~",
    "Group * progress_bin * medication + Sex + session +",
    random_rhs
  )
  model_formula <- as.formula(formula_text)

  model <- tryCatch(
    glmmTMB(
      model_formula,
      family = beta_family(),
      data = df_block,
      control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
    ),
    error = function(e) e
  )

  if (inherits(model, "error")) {
    return(list(
      summary = tibble(
        fix = fix_label,
        random_structure = random_structure,
        formula = formula_text,
        converged = FALSE,
        pdHess = FALSE,
        AIC = NA_real_,
        BIC = NA_real_,
        logLik = NA_real_,
        error = conditionMessage(model)
      ),
      anova = tibble()
    ))
  }

  sdr <- model$sdr
  converged <- model$fit$convergence == 0
  pdHess <- isTRUE(sdr$pdHess)

  anova_table <- tryCatch(
    car::Anova(model, type = 3) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("term") %>%
      as_tibble() %>%
      mutate(fix = fix_label, random_structure = random_structure),
    error = function(e) {
      tibble(
        term = NA_character_,
        Chisq = NA_real_,
        Df = NA_real_,
        `Pr(>Chisq)` = NA_real_,
        fix = fix_label,
        random_structure = random_structure,
        anova_error = conditionMessage(e)
      )
    }
  )

  list(
    summary = tibble(
      fix = fix_label,
      random_structure = random_structure,
      formula = formula_text,
      converged = converged,
      pdHess = pdHess,
      AIC = AIC(model),
      BIC = BIC(model),
      logLik = as.numeric(logLik(model)),
      error = NA_character_
    ),
    anova = anova_table
  )
}

results <- pmap(
  expand_grid(fix = c("Eyes", "Background"), candidate_random),
  ~ fit_candidate(..1, ..2, ..3)
)

model_comparison <- map_dfr(results, "summary") %>%
  group_by(fix) %>%
  mutate(delta_AIC = AIC - min(AIC, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(fix, delta_AIC)

anova_results <- map_dfr(results, "anova") %>%
  relocate(fix, random_structure, term)

key_terms <- anova_results %>%
  filter(term %in% c("medication", "session", "Group:medication", "Group:progress_bin:medication")) %>%
  select(fix, random_structure, term, Chisq, Df, `Pr(>Chisq)`) %>%
  arrange(fix, random_structure, term)

write_csv(model_comparison, file.path(output_dir, "01_random_structure_model_comparison.csv"))
write_csv(anova_results, file.path(output_dir, "02_random_structure_anova.csv"))
write_csv(key_terms, file.path(output_dir, "03_random_structure_key_terms.csv"))

write_lines(
  c(
    "Non-overlapping time-bin random-structure sensitivity check",
    "",
    "Input:",
    "duration_sliding_window_data.csv, restricted to 10%-wide windows starting at 0%, 10%, ..., 90%",
    "",
    "Fixed-effects part:",
    "window_fix_duration_prop_beta ~ Group * progress_bin * medication + Sex + session",
    "",
    "Random structures compared:",
    paste(candidate_random$random_structure, candidate_random$random_rhs, sep = ": "),
    "",
    "Interpretation note:",
    "This is a sensitivity check for whether medication/session effects and Group x medication effects depend on random-structure assumptions."
  ),
  file.path(output_dir, "00_model_notes.txt")
)

print(model_comparison)
print(key_terms)
