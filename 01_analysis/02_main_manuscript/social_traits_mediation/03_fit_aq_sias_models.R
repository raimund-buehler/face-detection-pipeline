# Manuscript section: Dimensional analyses
# Analysis family: AQ and SIAS models
# Original source path: scripts/publication/analysis/AQ_SIAS.r
# Primary input dataset(s): data prepared by social_traits_mediation/shared_aq_sias_plot_utils.R from manuscript analysis tables
# Primary output(s): model summaries and emtrends results for AQ, SIAS, social scales, and IRI
# Known TODOs: AQ/SIAS helper remains the provenance entry point for dataset preparation; figure generation moved to analysis/04_figures/main/social_traits_mediation/03_plot_aq_sias_panels.R
# Scientific logic note: scientific logic and model formulas are unchanged from source; plotting was separated for structural cleanup

library(here)
library(tidyverse)
library(lmerTest)
library(emmeans)
library(effectsize)
library(performance)
source(here("01_analysis", "02_main_manuscript", "social_traits_mediation", "shared_aq_sias_plot_utils.R"))
source(here("01_analysis", "shared", "model_output_utils.R"))

df <- prepare_aq_sias_data()

fit_model <- function(data, var_name, dep_name) {
  formula_str <- paste(dep_name, "~", var_name, "*fix*medication + Sex + session + (1 | ID)")
  formula <- as.formula(formula_str)

  model <- lmer(formula, data = data)
  output_dir <- here("02_outputs", "model_outputs", "main_manuscript", "social_traits_mediation", "aq_sias")

  print(summary(model))
  print(anova(model))
  write_model_anova(
    model,
    output_dir,
    paste0("anova_", tolower(var_name), "_", tolower(dep_name), ".csv"),
    paste0(var_name, "_", dep_name)
  )

  emtrends_result <- emtrends(model, ~fix, var = var_name, infer = c(TRUE, TRUE))
  print(emtrends_result)

  list(
    model = model,
    summary = summary(model),
    anova = anova(model),
    emtrends = emtrends_result
  )
}

df$AQ_scaled <- datawizard::standardise(df$AQ, center = TRUE, scale = FALSE)

result_aq_sqrt <- fit_model(df, "AQ_scaled", "sqrt_rate_all")

check_model(result_aq_sqrt$model)
check_collinearity(result_aq_sqrt$model)
car::vif(result_aq_sqrt$model)

model_vcov <- vcov(result_aq_sqrt$model)
condition_number <- kappa(model_vcov)
condition_number

eigenvalues <- eigen(model_vcov)$values
eigenvalues

# Eyes
t_to_r(t = -2.265, df_error = 234)

# Background
t_to_r(t = 2.954, df_error = 234)

result_sias_sqrt <- fit_model(df, "SIAS", "sqrt_rate_all")

# Background
t_to_r(t = 3.439, df_error = 235)

result_soc_sqrt <- fit_model(df, "Soz_Scales", "sqrt_rate_all")
result_iri_sqrt <- fit_model(df, "IRI", "sqrt_rate_all")

summary(result_aq_sqrt$model)
