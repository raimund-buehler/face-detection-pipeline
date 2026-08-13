Cortisol-reactivity to eye fixation-duration robustness checks

Outcome: percentage_fix_duration_beta for fix_on_eyes.
Base model: outcome ~ Group * predictor + medication + session + (1 | ID).

Robustness variants:
raw full sample
winsorized 5th/95th percentile predictor
rank-transformed predictor
trimmed 2.5th/97.5th percentile predictor
trimmed 5th/95th percentile predictor

See 06_interaction_robustness_summary.csv for the Group x cortisol robustness summary.
