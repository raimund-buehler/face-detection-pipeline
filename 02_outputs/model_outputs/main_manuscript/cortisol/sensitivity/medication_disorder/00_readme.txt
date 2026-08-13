Cortisol sensitivity analysis: ASD medication/disorder type and T2-T3 cortisol reactivity

Primary model:
T2T3 ~ Group * medication + (1 | anon_id)

Sensitivity question:
Do ASD medication/disorder type indicators moderate the study-medication effect on T2-T3?

Type indicators are ASD-specific dummy variables. Controls are coded 0 on all indicators.
This avoids treating absence of ASD medication/disorder records in controls as a comparable clinical category.

Primary reviewer-facing output:
05_type_moderation_lrt_tests.csv
