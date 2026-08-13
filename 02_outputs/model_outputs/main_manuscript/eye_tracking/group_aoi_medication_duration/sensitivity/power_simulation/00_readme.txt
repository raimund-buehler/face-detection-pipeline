Simulation-based sensitivity analysis for the main fixation-duration beta model

Goal:
Assess how often the omnibus Group x AOI x medication term is detected under plausible active-medication attenuation scenarios.

Model:
percentage_fix_duration ~ Group * fix * medication + Group * fix * Sex + session

Simulation design:
- Null baseline: fitted means from a reduced beta model without the Group x AOI x medication interaction.
- Reduced-model AIC = -1026.463; full-model AIC = -1014.595.
- For ASD rows under active medications (BOTH, NAL, OXT), the response mean was shifted upward for Eyes and downward for Background by the same absolute amount.
- This creates a pharmacological attenuation pattern relative to placebo while leaving the rest of the design unchanged.
- Simulated scenarios: 0%, 2%, 5%, 8%, 12%, 16%, 20%, 24%.
- Simulations per scenario: 100.

Interpretation:
- This is a sensitivity / detectability analysis, not observed power based on the realized p-value.
- The 0% scenario estimates the false-positive rate for the focal omnibus term under the reduced-model null.
- The 80% threshold is estimated by interpolation across the simulated grid and should be treated as approximate.
