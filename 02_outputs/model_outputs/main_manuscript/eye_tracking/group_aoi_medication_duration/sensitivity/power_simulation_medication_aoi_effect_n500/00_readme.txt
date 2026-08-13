Simulation-based sensitivity analysis for hypothesized AOI-specific medication effects

Goal:
Assess detectability of the preregistered-style medication hypothesis in the primary beta model.

Model:
percentage_fix_duration ~ Group * fix * medication + Group * fix * Sex + session

Simulation design:
- Null baseline: fitted means and precision from a reduced beta model without medication terms.
- Reduced-model AIC = -1050.689; full-model AIC = -1014.595.
- Medication effect was imposed equally across groups.
- In BOTH, eye-directed fixation was increased and background-directed fixation decreased by the full shift amount.
- In OXT and NAL, eye-directed fixation was increased and background-directed fixation decreased by half that amount.
- Placebo was left unchanged.
- Simulated scenarios: 0%, 8%, 12%, 16%.
- Simulations per scenario: 500.

Interpretation:
- The focal omnibus target is AOI x medication, because the hypothesized medication effect was AOI-specific rather than a global medication main effect.
- The Group x AOI x medication term is also tracked to confirm that the injected effect was not group-specific.
