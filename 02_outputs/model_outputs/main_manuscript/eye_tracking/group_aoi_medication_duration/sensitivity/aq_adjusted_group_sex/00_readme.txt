AQ-adjusted Group x AOI x Sex sensitivity

AQ imbalance model:
AQ ~ Group * Sex

Duration sensitivity models:
aq_complete_main: percentage_fix_duration_beta ~ Group * fix * medication + Group * fix * Sex + session
aq_adjusted: percentage_fix_duration_beta ~ Group * fix * medication + Group * fix * Sex + AQ_z * fix + session
aq_adjusted_by_sex: percentage_fix_duration_beta ~ Group * fix * medication + Group * fix * Sex + AQ_z * fix * Sex + session

AQ_z is standardized across the AQ-complete analysis sample.
Outputs 06-08 contain the omnibus model tests, group contrasts, and AQ simple slopes.
