AUCi deep sensitivity for eye fixation-duration association

Base model: percentage_fix_duration_beta ~ Group * AUCi_z + medication + session + (1 | ID).
Additional checks:
leave-one-session-out GLMMs
leave-one-subject-out GLMMs
Cook's distance from comparable logit-linear model
linear model after removing Cook's D > 4/n
robust linear model using MASS::rlm
rank-transformed AUCi model
1000 permutations of Group labels for the Group x AUCi interaction
