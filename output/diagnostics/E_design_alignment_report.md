# E. Design alignment audit

Generated: 2026-06-26 12:19:54 

| Status | Count |
|---|---|

   FAIL HANDLED    INFO    PASS 
      4       2       1      25 

## Status summary

### FAIL (4)

- **12c_morph_pigment_over_wound_blme.rds / fixed structure**: Restricted to wounded corals; wound is a stratification, not a covariate
- **12c_morph_polyps_out_blme.rds / fixed structure**: Restricted to wounded corals; wound is a stratification, not a covariate
- **12c_morph_tip_exist_blme.rds / fixed structure**: Restricted to wounded corals; wound is a stratification, not a covariate
- **12c_morph_tip_extension_blme.rds / fixed structure**: Restricted to wounded corals; wound is a stratification, not a covariate

### HANDLED (2)

- **12_pam_lmm / random slope on day?**: Random slope considered and intentionally omitted; expected singular/overfit with sparse per-id trajectories
- **12_color_lmm / ordinal data on Gaussian**: Addressed by 12b CLMM robustness check; KS-violation noted in Section 10

### INFO (1)

- **12_pam_lmm / n observations**: n_corals=48, days=-1,0,3,6,9,12,14

### PASS (25)

- **12_pam_lmm / fixed structure**: 4-way fixed structure required for genet x heat x wound x time
- **12_pam_lmm / random structure**: tank for tank-level confounds, id for repeated measures on same coral
- **12_pam_lmm / balance (tank x treatment)**: 8 unique tank x treatment cells expected (4 tanks per temp, fully nested)
- **12_pam_lmm / (1|tank) sufficient?**: tanks are uniquely IDed; no need for tank:treatment crossed term
- **12_color_lmm / fixed structure**: same 4-way as PAM
- **12_bw_lm / fixed structure**: no day term; tank retained as treatment-assignment block
- **12_bw_lm / random effects**: coral ID omitted because each coral has one endpoint growth observation
- **12_bw_lm / n**: Cells in design: 12 (2 trt x 2 wound x 3 genet = 12)
- **12_bw_lm / df residual**: Adequate residual df for 3-way fixed structure
- **12_zoox_lmm / fixed structure**: biopsy_day_c is centered at Day 1; 4-way as for PAM
- **12_zoox_lmm / drop (1|id) — destructive**: Correct — destructive sampling means each id has 1 obs
- **12_zoox_lmm / biopsy_day_c centering**: Centered correctly at Day 1 baseline
- **12c_morph_hole_in_center_blme.rds / fixed structure**: Restricted to wounded corals; wound is a stratification, not a covariate
- **12c_morph_hole_in_center_blme.rds / Cauchy(0,2.5) prior**: Prior addresses the morphology separation issue
- **12c_morph_new_corallites_on_tip_blme.rds / fixed structure**: Restricted to wounded corals; wound is a stratification, not a covariate
- **12c_morph_new_corallites_on_tip_blme.rds / Cauchy(0,2.5) prior**: Prior addresses the morphology separation issue
- **12c_morph_pigment_over_wound_blme.rds / Cauchy(0,2.5) prior**: Prior addresses the morphology separation issue
- **12c_morph_polyps_out_blme.rds / Cauchy(0,2.5) prior**: Prior addresses the morphology separation issue
- **12c_morph_tip_exist_blme.rds / Cauchy(0,2.5) prior**: Prior addresses the morphology separation issue
- **12c_morph_tip_extension_blme.rds / Cauchy(0,2.5) prior**: Prior addresses the morphology separation issue
- **12c_morph_wound_smoothed_blme.rds / fixed structure**: Restricted to wounded corals; wound is a stratification, not a covariate
- **12c_morph_wound_smoothed_blme.rds / Cauchy(0,2.5) prior**: Prior addresses the morphology separation issue
- **14_cox_* / stratification**: Per-genet HRs in separate per-genet models (low EPV caveat in Section 10)
- **14_cox_* / wounded-only at risk**: Confirmed in script 14
- **15_pca / centering + scaling**: Confirmed in script 15

