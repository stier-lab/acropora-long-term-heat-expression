# F. Model reproducibility check

Generated: 2026-06-29 13:33:47 

| Status | Count |
|---|---|
| PASS | 5 |
| DRIFT | 2 |
| HANDLED | 4 |

## Per-model verdicts

- **12_pam_lmm**: PASS (max coef drift = 0, logLik diff = 0)
- **12_color_lmm**: PASS (max coef drift = 0, logLik diff = 0)
- **12_bw_lm**: PASS (max coef drift = 0, logLik diff = 0)
- **12_zoox_lmm**: PASS (max coef drift = 0, logLik diff = 0)
- **12c_morph_polyps_out_blme**: DRIFT (max coef drift = 1.2, logLik diff = 2.31)
- **12c_morph_axial_polyp_formation_blme**: HANDLED (max coef drift = NA, logLik diff = NA)
- **12c_morph_wound_smoothed_blme**: HANDLED (max coef drift = NA, logLik diff = NA)
- **12c_morph_pigment_over_wound_blme**: DRIFT (max coef drift = 599, logLik diff = 13)
- **12c_morph_tip_exist_blme**: HANDLED (max coef drift = NA, logLik diff = NA)
- **12c_morph_tip_extension_blme**: HANDLED (max coef drift = NA, logLik diff = NA)
- **12c_morph_new_corallites_on_tip_blme**: PASS (max coef drift = 0, logLik diff = 0)
