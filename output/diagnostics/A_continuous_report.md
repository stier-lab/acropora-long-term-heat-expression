# Continuous-response model diagnostics
Generated: 2026-06-26 17:36:22.563342

Models reviewed: 12_pam_lmm, 12_color_lmm, 12_zoox_lmm, 12_bw_lm

## Summary
- Total checks: 36
- PASS: 27
- HANDLED: 6
- WARN: 3
- FAIL: 0

## 12_pam_lmm
- **isSingular** [PASS]: stat=0 Singular fit means a variance component is at/near zero
- **variance_components_near_zero** [PASS]: stat=0 VarCorr (grp:var): id:5.78e-05; tank:2.31e-04; Residual:6.31e-04
- **optimizer_convergence_messages** [PASS]: stat=0 
- **DHARMa_KS_uniformity** [PASS]: stat=0.06558 Kolmogorov-Smirnov on scaled residuals
- **DHARMa_dispersion** [PASS]: stat=0.8966 Ratio of obs/sim residual variance
- **DHARMa_outliers** [PASS]: stat=2 Excess outliers vs expected
- **cooks_distance_max** [WARN]: stat=1.077 Top-3 obs idx: 325,317,320; cd: 1.077,0.804,0.758
- **emmeans_direction_pam_fvfm** [PASS]: stat=0.08747 contrast 28C - 31C = 0.0875 (p=0.000112); observed heated decrease

## 12_color_lmm
- **isSingular** [PASS]: stat=0 Singular fit means a variance component is at/near zero
- **variance_components_near_zero** [PASS]: stat=0 VarCorr (grp:var): id:3.53e-02; tank:3.34e-02; Residual:4.86e-02
- **optimizer_convergence_messages** [PASS]: stat=0 
- **DHARMa_KS_uniformity** [HANDLED]: stat=0.2045 Kolmogorov-Smirnov on scaled residuals; handled by ordinal CLMM robustness model 12b_color_clmm
- **DHARMa_dispersion** [PASS]: stat=0.8499 Ratio of obs/sim residual variance
- **DHARMa_outliers** [PASS]: stat=3 Excess outliers vs expected
- **cooks_distance_max** [WARN]: stat=2.764 Top-3 obs idx: 296,274,324; cd: 2.764,2.761,2.739
- **emmeans_direction_color_dscale** [PASS]: stat=1.212 contrast 28C - 31C = 1.2122 (p=9.06e-05); observed heated decrease

## 12_zoox_lmm
- **isSingular** [PASS]: stat=0 Singular fit means a variance component is at/near zero
- **variance_components_near_zero** [PASS]: stat=0 VarCorr (grp:var): tank:3.31e-02; Residual:2.05e-01
- **optimizer_convergence_messages** [PASS]: stat=0 
- **DHARMa_KS_uniformity** [HANDLED]: stat=0.12 Kolmogorov-Smirnov on scaled residuals; handled by explicit top-four residual sensitivity check
- **DHARMa_dispersion** [PASS]: stat=0.8714 Ratio of obs/sim residual variance
- **DHARMa_outliers** [HANDLED]: stat=4 Excess outliers vs expected; handled by explicit top-four residual sensitivity check
- **cooks_distance** [HANDLED]: influence.ME::influence failed; handled with saved residual plots and direction checks
- **emmeans_direction_log_zoox** [PASS]: stat=1.888 contrast 28C - 31C = 1.8882 (p=1.56e-07); observed heated decrease
- **outlier_sensitivity_top4** [PASS]: stat=1.859 dropped largest |scaled residual| rows 187,130,186,189; full est=1.888, refit est=1.859

## 12_bw_lm
- **isSingular** [PASS]: stat=0 Singular fit means a variance component is at/near zero
- **variance_components_near_zero** [PASS]: stat=0 VarCorr (grp:var): tank:1.35e+00; Residual:2.88e+00
- **optimizer_convergence_messages** [PASS]: stat=0 
- **DHARMa_KS_uniformity** [PASS]: stat=0.1135 Kolmogorov-Smirnov on scaled residuals
- **DHARMa_dispersion** [PASS]: stat=0.7863 Ratio of obs/sim residual variance
- **DHARMa_outliers** [PASS]: stat=0 Excess outliers vs expected
- **shapiro_residual_normality** [WARN]: stat=0.9425 Shapiro-Wilk on LMM conditional residuals
- **VIF** [HANDLED]: Saturated 3-way fixed structure with tank random intercept; VIFs uninterpretable. Skipped.; handled by explicit full-factorial design statement and tank random intercept; no additive VIF interpretation
- **cooks_distance_max** [HANDLED]: stat=0.3013 Top-3 obs idx: 38,41,45; cd: 0.301,0.118,0.110; handled by top-three Cook's-distance sensitivity check
- **cooks_top3_sensitivity** [PASS]: stat=2.806 dropped top Cook's rows 38,41,45; full est=2.894, refit est=2.806
- **emmeans_direction_growth_pct** [PASS]: stat=2.894 contrast 28C - 31C = 2.8943 (p=0.0233); observed heated decrease

