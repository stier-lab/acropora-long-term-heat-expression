# Continuous-response model diagnostics
Generated: 2026-06-14 09:29:38.922942

Models reviewed: 12_pam_lmm, 12_color_lmm, 12_zoox_lmm, 12_bw_lm

## Summary
- Total checks: 32
- PASS: 24
- HANDLED: 8
- WARN: 0
- FAIL: 0

## 12_pam_lmm
- **isSingular** [PASS]: stat=0 Singular fit means a variance component is at/near zero
- **variance_components_near_zero** [PASS]: stat=0 VarCorr (grp:var): id:3.84e-05; tank:1.92e-04; Residual:6.60e-04
- **optimizer_convergence_messages** [PASS]: stat=0 
- **DHARMa_KS_uniformity** [PASS]: stat=0.0575 Kolmogorov-Smirnov on scaled residuals
- **DHARMa_dispersion** [PASS]: stat=0.9074 Ratio of obs/sim residual variance
- **DHARMa_outliers** [PASS]: stat=2 Excess outliers vs expected
- **cooks_distance** [HANDLED]: influence.ME::influence failed; handled with saved residual plots and direction checks
- **emmeans_direction_pam_fvfm** [PASS]: stat=0.08351 contrast 28C - 31C = 0.0835 (p=6.76e-05); observed heated decrease

## 12_color_lmm
- **isSingular** [PASS]: stat=0 Singular fit means a variance component is at/near zero
- **variance_components_near_zero** [PASS]: stat=0 VarCorr (grp:var): id:2.45e-02; tank:2.39e-02; Residual:5.78e-02
- **optimizer_convergence_messages** [PASS]: stat=0 
- **DHARMa_KS_uniformity** [HANDLED]: stat=0.1828 Kolmogorov-Smirnov on scaled residuals; handled by ordinal CLMM robustness model 12b_color_clmm
- **DHARMa_dispersion** [PASS]: stat=0.89 Ratio of obs/sim residual variance
- **DHARMa_outliers** [PASS]: stat=4 Excess outliers vs expected
- **cooks_distance** [HANDLED]: influence.ME::influence failed; handled with saved residual plots and direction checks
- **emmeans_direction_color_dscale** [PASS]: stat=1.153 contrast 28C - 31C = 1.1526 (p=3.33e-05); observed heated decrease

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
- **shapiro_residual_normality** [PASS]: stat=0.9587 Shapiro-Wilk on OLS residuals
- **breusch_pagan_heteroscedasticity** [PASS]: stat=16.5 BP test for non-constant variance
- **VIF** [HANDLED]: Saturated 3-way interaction lm; VIFs uninterpretable. Skipped.; handled by explicit full-factorial design statement; no additive VIF interpretation
- **DHARMa_KS_uniformity** [PASS]: stat=0.1148 DHARMa KS on simulated residuals
- **cooks_distance_max** [HANDLED]: stat=0.2856 Top-3 obs idx: 38,41,2; cd: 0.286,0.150,0.107; handled by top-three Cook's-distance sensitivity check
- **cooks_top3_sensitivity** [PASS]: stat=3.15 dropped top Cook's rows 38,41,2; full est=2.886, refit est=3.150
- **emmeans_direction_growth_pct** [PASS]: stat=2.886 contrast 28C - 31C = 2.8864 (p=2.43e-05); observed heated decrease

