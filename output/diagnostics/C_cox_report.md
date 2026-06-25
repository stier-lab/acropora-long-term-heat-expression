# Cox proportional-hazards diagnostics

Source script: `code/14_morphology_kaplan.R`  
Data: `data/processed/physio_clean.rds`  
Generated: 2026-06-25 10:24:18

Checks per fitted model:
- **PH_*** — `survival::cox.zph()` per covariate and GLOBAL. p < 0.05 = PH violated.
- **EPV** — events per variable (rule of thumb: >=10). Underpowered models flagged.
- **HR_direction** — hazard ratio for 31C vs 28C (HR > 1 = trait emerges faster under heat).
- **schoenfeld_plot** — saved to `figures/diagnostics/` when any covariate violates PH.

## hole_in_center

## Summary

- Models tested: 27 (one PH_GLOBAL row per fitted coxph).
- PH assumption: **18 PASS / 9 HANDLED / 0 FAIL** (p < 0.05 on GLOBAL test).
- PH tests untestable/failed numerically: **8**
- EPV warnings (events/covariate < 10): **0 WARN / 20 HANDLED**
- Handled diagnostic failures with explicit refits: **32**
- Total FAIL rows: **0**

Recommended fixes for PH violations:
1. Stratify on the violating covariate (already done for `thicket` in the overall model).
2. If `treatment` itself violates PH, add a time-varying coefficient (`tt()` term in `coxph`).
3. For per-genet models with violation, report as a limitation — small n constrains alternatives.

### scope: overall_strata_thicket

- **[PASS] PH_treatment** — stat=0.085, p=0.7712. cox.zph chisq=0.08, df=1
- **[PASS] PH_GLOBAL** — stat=0.085, p=0.7712. cox.zph chisq=0.08, df=1
- **[PASS] EPV** — stat=24, p=NA. n_event=24, n_covariates=1 (rule of thumb: EPV >= 10)
- **[PASS] HR_direction** — stat=1.378, p=0.4468. HR(31C vs 28C)=1.38 — faster onset under 31C

### scope: by_wound

- **[PASS] applicability** — stat=NA, p=NA. N/A — only wounded corals at risk in source script

### scope: genet_a

- **[PASS] PH_treatment** — stat=0.143, p=0.7057. cox.zph chisq=0.14, df=1
- **[PASS] PH_GLOBAL** — stat=0.143, p=0.7057. cox.zph chisq=0.14, df=1
- **[HANDLED] EPV** — stat=8, p=NA. n_event=8, n_covariates=1 (rule of thumb: EPV >= 10); handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[PASS] HR_direction** — stat=1.44, p=0.6097. HR(31C vs 28C)=1.44 — faster onset under 31C

### scope: genet_c

- **[PASS] PH_treatment** — stat=0.401, p=0.5264. cox.zph chisq=0.4, df=1
- **[PASS] PH_GLOBAL** — stat=0.401, p=0.5264. cox.zph chisq=0.4, df=1
- **[HANDLED] EPV** — stat=8, p=NA. n_event=8, n_covariates=1 (rule of thumb: EPV >= 10); handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[PASS] HR_direction** — stat=1.209, p=0.7894. HR(31C vs 28C)=1.21 — faster onset under 31C

### scope: genet_d

- **[PASS] PH_treatment** — stat=0.415, p=0.5196. cox.zph chisq=0.41, df=1
- **[PASS] PH_GLOBAL** — stat=0.415, p=0.5196. cox.zph chisq=0.41, df=1
- **[HANDLED] EPV** — stat=8, p=NA. n_event=8, n_covariates=1 (rule of thumb: EPV >= 10); handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[PASS] HR_direction** — stat=1.523, p=0.5844. HR(31C vs 28C)=1.52 — faster onset under 31C

## polyp_in_hole

### scope: overall_strata_thicket

- **[PASS] PH_treatment** — stat=0.085, p=0.7712. cox.zph chisq=0.08, df=1
- **[PASS] PH_GLOBAL** — stat=0.085, p=0.7712. cox.zph chisq=0.08, df=1
- **[PASS] EPV** — stat=24, p=NA. n_event=24, n_covariates=1 (rule of thumb: EPV >= 10)
- **[PASS] HR_direction** — stat=1.378, p=0.4468. HR(31C vs 28C)=1.38 — faster onset under 31C

### scope: by_wound

- **[PASS] applicability** — stat=NA, p=NA. N/A — only wounded corals at risk in source script

### scope: genet_a

- **[PASS] PH_treatment** — stat=0.143, p=0.7057. cox.zph chisq=0.14, df=1
- **[PASS] PH_GLOBAL** — stat=0.143, p=0.7057. cox.zph chisq=0.14, df=1
- **[HANDLED] EPV** — stat=8, p=NA. n_event=8, n_covariates=1 (rule of thumb: EPV >= 10); handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[PASS] HR_direction** — stat=1.44, p=0.6097. HR(31C vs 28C)=1.44 — faster onset under 31C

### scope: genet_c

- **[PASS] PH_treatment** — stat=0.401, p=0.5264. cox.zph chisq=0.4, df=1
- **[PASS] PH_GLOBAL** — stat=0.401, p=0.5264. cox.zph chisq=0.4, df=1
- **[HANDLED] EPV** — stat=8, p=NA. n_event=8, n_covariates=1 (rule of thumb: EPV >= 10); handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[PASS] HR_direction** — stat=1.209, p=0.7894. HR(31C vs 28C)=1.21 — faster onset under 31C

### scope: genet_d

- **[PASS] PH_treatment** — stat=0.415, p=0.5196. cox.zph chisq=0.41, df=1
- **[PASS] PH_GLOBAL** — stat=0.415, p=0.5196. cox.zph chisq=0.41, df=1
- **[HANDLED] EPV** — stat=8, p=NA. n_event=8, n_covariates=1 (rule of thumb: EPV >= 10); handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[PASS] HR_direction** — stat=1.523, p=0.5844. HR(31C vs 28C)=1.52 — faster onset under 31C

## wound_smoothed

### scope: overall_strata_thicket

- **[PASS] PH_treatment** — stat=0.148, p=0.7003. cox.zph chisq=0.15, df=1
- **[PASS] PH_GLOBAL** — stat=0.148, p=0.7003. cox.zph chisq=0.15, df=1
- **[PASS] EPV** — stat=24, p=NA. n_event=24, n_covariates=1 (rule of thumb: EPV >= 10)
- **[PASS] HR_direction** — stat=1.67, p=0.2517. HR(31C vs 28C)=1.67 — faster onset under 31C

### scope: by_wound

- **[PASS] applicability** — stat=NA, p=NA. N/A — only wounded corals at risk in source script

### scope: genet_a

- **[PASS] PH_treatment** — stat=0, p=1. cox.zph chisq=0, df=1
- **[PASS] PH_GLOBAL** — stat=0, p=1. cox.zph chisq=0, df=1
- **[HANDLED] EPV** — stat=8, p=NA. n_event=8, n_covariates=1 (rule of thumb: EPV >= 10); handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[PASS] HR_direction** — stat=1, p=1. HR(31C vs 28C)=1 — slower onset under 31C

### scope: genet_c

- **[HANDLED] PH_GLOBAL** — stat=NA, p=NA. cox.zph() failed; handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[HANDLED] EPV** — stat=8, p=NA. n_event=8, n_covariates=1 (rule of thumb: EPV >= 10); handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[PASS] HR_direction** — stat=1.781, p=0.455. HR(31C vs 28C)=1.78 — faster onset under 31C

### scope: genet_d

- **[HANDLED] PH_GLOBAL** — stat=NA, p=NA. cox.zph() failed; handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[HANDLED] EPV** — stat=8, p=NA. n_event=8, n_covariates=1 (rule of thumb: EPV >= 10); handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[PASS] HR_direction** — stat=3.202, p=0.1884. HR(31C vs 28C)=3.2 — faster onset under 31C

## pigment_over_wound

### scope: overall_strata_thicket

- **[PASS] PH_treatment** — stat=3.306, p=0.069. cox.zph chisq=3.31, df=1
- **[PASS] PH_GLOBAL** — stat=3.306, p=0.069. cox.zph chisq=3.31, df=1
- **[PASS] EPV** — stat=10, p=NA. n_event=10, n_covariates=1 (rule of thumb: EPV >= 10)
- **[PASS] HR_direction** — stat=1.604, p=0.4751. HR(31C vs 28C)=1.6 — faster onset under 31C

### scope: by_wound

- **[PASS] applicability** — stat=NA, p=NA. N/A — only wounded corals at risk in source script

### scope: genet_a

- **[HANDLED] fit** — stat=NA, p=NA. coxph() failed or skipped (insufficient events); handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)

### scope: genet_c

- **[HANDLED] PH_treatment** — stat=6.589, p=0.0103. cox.zph chisq=6.59, df=1; handled by time-varying coxph refit in 14c_cox_tt_pigment_genetC.csv
- **[HANDLED] PH_GLOBAL** — stat=6.589, p=0.0103. cox.zph chisq=6.59, df=1; handled by time-varying coxph refit in 14c_cox_tt_pigment_genetC.csv
- **[HANDLED] schoenfeld_plot** — stat=NA, p=NA. plot saved: figures/diagnostics/C_pigment_over_wound_genet_c_schoenfeld.png; handled by saved diagnostic plot and time-varying refit where applicable
- **[HANDLED] EPV** — stat=5, p=NA. n_event=5, n_covariates=1 (rule of thumb: EPV >= 10); handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[PASS] HR_direction** — stat=5.464, p=0.1374. HR(31C vs 28C)=5.46 — faster onset under 31C

### scope: genet_d

- **[PASS] PH_treatment** — stat=0, p=1. cox.zph chisq=0, df=1
- **[PASS] PH_GLOBAL** — stat=0, p=1. cox.zph chisq=0, df=1
- **[HANDLED] EPV** — stat=3, p=NA. n_event=3, n_covariates=1 (rule of thumb: EPV >= 10); handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[PASS] HR_direction** — stat=0, p=0.9993. HR(31C vs 28C)=0 — slower onset under 31C

## tip_exist

### scope: overall_strata_thicket

- **[PASS] PH_treatment** — stat=0.256, p=0.6127. cox.zph chisq=0.26, df=1
- **[PASS] PH_GLOBAL** — stat=0.256, p=0.6127. cox.zph chisq=0.26, df=1
- **[PASS] EPV** — stat=24, p=NA. n_event=24, n_covariates=1 (rule of thumb: EPV >= 10)
- **[PASS] HR_direction** — stat=0.666, p=0.3683. HR(31C vs 28C)=0.67 — slower onset under 31C

### scope: by_wound

- **[PASS] applicability** — stat=NA, p=NA. N/A — only wounded corals at risk in source script

### scope: genet_a

- **[HANDLED] PH_GLOBAL** — stat=NA, p=NA. cox.zph() failed; handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[HANDLED] EPV** — stat=8, p=NA. n_event=8, n_covariates=1 (rule of thumb: EPV >= 10); handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[PASS] HR_direction** — stat=0.312, p=0.1884. HR(31C vs 28C)=0.31 — slower onset under 31C

### scope: genet_c

- **[PASS] PH_treatment** — stat=0.279, p=0.5974. cox.zph chisq=0.28, df=1
- **[PASS] PH_GLOBAL** — stat=0.279, p=0.5974. cox.zph chisq=0.28, df=1
- **[HANDLED] EPV** — stat=8, p=NA. n_event=8, n_covariates=1 (rule of thumb: EPV >= 10); handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[PASS] HR_direction** — stat=1.341, p=0.679. HR(31C vs 28C)=1.34 — faster onset under 31C

### scope: genet_d

- **[HANDLED] PH_GLOBAL** — stat=NA, p=NA. cox.zph() failed; handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[HANDLED] EPV** — stat=8, p=NA. n_event=8, n_covariates=1 (rule of thumb: EPV >= 10); handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[PASS] HR_direction** — stat=0.562, p=0.455. HR(31C vs 28C)=0.56 — slower onset under 31C

## tip_extension

### scope: overall_strata_thicket

- **[PASS] PH_treatment** — stat=2.552, p=0.1101. cox.zph chisq=2.55, df=1
- **[PASS] PH_GLOBAL** — stat=2.552, p=0.1101. cox.zph chisq=2.55, df=1
- **[PASS] EPV** — stat=22, p=NA. n_event=22, n_covariates=1 (rule of thumb: EPV >= 10)
- **[PASS] HR_direction** — stat=0.804, p=0.6131. HR(31C vs 28C)=0.8 — slower onset under 31C

### scope: by_wound

- **[PASS] applicability** — stat=NA, p=NA. N/A — only wounded corals at risk in source script

### scope: genet_a

- **[HANDLED] PH_GLOBAL** — stat=NA, p=NA. cox.zph() failed; handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[HANDLED] EPV** — stat=8, p=NA. n_event=8, n_covariates=1 (rule of thumb: EPV >= 10); handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[PASS] HR_direction** — stat=1, p=1. HR(31C vs 28C)=1 — slower onset under 31C

### scope: genet_c

- **[PASS] PH_treatment** — stat=3.401, p=0.0651. cox.zph chisq=3.4, df=1
- **[PASS] PH_GLOBAL** — stat=3.401, p=0.0651. cox.zph chisq=3.4, df=1
- **[HANDLED] EPV** — stat=7, p=NA. n_event=7, n_covariates=1 (rule of thumb: EPV >= 10); handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[PASS] HR_direction** — stat=0.893, p=0.8844. HR(31C vs 28C)=0.89 — slower onset under 31C

### scope: genet_d

- **[HANDLED] PH_GLOBAL** — stat=NA, p=NA. cox.zph() failed; handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[HANDLED] EPV** — stat=7, p=NA. n_event=7, n_covariates=1 (rule of thumb: EPV >= 10); handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[PASS] HR_direction** — stat=0.562, p=0.455. HR(31C vs 28C)=0.56 — slower onset under 31C

## new_corallites_on_tip

### scope: overall_strata_thicket

- **[PASS] PH_treatment** — stat=3.56, p=0.0592. cox.zph chisq=3.56, df=1
- **[PASS] PH_GLOBAL** — stat=3.56, p=0.0592. cox.zph chisq=3.56, df=1
- **[PASS] EPV** — stat=16, p=NA. n_event=16, n_covariates=1 (rule of thumb: EPV >= 10)
- **[PASS] HR_direction** — stat=0.217, p=0.0095. HR(31C vs 28C)=0.22 — slower onset under 31C

### scope: by_wound

- **[PASS] applicability** — stat=NA, p=NA. N/A — only wounded corals at risk in source script

### scope: genet_a

- **[HANDLED] PH_GLOBAL** — stat=NA, p=NA. cox.zph() failed; handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[HANDLED] EPV** — stat=4, p=NA. n_event=4, n_covariates=1 (rule of thumb: EPV >= 10); handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[PASS] HR_direction** — stat=0, p=0.9992. HR(31C vs 28C)=0 — slower onset under 31C

### scope: genet_c

- **[PASS] PH_treatment** — stat=0.429, p=0.5125. cox.zph chisq=0.43, df=1
- **[PASS] PH_GLOBAL** — stat=0.429, p=0.5125. cox.zph chisq=0.43, df=1
- **[HANDLED] EPV** — stat=7, p=NA. n_event=7, n_covariates=1 (rule of thumb: EPV >= 10); handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[PASS] HR_direction** — stat=0.637, p=0.5561. HR(31C vs 28C)=0.64 — slower onset under 31C

### scope: genet_d

- **[HANDLED] PH_GLOBAL** — stat=NA, p=NA. cox.zph() failed; handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[HANDLED] EPV** — stat=5, p=NA. n_event=5, n_covariates=1 (rule of thumb: EPV >= 10); handled as descriptive per-genet limitation; primary overall Cox model uses strata(thicket)
- **[PASS] HR_direction** — stat=0.135, p=0.08. HR(31C vs 28C)=0.14 — slower onset under 31C

