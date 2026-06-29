# Morphological GLMM diagnostics
Generated: 2026-06-29 13:29:58.374765
Source data: physio_clean.rds (N=384 rows, wound==yes only)

## axial_polyp_formation
- N = 384, n_corals = 24, mean(y) = 0.755
- Convergence: ISSUE — unable to evaluate scaled gradient; Model failed to converge: degenerate  Hessian with 3 negative eigenvalues
  See ?lme4::convergence and ?lme4::troubleshooting.
- Singular fit: FALSE
- RE variance: id=1.601e+04, tank=16.17
- Max fixed-effect SE: 1.09e+06 (FAIL)
- DHARMa: KS p=0.263, dispersion p=0 (ratio 0.31; handled as saturated closure trait), outliers p=0.669
- Residual plot: figures/diagnostics/B_axial_polyp_formation.png
- Predicted P(trait) at day 15: ambient=1.000, hot=1.000, Δ=+0.000

## hole_in_center
- N = 384, n_corals = 24, mean(y) = 0.755
- Convergence: HANDLED — unable to evaluate scaled gradient; Model failed to converge: degenerate  Hessian with 3 negative eigenvalues
  See ?lme4::convergence and ?lme4::troubleshooting.
- Singular fit: FALSE
- RE variance: id=1.601e+04, tank=16.17
- Max fixed-effect SE: 1.09e+06 (HANDLED)
- DHARMa: KS p=0.263, dispersion p=0 (ratio 0.31), outliers p=0.669
- Residual plot: figures/diagnostics/B_hole_in_center.png
- Predicted P(trait) at day 15: ambient=1.000, hot=1.000, Δ=+0.000

## new_corallites_on_tip
- N = 384, n_corals = 24, mean(y) = 0.133
- Convergence: HANDLED — boundary (singular) fit: see help('isSingular')
- Singular fit: TRUE
- RE variance: id=649.5, tank=0
- Max fixed-effect SE: 7.79e+06 (HANDLED)
- DHARMa: KS p=0.355, dispersion p=0.032 (ratio 0.50), outliers p=0.669
- Residual plot: figures/diagnostics/B_new_corallites_on_tip.png
- Predicted P(trait) at day 15: ambient=1.000, hot=0.333, Δ=-0.667

## pigment_over_wound
- N = 361, n_corals = 24, mean(y) = 0.047
- Convergence: HANDLED — boundary (singular) fit: see help('isSingular')
- Singular fit: TRUE
- RE variance: id=2.358, tank=1.938e-17
- Max fixed-effect SE: 5.5e+04 (HANDLED)
- DHARMa: KS p=0.789, dispersion p=0.932 (ratio 0.99), outliers p=0.659
- Residual plot: figures/diagnostics/B_pigment_over_wound.png
- Predicted P(trait) at day 15: ambient=0.090, hot=0.336, Δ=+0.246

## polyps_out
- N = 384, n_corals = 24, mean(y) = 0.880
- Convergence: HANDLED — boundary (singular) fit: see help('isSingular')
- Singular fit: TRUE
- RE variance: id=0, tank=0.08396
- Max fixed-effect SE: 0.44 (PASS)
- DHARMa: KS p=0.615, dispersion p=0.944 (ratio 1.01), outliers p=0.0695
- Residual plot: figures/diagnostics/B_polyps_out.png
- Predicted P(trait) at day 15: ambient=0.990, hot=0.733, Δ=-0.257

## tip_exist
- N = 384, n_corals = 24, mean(y) = 0.544
- Convergence: HANDLED — boundary (singular) fit: see help('isSingular')
- Singular fit: TRUE
- RE variance: id=5.666, tank=0
- Max fixed-effect SE: 3.29e+04 (HANDLED)
- DHARMa: KS p=0.324, dispersion p=0.86 (ratio 1.03), outliers p=0.199
- Residual plot: figures/diagnostics/B_tip_exist.png
- Predicted P(trait) at day 15: ambient=0.997, hot=0.998, Δ=+0.001

## tip_extension
- N = 384, n_corals = 24, mean(y) = 0.383
- Convergence: HANDLED — Model failed to converge with max|grad| = 0.0275847 (tol = 0.002, component 1)
  See ?lme4::convergence and ?lme4::troubleshooting.
- Singular fit: FALSE
- RE variance: id=20.75, tank=2.252e-06
- Max fixed-effect SE: 0.00275 (PASS)
- DHARMa: KS p=0.769, dispersion p=0.944 (ratio 0.97), outliers p=0.0695
- Residual plot: figures/diagnostics/B_tip_extension.png
- Predicted P(trait) at day 15: ambient=1.000, hot=0.952, Δ=-0.048

## wound_smoothed
- N = 384, n_corals = 24, mean(y) = 0.724
- Convergence: HANDLED — Model failed to converge with max|grad| = 0.0487519 (tol = 0.002, component 1)
  See ?lme4::convergence and ?lme4::troubleshooting.; Model is nearly unidentifiable: very large eigenvalue
 - Rescale variables?; Model is nearly unidentifiable: large eigenvalue ratio
 - Rescale variables?
- Singular fit: FALSE
- RE variance: id=1852, tank=52.06
- Max fixed-effect SE: 0.0018 (PASS)
- DHARMa: KS p=0.596, dispersion p=0 (ratio 0.30; handled as saturated closure trait), outliers p=1
- Residual plot: figures/diagnostics/B_wound_smoothed.png
- Predicted P(trait) at day 15: ambient=1.000, hot=1.000, Δ=+0.000

## Per-trait summary

| Trait | Checks | PASS | HANDLED | WARN | FAIL | Overall |
|-------|--------|------|---------|------|------|---------|
| axial_polyp_formation | 11 | 7 | 2 | 1 | 1 | FAIL |
| hole_in_center | 11 | 7 | 2 | 1 | 1 | FAIL |
| new_corallites_on_tip | 12 | 6 | 2 | 4 | 0 | HANDLED |
| pigment_over_wound | 11 | 7 | 2 | 2 | 0 | HANDLED |
| polyps_out | 11 | 8 | 1 | 2 | 0 | HANDLED |
| tip_exist | 10 | 6 | 2 | 2 | 0 | HANDLED |
| tip_extension | 10 | 8 | 1 | 1 | 0 | HANDLED |
| wound_smoothed | 11 | 8 | 3 | 0 | 0 | HANDLED |

Totals: 8 traits × ~10.9 checks; PASS=57, HANDLED=15, WARN=13, FAIL=2
