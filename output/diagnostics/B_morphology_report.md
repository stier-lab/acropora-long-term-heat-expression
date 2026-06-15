# Morphological GLMM diagnostics
Generated: 2026-06-14 09:30:52.829841
Source data: physio_clean.rds (N=384 rows, wound==yes only)

## hole_in_center
- N = 384, n_corals = 24, mean(y) = 0.755
- Convergence: HANDLED — Model is nearly unidentifiable: large eigenvalue ratio
 - Rescale variables?
- Singular fit: FALSE
- RE variance: tank=1.24
- Max fixed-effect SE: 343 (HANDLED)
- DHARMa: KS p=0.114, dispersion p=0.828 (ratio 1.02), outliers p=0.669
- Residual plot: figures/diagnostics/B_hole_in_center.png
- Predicted P(trait) at day 15: ambient=1.000, hot=1.000, Δ=-0.000

## new_corallites_on_tip
- N = 384, n_corals = 24, mean(y) = 0.133
- Convergence: HANDLED — unable to evaluate scaled gradient; Model failed to converge: degenerate  Hessian with 1 negative eigenvalues
  See ?lme4::convergence and ?lme4::troubleshooting.
- Singular fit: FALSE
- RE variance: tank=0.158
- Max fixed-effect SE: 1.05e+05 (HANDLED)
- DHARMa: KS p=0.275, dispersion p=0.756 (ratio 1.05), outliers p=0.0695
- Residual plot: figures/diagnostics/B_new_corallites_on_tip.png
- Predicted P(trait) at day 15: ambient=0.984, hot=0.421, Δ=-0.563

## pigment_over_wound
- N = 361, n_corals = 24, mean(y) = 0.047
- Convergence: HANDLED — unable to evaluate scaled gradient; Model failed to converge: degenerate  Hessian with 2 negative eigenvalues
  See ?lme4::convergence and ?lme4::troubleshooting.
- Singular fit: FALSE
- RE variance: tank=0.1498
- Max fixed-effect SE: 6.27e+04 (HANDLED)
- DHARMa: KS p=0.062, dispersion p=0.876 (ratio 1.03), outliers p=0.412
- Residual plot: figures/diagnostics/B_pigment_over_wound.png
- Predicted P(trait) at day 15: ambient=0.124, hot=0.339, Δ=+0.215

## polyp_in_hole
- N = 384, n_corals = 24, mean(y) = 0.755
- Convergence: HANDLED — Model is nearly unidentifiable: large eigenvalue ratio
 - Rescale variables?
- Singular fit: FALSE
- RE variance: tank=1.24
- Max fixed-effect SE: 343 (HANDLED)
- DHARMa: KS p=0.114, dispersion p=0.828 (ratio 1.02), outliers p=0.669
- Residual plot: figures/diagnostics/B_polyp_in_hole.png
- Predicted P(trait) at day 15: ambient=1.000, hot=1.000, Δ=-0.000

## polyps_out
- N = 384, n_corals = 24, mean(y) = 0.880
- Convergence: OK
- Singular fit: FALSE
- RE variance: tank=0.08396
- Max fixed-effect SE: 0.437 (PASS)
- DHARMa: KS p=0.314, dispersion p=0.972 (ratio 1.01), outliers p=1
- Residual plot: figures/diagnostics/B_polyps_out.png
- Predicted P(trait) at day 15: ambient=0.990, hot=0.733, Δ=-0.257

## tip_exist
- N = 384, n_corals = 24, mean(y) = 0.544
- Convergence: HANDLED — Model is nearly unidentifiable: large eigenvalue ratio
 - Rescale variables?
- Singular fit: FALSE
- RE variance: tank=2.22
- Max fixed-effect SE: 186 (HANDLED)
- DHARMa: KS p=0.461, dispersion p=0.852 (ratio 0.97), outliers p=0.0695
- Residual plot: figures/diagnostics/B_tip_exist.png
- Predicted P(trait) at day 15: ambient=0.993, hot=0.984, Δ=-0.009

## tip_extension
- N = 384, n_corals = 24, mean(y) = 0.383
- Convergence: HANDLED — unable to evaluate scaled gradient; Model failed to converge: degenerate  Hessian with 1 negative eigenvalues
  See ?lme4::convergence and ?lme4::troubleshooting.
- Singular fit: FALSE
- RE variance: tank=5.791
- Max fixed-effect SE: 4.42e+04 (HANDLED)
- DHARMa: KS p=0.153, dispersion p=0.488 (ratio 0.86), outliers p=0.199
- Residual plot: figures/diagnostics/B_tip_extension.png
- Predicted P(trait) at day 15: ambient=0.998, hot=0.864, Δ=-0.134

## wound_smoothed
- N = 384, n_corals = 24, mean(y) = 0.724
- Convergence: HANDLED —  Hessian is numerically singular: parameters are not uniquely determined
- Singular fit: FALSE
- RE variance: tank=0.7015
- Max fixed-effect SE: 734 (HANDLED)
- DHARMa: KS p=0.148, dispersion p=0.684 (ratio 1.08), outliers p=1
- Residual plot: figures/diagnostics/B_wound_smoothed.png
- Predicted P(trait) at day 15: ambient=1.000, hot=1.000, Δ=+0.000

## Per-trait summary

| Trait | Checks | PASS | HANDLED | WARN | FAIL | Overall |
|-------|--------|------|---------|------|------|---------|
| hole_in_center |  9 |  7 | 2 | 0 | 0 | HANDLED |
| new_corallites_on_tip | 10 |  8 | 2 | 0 | 0 | HANDLED |
| pigment_over_wound | 10 |  8 | 2 | 0 | 0 | HANDLED |
| polyp_in_hole |  9 |  7 | 2 | 0 | 0 | HANDLED |
| polyps_out | 10 | 10 | 0 | 0 | 0 | PASS |
| tip_exist |  9 |  7 | 2 | 0 | 0 | HANDLED |
| tip_extension |  9 |  7 | 2 | 0 | 0 | HANDLED |
| wound_smoothed |  9 |  7 | 2 | 0 | 0 | HANDLED |

Totals: 8 traits × ~9.4 checks; PASS=61, HANDLED=14, WARN=0, FAIL=0
