# PCA + Genet LRT Diagnostics
_Run: 2026-06-25 09:05:14_

## PCA (15_multivariate.R)

- n = 48 corals on 4 vars (pam_end, color_end, growth_pct, zoox_end)
- Variance explained: PC1=82.4%, PC2=13.2%, PC3=2.8%, PC4=1.6%
- PC1+PC2 = **95.5%** (PASS)
- Kaiser: 1 PCs with eigenvalue > 1 (eigenvals 3.29, 0.53, 0.11, 0.07)
- Center = TRUE, Scale = TRUE (PASS)

### Loadings

- **PC1** (82.4%): color_end=+0.53, zoox_end=+0.53, pam_end=+0.52, growth_pct=+0.41
- **PC2** (13.2%): growth_pct=+0.91, pam_end=-0.29, color_end=-0.23, zoox_end=-0.20
- **PC3** (2.8%): pam_end=-0.78, zoox_end=+0.58, color_end=+0.23, growth_pct=-0.06
- **PC4** (1.6%): color_end=-0.78, zoox_end=+0.58, pam_end=+0.21, growth_pct=-0.00

## Genet LRTs (13_genet_interaction.R)

### pam_fvfm
- n_obs = 336; fixed params: null=10, full=24 (Δfx=14). Total npar Δ=14 (reported df=14).
- REML status: null=FALSE, full=FALSE, saved=FALSE (all should be FALSE)
- Convergence: null singular=FALSE, full singular=FALSE, warnings=0
- LRT χ²(14) = 79.03, p = 4.28e-11 (recomputed 4.28e-11)
### color_dscale
- n_obs = 336; fixed params: null=10, full=24 (Δfx=14). Total npar Δ=14 (reported df=14).
- REML status: null=FALSE, full=FALSE, saved=FALSE (all should be FALSE)
- Convergence: null singular=FALSE, full singular=FALSE, warnings=0
- LRT χ²(14) = 169.31, p = 9.43e-29 (recomputed 9.43e-29)
### log_zoox
- n_obs = 192; fixed params: null=10, full=24 (Δfx=14). Total npar Δ=14 (reported df=14).
- REML status: null=FALSE, full=FALSE, saved=FALSE (all should be FALSE)
- Convergence: null singular=FALSE, full singular=FALSE, warnings=0
- LRT χ²(14) = 64.30, p = 2.04e-08 (recomputed 2.04e-08)

### growth_areal
- ML LMM LRT (tank random intercept); χ²(6) = 6.50, p = 0.37
- Fixed params: null=6, full=12 (Δ=6, reported df=6)

## Effect direction vs. emmeans

- **PAM Fv/Fm**: Δ(31C-28C) by genet = a=-0.148, c=-0.037, d=-0.131 — additive (parallel)
- **Color (D-scale)**: Δ(31C-28C) by genet = a=-1.875, c=-0.312, d=-1.562 — additive (parallel)
- **Calcification (mg cm⁻² d⁻¹)**: Δ(31C-28C) by genet = a=-3.685, c=-1.381, d=-3.549 — additive (parallel)
- **log10 symbionts cm⁻² (raw summary)**: Δ(31C-28C) by genet = a=-0.583, c=-0.214, d=-0.517 — additive (parallel)

## Summary verdicts

- PCA: **PASS**
- pam_fvfm LRT: **PASS**
- color_dscale LRT: **PASS**
- log_zoox LRT: **PASS**
- growth_areal LMM LRT: **PASS**
