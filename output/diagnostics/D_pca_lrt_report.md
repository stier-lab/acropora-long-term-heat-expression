# PCA + Genet LRT Diagnostics
_Run: 2026-06-13 07:04:05_

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
- n_obs = 336; fixed params: null=8, full=24 (Δfx=16). Total npar Δ=15 (reported df=15).
- REML status: null=FALSE, full=FALSE, saved=FALSE (all should be FALSE)
- Convergence: null singular=FALSE, full singular=FALSE, warnings=0
- LRT χ²(15) = 90.49, p = 8.02e-13 (recomputed 8.02e-13)
### color_dscale
- n_obs = 336; fixed params: null=8, full=24 (Δfx=16). Total npar Δ=15 (reported df=15).
- REML status: null=FALSE, full=FALSE, saved=FALSE (all should be FALSE)
- Convergence: null singular=FALSE, full singular=FALSE, warnings=0
- LRT χ²(15) = 177.72, p = 6.85e-30 (recomputed 6.85e-30)
### log_zoox
- n_obs = 192; fixed params: null=8, full=24 (Δfx=16). Total npar Δ=15 (reported df=15).
- REML status: null=FALSE, full=FALSE, saved=FALSE (all should be FALSE)
- Convergence: null singular=FALSE, full singular=FALSE, warnings=0
- LRT χ²(15) = 73.57, p = 1.03e-09 (recomputed 1.03e-09)

### growth_areal
- lm-based F-test (no time dim); F(6, 36) = 0.62, p = 0.716
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
- growth_areal F-test: **PASS**
