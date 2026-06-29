# LTH Results Summary

> **Authoritative phenotype results narrative** · Updated 2026-06-12 · Index: [`README.md`](README.md) · every number traces to `output/tables/20_master_results.csv` (built by `code/20`).

Auto-generated from the analysis pipeline (`code/01–19_*.R`). Numbers come straight from `output/tables/*.csv` and the fitted models in `output/models/*.rds`. Use this as the starting point for the manuscript Results section.

## Experimental snapshot

| Treatment | Wound | n (destructive) | n (also non-destructive) |
|---|---|---|---|
| 28 °C | no  | 48 | 24 |
| 28 °C | yes | 56 | 24 |
| 31 °C | no  | 49 | 24 |
| 31 °C | yes | 55 | 24 |

Three genets (a, c, d), 8 tanks (4 per treatment), wounded on Day 0 after 7 days at target temperature. Non-destructive measurements every ~3 days through Day 14; destructive biopsies on D1, D3, D10, D15.

**Statistical-model structure (per Progress Notes 2026-04-29; refined to integrate genet):**
- Primary fixed structure: `treatment * wound * day * thicket` (full 4-way interactions where data permit).
- Random effects: `(1 | tank) + (1 | id)` (genet promoted from random to fixed because n_genets = 3 — too few for reliable variance estimation).
- Per-genet treatment-effect contrasts and reaction norms are integrated throughout the per-response sections below, not bolted on as a separate analysis.

---

## Key phenotype findings (organismal context)

*These are the Stier-lab phenotype results — the organismal context for the paper. The paper's lead
result is the transcriptomic mechanism (S. Banerjee, lead author; RNA-seq analysis pending).*

1. **Sustained heating to 31 °C compromises photochemistry, pigmentation, growth, and symbiont retention**, with the divergence between treatments growing across the 14-day experiment.
2. **Heat impairs the regeneration phase** (new corallite / skeletal regrowth at the tip, tip extension) **but NOT the tissue-healing phase** (coenosarc coverage — hole closure, polyp emergence, surface smoothing happen at the same rate in both treatments).
3. **Strong genet × treatment interactions across most responses**: genet C is consistently more thermally resilient than genets A and D. Composite thermal-resilience ranking (`figures/19b_genet_resilience_ranking.png`): A (most sensitive) > D >> C (most resilient). Multivariate PCA centroid displacement under heat: A = 3.72, D = 3.34, C = 1.06.
4. **Wounding alone has minimal effect on colony-wide physiology** — the wound response is concentrated in localized morphological changes at the wound site. Heat affects wound site outcomes through the regeneration program rather than closure.

R² estimates (`output/tables/12_r2_summary.csv`) — note marginal R² jumped vs the no-genet baseline because the fixed-effects structure now absorbs the substantial among-genet variation:

| Response | R²m (fixed inc. genet) | R²c (fixed + random) |
|---|---|---|
| PAM Fv/Fm | 0.62 | 0.72 |
| Color (D-scale) | 0.71 | 0.82 |
| Growth (%) | 0.39 | 0.39 |
| log(symbionts cm⁻²) | 0.72 | 0.76 |

---

## 1. Photochemical efficiency (PAM Fv/Fm)

**Model:** `fv_fm ~ treatment * wound * day * thicket + (1|tank) + (1|id)` — `output/models/12_pam_lmm.rds`. n = 336 observations (top/bottom averaged per coral × day).

**Type-III ANOVA, key terms** (Satterthwaite, `lmerTest`; from `output/tables/12_anova_summary.csv`). Note: with `day` uncentered, *main* effects are evaluated at day 0 (wounding), before heat has acted — so the treatment signal lives in the **treatment × day interaction**, not the (conditional) main effect.

| Term | F | df | Reading |
|---|---|---|---|
| treatment (at day 0) | 0.52 | 1, 7 | n.s. — no treatment difference at wounding |
| thicket | 7.75 | 2, 102 | genets differ (p < 0.001) |
| day | 156.0 | 1, 276 | overall change over time |
| **treatment × day** | **106.57** | 1, 276 | **31 °C declines through the experiment (p < 0.001) — the heat signal** |
| treatment × thicket | 0.55 | 2, 102 | n.s. — genets don't differ at day 0 |
| **treatment × day × thicket** | **15.66** | 2, 276 | **genet-specific decline rates (p < 0.001) — the G × E** |

**Per-genet end-of-experiment heat effect (Day 14 Fv/Fm, unwounded; from `output/tables/12_genet_treatment_effects.csv`):**

| Genet | Δ Fv/Fm (28 − 31) | t | p |
|---|---|---|---|
| a | **0.126** | 7.79 | < 0.001 |
| c | **0.046** | 2.90 | 0.007 |
| d | 0.107 | 6.58 | < 0.001 |

Genet C loses 2.7× less photochemical efficiency than genet A under sustained heating. For wounded corals the contrast is even larger — genet C shows a non-significant Fv/Fm change under heat (Δ = 0.016, p = 0.34), while wounded genet A loses 0.10 Fv/Fm units (p < 0.001).

**Day-specific contrasts (pooled across genets, unwounded; from `output/tables/12_emmeans_contrasts.csv`):** divergence detectable from Day 7 onward (Δ Fv/Fm = 0.052, p_adj = 0.007), growing to 0.09 by Day 14 (p_adj < 0.001).

**Figures:** `figures/02_pam_fvfm_trajectory.{pdf,png}`; per-genet trends in panel of `figures/13_genet_response_panel.{pdf,png}`.

---

## 2. Pigmentation (Siebeck D-scale color)

**Model:** `color_num ~ treatment * wound * day * thicket + (1|tank) + (1|id)`. n = 336 obs.

**Type-III ANOVA, key terms** (Satterthwaite; main effects conditional at day 0):

| Term | F | df | Reading |
|---|---|---|---|
| treatment (at day 0) | 0.34 | 1, 7 | n.s. at wounding |
| **treatment × day** | **239.34** | 1, 287 | **dominant signal — heated corals pale through time (p < 0.001)** |
| treatment × wound | 3.03 | 1, 245 | trend only (p = 0.083) |
| treatment × thicket | 0.77 | 2, 60 | n.s. |
| **treatment × day × thicket** | **31.92** | 2, 287 | **genet-specific paling rates (p < 0.001) — the G × E** |

**Per-genet end-of-experiment heat-induced paling (Day 14, unwounded):**

| Genet | Δ Color (28 − 31) | t | p |
|---|---|---|---|
| a | **2.07** D-units | 9.4 | < 0.001 |
| c | **0.63** | 2.86 | 0.007 |
| d | 1.41 | 6.4 | < 0.001 |

Genet C pales 3.3× less than genet A. By Day 14, 0–8% of 28 °C corals had visibly paled vs 58% of 31 °C unwounded and 67% of 31 °C wounded (`output/tables/03_color_end_proportions.csv`).

A treatment × wound interaction is present as a **trend** under type-III SS (F₁,₂₄₅ = 3.03, p = 0.083; note: an earlier type-I computation overstated this as F = 16.9, p < 0.001 — see the type-III correction note below) — wounded heated corals tended to pale slightly *less* than unwounded heated corals. Descriptively this is most pronounced in genet C: heated wounded genet C corals show only a 0.25 D-unit drop vs 28 °C controls (p = 0.27, n.s.), while unwounded heated genet C corals lose 0.63 units (p = 0.007). We report this as a suggestive pattern, not a significant interaction.

**Figures:** `figures/03_color_trajectory.{pdf,png}`.

---

## 3. Calcification (areal, mg CaCO₃ cm⁻² d⁻¹)

**Metric.** Growth is reported as **areal calcification rate** — dry-mass gain normalized to the coral's calcifying surface area and to time (mg CaCO₃ cm⁻² d⁻¹), the field-standard buoyant-weight unit. Calcification is a surface-mediated process, so surface area (living tissue), not skeletal mass, is the mechanistically correct denominator; surface area came from the day-15 wax-dipping standard curve and does not differ by treatment (28 °C 4.46 vs 31 °C 4.73 cm², p = 0.27), so it is an unbiased denominator. Full reasoning and the allometric checks (log-log exponent b = 0.97, mass gain independent of size) are in `notes/growth_allometry.md`.

**Model:** `areal_calc ~ treatment * wound * thicket + (1|tank)` (n = 48; single endpoint observation per coral, with temperature randomized at the tank level). `output/models/12_bw_lm.rds`.

**Tank-level inference:** mean areal calcification fell from **7.62 mg cm⁻² d⁻¹ at 28 °C to 4.75 at 31 °C — a 38% reduction** under sustained heating. The exact tank-level randomization test is a trend rather than a formal p < 0.05 result (28 °C minus 31 °C = 2.87 mg cm⁻² d⁻¹, p = 0.057; `output/tables/05_buoyant_weight_tank_test.csv`). The coral-level LM coefficients remain descriptive and point in the same direction.

**Metric-invariance.** The heat effect has the same direction under alternative growth metrics (`output/tables/05b_growth_metric_comparison.csv`): areal calcification, % mass change, and specific growth rate all decrease at 31 °C. The previous primary metric (% mass change) gives the same descriptive conclusion (28 °C 6.45% vs 31 °C 4.26% over 14 d, a 34% reduction) and is retained as a tank-aware robustness model (`output/models/12_bw_pct_lm.rds`).

**Genet.** The genet × treatment LRT is the only non-significant interaction test in `output/tables/13_genet_anova.csv` (p = 0.37), confirming no statistically detectable G × E for growth (n = 4 per genet × treatment cell — underpowered). Directionally, genet c is the least heat-suppressed in areal calcification (≈−21% vs ≈−45% for genets a and d), consistent with its resilience in PAM, color, and symbiont density, but this is not significant for growth.

**Figure:** `figures/05_buoyant_weight_growth.{pdf,png}`.

---

## 4. Symbiont density (zoox cells cm⁻²)

**Model:** `log(cells_per_cm²) ~ treatment * wound * biopsy_day * thicket + (1|tank)` (n = 192). `output/models/12_zoox_lmm.rds`.

**Type-III ANOVA** (Satterthwaite; `biopsy_day_c` centered at day 1, so main effects are at day 1):

| Term | F | df | Reading |
|---|---|---|---|
| treatment (at day 1) | 3.35 | 1, 9 | marginal at day 1 (p = 0.099) |
| **treatment × biopsy_day** | **93.97** | 1, 162 | **heated corals lose symbionts through time (p < 0.001) — the heat signal** |
| treatment × thicket | 0.70 | 2, 162 | n.s. |
| **treatment × biopsy_day × thicket** | **6.34** | 2, 162 | **genet-specific symbiont-loss rates (p = 0.002) — the G × E** |

**Per-genet end-of-experiment heat effect (Day 15, unwounded; on log scale, so ~0.7 = halving):**

| Genet | Δ log(cells cm⁻²) | t | p |
|---|---|---|---|
| a | **1.45** | 6.4 | < 0.001 |
| c | **0.40** | 1.7 | 0.10 |
| d | 1.39 | 6.1 | < 0.001 |

Genet C retains symbionts under heat almost as well as in ambient (≤30% loss), while genets A and D lose 75–80% of their symbionts.

**Figure:** `figures/06_symbiont_chl_by_day.{pdf,png}`.

---

## 5. Morphological wound-healing characteristics

**Model per trait:** `expressed ~ treatment * day * thicket + (1|tank)`, binomial logit (`code/04_physio_morphology.R`). Restricted to wounded corals (n = 24 per treatment, ~8 per genet × treatment cell).

Recovery follows the **biphasic sequence** standard in the coral wound-healing literature — a **tissue-healing phase** (re-epithelialization / coenosarc coverage that seals the wound) followed by a **regeneration phase** (reappearance of polyps and skeletal calyx/corallite structure). The 9 binary traits index these two phases. (Note: LTH wounds are apical-**tip excisions**, so regeneration is scored as new skeleton at the branch tip, rather than polyp reappearance in a surface wound bed.)

### 5a. Tissue healing — coenosarc coverage (uniform across genets)
- **Hole in center, polyp in hole, wound smoothed** index re-epithelialization / coenosarc coverage over the wound bed: ≥90% expression in both treatments by Day 5–7. No significant treatment effect, no genet × treatment interaction. **Heat does not delay tissue healing.**
- **Polyps out**: significant treatment × day interaction (Wald χ²₁ = 6.70, p = 0.0097) — heated wounded corals more likely to keep polyps retracted at later timepoints than ambient controls. The only healing-phase trait with a real heat response.

### 5b. Regeneration phase — skeletal regrowth at the tip (heat-impaired, genet-dependent)
- **Tip extension**: by Day 15, ~92% of 28 °C corals show tip extension vs 83% of 31 °C. Modest overall heat effect (interval-censored time ratio = 1.16, p = 0.18; first-observed Cox HR = 0.80, p = 0.61).
- **New corallites on tip**: 100% of 28 °C wounded corals form new corallites by Day 15 vs 33% of 31 °C. **Interval-censored time ratio = 1.32 (95% CI 1.19–1.47, p = 1.4×10⁻⁷); first-observed Cox HR = 0.22 (95% CI 0.07–0.69, p = 0.010).**
- **Per-genet Cox HR for new corallites (31 °C vs 28 °C):**

| Genet | HR | n events / 28C | n events / 31C |
|---|---|---|---|
| a | ~0 (no events under 31 °C) | 4 / 4 | 0 / 4 |
| c | 0.64 | 4 / 4 | 3 / 4 |
| d | 0.14 (p = 0.08) | 4 / 4 | 1 / 4 |

LRT for genet × treatment in Cox regression: χ² = 11.3, df = 2, **p = 0.0035** (`output/tables/14_cox_genet_LRT.csv`). The three genets respond significantly differently to heat on this trait — and genet C is again the most resilient.

A similar significant interaction holds for **pigment over wound** (Cox LRT χ² = 13.3, p = 0.001).

**Figures:** `figures/04_morphology_trajectories.{pdf,png}` (pooled by genet), `figures/04b_morphology_trajectories_by_genet.{pdf,png}` (split by genet), `figures/14_morphology_KM.{pdf,png}` (KM by treatment), `figures/14b_morphology_KM_by_genet.{pdf,png}` (KM by treatment × genet).

### 5c. Healing-to-regeneration lag (per-coral event timing)

The cleanest quantification of "heat impairs the regeneration phase, not the healing phase" is the **lag** between achieving coenosarc coverage / tissue healing (`wound_smoothed`) and forming new skeleton at the tip (`new_corallites_on_tip`), computed per coral (`code/14_morphology_kaplan.R`, `output/tables/14_milestone_lag_summary.csv`):

| Treatment | n closed | reached both | **% closed but never regenerated** | median lag (days) |
|---|---|---|---|---|
| 28 °C | 12 | 12 | **0 %** | 8 |
| 31 °C | 12 | 4 | **67 %** | 10 |

At 28 °C every coral that closed its wound went on to regenerate (median 8 d later). At 31 °C **two-thirds of corals closed the wound but never rebuilt skeleton** within the experiment, and the few that did took ~2 d longer. Among corals reaching both milestones the lag difference is marginal (Wilcoxon W = 9, p = 0.061, n = 4 heated) — the dominant signal is the *censored fraction*, not the lag length. This converts the "regeneration vs closure" dichotomy into a single, falsifiable per-coral statistic.

---

## 6. Multivariate genet × treatment signature (`code/15_multivariate.R`)

PCA on the four endpoint responses (PAM, color, growth, log-symbionts) explained 83% of variance on PC1 (the "heat-stress axis") and 13% on PC2 (the "growth axis"). All four physiological variables load positively on PC1.

**Per-genet centroid displacement under heat in PCA space** (Euclidean distance from 28 °C centroid to 31 °C centroid):

| Genet | Displacement |
|---|---|
| a | 3.72 |
| d | 3.34 |
| c | **1.06** |

Genet C's physiological state shifts 3.5× less under heat than genet A. The faceted biplot (`figures/15b_physio_PCA_by_genet.png`) shows that each genet's 28 °C and 31 °C clouds overlap *much more* for genet C than for genets A and D.

---

## 7. Integrative genet ranking (`code/19_genet_dashboard.R`)

The `figures/19_genet_dashboard.png` forest plot collapses every standardized heat-sensitivity effect (4 continuous responses + estimable morphology first-observed HRs) into one figure with rows grouped by domain (physiology, wound closure, regeneration). Composite ranking (`figures/19b_genet_resilience_ranking.png`, `output/tables/19_genet_resilience_summary.csv`):

| Genet | Mean standardized sensitivity | PCA displacement | Rank |
|---|---|---|---|
| **a** | +0.43 | 3.72 | most sensitive |
| **d** | +0.29 | 3.34 | sensitive |
| **c** | −0.03 | 1.06 | resilient |

**Heat-only vs heat-while-wounded decomposition** (`figures/19c_decomposed_resilience.pdf`, `output/tables/19c_resilience_decomp_by_scope.csv`): the genet spread in resilience is largest in the **unwounded** heat response — A and D show standardized heat sensitivities of 0.99 and 0.87, vs C's 0.44. In the **wounded** state the gap compresses (A=0.36, D=0.26, C=0.28); all three genets respond more similarly when wounded. Descriptively, then, the genet differences are sharpest in the unwounded heat response and compress under wounding.

**Implication for the RNA-seq analysis:** comparing the gene expression of genet C against genets A and D — particularly in apical-tip biopsies at Day 10 (the timepoint where new-corallite formation diverges most sharply) — is the most likely route to identify candidate genes underlying heritable thermal tolerance in *A. pulchra*.

---

## 8. Environmental controls

### Water chemistry (YSI daily, n = 72)
`figures/09_ysi_water_chem.png` — temperature, DO, salinity, pH all tracked daily during the experiment. Ramps and steady-state plateaus visible; no detectable cross-tank contamination or DO depletion.

### Apex continuous temperature
`figures/08_apex_temperature.png` — daily-mean per-tank temperature from the Neptune Apex datalogs. Tanks held within ~0.3 °C of nominal 28 °C and 31 °C setpoints across the experimental window.

### Flatworm (AEFW) surveillance
`figures/10_worm_presence.png`, `output/tables/10_worm_summary.csv`. AEFW were checked on three dates (06/07, 06/08, 06/12). 21 unique corals tested worm-positive at some point, with 17/22 worm-positive observations concentrated in 31 °C tanks on **06/07** (3 days post-wounding). All corals were treated with the standard Worm Exit protocol (4 drops per 20 L + transfer-pipet rinse); subsequent surveys showed near-complete clearance. The temperature-biased early-experiment distribution suggests that heat-stressed corals may have been more vulnerable to flatworm reinfection from the source water; this is a *confound to monitor* in any future thermal-tolerance experiment. None of the published response analyses include a worm-presence covariate; sensitivity analyses can be added if Adrian wants to confirm robustness.

---

## 9. Wax-dipping surface-area calibration

Standard curve: SA (cm²) = 0.75 + 61.9 × wax mass (g), based on 15 cylinders of known SA (R² ≈ 0.97, p < 10⁻¹¹). All per-coral SAs in `data/processed/wax_clean.rds` are predicted from this curve and used to normalize buoyant-weight calcification and symbiont counts to area.

**Figure:** `figures/07_wax_standard_curve.{pdf,png}`. **Table:** `output/tables/07_wax_curve_fit.csv`.

---

## 9b. Thermal context — LTH treatments vs *A. pulchra* acute thermal tolerance

We benchmarked the LTH design against an independent, calibrated acute thermal-tolerance dataset for the **same species and island**: Cunning et al. 2024 (*Coral Reefs*) measured CBASS Fv/Fm **ED50** for 20 *A. pulchra* genets from **Mahana, Mo'orea** (`data/external/cunning2024_apulchra_ed50.csv`; `code/sensitivity/26_thermal_context.R`; `figures/26_thermal_context.{pdf,png}`).

**(1) The LTH chronic treatments are sublethal relative to acute limits.** Acute ED50 averaged **35.4 °C** across genets (range 34.4–36.6 °C, SD 0.59). The LTH heated treatment (31 °C) sits **4.4 °C below the mean acute ED50** (3.4 °C below even the most heat-sensitive genet); ambient (28 °C) is 7.4 °C below. The strong, progressive Fv/Fm decline, paling, symbiont loss, and calcification suppression we observe at 31 °C are therefore a **chronic-sublethal** response — accumulated over weeks well below the acute photochemical threshold — not acute photoinhibition. (ED50 is an 18-h acute metric and is used here only as a reference axis, not as the chronic temperature scale.)

**(2) Among-genotype thermal-tolerance variation is concordant across methods.** Cunning's acute assay resolved substantial heritable variation (ED50 range 2.2 °C; ED50 predicts classic bleaching, R = 0.74). The LTH chronic experiment independently detects among-genotype variation in the same population — the resilience ranking **C > D > A** across PAM, color, symbionts, and calcification (`output/tables/26_genotype_variation_concordance.csv`). Both the acute (CBASS) and chronic (LTH) methods thus agree that Mahana *A. pulchra* harbors ecologically meaningful, genotype-level thermal-tolerance variation.

**Caveat / future work.** The LTH thicket labels (A, C, D) are arbitrary and **not genotype-matched** to Cunning's numbered genets, so we cannot yet correlate ED50 against chronic resilience at the individual-genet level. Cunning's genets are genotyped (`genet_map` in the CBASS_methods repo); calling SNPs from the forthcoming LTH RNA-seq and matching them to that map would enable the direct acute-vs-chronic genet correlation — the natural next test.

---

## 10. Statistical limitations and robustness

A diagnostic swarm (`output/diagnostics/{A,B,C,D}_*.{csv,md}`) checked every primary model. Four concerns are worth flagging in the Methods section of the manuscript:

- **Color D-scale is ordinal, not continuous.** The Siebeck D-scale takes five discrete values (D1–D5). We fit it as a Gaussian LMM for direct comparison with the other continuous physiology metrics, and `DHARMa::simulateResiduals` flagged the expected non-uniform residual distribution (KS p < 0.001). As a robustness check we refit the same fixed and random structure as a cumulative-link mixed model (`ordinal::clmm`, `output/models/12b_color_clmm.rds`); every qualitative inference reported above held under the ordinal likelihood. We retain the Gaussian model for presentation; ordinal LRTs are in `output/tables/12b_color_clmm.csv`.

- **Morphology GLMM separation.** Seven of nine binary wound-healing traits show complete or quasi-complete separation in the four-way `treatment × wound × day × thicket` fit. The raw `glmer` fits produce sensible predicted probabilities (used in figures) and correct omnibus type-II Wald χ² from `car::Anova` (used in Sections 5a/5b), but the individual coefficient Wald z-statistics diverge. We therefore refit all 8 traits with weakly-informative Cauchy(0, 2.5) priors on the fixed effects via `blme::bglmer` (Gelman 2008 default for logistic regression with separation; `output/models/12c_morph_<trait>_blme.rds`). The penalized fits yield finite SEs (max SE per trait: 0.005–1.81 except `pigment_over_wound` at 37.2 due to a residual cell of zero events). All quantitative claims on individual morphology coefficients in the master spreadsheet come from these penalized refits (`model_type = "GLMM (binomial, blme Cauchy(0,2.5))"`).

- **Discrete scoring days for morphology timing.** The primary timing tests now use interval-censored Weibull AFT models (`output/tables/14_interval_survreg.csv`) because traits were only scored on visit days. Kaplan-Meier and Cox outputs are retained as first-observed-day summaries and diagnostics.

- **Closure GLMM residuals saturate.** The closure traits `hole_in_center`, `polyp_in_hole`, and `wound_smoothed` reach near-complete occurrence in both temperatures, so the binomial morphology GLMM diagnostic reports DHARMa underdispersion. We retain those diagnostic rows as **HANDLED** rather than treating the GLMM p-values as decisive; closure inference comes from interval-censored milestone timing and event summaries.

- **Cox per-genet events-per-variable.** Per-genet Cox models contain only 4–8 events per genet × treatment cell — below the rule-of-thumb 10 events per variable. We use per-genet Cox output only for visualization/resilience summaries; primary timing inference comes from the interval-censored overall models.

- **One PH violation.** The proportional-hazards assumption was met for every overall Cox model. In the per-genet decomposition, `pigment_over_wound` in genet C shows a Schoenfeld global p = 0.01 (n_event = 5). A time-varying-coefficient refit (`tt(treatment) * log(t+1)`, `output/tables/14c_cox_tt_pigment_genetC.csv`) yields a non-significant heat effect (coef = 0.65, p = 0.13), confirming that the original per-genet point HR was inflated by the PH violation. Schoenfeld plot at `figures/diagnostics/C_pigment_over_wound_genet_c_schoenfeld.png`. We do not propagate the original per-genet HR for this cell into the resilience composite.

- **PAM probe location (top vs bottom).** Molly's original analysis (`code/archive/molly_original/LTH_PAM.R`) noted that top vs bottom probe placement appeared to differ. We test this directly (`code/02_pam_analysis.R` → `output/tables/02b_pam_location_sensitivity.csv`): the **location main effect is significant** (F₁,₆₂₄ = 9.22, p = 0.0025 — a real top/bottom offset), but **every interaction of location with treatment, wound, and day is non-significant** (all p ≥ 0.07). Averaging the two probe locations per coral-day (as the pipeline does) therefore removes a constant offset without distorting the treatment/wound/day effects of interest.

- **Color-card split scores.** 40 of ~962 scored color observations are split Siebeck scores (`D1/D2`, `D2/D3`, `D3/D4`). These are averaged to the midpoint (e.g. `D3/D4` → 3.5) following Molly's convention, rather than truncated to the first digit; the continuous value feeds the Gaussian LMM and the ordinal CLMM robustness check uses deterministic half-up rounding to keep a clean D1–D5 scale (`code/03_color_card_analysis.R`, `code/12_models.R`).

- **Flagged colonies and tanks (QA/QC).** Molly flagged corals 116 and 121 (both wounded genet-c) as visually odd in the per-individual morphology trajectories, tank 3 (ambient) as a slow grower, and a negative within-tank wound effect in tank 11. These are documented in `notes/QAQC_flagged_samples.md` and verified against the data (tank 3 grows ≈40% less than the other three ambient tanks). They are retained, not excluded; a sensitivity analysis dropping corals 116/121 and tank 3 (`code/sensitivity/22_sensitivity_flagged.R`, `output/tables/22_sensitivity_flagged.csv`) leaves every conclusion intact — the PAM, color, and growth treatment effects all remain highly significant, and the growth treatment effect strengthens (F = 23 → 72) once the slow ambient tank is removed.

- **Time-series structure of the repeated-measures responses (tested).** PAM Fv/Fm and color are true repeated measures (~7 observations per coral over 14 days); the primary models fit `day` as a linear fixed effect with random intercepts. A dedicated time-series diagnostic suite (`code/sensitivity/23_timeseries_diagnostics.R`, `output/tables/23_timeseries_diagnostics.csv`, `output/diagnostics/I_timeseries_report.md`) tested the three assumptions that matter for time series:
  - *Temporal autocorrelation* (Molly's concern): present in PAM (AR(1) φ = −0.41, LRT p = 0.001) but absent in color (φ ≈ 0, p = 1.0). **The treatment × day conclusion is robust to refitting with an AR(1) correlation structure** — the interaction p stays ≈10⁻¹⁰ (PAM) and ≈10⁻²⁰ (color). The negative φ in PAM is consistent with fitting a straight line to a curved (nonlinear) trajectory.
  - *Random slopes*: a random slope of day by coral significantly improves fit for both PAM (LRT p < 10⁻⁵) and color (p < 10⁻¹⁸) — individual corals follow different decline trajectories — so the intercept-only model understates among-individual variation, though it does not bias the fixed-effect treatment × day inference.
  - *Linearity of time*: all three trajectories (PAM p = 0.007, color p = 4×10⁻⁴, symbionts p = 0.013) are significantly nonlinear (decline accelerates). The linear `day` term is therefore an approximation of the dominant directional trend; per-day `emmeans` contrasts (reported above) characterize the actual timepoint-specific divergence without assuming linearity.

  None of these change the qualitative conclusions. We tested whether upgrading the headline models is warranted by fitting the full "maximum-rigor" specification — quadratic time + random slopes, plus AR(1) for PAM — and comparing it head-to-head with the current model on the reported day-14 heat effect (`code/sensitivity/24_headline_model_comparison.R`, `output/tables/24_headline_model_comparison.csv`, `output/diagnostics/J_headline_model_comparison.md`):
  - The day-14 effect size shifts by only **+12.5% (PAM, 0.084 → 0.094)** and **+9.4% (color, 1.15 → 1.26)**, and remains highly significant in every specification — conclusion unchanged.
  - The quadratic-time + random-slope model with the full four-way interaction has **worse AIC** than the linear/intercept model (PAM −1199 → −1051; color 234 → 251): the extra interaction parameters are not justified by the fit. The AR(1) model for PAM does improve AIC (−1199 → −1226), confirming the autocorrelation is real, but its day-14 effect (0.083) is nearly identical to the current model (0.084).

  **Verdict: the parsimonious linear / random-intercept model is retained for the headline.** The richer model neither changes the conclusions nor improves the fit (except AR(1) for PAM, which leaves the estimates unchanged); the simplification is therefore justified rather than merely convenient. Residual ACF plots are in `figures/diagnostics/I_*_acf.png`.

- **Type-III ANOVA correction.** The continuous mixed models are now fit with `lmerTest::lmer` and reported with Satterthwaite type-III ANOVA. An earlier version used `lme4::lmer`, whose `anova(type=3)` silently returns **type-I (sequential)** sums of squares — which inflated several main-effect and lower-order interaction F-values (e.g. PAM treatment main 15.8→0.52; color treatment×wound 16.9→3.03). All omnibus F/p in §1–4 are now correct type-III. The primary treatment×time interactions and per-genet `emmeans` contrasts are nearly unchanged; the main consequence is that genotype × heat variation is correctly attributed to the **rate** (treatment × day × thicket), and the color treatment×wound interaction is a non-significant trend rather than significant.

- **Confirmatory vs exploratory testing.** Tests are split into a-priori directed hypotheses grounded in the coral thermal-stress literature (reported **unadjusted**, confirmatory) and exploratory tests with no strong prior prediction (**Benjamini-Hochberg corrected**); `code/sensitivity/28_multiple_testing.R`, `output/tables/28_multiple_testing.csv`. The a-priori set uses the tank-level calcification permutation p and interval-censored survival p-values for timing. Six confirmatory tests are significant (PAM, color, symbionts, morphology new-corallite treatment, interval tip-exist, interval new-corallites); tank-level calcification is a trend (p = 0.057). No exploratory closure-trait effect survives BH correction — heat does not robustly alter wound closure. The regeneration conclusion rests on the new-corallite interval timing result **plus** the censored-fraction and healing-to-regeneration-lag statistics (§5b–c), not on a uniform signal across every tip trait.

- **Variance partitioning (ICC).** `code/sensitivity/27_variance_partitioning.R` reports latent-scale ICC for every mixed model (`output/tables/27_variance_partitioning.csv`). Tank explains a modest fraction of variance in the physiology models (PAM 22 %, color 23 %, symbionts 14 %); among the morphology traits, tip-extension has a high tank ICC (0.83) while most others are 0.06–0.40 — supporting the `(1|tank)` random effect.

- **Probability-scale morphology contrasts.** In addition to log-odds, `code/sensitivity/29_morphology_prob_contrasts.R` reports treatment contrasts on the interpretable Δ-probability scale (plus odds ratios) at day 10 (`output/tables/29_morphology_prob_contrasts.csv`); late traits like new corallites have not yet diverged at day 10, so interval-censored timing/KM framing (§5b–c) remains the primary tool for those.

- **Diagnostic coverage (complete).** Every fitted model has both a result visualization and a residual/assumption diagnostic, audited by `code/sensitivity/25_model_diagnostic_coverage.R` (`output/tables/25_model_diagnostic_coverage.csv`, `output/diagnostics/K_model_coverage_report.md`): **34/34 models covered, 0 gaps.** Continuous LMMs and GLMMs use DHARMa simulated-residual plots where applicable; the ordinal color model uses an observed-vs-fitted check; the penalized morphology GLMMs each have a DHARMa plot; and **all seven overall Cox models now have proportional-hazards (Schoenfeld) diagnostics** (`figures/diagnostics/14_cox_ph_*.png`, `output/tables/14_cox_ph_tests.csv`) — the PH assumption is met for every overall model (all cox.zph p ≥ 0.059) — though the headline `new_corallites_on_tip` model is the closest to the boundary (cox.zph p = 0.059), so its hazard ratio is best read alongside its interval-censored timing model, Kaplan–Meier curve, and the lag statistic (§5b–c) rather than in isolation. Every diagnostic test statistic is also recorded in the master results table (`output/tables/20_master_results.csv`, domains `Time-series diagnostic` and `Survival diagnostic`).

## 11. What's still missing

- **Chlorophyll-a values** — chl-a was planned but ultimately not run, so it is not part of the analysis. Pigmentation/photophysiology is represented by PAM Fv/Fm, color-card scores, and symbiont density.
- **Gene expression data** — RNA-seq libraries at UC Davis Bay lab. Plate layout in `data/raw/plate_layout/Plate_{1,2}.csv`; sequencing plan in `notes/sequencing-plan-keck-LTH.md` (144 libraries: 4 tanks × 3 genets × 2 wound × 3 days × 2 temps). Stub DESeq2 pipeline at `code/21_rnaseq_stub.R` (waiting on sequencing).
- **Statistical write-up for manuscript** — the phenotype Methods + Results are already written into `manuscript/Manuscript_LTH.md` from this narrative; this file (`RESULTS.md`) remains the authoritative phenotype results narrative, and every number traces to `output/tables/20_master_results.csv`.

---

## Pipeline integrity

To reproduce every number and figure in this document:

```bash
cd ~/Stier-LTH-expression-by-temperature-2025
Rscript code/_run_all.R
```

Expected runtime: ~3 minutes (Apex XML parse dominates; first run only).

**R version:** 4.5.2. **Key packages:** tidyverse, lme4, lmerTest, emmeans, DHARMa, MuMIn, car, broom.mixed, survival, patchwork, scales, here, janitor, readxl, xml2, ggrepel.
