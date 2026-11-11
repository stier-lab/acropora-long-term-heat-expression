# Model diagnostics audit: LTH core analyses

2026-11-11. Run with the `stats-diagnostics` skill on every statistical result in `code/01`–`08`. Full table: `model_diagnostics_2026-11-11.csv` (18 analysis families). Call scan: `stats_call_scan_2026-11-11.csv`. Scripts to reproduce: `scripts/` (run from the repo root after `code/run_all.R`).

Status labels: **yes** = adequate for the current claim; **handled** = a concern exists but a specific check shows it doesn't change the conclusion; **partial** = a specific gap remains; **not applicable** = no fitted model (counts, summaries).

## Fixes applied (2026-11-11, commit after 785a06b)

All four fixes and the three "state in the text" items below were applied the same day. Numbers now in `output/paper/`:

| What | Before | After |
|---|---|---|
| Fv/Fm temperature × day | F₁,₂₂₈ = 120 | F₁,₇ = 18.2, p = 0.004 |
| Color temperature × day | F₁,₂₂₈ = 288 | F₁,₆ = 34.3, p < 0.001 |
| Symbionts temperature × day | F₁,₁₆₂ = 95 | F₁,₇ = 19.7, p = 0.003 |
| Symbionts temperature × wound × day | p = 0.11 | p = 0.079; one tank left out: 0.009–0.29 (`22b_sensitivity_wound_growth.csv`) |
| Patch × temperature (Fv/Fm, color, symbionts, growth) | χ² on 14, 14, 14, 6 df | χ²₈ = 60.1, 108.0, 54.6; χ²₄ = 13.0 (all p ≤ 0.011) |
| Tip extension time ratio | 1.36, p = 0.023 (Weibull) | 1.16, p = 0.135 (log-logistic) |
| Budding time ratio | 1.32, p = 1e-7 (Weibull) | 1.28, p = 3e-7 (log-logistic) |
| Wound smoothing time ratio | 0.83, p = 0.022 | 0.89, p = 0.14 |

What changed in the code:

- `03_physiology.R`: Fv/Fm and color use `(1 + day || tank) + (1 | id)`; symbionts use `(1 + biopsy_day_c || tank)`. Per-patch heat contrasts gain a Holm-adjusted `p_holm` column (6 per response).
- `05_source_patches.R`: the reduced model drops only the terms containing source patch × temperature (8 df; 4 for growth), with the same random effects as 03. This now matches the methods text.
- `04_regeneration.R`: interval-censored models use the log-logistic distribution (best or within 0.6 AIC for every reported stage; pigment over the wound falls back to lognormal). `14_interval_survreg_distributions.csv` records the AIC comparison. One correction to the audit below: lognormal is *not* within 2.3 AIC of the best for every stage (it is 7.8 behind for tip formation), so log-logistic was used instead.
- `07_supplement.R`: the flagged-fragment refits use the 03 structure (Fv/Fm p = 0.028, color p = 0.007 without fragments 116, 121 and tank 3). The new `22b_sensitivity_wound_growth.csv` covers the symbiont wound term (with the new model, two fragments are flagged as outliers, 163 and 212; without them p = 0.009) and growth's dependence on tank 3 (without it p = 0.002; leaving out tank 12 gives p = 0.108; 7 of 48 single-fragment drops give p > 0.05).
- Team summary: methods and results text updated (tank slopes, log-logistic, tip extension "no clear difference", the symbiont wound caveat, growth and tank 3), plus "out of date" notes listing the superseded numbers.

Status after the fixes: A, C, F and I move from partial to handled; E becomes yes (Holm reported); N becomes handled. Remaining partial: J (Cox p-values are fragment-level; cite the tank permutation) and O (exploratory genetics).

## Bottom line (audit as run, before the fixes)

No headline conclusion changes. Heat lowered Fv/Fm, color and symbiont density over time and stopped most wounded fragments from budding by Day 15, under every check. Four things need fixing before the numbers go into the manuscript:

1. **The physiology F tests overstate the evidence (A, B, C).** Temperature was applied to tanks, but the models let every observation inform the temperature × day slope separately. Models with a day slope per tank fit much better (ΔAIC −28 to −36). In them, temperature × day is still clear, but on about 6 denominator df instead of 228:

   | Response | Reported | With a day slope per tank |
   |---|---|---|
   | Fv/Fm | F = 120, p ≈ 0 | F = 18–20, p ≈ 0.004 |
   | Color | F = 288 | F = 34–38, p < 0.001 |
   | Symbionts | F = 95 | F = 17–20, p ≈ 0.003–0.006 |

   Fix: refit with `(1 + day || tank) + (1 | id)` (symbionts: `(1 + day || tank)`) and report those tests.

2. **"Clipping added nothing" is fragile for symbionts (C).** The wound × heat × day term is p = 0.11 as reported, 0.08 with tank slopes, and 0.003 without five heated Day-15 patch-A fragments. Three of those were nearly bleached out (13,000–19,000 cells/cm²), and they look like real bleaching, not errors. Its direction is that heated *wounded* fragments kept more symbionts (Day-15 geometric means 135,000 vs 105,000 cells/cm²). Say that the wound result for symbionts depends on the most-bleached fragments. For color the same term is unstable too (p = 0.02–0.90 when any one tank is dropped).

3. **The source-patch test isn't the one named (F).** "Patch × temperature" compares 14 df: every patch interaction, including patch × wound and patch × day. A targeted patch × temperature test with day slopes still gives p < 10⁻⁵ for Fv/Fm, color and symbionts, and p = 0.014 for growth. Fix: switch to the targeted test or rename it.

4. **Tip extension "took somewhat longer under heat" depends on the Weibull distribution (I).** Weibull fits these intervals worst (ΔAIC 12). Better-fitting distributions give p = 0.09 (lognormal) and 0.14 (log-logistic). Wound smoothing behaves the same way. Budding's time ratio (1.26–1.32) is stable across distributions. Fix: choose the distribution by AIC (log-logistic is best or near-best for every stage) and drop or soften the tip-extension timing claim.

Also state in the text:

- **Growth depends on tank 3 (D).** Temperature p = 0.046 as reported. Dropping any one fragment gives 0.023–0.069, and dropping any one tank 0.002–0.11. Ambient tank 3, flagged as slow to equilibrate, grew like the heated tanks. The tank permutation gives p = 0.057 with tank 3 and 0.029 without it.
- **Budding p-values below 10⁻⁴ (G, I, J, L) treat 24 fragments in 8 tanks as independent.** Clustered errors aren't usable with 8 tanks. The tank permutation (H: every ambient tank 3/3, every heated tank 1/3, p = 0.029) is the tank-level evidence.
- **The per-patch heat contrasts (E) are unadjusted.** With Holm adjustment only one call changes (color, patch C unwounded, 0.025 → 0.050).

Already handled: the budding patch × temperature Cox test (K, now reported as counts) and color's non-normal residuals (B, tank permutation p = 0.029). Supporting analyses were fine: per-patch contrasts, the Fisher tests, the permutation designs, the genetics model (exploratory), the wax standard curve, and the Figure 2 error bars (tank means).

## Analysis groups and the question each answers

| Group | Question | Families |
|---|---|---|
| Tank conditions | Did the tanks hold their temperatures and water quality? | Q |
| Physiology | Did heat, clipping, or both change Fv/Fm, color, symbionts and growth? | A, B, C, D, R |
| Source patches | Did the heat effect differ among patches A, C and D? | E, F, K |
| Regeneration | Did heat change whether and when wounded fragments healed and budded? | G, H, I, J, L, M |
| Supplement | Do results hold without flagged fragments? | N |
| Genetics | Does heat affect symbionts within genetic clusters? | O |
| Calibration | Are surface areas consistent with the wax curve? | P |

## Status by family

| ID | Family | Role | Status (as audited) | Changes the conclusion? |
|---|---|---|---|---|
| A | Fv/Fm over time | primary | partial | No; reported F and df overstated |
| B | Color over time | primary | handled | No; F overstated, wound term unstable |
| C | Symbiont density | primary | partial | No for heat; wound result fragile |
| D | Growth | primary | handled | No; depends on tank 3 |
| E | Heat effect per patch | support | partial | No; state family and Holm p |
| F | Patch × temperature | primary | partial | No; test mislabeled |
| G | Stage present Day 15 (Fisher) | primary | handled | No |
| H | Budding by tank (permutation) | primary | yes | No |
| I | Time to stage (AFT) | primary | partial | Yes for tip extension timing |
| J | Budding hazard (Cox) | support | partial | No |
| K | Budding patch × temperature | exploratory | handled | Fixed 2026-11-11 |
| L | Log-rank | support | handled | No |
| M | Reached by Day 15; tissue cover | derived | not applicable | No |
| N | Without flagged fragments | sensitivity | partial | No; refit with tank slopes later |
| O | Symbionts within clusters | exploratory | partial | No claim yet |
| P | Wax standard curve | derived | yes | No |
| Q | Tank conditions | derived | not applicable | No |
| R | Figure 2 error bars | derived | yes | No |

## Diagnostic artifacts

| Artifact | Made by | Reports |
|---|---|---|
| `output/tables/18_model_checks.csv` | `08_checks.R` | Convergence, singular fits, DHARMa uniformity / dispersion / outliers for A–D |
| `figures/diagnostics/08_dharma_*.png` | `08_checks.R` | DHARMa residual plots for A–D (not tracked) |
| `output/tables/14_cox_ph_tests.csv`, `figures/diagnostics/14_cox_ph_*.png` | `08_checks.R` | Schoenfeld PH tests for J (all p > 0.05) |
| `output/tables/18_validation_summary.csv` | `08_checks.R` | Data validation (0 FAIL) |
| `output/tables/12b_color_tank_test.csv` | `07_supplement.R` | Tank permutation for B |
| `output/tables/05_buoyant_weight_tank_test.csv` | `03_physiology.R` | Tank permutation for D |
| `output/tables/22_sensitivity_flagged.csv` | `07_supplement.R` | N |
| `output/tables/43_final_snp_model_diagnostics.csv`, `43_final_snp_design_balance.csv` | `06_genetics.R` | O |
| `scripts/diag_lmm.R` | this audit | A–D: variance components, heteroscedasticity, autocorrelation, curvature, random slopes, leave-one-tank-out |
| `scripts/diag_lmm2.R`, `diag_lmm3.R` | this audit | A–C: slopes per tank and fragment, unequal variance (nlme), symbiont outliers |
| `scripts/diag_patch.R` | this audit | F: reported vs targeted patch × temperature tests |
| `scripts/diag_surv.R` | this audit | G–J: events and censoring, separation, AFT distributions, clustered SEs, ties |
| `scripts/diag_rest.R` | this audit | D, E, G, H, N, O: influence, contrast families, Fisher tables, permutations |

Not yet in the pipeline: the tank-slope refits, leave-one-tank-out and AFT distribution checks exist only in `scripts/`. Once the models are updated, the ones that support reported numbers should move into `08_checks.R`.
