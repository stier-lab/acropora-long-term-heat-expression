# Figure Index — LTH (Heat × Wound × *A. pulchra*, Mo'orea 2025)

> **Catalog of every figure → its script** · Updated 2026-07-02 · Index: [`README.md`](../README.md) · figures built by `code/13–16,19,26`.

Lists all 21 PNG figures in `figures/` (each has a matching PDF), the
`code/NN_*.R` script that makes each, and where it appears in the manuscript
(`manuscript/Manuscript_LTH.md`).

---

## Recommended reading order for a first-time reader

Five figures that tell the core story from start to finish:

1. **`08_apex_temperature.png`** — *(a) treatment check:* tank temperatures held at 28 vs 31 °C (the heat ramp and steady state).
2. **`14_morphology_KM.png`** — *(b) healing vs regeneration:* Kaplan–Meier curves showing that tissue still heals under heat but skeleton fails to regrow.
3. **`16_manuscript_fig1.png`** — *(b+c summary:* the 4-panel publication figure — temperature, Fv/Fm, KM healing-vs-regeneration, PCA — the whole story in one figure.
4. **`19_genet_dashboard.png`** — *(c) genet variation:* forest plot of standardized heat sensitivity across all responses, showing the tough genet (c) versus the sensitive ones (a, d).
5. **`26_thermal_context.png`** — *(d) chronic-sublethal thermal context:* puts 28/31 °C next to Cunning et al. 2024 acute CBASS ED50s for the same species and island — 31 °C is chronic but sublethal.

---

## Manuscript reference check

How the manuscript's figure references map to files:

| Manuscript ref | Refers to | Real file? |
|---|---|---|
| Fig 1 (A–D) | Concept plate: colony photo, mesocosm, timeline, wound micrograph | **Not a code-generated stats figure.** Made by `plate_fig.R` (none of the 21 files is named `figures/*1*plate*.png`). Note that `16_manuscript_fig1.png` is a *different* 4-panel figure, NOT the manuscript's "Fig 1." |
| Fig. 2A | Apex temperature trace | maps to `08_apex_temperature.png` |
| Fig. 2B | PAM Fv/Fm trajectory | maps to `02_pam_fvfm_trajectory.png` |
| Fig. 2C | Color (Siebeck D-scale) trajectory | maps to `03_color_trajectory.png` |
| Fig. 2D | Growth (% skeletal mass change) | maps to `05_buoyant_weight_growth.png` |
| Fig. 2E | Symbiont density by biopsy day | maps to `06_symbiont_density_by_day.png` |
| Fig. 3A | PCA biplot of end physiology | maps to `15_physio_PCA_biplot.png` |
| Fig. 3B | Kaplan–Meier morphology (healing vs regeneration) | maps to `14_morphology_KM.png` |
| Fig. 3C | Genet resilience across responses | maps to `19_genet_dashboard.png` / `19b_genet_resilience_ranking.png` |

**Important:** In the manuscript, Fig. 2 and Fig. 3 are *composite* figures. Each lettered panel comes from one of the scripts above, assembled at publication. `code/16` also builds a standalone 4-panel overview (`16_manuscript_fig1`). With this panel-to-script mapping applied, every manuscript figure reference points to a real file.

**Missing-reference flags:** none. Every Fig 1/2/3 reference matches a generated panel.

**Orphan figures (never cited as a numbered manuscript Fig):**
`04_morphology_trajectories.png`, `04b_morphology_trajectories_by_genet.png`,
`07_wax_standard_curve.png`,
`08b_apex_temperature_full.png`, `09_ysi_water_chem.png`,
`10_worm_presence.png`, `13_genet_response_panel.png`,
`14b_morphology_KM_by_genet.png`, `15b_physio_PCA_by_genet.png`,
`16_manuscript_fig1.png` (standalone 4-panel build, not cited by number),
`19b_genet_resilience_ranking.png`, `19c_decomposed_resilience.png`,
`26_thermal_context.png`. These are supporting figures, diagnostics, by-genet
variants, and extra context — not numbered main-text figures.

---

## Full catalog

| File | What it shows | Generating script | Manuscript fig |
|---|---|---|---|
| `02_pam_fvfm_trajectory.png` | PAM Fv/Fm trajectories by treatment × wound × day, 95% CI | `code/02_pam_analysis.R` | Fig. 2B |
| `03_color_trajectory.png` | Color-card pigmentation (Siebeck D-scale) trajectory + paling tallies | `code/03_color_card_analysis.R` | Fig. 2C |
| `04_morphology_trajectories.png` | 8 binary wound-healing traits, cumulative proportion by day × treatment | `code/04_physio_morphology.R` | — |
| `04b_morphology_trajectories_by_genet.png` | Same morphology trajectories, split by genet | `code/04_physio_morphology.R` | — |
| `05_buoyant_weight_growth.png` | Growth (% skeletal mass change) from buoyant weight | `code/05_buoyant_weight.R` | Fig. 2D |
| `06_symbiont_density_by_day.png` | Symbiont density (cells/cm²) by biopsy day | `code/06_symbiont_chl.R` | Fig. 2E |
| `07_wax_standard_curve.png` | Wax-dipping surface-area calibration standard curve | `code/07_wax_dipping.R` | — |
| `08_apex_temperature.png` | Apex per-tank temperature: heat ramp + steady-state (treatment validation) | `code/08_apex_temperature.R` | Fig. 2A |
| `08b_apex_temperature_full.png` | Full-datalog Apex temperature incl. ramps and cooldowns | `code/08_apex_temperature.R` | — |
| `09_ysi_water_chem.png` | Daily YSI spot checks: temp, DO%, DO mg/L, salinity, pH | `code/09_ysi_water_chem.R` | — |
| `10_worm_presence.png` | Acropora-eating flatworm presence tally per tank/treatment (contamination check) | `code/10_worms.R` | — |
| `13_genet_response_panel.png` | 4-panel reaction norms: each genet's mean response shift under heating | `code/13_genet_interaction.R` | — |
| `14_morphology_KM.png` | Kaplan–Meier time-to-onset curves by treatment (healing vs regeneration traits) | `code/14_morphology_kaplan.R` | Fig. 3B |
| `14b_morphology_KM_by_genet.png` | Kaplan–Meier curves by treatment × genet | `code/14_morphology_kaplan.R` | — |
| `15_physio_PCA_biplot.png` | PCA biplot of end-of-experiment physiology (thermal-stress axis) | `code/15_multivariate.R` | Fig. 3A |
| `15b_physio_PCA_by_genet.png` | PCA biplot faceted by genet | `code/15_multivariate.R` | — |
| `16_manuscript_fig1.png` | Standalone 4-panel publication figure: temp / Fv/Fm / KM healing-vs-regen / PCA | `code/16_main_figure.R` | — (integrative build) |
| `19_genet_dashboard.png` | Forest plot of standardized heat sensitivity across all responses, per genet | `code/19_genet_dashboard.R` | Fig. 3C |
| `19b_genet_resilience_ranking.png` | Composite thermal-resilience score ranking per genet | `code/19_genet_dashboard.R` | Fig. 3C (companion) |
| `19c_decomposed_resilience.png` | Heat sensitivity decomposed by wound state (wounded vs unwounded) | `code/19_genet_dashboard.R` | — |
| `26_thermal_context.png` | LTH 28/31 °C vs Cunning et al. 2024 acute CBASS ED50s — chronic-sublethal context | `code/sensitivity/26_thermal_context.R` | — |
