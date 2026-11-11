# Which outputs are current

Updated 2026-11-11. Use this page to tell current results from historical files before reading figures or tables.

## Current

| What | Where | Notes |
|---|---|---|
| Every reported number | `output/paper/*_numbers.csv` | Single source of truth since 2026-11-11. One file per results section, written by the last section of `code/02`–`05` (run `code/run_all.R`). |
| Plain-language results | `docs/team_summary/LTH_results_summary.html` | Rebuilt 2026-11-11 with a statistics drop-down (methods and results text) under each finding, read from `output/paper/`. |
| Physiology figure | `figures/45_physiology_figure2.*` | Publication candidate. Legend: `manuscript/FIGURE_2_LEGEND.md`. |
| Stage-timing figure | `figures/39b_stage_timing_cumulative_option.*` | Paper figure and team-summary display (decided 2026-11-11). `41_stage_progression` is the alternative. Both keep fragments that never reached a stage in the denominator. |
| Genetics | `figures/43_final_snp_*`, `output/tables/43_final_snp_*` | Final host clusters from Rachael (2026-10-13), built by `code/06_genetics.R`. Symbiont read shares are rough. |
| Morphology scores | `data/raw/physio_morphology/data.csv`, `data/raw/microscope_physio/data.csv` | October 16 corrected exports. Cell-level changes: `extras/output/tables/44_scoring_changes.csv`. |
| Water chemistry | `figures/09_ysi_water_chem.*`, `output/tables/09_ysi_*` | October 27 refresh from the corrected daily log. |
| Figure catalog | `figures/FIGURE_INDEX.md` | Rewritten 2026-11-11 for the figures the core scripts make. |
| Model diagnostics audit | `docs/diagnostics/MODEL_DIAGNOSTICS_2026-11-11.md` | All 18 analysis families checked 2026-11-11; the fixes were applied the same day (see below and the report's "Fixes applied" section). No headline conclusion changed. |
| Everything else | `extras/` | Frozen 2026-11-11: exploratory analyses, older audits, the daily yes/no models, the master results table, and their outputs. Not rerun. See `extras/README.md`. |

## Historical: keep for provenance, don't cite

| What | Superseded by |
|---|---|
| `extras/figures/32_prelim_snp_*`, `extras/output/tables/32_prelim_snp_*` | `code/06_genetics.R` (final SNP set) |
| `extras/figures/39a_*`, `39c_*`, `40_stage_mean_review*` | `39b`. These options can misrepresent the 12/12 vs 4/12 budding result (commit 2ffa570). Removed from the team summary 2026-11-11. |
| `data/raw/*/data_before_2026-10-16.csv` | The October 16 corrected `data.csv` in the same folder |
| `data/metadata/119_before_correction_anova_2026-10-16.csv` | Current models with 119 as wounded |
| `extras/output/tables/20_master_results*.csv` | `output/paper/`. The master table also held exploratory and daily yes/no model results. |
| `RESULTS.docx`, `manuscript/Manuscript_LTH.docx`, `extras/figures/FIGURE_INDEX_before_merge.*` | Not reconciled with the October refreshes. Use the team summary and `output/paper/`. |
| `notes/archive/` | Old plans and meeting notes |


## Why files show as changed after running the pipeline

`code/run_all.R` rewrites every table, figure and processed `.rds` file, and all of those are tracked in git. Two things follow.

1. **Real changes.** When code or raw data change, most outputs genuinely change. Between 2026-10-13 and 2026-10-27, three rounds of updates (final SNPs, October 16 scoring corrections, the fragment 119 label and the YSI refresh) were run but only partly committed. That is why about 170 files showed as changed. They were committed on 2026-11-10.
2. **Noise.** Some files change on every run even when nothing real did:
   - PDFs embed the time they were written and a random font tag (the checker treats a PDF as unchanged when its PNG twin is unchanged)
   - `output/session_info.txt` depends on which packages were loaded
   - text outputs can gain or lose trailing spaces

To sort the two after a run:

```bash
Rscript code/run_all.R
python3 tools/check_output_changes.py            # lists real changes vs noise
python3 tools/check_output_changes.py --restore  # undoes the noise-only files
git status                                            # what's left is real
```

Commit code, raw-data and output changes together, with a message that says what changed and why, so the outputs in git always match the code that produced them.

## Reporting decisions

- **Growth temperature test (2026-11-10, Adrian):** report the model test (temperature, F1,6 = 6.33, p = 0.046, from `code/03_physiology.R`). The tank permutation (p = 0.057, `05_buoyant_weight_tank_test.csv`) goes in the supplement as a check.
- **Stage-timing figure (2026-11-11, Adrian):** `39b` (cumulative percent reached by day) is Figure 1B. `40_stage_mean_review` is dropped from the summary.
- **Daily yes/no morphology models (2026-11-11):** exploratory only; don't report their p-values. They change with the computer and package versions (`docs/provenance/rerun_check_2026-11-10.md`). Regeneration results use stage counts, tank counts and the timing models.
- **Code layout (2026-11-11):** the paper's code is ten files in `code/` (`00_setup.R`, `01`–`08`, `run_all.R`); everything else is frozen in `extras/`. Earlier layouts are tagged `pre-refactor-2026-11-11` and `pre-merge-2026-11-11`. Name map in `code/README.md`.
- **Color model check (2026-11-11):** color residuals are not normal (a 1–6 score where almost every ambient fragment stayed at 4; `output/tables/18_model_checks.csv`). An ordered-category model cannot be estimated on these data (its old output is frozen in `extras/`). The check is a tank-level exact permutation on final color (`output/tables/12b_color_tank_test.csv`, `code/07_supplement.R`): every heated tank had 4–6 of 6 fragments paled versus 0–1 in ambient tanks, p = 0.029, the smallest possible with 4 tanks per side.
- **Physiology models (2026-11-11):** Fv/Fm and color use a random intercept and day slope for each tank plus an intercept for each fragment; symbionts use an intercept and day slope for each tank. Temperature × day is now tested on about 6–7 denominator df (Fv/Fm F = 18.2, p = 0.004; color F = 34.3, p < 0.001; symbionts F = 19.7, p = 0.003). Earlier F values (95–288 on 160–230 df) are out of date.
- **Symbiont wound result (2026-11-11):** temperature × wound × day p = 0.079, but 0.009–0.29 when one tank is left out (`output/tables/22b_sensitivity_wound_growth.csv`); where there is a difference, heated wounded fragments kept more symbionts. Report it with that caveat.
- **Growth and tank 3 (2026-11-11):** the temperature effect (p = 0.046) depends on ambient tank 3: without it p = 0.002; leaving out tank 12 instead gives p = 0.108. Say so in the text.
- **Patch × temperature test (2026-11-11):** now drops only the terms containing source patch × temperature (8 df; 4 for growth), matching the methods text. Per-patch heat contrasts carry Holm-adjusted p (`p_holm`).
- **Timing models (2026-11-11):** interval-censored log-logistic (best or within 0.6 AIC for every reported stage; `14_interval_survreg_distributions.csv`), replacing Weibull. Tip extension no longer differs clearly (time ratio 1.16, p = 0.135); budding 1.28, p < 0.001. For significance, cite the tank permutation, not fragment-level p-values.
- **Budding by source patch (2026-11-11):** don't report the Cox source patch × temperature p-value (χ² = 11.3, p = 0.0035). No heated patch-A fragment budded, so the model's coefficients are infinite. Report the counts instead (heated: A 0/4, C 3/4, D 1/4; ambient 4/4 each; `output/paper/regeneration_numbers.csv`).
- **Color baseline labels (2026-11-11):** the color sheet records every fragment as unwounded at the Day −4 baseline (before clipping). `code/01_clean_data.R` now takes wound groups from the metadata, as the Fv/Fm data does. Models use Day 0 onward, so no reported number changed; the Day −4 points in `03_color_trajectory` now show all four groups.
