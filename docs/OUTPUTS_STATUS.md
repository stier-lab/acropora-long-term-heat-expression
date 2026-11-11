# Which outputs are current

Updated 2026-11-10. Use this page to tell current results from historical files before reading figures or tables.

## Current

| What | Where | Notes |
|---|---|---|
| Every reported number | `output/tables/20_master_results.csv` | Single source of truth. `code/30_manuscript_audit.R` flags manuscript text that drifts from it. |
| Plain-language results | `docs/team_summary/LTH_results_summary.html` | Rebuilt 2026-10-27 with fragment 119 confirmed wounded. |
| Physiology figure | `figures/45_physiology_figure2.*` | Publication candidate. Legend: `manuscript/FIGURE_2_LEGEND.md`. |
| Stage-timing figures | `figures/39b_stage_timing_cumulative_option.*` (main display in the team summary), `figures/41_stage_progression.*` (recommended in `FIGURE_INDEX.md`) | Both keep fragments that never reached a stage in the denominator. Pick one for the paper. |
| Genetics | `figures/43_final_snp_*`, `output/tables/43_final_snp_*` | Final host clusters from Rachael (2026-10-13). Symbiont read shares are rough. |
| Morphology scores | `data/raw/physio_morphology/data.csv`, `data/raw/microscope_physio/data.csv` | October 16 corrected exports. Cell-level changes: `output/tables/44_scoring_changes.csv`. |
| Water chemistry | `figures/09_ysi_water_chem.*`, `output/tables/09_ysi_*` | October 27 refresh from the corrected daily log. |
| Figure catalog | `figures/FIGURE_INDEX.md` | See the caveat below about the "Start Here" table. |

## Historical: keep for provenance, don't cite

| What | Superseded by |
|---|---|
| `figures/32_prelim_snp_*`, `output/tables/32_prelim_snp_*`, `code/32_prelim_snp_phenotype_integration.R` | Script 43 (final SNP set) |
| `figures/39a_*`, `39c_*`, `40_stage_mean_review*` | `39b` or `41`. These options can misrepresent the 12/12 vs 4/12 budding result (commit 2ffa570). The team summary still shows `40` with caveats. |
| `data/raw/*/data_before_2026-10-16.csv` | The October 16 corrected `data.csv` in the same folder |
| `data/metadata/119_before_correction_anova_2026-10-16.csv` | Current models with 119 as wounded |
| `output/models/superseded/` | Failed or stale robustness fits, archived automatically by `code/12_models.R` |
| `RESULTS.docx`, `manuscript/Manuscript_LTH.docx`, `figures/FIGURE_INDEX.docx` | Not reconciled with the October refreshes. Use the team summary and the master results table. |
| `notes/archive/` | Old plans and meeting notes |

Open caveat: the "Start Here" table in `figures/FIGURE_INDEX.md` still points to `36_molly_stage_timing` and `36b_molly_stage_violin`, while the same file later recommends `41_stage_progression`. Treat 41 as the current stage-timing display until the index is updated.

## Why files show as changed after running the pipeline

`code/_run_all.R` rewrites every table, figure and processed `.rds` file, and all of those are tracked in git. Two things follow.

1. **Real changes.** When code or raw data change, most outputs genuinely change. Between 2026-10-13 and 2026-10-27, three rounds of updates (final SNPs, October 16 scoring corrections, the fragment 119 label and the YSI refresh) were run but only partly committed. That is why about 170 files showed as changed. They were committed on 2026-11-10.
2. **Noise.** Some files change on every run even when nothing real did:
   - PDFs embed the time they were written
   - `output/session_info.txt` depends on which packages were loaded
   - text outputs can gain or lose trailing spaces

To sort the two after a run:

```bash
Rscript code/_run_all.R
python3 code/tools/check_output_changes.py            # lists real changes vs noise
python3 code/tools/check_output_changes.py --restore  # undoes the noise-only files
git status                                            # what's left is real
```

Commit code, raw-data and output changes together, with a message that says what changed and why, so the outputs in git always match the code that produced them.

## Reporting decisions

- **Growth temperature test (2026-11-10, Adrian):** report the model test (temperature, F1,6 = 6.33, p = 0.046, from `12_models.R`). The tank permutation (p = 0.057, `05_buoyant_weight_tank_test.csv`) goes in the supplement as a check.
