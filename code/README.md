# Code guide

Updated 2026-11-11. The code is grouped by purpose. To write the paper you only need the `paper/` folder and its outputs.

## Start here

```bash
Rscript code/run_paper.R    # about 2 minutes: rebuilds everything the results section uses
```

That writes four files to `output/paper/`, one per results section. Each row is one number used in the text, with its test, degrees of freedom, p-value, sample size, and the table it came from:

| Results section | Script | Output |
|---|---|---|
| Tank temperature and water quality | `paper/tanks_water.R` | `output/paper/tanks_water_numbers.csv` |
| Photosynthesis, color, symbionts, growth (Figure 2) | `paper/physiology.R` | `output/paper/physiology_numbers.csv` |
| Healing and regeneration (Figure 1) | `paper/regeneration.R` | `output/paper/regeneration_numbers.csv` |
| Differences among source patches | `paper/source_patches.R` | `output/paper/source_patches_numbers.csv` |

Each paper script starts with a short methods description for its section. The statistics boxes in `docs/team_summary/LTH_results_summary.html` read these same files and turn them into draft methods and results text.

To rebuild everything, including exploratory analyses, model checks and the full results table (`output/tables/20_master_results.csv`), run `Rscript code/run_all.R` (about 5 minutes). Do this before sharing results or rebuilding the team summary.

## Folders

| Folder | Scripts | What's in it | Run by |
|---|---:|---|---|
| `paper/` | 4 | One script per results section; writes `output/paper/` | both |
| `analysis/` | 19 | Data cleaning and the analyses the paper uses | both (`20` and `43` only in `run_all.R`) |
| `exploratory/` | 11 | Analyses we ran to look; not planned for the main paper | `run_all.R` |
| `checks/` | 23 | Data validation, model diagnostics, sensitivity analyses, audits | `run_all.R` |
| `rnaseq/` | 4 | Handoff material for the RNA-seq team | `run_all.R` (31, 32) |
| `archive/` | 5 | Earlier drafts and Molly's original scripts | not run |
| `tools/` | 1 (Python) | `check_output_changes.py` | by hand |

Plus `00_setup.R` (packages, paths, colors; loaded by every script), `run_paper.R` and `run_all.R`.

In `analysis/`, the number at the start of each script name matches the prefix of the files it writes (for example, `12_models.R` writes `output/tables/12_*.csv`), so you can trace any output back to its script.

## `analysis/`

| Script | What it does |
|---|---|
| `01_metadata.R` | Cleans the master metadata: one row per fragment, with treatment, tank, source patch and wound status |
| `02_photosynthesis.R` | Cleans Fv/Fm readings; plots Fv/Fm over time |
| `03_color.R` | Converts color-card scores to numbers; plots color over time; final color table |
| `04_regeneration_scores.R` | Cleans the daily regeneration scores; plots the share of fragments at each stage. Its yes/no models are exploratory (see below) |
| `05_growth.R` | Calculates growth; tank permutation test for temperature |
| `06_symbionts.R` | Symbiont density per cm²; means by group and day |
| `07_surface_area.R` | Converts wax mass to surface area for the symbiont samples |
| `08_temperature.R` | Rebuilds tank temperatures from the Apex logs |
| `09_water_quality.R` | Salinity, oxygen and pH summaries |
| `11_microscope_experiment.R` | Separate 16-fragment microscope experiment |
| `11c_stage_counts.R` | Day 15 stage counts, Fisher and log-rank tests, for both regeneration datasets |
| `12_models.R` | Mixed models for Fv/Fm, color, symbionts and growth; also the daily yes/no morphology models (exploratory) |
| `13_source_patch_tests.R` | Tests whether the heat effect differs among source patches |
| `14_stage_timing_models.R` | Time to each regeneration stage: interval-censored, Cox and Kaplan–Meier |
| `20_master_results_table.R` | Collects every statistic in the repo into `output/tables/20_master_results.csv` |
| `36_stage_timing_tables.R` | Stage-by-stage counts and first-observation days; temperature coverage checks |
| `39_stage_timing_figures.R` | Stage-timing figures; `39b` (cumulative percent by day) is the paper figure |
| `43_genetic_clusters.R` | Final genetic clusters from the Bay lab joined to the phenotype data |
| `45_figure2_physiology.R` | Figure 2: Fv/Fm, symbionts and growth |

**Don't report** p-values from the daily yes/no morphology models (`12_models.R` morphology section; `12_morph_*` and `12c_*` tables; `04_regeneration_scores.R` model tables). They change with the computer and package versions; see `docs/provenance/rerun_check_2026-11-10.md`. The regeneration results use stage counts, tank counts and the timing models instead.

## `exploratory/`

`15_physiology_pca.R`, `19_source_patch_dashboard.R`, `26_thermal_context.R` (comparison with published heat tolerance, Cunning et al. 2024), `26b_trinity_hauru_context.R`, `29_morphology_prob_contrasts.R`, `33`–`35` and `37` (growth and regeneration trade-offs), `40_stage_mean_among_reachers.R` (average day among fragments that reached each stage; misleading for budding under heat), `41_stage_progression_figure.R` (an alternative stage-timing figure).

## `checks/`

`10_flatworms.R`, `17_figure_audit.R`, `18_data_validation.R`, `22_sensitivity_flagged.R` (results without the QA-flagged fragments and tank), `23`–`25` and `27` (model adequacy and coverage), `28_multiple_testing.R`, `30_manuscript_audit.R` (compares the older manuscript draft with the tables), `38_analysis_inventory.R`, `40b_stage_source_audit.R` (not run), `42_pam_location_and_score_audit.R`, `44_corrected_scoring_audit.R` (Molly's Oct 16 score corrections), `46_wound_label_validation.R` (fragment 119), and `diagnostics/A`–`H` (residuals, convergence, reproducibility and coverage for every model). Their outputs in `output/diagnostics/` and `figures/diagnostics/` are not tracked in git.

## `rnaseq/` and `archive/`

- `rnaseq/`: `21_rnaseq_import_stub.R` (optional starting code), `31_rnaseq_covariate_table.R` (phenotype data per library), `32_prelim_snp_integration.R` (September preliminary clusters, superseded by `analysis/43`; still run because two later scripts read its join table), `plate_layout_figure.R`.
- `archive/`: `16_main_figure_draft.R` (an earlier Figure 1 draft; its figure is no longer rebuilt) and `molly_original/` (Molly's first PAM, growth, color and morphology scripts).

## After a run

Rerunning rewrites tracked files. To see which changes are real:

```bash
python3 code/tools/check_output_changes.py            # real changes vs noise
python3 code/tools/check_output_changes.py --restore  # undo the noise-only files
```

PDFs change on every run even when nothing real did (R embeds the time and a random font tag), so the tool treats a PDF as unchanged when its PNG twin is unchanged. Run R in a UTF-8 locale (the Mac default); otherwise symbols like ° and × are written as byte codes.

## Where things moved (Nov 11, 2026)

The old flat `code/` layout is tagged `pre-refactor-2026-11-11`. Old → new names: `01_load_clean_metadata` → `analysis/01_metadata`; `02_pam_analysis` → `analysis/02_photosynthesis`; `03_color_card_analysis` → `analysis/03_color`; `04_physio_morphology` → `analysis/04_regeneration_scores`; `05_buoyant_weight` → `analysis/05_growth`; `06_symbiont_chl` → `analysis/06_symbionts`; `07_wax_dipping` → `analysis/07_surface_area`; `08_apex_temperature` → `analysis/08_temperature`; `09_ysi_water_chem` → `analysis/09_water_quality`; `11_microscope_physio` → `analysis/11_microscope_experiment`; `11c_morphology_dataset_plots` → `analysis/11c_stage_counts`; `13_genet_interaction` → `analysis/13_source_patch_tests`; `14_morphology_kaplan` → `analysis/14_stage_timing_models`; `36_molly_followup_checks` → `analysis/36_stage_timing_tables`; `39_stage_timing_visual_options` → `analysis/39_stage_timing_figures`; `43_final_snp_integration` → `analysis/43_genetic_clusters`; `45_physiology_figure2` → `analysis/45_figure2_physiology`; `_run_all.R` → `run_all.R`. Files named `36_molly_*` and `36b_molly_stage_violin` are now `36_*` and `36b_stage_violin`. The `sensitivity/` and `diagnostics/` folders are now inside `checks/` (three context and contrast scripts went to `exploratory/`). Dated notes in `docs/provenance/` keep the old names.
