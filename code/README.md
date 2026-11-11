# Code guide

Updated 2026-11-11. Ten files make every number, table and figure the paper uses.

```bash
Rscript code/run_all.R    # under a minute; rebuilds everything below
```

| File | What it does | Paper section |
|---|---|---|
| `00_setup.R` | Packages, folder paths, colors, plot theme, random seed. Loaded by every script. | none |
| `01_clean_data.R` | Turns each raw coral spreadsheet into a clean table in `data/processed/`. No statistics. | Methods |
| `02_tanks_water.R` | Tank temperature from the Apex logs; salinity, oxygen and pH from the YSI meter. | Methods, tank conditions |
| `03_physiology.R` | Mixed models for Fv/Fm, color, symbionts and growth; Figure 2. | Results 1 |
| `04_regeneration.R` | Day 15 stage counts, timing models, stage table; Figure 1B. | Results 2 |
| `05_source_patches.R` | Does the heat effect differ among source patches A, C and D? | Results 3 |
| `06_genetics.R` | Final genetic clusters; one-row-per-library table for the RNA-seq team. | Results 4 |
| `07_supplement.R` | Microscope-experiment figure, results without flagged fragments, tank-level color check. | Supplement |
| `08_checks.R` | Data validation, model residual checks, Cox proportional-hazards checks. | none (QA) |
| `run_all.R` | Runs 01–08 in order. | none |

Each script starts with a header that says what it reads and writes, and is split into numbered sections.

## Where the reported numbers are

Scripts 02–05 end with a section that writes one table to `output/paper/`:

| Table | Written by |
|---|---|
| `output/paper/tanks_water_numbers.csv` | `02_tanks_water.R`, section 3 |
| `output/paper/physiology_numbers.csv` | `03_physiology.R`, section 6 |
| `output/paper/regeneration_numbers.csv` | `04_regeneration.R`, section 5 |
| `output/paper/source_patches_numbers.csv` | `05_source_patches.R`, section 3 |

Each row is one number used in the text, with its test, degrees of freedom, p-value, sample size and the table it came from. The statistics boxes in `docs/team_summary/LTH_results_summary.html` read these files. Don't type numbers into the draft by hand. Copy them from here, or rerun and recopy.

## Output names

Output files kept the names they had before the scripts were merged, so links in the draft, the team summary and the Zenodo plan still work. The number at the front of an output name tells you which script writes it:

| Output prefix | Written by |
|---|---|
| `01_`, `07_`, `10_`, `11_microscope_design_crosswalk` | `01_clean_data.R` |
| `08_`, `09_` | `02_tanks_water.R` |
| `03_`, `05_`, `06_`, `12_` (except `12b_`), `45_` | `03_physiology.R` |
| `11f_`, `14_` (except `14_cox_ph_tests`), `36_`, `36b_`, `39_`, `39b_` | `04_regeneration.R` |
| `13_` | `05_source_patches.R` |
| `31_`, `43_` | `06_genetics.R` |
| `11_microscope_trait_*`, `12b_`, `22_` | `07_supplement.R` |
| `14_cox_ph_tests`, `18_` | `08_checks.R` |

## Not reported

- **Daily yes/no regeneration models.** Their p-values change with the computer and package versions (`docs/provenance/rerun_check_2026-11-10.md`). They are frozen in `extras/`. The regeneration results use stage counts, tank counts and the timing models instead.
- **The average day a stage was reached** for stages some fragments never reached. For budding under heat, that average describes only 4 of 12 fragments.

## After a run

Rerunning rewrites tracked files. To see which changes are real:

```bash
python3 tools/check_output_changes.py            # real changes vs noise
python3 tools/check_output_changes.py --restore  # undo the noise-only files
```

PDFs change on every run even when nothing real did (R embeds the time and a random font tag), so the tool treats a PDF as unchanged when its PNG twin is unchanged. Run R in a UTF-8 locale (the Mac default); otherwise symbols such as ° and × are written as byte codes.

## `extras/`

Everything else is in `extras/`: exploratory analyses, older audits and diagnostics, the RNA-seq starter code, Molly's original scripts, and the scripts as they were before the merge. Those scripts are frozen and `run_all.R` does not run them. See `extras/README.md`.

## Where things went (Nov 11, 2026)

The 26 scripts the paper used were merged into 01–08. Two tags mark the earlier layouts: `pre-refactor-2026-11-11` (one flat folder) and `pre-merge-2026-11-11` (grouped folders, just before the merge).

| Old script | Now |
|---|---|
| `analysis/01_metadata`, `02_photosynthesis` (cleaning), `03_color` (cleaning), `04_regeneration_scores` (cleaning), `05_growth` (calculation), `06_symbionts` (cleaning), `07_surface_area`, `11_microscope_experiment` (cleaning), `checks/10_flatworms` | `01_clean_data.R` |
| `analysis/08_temperature`, `09_water_quality`, `paper/tanks_water` | `02_tanks_water.R` |
| `analysis/12_models` (continuous models), `03_color`, `05_growth` (tank test, figure), `06_symbionts` (summary, figure), `45_figure2_physiology`, `paper/physiology` | `03_physiology.R` |
| `analysis/11c_stage_counts` (tests), `14_stage_timing_models`, `36_stage_timing_tables`, `39_stage_timing_figures` (39b), `paper/regeneration` | `04_regeneration.R` |
| `analysis/13_source_patch_tests`, `paper/source_patches` | `05_source_patches.R` |
| `rnaseq/31_rnaseq_covariate_table`, `analysis/43_genetic_clusters` | `06_genetics.R` |
| `analysis/11_microscope_experiment` (figure), `checks/22_sensitivity_flagged` | `07_supplement.R`, which also adds a tank-level color check in place of the ordinal color model (that model could not be estimated) |
| `checks/18_data_validation`, `analysis/14_stage_timing_models` (PH tests), `checks/diagnostics/A` (residual checks, rewritten shorter) | `08_checks.R` |
| `analysis/12_models` and `04_regeneration_scores` (daily yes/no models), `20_master_results_table`, everything in `exploratory/`, other `checks/`, `rnaseq/21`, `32`, plate figure, `archive/` | `extras/` (frozen) |
| `run_paper.R`, `tools/check_output_changes.py` | `run_all.R` now takes under a minute; the checker moved to `tools/` at the repo root |
