# Code map

Updated 2026-11-10. This page lists every script, what it does, and whether its results are current. Which outputs are current, and why files show as changed after a run: `docs/OUTPUTS_STATUS.md`.

## At a glance

| Folder | R scripts | Run by `_run_all.R` |
|---|---:|---|
| `code/` | 42 | 38 (not `00_setup.R`, which every script loads itself; not `_run_all.R`; not `21_rnaseq_stub.R` or `plate_fig.R`) |
| `code/sensitivity/` | 10 | 9 (not `40b_stage_source_audit.R`) |
| `code/diagnostics/` | 8 | all 8 |
| `code/archive/molly_original/` | 4 | none (Molly's original scripts, kept for reference) |
| `code/tools/` | 1 Python | none (run by hand after the pipeline) |
| **Total** | **64 R + 1 Python**, about 17,300 lines | |

Run everything from the repo root with `Rscript code/_run_all.R` (about 5 minutes). Every script starts by loading `code/00_setup.R` (packages, paths, colors), so most scripts can also be run on their own once `data/processed/` exists.

## Where to find the numbers for each result

All reported statistics are collected in one table, `output/tables/20_master_results.csv`, built by `20_master_results_table.R`. The statistics boxes in `docs/team_summary/LTH_results_summary.html` give draft methods and results text for each finding, with its source.

| Result | Script that runs the test | Table to read |
|---|---|---|
| Tank temperatures | `08_apex_temperature.R` | `08_apex_treatment_summary.csv` |
| Salinity, oxygen, pH | `09_ysi_water_chem.R` | `09_ysi_study_summary.csv`, `09_ysi_treatment_summary.csv` |
| Fv/Fm, color, symbionts, growth models | `12_models.R` | `20_master_results.csv` (model type "LMM"; growth "LM"), `12_anova_summary.csv` |
| Growth tank test | `05_buoyant_weight.R` | `05_buoyant_weight_tank_test.csv` |
| Figure 2 means (Fv/Fm, symbionts, growth) | `45_physiology_figure2.R` | `45_figure2_summary.csv` |
| Final color | `03_color_card_analysis.R` | `03_color_end_proportions.csv` |
| Day 15 symbionts | `06_symbiont_chl.R` | `06_symbiont_chl_summary.csv` |
| Regeneration stage counts (Fisher, log-rank) | `11c_morphology_dataset_plots.R` | `11f_morphology_dataset_endpoint_tests.csv`, `11f_morphology_dataset_logrank_tests.csv` |
| Time to each stage | `14_morphology_kaplan.R` | `14_interval_survreg.csv`, `14_cox_hazard_ratios.csv` |
| Stage percentages by day | `39_stage_timing_visual_options.R` | `39_stage_cumulative_percent.csv` |
| Tissue cover (microscope experiment) | `11_microscope_physio.R` | `11_microscope_event_summary.csv` |
| Source-patch differences | `13_genet_interaction.R` | `13_genet_anova.csv`, LRT rows of `20_master_results.csv` |
| Final genetic clusters | `43_final_snp_integration.R` | `43_final_snp_*.csv` |

**Don't report** the p-values from the daily yes/no morphology models (`12_models.R` morphology section, `04_physio_morphology.R` tables, `12_morph_*` and `12c_*` tables). They change with computer and package versions; see `docs/provenance/rerun_check_2026-11-10.md`.

## Status labels

- **Current:** feeds a result or figure in the team summary or manuscript.
- **Data prep:** cleans raw data into `data/processed/`.
- **Check:** validates data, models or outputs; doesn't produce a reported result.
- **Exploratory:** analyses we ran to look; not planned for the main paper.
- **Superseded:** kept so older outputs can be traced; don't use for new text.
- **Handoff:** for the RNA-seq team.

## Main scripts (`code/`)

| Script | What it does | Status |
|---|---|---|
| `00_setup.R` | Packages, file paths, plot theme and colors. Loaded by every script. | Utility |
| `_run_all.R` | Runs the whole pipeline in order and records package versions. | Utility |
| `01_load_clean_metadata.R` | Cleans the master metadata: one row per fragment, with treatment, tank, source patch and wound status. | Data prep |
| `02_pam_analysis.R` | Cleans Fv/Fm readings; plots Fv/Fm over time. Its contrasts are earlier versions of the `12_models.R` results. | Data prep + figure |
| `03_color_card_analysis.R` | Converts color-card scores to numbers; plots color over time; final color table. | Current |
| `04_physio_morphology.R` | Cleans the daily regeneration scores; plots the share of fragments at each stage. Its model tables are not reportable (see above). | Data prep + figure |
| `05_buoyant_weight.R` | Calculates growth; tank-level permutation test for temperature. | Current |
| `06_symbiont_chl.R` | Symbiont density per cm²; means by group and day. | Current |
| `07_wax_dipping.R` | Converts wax mass to surface area for the symbiont samples. | Data prep |
| `08_apex_temperature.R` | Rebuilds tank temperatures from the Apex logs; treatment means. | Current |
| `09_ysi_water_chem.R` | Salinity, oxygen and pH summaries and methods wording. | Current |
| `10_worms.R` | Tallies flatworm-positive fragments by tank. | Check |
| `11_microscope_physio.R` | Separate 16-fragment microscope experiment: stage timing and tissue cover. | Current |
| `11c_morphology_dataset_plots.R` | Compares the two regeneration datasets; Day 15 stage counts, Fisher and log-rank tests. | Current |
| `12_models.R` | Main models for Fv/Fm, color, symbionts and growth (current); daily yes/no morphology models (not reportable). Saves models to `output/models/`. | Current (physiology) |
| `13_genet_interaction.R` | Tests whether the heat effect differs among source patches. | Current |
| `14_morphology_kaplan.R` | Time to each regeneration stage: interval-censored, Cox and Kaplan–Meier. | Current |
| `15_multivariate.R` | PCA of the four physiology responses. | Exploratory |
| `16_main_figure.R` | Earlier four-panel Figure 1 draft. The current physiology figure is `45`. | Superseded |
| `17_figure_audit.R` | Flags figures that failed to render or are oversized. | Check |
| `18_data_validation.R` | Checks sample sizes, missing values and ID matches before analysis. | Check |
| `19_genet_dashboard.R` | Standardized heat penalties by source patch across responses. | Exploratory |
| `20_master_results_table.R` | Collects every statistic into `20_master_results.csv`. | Current |
| `21_rnaseq_stub.R` | Optional starting code for importing RNA-seq counts. Not run. | Handoff |
| `30_manuscript_audit.R` | Checks that numbers in `manuscript/Manuscript_LTH.docx` match the tables. That draft is older, so mismatches are expected. | Check |
| `31_rnaseq_covariate_table.R` | One row per RNA-seq library with its phenotype data. | Handoff |
| `32_prelim_snp_phenotype_integration.R` | September preliminary genetic clusters. Replaced by `43`. | Superseded |
| `33_growth_tradeoff_screen.R` | Growth vs heat-tolerance screen. | Exploratory |
| `34_regeneration_tradeoff_screen.R` | Growth vs regeneration screens. | Exploratory |
| `35_tradeoff_summary.R` | Four-panel trade-off figure for the team summary. | Exploratory |
| `36_molly_followup_checks.R` | Stage-timing tables and Molly's September follow-up checks. Its tables are current; its figures `36` and `36b` are older displays. | Current (tables) |
| `37_tradeoff_formal_tests.R` | Formal trade-off tests and diagnostics. | Exploratory |
| `38_analysis_diagnostic_inventory.R` | Lists every analysis and which checks cover it. | Check |
| `39_stage_timing_visual_options.R` | Three stage-timing figure designs. `39b` (cumulative) is the main display in the team summary; `39a` and `39c` are not used. | Current (`39b`) |
| `40_molly_stage_mean_review.R` | Average day each stage was reached, among fragments that reached it. Shown in the summary with caveats; misleading for budding under heat. | Use with care |
| `41_stage_progression_figure.R` | Stage timing with fragments that never reached a stage kept visible. Recommended in `figures/FIGURE_INDEX.md`. | Current (alternative to `39b`) |
| `42_pam_location_and_score_audit.R` | Tip vs base Fv/Fm comparison; review of individual morphology scores. | Check |
| `43_final_snp_integration.R` | Final genetic clusters from the Bay lab (Oct 2026) joined to the phenotype data. | Current |
| `44_corrected_scoring_audit.R` | Lists every score changed in Molly's Oct 16 corrections. | Check |
| `45_physiology_figure2.R` | Figure 2 candidate: Fv/Fm, symbionts and growth. | Current |
| `46_wound_label_validation.R` | Confirms fragment 119 is wounded and compares results before and after the fix. | Check |
| `plate_fig.R` | Molly's RNA-seq plate-layout figure. Not run. | Handoff |

## Sensitivity analyses (`code/sensitivity/`)

| Script | What it does | Status |
|---|---|---|
| `22_sensitivity_flagged.R` | Reruns the temperature tests without the fragments and tank Molly flagged in QA. | Check |
| `23_timeseries_diagnostics.R` | Checks whether straight-line time trends and fragment intercepts are adequate for the repeated measures. | Check |
| `24_headline_model_comparison.R` | Compares the main Fv/Fm and color models with more flexible versions. | Check |
| `25_model_diagnostic_coverage.R` | Makes sure every saved model has a residual plot. | Check |
| `26_thermal_context.R` | Compares our temperatures with published heat-tolerance values for *A. pulchra* in Mo'orea (Cunning et al. 2024). | Context |
| `26b_trinity_hauru_context.R` | Matches our source patches to Trinity Conn's Hauru colonies by location. | Context |
| `27_variance_partitioning.R` | How much variation sits among tanks, fragments and residuals. | Check |
| `28_multiple_testing.R` | Separates planned tests from exploratory ones and adjusts the exploratory p-values. | Check |
| `29_morphology_prob_contrasts.R` | Morphology differences as changes in probability. Built on the non-reportable daily models. | Check |
| `40b_stage_source_audit.R` | Checks disputed scores against archived sheets. Not run by the pipeline. | Check |

## Model checks (`code/diagnostics/`)

All run by `_run_all.R`. They write to `output/diagnostics/` and `figures/diagnostics/`, which are not tracked in git.

| Script | What it checks |
|---|---|
| `A_continuous_lmm.R` | Residuals, convergence and influential points for the four physiology models |
| `B_morphology_glmm.R` | The same for the daily morphology models, plus separation |
| `C_cox_diagnostics.R` | Proportional-hazards assumption and events per variable for the timing models |
| `D_pca_lrt.R` | The PCA and the source-patch likelihood-ratio tests |
| `E_design_alignment.R` | That each model's structure matches the experimental design |
| `F_model_reproducibility.R` | That refitting each saved model gives the same estimates |
| `G_diagnostic_plots.R` | Fills in any missing diagnostic plots |
| `H_spreadsheet_coverage.R` | That every statistic reached `20_master_results.csv` |

## Other

- `code/archive/molly_original/`: Molly's original PAM, buoyant weight, color-card and morphology scripts. See the README in that folder.
- `code/tools/check_output_changes.py`: after a pipeline run, lists real changes vs timestamp and whitespace noise. Usage is in `docs/OUTPUTS_STATUS.md`.

## Known loose ends

- Nine scripts start with a short comment instead of the standard Purpose/Input/Output header: `08`, `09`, `40`, `42`, `43`, `44`, `45`, `46` and `sensitivity/40b`.
- Several `37_tradeoff_*` and `42_pam_csv_parse_audit.csv` tables store full file paths from Adrian's Mac, so they show as changed on any other computer.
- Script and figure names containing "molly" (`36`, `40`) record who asked for the analysis; they are not a separate dataset.
