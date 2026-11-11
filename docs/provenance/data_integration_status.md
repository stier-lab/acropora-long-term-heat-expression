# Data Integration Status

Last checked: 2026-09-06

Status categories:

- `full_analysis`: cleaned from raw data, written to processed/output artifacts,
  and included in the master pipeline.
- `supporting_analysis`: cleaned and included in the master pipeline, but mainly
  as design, calibration, environmental, or quality-control context rather than
  as a primary response model.
- `handoff_only`: used to build a downstream handoff table; no statistical
  analysis is run in this repo.
- `preliminary_analysis`: imported and analyzed only as an explicitly
  preliminary/exploratory layer; not final phenotype evidence.
- `provenance_only`: imported and characterized, but not read by the current
  analysis pipeline.

| source group | status | current use |
|---|---|---|
| `metadata` | `supporting_analysis` | Cleaned by `code/01_clean_data.R`; writes `data/processed/coral_metadata.rds` and the metadata summary. This is the design spine for downstream joins and validation. |
| `pam` | `full_analysis` | Cleaned/analyzed by `code/01_clean_data.R` (cleaning) and `code/03_physiology.R`; writes `data/processed/pam_clean.rds`, PAM figures, and treatment contrasts. |
| `color_card` | `full_analysis` | Cleaned/analyzed by `code/01_clean_data.R` (cleaning) and `code/03_physiology.R`; writes `data/processed/color_clean.rds`, color figures, end proportions, and contributes to the main models. |
| `physio_morphology` | `full_analysis` | Cleaned/analyzed by `code/01_clean_data.R` (cleaning) and `code/04_regeneration.R`; writes `data/processed/physio_clean.rds`, morphology trajectory figures, GLMM summaries, and feeds survival/KM/model scripts. Also feeds the cross-dataset morphology plots, endpoint/timing tests, and diagnostic coverage checks in `code/04_regeneration.R`. |
| `microscope_physio` | `full_analysis` | Cleaned/analyzed by `code/01_clean_data.R` (cleaning), `code/04_regeneration.R` (tests) and `code/07_supplement.R` as a separate photo-only data set; writes cleaned data, design/event/trajectory tables, and microscope figures. Also feeds the cross-dataset morphology plots, endpoint/timing tests, and diagnostic coverage checks in `code/04_regeneration.R`. |
| `wax_dipping` | `supporting_analysis` | Cleaned/analyzed by `code/01_clean_data.R`; writes `data/processed/wax_clean.rds`, the wax standard-curve figure, and surface areas used by growth/symbiont calculations. |
| `buoyant_weight` | `full_analysis` | Cleaned/analyzed by `code/01_clean_data.R` (calculation) and `code/03_physiology.R`; writes `data/processed/buoyant_weight_clean.rds`, growth figures, growth model tables, tank tests, and metric comparisons. |
| `symbiont_counts` | `full_analysis` | Cleaned/analyzed by `code/01_clean_data.R` (cleaning) and `code/03_physiology.R`; writes `data/processed/symbiont_chl_clean.rds`, symbiont-density figure, summaries, and RNA-seq covariates. Chlorophyll-a is explicitly handled as not run. |
| `apex` | `supporting_analysis` | Updated 2026-09-11: Molly's recovered XML and compiled sheet fill the late gap. All eight tanks have 144 ten-minute readings per day through Days 0-15 (and Day 16). `code/02_tanks_water.R` verifies duplicate agreement, checks the compiled sheet, uses dated tank assignments, and exports hourly/daily data plus coverage and clock-sensitivity tables. See `docs/provenance/molly_followup_2026-09-11.md`. |
| `ysi` | `supporting_analysis` | Updated October 27 from Molly's daily-log spreadsheet. Script 09 imports the dated raw snapshot, checks units and legacy overlap, and exports cleaned data, coverage, descriptive summaries, and water-chemistry figures. June 4–19 includes 112 readings from eight tanks on 14 dates; June 4 and 7 are missing. SD across readings and SE among tank averages are reported separately. See `docs/provenance/ysi_refresh_2026-10-27.md`. |
| `worm_presence` | `supporting_analysis` | Cleaned/summarized by `code/01_clean_data.R`; writes `data/processed/worm_clean.rds`, worm-presence figure, and summary table. Used as contamination/QC context. |
| `trinity_hauru_colonies` | `supporting_analysis` | Preliminary external Hauru colony table attached to Trinity Conn's 2026-08-10 email and forwarded by Molly Brzezinski on 2026-09-03. Read by `extras/code/exploratory/26b_trinity_hauru_context.R`, which writes a cleaned CSV, nearest-coordinate table, and map. This is context for future DNA-marker matching only, not a confirmed genotype match to LTH source patches. |
| `plate_layout` | `handoff_only` | Read by `code/06_genetics.R`; writes RNA-seq library lookup and phenotype-covariate handoff tables. No expression model is fit in this repo. |
| `rnaseq` | final host calls; exploratory phenotype analysis | October 12 `FINAL_LTH_genoclusters_syms.csv` supersedes the September file. Script 43 joins all 144 final host assignments and rough symbiont read shares to the Days 1/3/10 covariates. Memberships are unchanged; PCs are updated. Density-model diagnostics remain partial. No Day 15 genotypes or expression counts were delivered. |
| `shipping` | `provenance_only` | Imported/codebooked sample handling and storage metadata. The current pipeline does not read these sheets directly. |
| `daily_health_log` | `provenance_only` | Imported/codebooked for completeness; mostly a filled design/template log and not used by the current pipeline. |
| `printable_data_sheets` | `provenance_only` | Imported/codebooked collection templates for YSI, PAM, color-card, and daily-health logs; not analysis inputs. |
| `wax_dipping_standard_curve` | `provenance_only` | Standalone standard-curve copy retained for traceability. The current pipeline uses `data/raw/wax_dipping/Standard_curve.csv`. |
| `buoyant_weight_archive` | `provenance_only` | Duplicate/working buoyant-weight calculation workbook retained for traceability. The current pipeline uses `data/raw/buoyant_weight/` tab CSVs. |
| `chlorophyll` | `provenance_only` | Example calculation workbook only. Chlorophyll-a was planned but not populated/run for this project. |
| Drive method docs | `provenance_only` | Exported to `docs/provenance/drive_methods/`; used to document historical design decisions, methods, photo-source paths, and interpretation. The root README, `docs/README.md`, and current scripts override stale planning language in those exports. |
| NAS image sources | `provenance_only` | Local NAS paths and file counts are documented in `docs/provenance/nas_image_sources.md`; raw images are not imported or pixel-analyzed in this repo. |

Phenotype validation status after this classification:

`39 PASS, 3 HANDLED, 0 WARN, 0 FAIL`

Molly follow-up outputs added 2026-09-05 and updated 2026-09-06:

- `output/tables/36_temperature_coverage_days0_16.csv`
- `output/tables/36_stage_timing.csv`
- `output/tables/36_stage_timing_summary.csv`
- `output/tables/36_followup_action_items.csv`
- `figures/36_stage_timing.png`
