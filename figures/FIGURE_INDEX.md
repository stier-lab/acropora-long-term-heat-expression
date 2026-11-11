# Figure Index

Updated 2026-11-11. Every figure in this folder is rebuilt by `code/run_all.R` (each as `.pdf` and `.png`). Frozen figures from exploratory and older analyses are in `extras/figures/`, with the previous index at `extras/figures/FIGURE_INDEX_before_merge.md`. Reader-facing text says **source patch** for A, C and D; some file names still say `genet`.

## Paper

| Figure | Content | Script |
|---|---|---|
| `39b_stage_timing_cumulative_option` | **Figure 1B.** Percent of wounded fragments that had reached each healing or regeneration stage by each day; fragments that never reached a stage stay in the denominator. Decided 2026-11-11. | `04_regeneration.R` |
| `45_physiology_figure2` | **Figure 2.** Fv/Fm (A), symbiont density (B) and 15-day skeletal growth (C) by temperature and wound status; error bars use the four tank means per temperature. Legend: `manuscript/FIGURE_2_LEGEND.md`. | `03_physiology.R` |

## Supplement candidates

| Figure | Content | Script |
|---|---|---|
| `08_apex_temperature` | Daily mean temperature in each tank, Day −7 to 15. | `02_tanks_water.R` |
| `09_ysi_water_chem` | Salinity, dissolved oxygen, pH and temperature spot checks. | `02_tanks_water.R` |
| `07_wax_standard_curve` | Wax-dipping standard curve used for surface area. | `01_clean_data.R` |
| `14_morphology_KM` | Kaplan–Meier curves for each regeneration stage, by temperature. | `04_regeneration.R` |
| `11_microscope_trait_trajectories` | Microscope experiment: share of the 16 fragments showing each stage, by day. | `07_supplement.R` |
| `13_genet_response_panel` | End-of-experiment means by source patch and temperature (Fv/Fm and color Day 14, symbionts Day 15, growth Days 0–15). Used in the team summary, section 4. | `05_source_patches.R` |
| `43_final_snp_structure`, `43_final_snp_design`, `43_final_snp_symbiont_density`, `43_final_snp_symbiont_typing` | Final host clusters, how they spread across treatments, symbiont density within clusters, and rough symbiont read shares. | `06_genetics.R` |

## Team summary and checks

| Figure | Content | Script |
|---|---|---|
| `03_color_trajectory`, `05_buoyant_weight_growth`, `06_symbiont_density_by_day` | Single-response displays used in the team summary's supporting section. | `03_physiology.R` |
| `36_stage_timing`, `36b_stage_violin` | First day each fragment reached each stage, with X marks for fragments that never did. Supporting displays; don't use averages for stages some fragments never reached. | `04_regeneration.R` |
| `43_final_snp_diagnostics` | Residual checks for the cluster model. | `06_genetics.R` |

Residual and proportional-hazards plots from `08_checks.R` go to `figures/diagnostics/`, which is not tracked in git.
