# Figure Index

Updated 2026-11-12. Every figure in this folder is rebuilt by `code/run_all.R` (each as `.pdf` and `.png`). The few kept exploratory figures are in `extras/figures/`; older figures were deleted on 2026-11-13 and can be recovered from the `pre-prune-2026-11-13` tag. Reader-facing text says **source patch** for A, C and D; some file names still say `genet`.

## Key (no figure carries an embedded legend)

Since 2026-11-12 the artwork has no titles, subtitles, captions or legends; keys and sample sizes belong in the caption text. Mapping used by every script (`code/00_setup.R` palettes):

- **Temperature:** blue `#56B4E9` = 28 °C, orange `#D55E00` = 31 °C. In Figure 1B and `39d`, 28 °C is also solid with circles and 31 °C dashed with triangles.
- **Wounding:** Figure 2 and its panels: open symbols, white boxes, dashed lines = unwounded; filled symbols, gray boxes, solid lines = wounded. `47`: dashed = unwounded, solid = wounded. `03_color_trajectory`: gray = unwounded, black = wounded. `43_final_snp_symbiont_density`: circles = unwounded, triangles = clipped.
- **Source patch:** gold `#E69F00` = A, green `#009E73` = C, dark blue `#0072B2` = D (`13`, `05`, `43_final_snp_structure`).
- **Symbiont read group** (`43_final_snp_symbiont_typing`): green = A, pink = B, gold = C, dark blue = D.
- **Stage figures** (`36`, `36b`): X = not reached by Day 15; right-side labels in `36` are reached/total.
- **Multi-panel files** use letter tags (A, B, C): `45_physiology_figure2`, `43_final_snp_structure`, `43_final_snp_diagnostics`.

## Paper

| Figure | Content | Script |
|---|---|---|
| `39b_stage_timing_cumulative_option` | **Figure 1B.** Percent of wounded fragments (12 per temperature, main experiment) that had reached each of the four regeneration milestones by each day: 1 axial polyp forms, 3 axial corallite forms, 4 axial corallite extends, 5 radial corallites form. Fragments that never reached a milestone stay in the denominator. | `04_regeneration.R` |
| `45_physiology_figure2` | **Figure 2.** Fv/Fm (A), symbiont density (B) and 15-day skeletal growth (C) by temperature and wound status; error bars use the four tank means per temperature. Legend: `manuscript/FIGURE_2_LEGEND.md`. | `03_physiology.R` |

## Supplement

Supplementary figures 1 (phases schematic) and the subsampling schematic are drawings, not code output. Figure 1A and the photo parts of Supplementary fig. 3 are photographs. See `manuscript/MAIN_FIGURE_PLAN.md`.

| Figure | Content | Script |
|---|---|---|
| `39d_stage_cumulative_all_stages` | All six stages, including healing (microscope experiment) and stage 2, wound surface remodels. Same layout as Figure 1B. | `04_regeneration.R` |
| `47_other_daily_observations` | **Supplementary fig. 2.** Polyps extended, pigment over wound, algae on wound, and flatworm (AEF) signs, by day. All fragments were alive at every check. | `07_supplement.R` |
| `03_color_trajectory` | **Supplementary fig. 4** (color card). Color score over time by temperature and wound status. | `03_physiology.R` |
| `08_apex_temperature` | Daily mean temperature in each tank, Day −7 to 15. | `02_tanks_water.R` |
| `09_ysi_water_chem` | Salinity, dissolved oxygen, pH and temperature spot checks. | `02_tanks_water.R` |
| `07_wax_standard_curve` | Wax-dipping standard curve used for surface area. | `01_clean_data.R` |
| `14_morphology_KM` | Kaplan–Meier curves for each regeneration stage, by temperature. | `04_regeneration.R` |
| `11_microscope_trait_trajectories` | **Supplementary fig. 3** (data panel). Microscope experiment: share of the 16 fragments showing healing and stages 1-5, by day. The microscope photos and stage diagrams are separate artwork. | `07_supplement.R` |
| `13_genet_response_panel` | End-of-experiment means by source patch and temperature. | `05_source_patches.R` |
| `43_final_snp_structure`, `43_final_snp_design`, `43_final_snp_symbiont_density`, `43_final_snp_symbiont_typing` | Final host clusters, how they spread across treatments, symbiont density within clusters, and rough symbiont read shares. | `06_genetics.R` |

## Team summary and checks

| Figure | Content | Script |
|---|---|---|
| `05_buoyant_weight_growth`, `06_symbiont_density_by_day` | Single-response displays used in the team summary's supporting section. | `03_physiology.R` |
| `36_stage_timing`, `36b_stage_violin` | First day each fragment reached each stage, with X marks for fragments that never did. Supporting displays; don't use averages for stages some fragments never reached. | `04_regeneration.R` |
| `43_final_snp_diagnostics` | Residual checks for the cluster model. | `06_genetics.R` |

Residual and proportional-hazards plots from `08_checks.R` go to `figures/diagnostics/`, which is not tracked in git.
