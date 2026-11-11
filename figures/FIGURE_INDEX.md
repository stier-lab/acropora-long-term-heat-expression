# Figure Index

## Publication Candidate (October 16)

`45_physiology_figure2.pdf` and `.png` combine photosynthesis efficiency,
symbiont density, and 15-day skeletal growth. Temperature columns and wound
symbols are consistent; descriptive error bars use four tank means. See
`manuscript/FIGURE_2_LEGEND.md`. Confirm fragment 119's wound label before
finalizing this candidate. The summary uses this figure as its main physiology
display; previous distributions remain in an expandable supporting section.

Morphology and trade-off figures were regenerated from the October 16 corrected
scoring sheets. Details: `docs/provenance/scoring_refresh_2026-10-16.md`.

## Current Genetic Figures (October 13)

The `43_final_snp_*` figures supersede the historical `32_prelim_snp_*` figures
for current interpretation. Script 43 uses final host clusters and separately
labels the new symbiont read-mapping shares as rough.

| Figure | Content |
|---|---|
| `43_final_snp_structure` | Final genetic coordinates and source-patch memberships |
| `43_final_snp_design` | Full cluster x temperature x wound x day cell counts, including zeros |
| `43_final_snp_symbiont_density` | Same-fragment density observations by cluster, day and treatment |
| `43_final_snp_symbiont_typing` | Rough relative mapped-read shares; not cell proportions |
| `43_final_snp_diagnostics` | Residual and Q-Q checks; inference remains exploratory |

Updated: 2026-09-09

This is the current figure catalog for the LTH heat x wound repo. The older
`FIGURE_INDEX.docx` is retained for provenance but may lag the current pipeline
and terminology. Reader-facing text should use **source patch** for A, C, and D;
some filenames still use legacy `genet` wording.

## Start Here

| Figure | Use |
|---|---|
| `08_apex_temperature.png` | Confirms the tanks held the 28 C and 31 C temperature treatments. |
| `14_morphology_KM.png` | Shows that heat mostly blocked new-corallite budding rather than wound sealing. |
| `39b_stage_timing_cumulative_option.png` | **Paper stage-timing figure (decided 2026-11-11).** Cumulative percent of fragments that had reached each stage by each day; fragments that never reached a stage stay in the denominator. Also the main display in the team summary. |
| `36_stage_timing.png` | Supporting exact-day healing/regeneration stage plot with stacked points and X marks for stages not reached. |
| `19e_source_heat_penalty_summary.png` | Plain-language summary of source-patch heat penalties across response types. |
| `19f_source_physiology_heat_penalties.png` | Easier view of the whole-coral physiology and growth heat penalties. |
| `19_genet_dashboard.png` | Legacy-name detailed source-patch heat-penalty plot across measurements. |
| `32_prelim_snp_structure.png` | Preliminary RNA-seq SNP structure: C is coherent, A and D are mixed. |
| `26b_trinity_hauru_context.png` | Preliminary Trinity/Ross Hauru coordinate context; candidate matches only, not confirmed identities. |
| `26_thermal_context.png` | Shows that 31 C is below the acute CBASS ED50 range. |
| `35_tradeoff_summary.png` | Exploratory growth, heat-tolerance, and regeneration trade-off screen. |

## Current Pipeline Figures

| File stem | What it shows | Main script | Current use |
|---|---|---|---|
| `02_pam_fvfm_trajectory` | Photosynthesis-efficiency score over time | `code/analysis/02_photosynthesis.R` | Physiology result |
| `03_color_trajectory` | Color-card paling over time | `code/analysis/03_color.R` | Physiology result |
| `04_morphology_trajectories` | Eight wound-healing and regrowth traits over time | `code/analysis/04_regeneration_scores.R` | Supporting morphology |
| `04b_morphology_trajectories_by_genet` | Same traits by legacy source-patch/genet label | `code/analysis/04_regeneration_scores.R` | Supporting source-patch view |
| `05_buoyant_weight_growth` | Skeletal mass change | `code/analysis/05_growth.R` | Growth result |
| `06_symbiont_density_by_day` | Algal symbiont density by day | `code/analysis/06_symbionts.R` | Symbiont result |
| `07_wax_standard_curve` | Wax-dipping surface-area calibration | `code/analysis/07_surface_area.R` | Calibration |
| `08_apex_temperature` | Main Apex tank-temperature record | `code/analysis/08_temperature.R` | Treatment validation |
| `08b_apex_temperature_full` | Full Apex temperature record | `code/analysis/08_temperature.R` | Supporting treatment validation |
| `09_ysi_water_chem` | Daily YSI water chemistry checks | `code/analysis/09_water_quality.R` | Husbandry context |
| `10_worm_presence` | Flatworm presence by treatment/tank | `code/checks/10_flatworms.R` | Contamination context |
| `11_microscope_trait_trajectories` | Microscope-photo healing and regrowth trajectories | `code/analysis/11_microscope_experiment.R` | Supporting photo cohort |
| `11b_microscope_trait_trajectories_by_genet` | Microscope-photo trajectories by source patch | `code/analysis/11_microscope_experiment.R` | Supporting photo cohort |
| `11c_morphology_dataset_all_trait_trajectories` | Main and microscope morphology datasets kept separate | `code/analysis/11c_stage_counts.R` | Cross-dataset check |
| `11d_morphology_dataset_shared_trait_comparison` | Shared morphology traits across datasets | `code/analysis/11c_stage_counts.R` | Cross-dataset check |
| `11e_morphology_dataset_event_summary` | Morphology event summary across datasets | `code/analysis/11c_stage_counts.R` | Cross-dataset check |
| `13_genet_response_panel` | Legacy-name response panel by source-patch label | `code/analysis/13_source_patch_tests.R` | Supporting source-patch view |
| `14_morphology_KM` | Kaplan-Meier time-to-onset curves | `code/analysis/14_stage_timing_models.R` | Main morphology result |
| `14b_morphology_KM_by_genet` | Kaplan-Meier curves by source-patch label | `code/analysis/14_stage_timing_models.R` | Supporting source-patch view |
| `15_physio_PCA_biplot` | End-of-study physiology PCA | `code/exploratory/15_physiology_pca.R` | Multivariate condition |
| `15b_physio_PCA_by_genet` | Physiology PCA by source-patch label | `code/exploratory/15_physiology_pca.R` | Supporting source-patch view |
| `16_manuscript_fig1` | Standalone four-panel overview | `code/archive/16_main_figure_draft.R` | Overview figure |
| `19_genet_dashboard` | Legacy-name detailed heat-penalty dashboard | `code/exploratory/19_source_patch_dashboard.R` | Detailed source-patch summary |
| `19b_genet_resilience_ranking` | Composite heat-penalty ranking | `code/exploratory/19_source_patch_dashboard.R` | Supporting source-patch summary |
| `19c_decomposed_resilience` | Heat penalties split by wounded/unwounded scope | `code/exploratory/19_source_patch_dashboard.R` | Supporting source-patch summary |
| `19d_wound_healing_heat_penalties` | Healing and regrowth heat penalties | `code/exploratory/19_source_patch_dashboard.R` | Wound-specific source-patch summary |
| `19e_source_heat_penalty_summary` | Plain-language heat-penalty summary | `code/exploratory/19_source_patch_dashboard.R` | Team-summary Section 4 |
| `19f_source_physiology_heat_penalties` | Whole-coral heat penalties by response | `code/exploratory/19_source_patch_dashboard.R` | Team-summary Section 4 |
| `26_thermal_context` | LTH temperatures versus acute CBASS ED50 | `code/exploratory/26_thermal_context.R` | Thermal context |
| `26b_trinity_hauru_context` | Preliminary Trinity/Ross Hauru coordinate screen | `code/exploratory/26b_trinity_hauru_context.R` | Candidate-match context |
| `32_prelim_snp_structure` | Preliminary SNP clusters and genetic PCs | `code/rnaseq/32_prelim_snp_integration.R` | Preliminary genetics |
| `32_prelim_snp_design_balance` | Preliminary SNP clusters by treatment design | `code/rnaseq/32_prelim_snp_integration.R` | Preliminary genetics |
| `32_prelim_snp_symbiont_heat_effects` | Direct SNP-by-symbiont heat-penalty check | `code/rnaseq/32_prelim_snp_integration.R` | Preliminary genetics |
| `33_growth_heat_tradeoff` | Growth versus heat-penalty screen | `code/33_growth_heat_tradeoff.R` | Exploratory trade-off |
| `34_regeneration_tradeoff_screens` | Regeneration trade-off screens | `code/34_regeneration_tradeoff_screens.R` | Exploratory trade-off |
| `35_tradeoff_summary` | Four-panel trade-off summary | `code/exploratory/35_tradeoff_summary.R` | Team-summary trade-off check |
| `36_stage_timing` | Exact-day healing/regeneration stage timing, with unreached late stages shown explicitly | `code/analysis/36_stage_timing_tables.R` | Supporting output for Molly follow-up |
| `36b_stage_violin` | Violin-style healing/regeneration stage timing, with individual points, non-reached stages, and reached/scored labels | `code/analysis/36_stage_timing_tables.R` | Team-summary Section 3 and Molly follow-up |
| `39b_stage_timing_cumulative_option` | **Paper stage-timing figure (2026-11-11).** Cumulative percent reached by day, never-reached fragments kept in the denominator | `code/analysis/39_stage_timing_figures.R` | Figure 1B and team-summary Section 3 |
| `41_stage_progression` | Alternative stage-timing display. Per-fragment first-observation day, an explicit never-reached zone on the same axis, and the proportion reached with binomial intervals. Suppresses a median wherever <80% of a group reached the milestone. | `code/exploratory/41_stage_progression_figure.R` | Alternative to 39b; supplement candidate |

## Notes

- **Stage-timing figures:** use `39b_stage_timing_cumulative_option` (decided 2026-11-11);
  `41_stage_progression` is the alternative. The other options are kept
  as the design comparison they were commissioned as, but `40_stage_mean_review`
  plots a mean first-observation day among fragments that *reached* the milestone,
  which renders the headline result as a 0.9-day gap (12.6 vs 13.5 d) when the
  finding is 12/12 versus 4/12; and `39a_stage_timing_violin_option` draws kernel
  densities from at most 12 reachers. Do not put either in the manuscript.
- Most figure stems have both `.png` and `.pdf` outputs.
- `figures/12_diagnostics/` contains diagnostic plots for model checks and is
  not part of the reader-facing figure set.
- The Trinity/Ross Hauru figure is preliminary context only. It should not be used to assign
  ED50 values to LTH source patches until DNA-marker matching is complete.
# October 27 physiology refresh

`45_physiology_figure2.{pdf,png}` now uses the confirmed wound assignment for fragment 119. PAM remains a repeated-measurement trajectory; symbiont density and skeletal growth use side-by-side wound-group boxplots. Standalone panels are `45a_photosynthesis`, `45b_symbiont_boxplots`, and `45c_growth_boxplots`, each in PDF and PNG. Driver: `code/analysis/45_figure2_physiology.R`; legend: `manuscript/FIGURE_2_LEGEND.md`. These files are generated locally; availability on GitHub requires a separate push.
