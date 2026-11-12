# Manuscript Figure Map

Updated 2026-11-12. This page maps each figure in Molly's working draft and supplement (Google Docs `working_draft` and `Supplemental`) to the file and script that produce it. Every data figure is rebuilt by `Rscript code/run_all.R`. The September plan this replaces is in `docs/provenance/MAIN_FIGURE_PLAN_2026-09-16.md`; its numbers are out of date.

## Main figures

| Draft figure | What the draft asks for | File | Script | Status |
|---|---|---|---|---|
| 1A | Photo panel of the regeneration stages | none (photos) | none | Artwork in the Drive Figure 1 folder. |
| 1B | Daily assessment of four regeneration milestones over 15 days: axial polyp regeneration, axial corallite formation, axial corallite extension, radial corallite formation | `figures/39b_stage_timing_cumulative_option` | `04_regeneration.R`, section 4 | Done. Four panels, 12 wounded fragments per temperature, main experiment. Stage numbers match Supplementary fig. 3. |
| 2A | Fv/Fm over time | `figures/45_physiology_figure2` (panel A; alone as `45a_photosynthesis`) | `03_physiology.R` | Done. Legend: `FIGURE_2_LEGEND.md`. |
| 2B | Symbiont density | `45_physiology_figure2` (panel B; `45b_symbiont_boxplots`) | `03_physiology.R` | Done. |
| 2C | Growth | `45_physiology_figure2` (panel C; `45c_growth_boxplots`) | `03_physiology.R` | Done. |
| 3, 4 | Transcriptomics | none yet | not in the core pipeline | Waiting on the expression analysis. |

## Supplementary figures

| Draft figure | What the draft asks for | File | Script | Status |
|---|---|---|---|---|
| S1 | Phases schematic | none (drawing) | none | Artwork. |
| S2 | Other macroscopic observations: pigment and algae on the wound, polyp behavior, disease state | `figures/47_other_daily_observations`; counts in `output/tables/47_other_daily_observations.csv` | `07_supplement.R`, section 1b | Done. Disease state is the flatworm (AEF) record; every fragment was recorded alive at every check, so there is no mortality panel. |
| S3 | Microscope fragments, healing and the five stages | `figures/11_microscope_trait_trajectories` (data panel) | `07_supplement.R`, section 1 | Data panel done; the microscope photos are artwork. |
| S4 | Color card data | `figures/03_color_trajectory` | `03_physiology.R` | Done. |
| S4 (second) | Subsampling schematic | none (drawing) | none | Artwork. The draft numbers two figures S4; one needs renumbering. |
| not in draft | All six stages, including healing and stage 2 | `figures/39d_stage_cumulative_all_stages` | `04_regeneration.R` | Optional. Useful if a reviewer asks where stage 2 went. |

## Corrections to the draft text and legends

These come from the current outputs. The numbers are in `output/paper/regeneration_numbers.csv` and `physiology_numbers.csv`.

1. **Tip extension on Day 8.** The draft says 100% of ambient fragments compared with none of the heated fragments. The data give 11 of 12 ambient (92%) and 7 of 12 heated (58%) (`tip_extension_pct_day8_*`).
2. **Radial corallites ("11ish and 12ish").** Ambient fragments first showed radial corallites on Days 12-15; the four heated fragments that budded did so on Days 12-14 (`budding_first_day_range_*`). By Day 12, 9 of 12 ambient and 1 of 12 heated fragments had budded. By Day 15, 12 of 12 and 4 of 12.
3. **"X% in heated."** Use the counts above. Write counts as well as percentages, since each is out of 12.
4. **Figure 2B units.** Symbiont density is 10^6 cells cm^-2 (per square centimeter of surface), not cm^-1.
5. **Figure 2 sample sizes.** Fv/Fm: the same 48 fragments at each visit, 24 per temperature and 12 per temperature-by-wound group. Symbionts: 48 distinct fragments per sampling day, 24 per temperature and 12 per group. Growth: 48 fragments, 12 per group. State the per-group n in the legend, since each plotted point or box is one group.
6. **Day 4 polyps.** All polyps were retracted on Day 4, the day every fragment was dipped to treat flatworms. Say so in the S2 legend, so readers don't read it as a heat effect.
