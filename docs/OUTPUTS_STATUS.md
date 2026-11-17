# Which outputs are current

Updated 2026-11-15. Use this page to tell current results from historical files before reading figures or tables.

## Current

| What | Where | Notes |
|---|---|---|
| Every reported number | `output/paper/*_numbers.csv` | Single source of truth, written by the core scripts. The team summary reads its numbers from here. |
| Plain-language results | `docs/team_summary/LTH_results_summary.html` | Rebuilt 2026-11-14. |
| Manuscript figure map | `manuscript/MAIN_FIGURE_PLAN.md` | Each figure in Molly's working draft and supplement, with the file and script that make it. |
| Figure 1B | `figures/39b_stage_timing_cumulative_option.*` | The four regeneration milestones, main experiment. All six stages: `39d_stage_cumulative_all_stages`. |
| Figure 2 | `figures/45_physiology_figure2.*` | Legend: `manuscript/FIGURE_2_LEGEND.md`. |
| Figure artwork | `figures/*` | Restyled 2026-11-15 ("Information is Beautiful" look; see `figures/FIGURE_INDEX.md`). No embedded titles, subtitles or captions (since 2026-11-12); since 2026-11-17 every figure carries its own legend for colour, shape and line type (key also in `figures/FIGURE_INDEX.md`). All outputs regenerated the same day under the package versions pinned in `renv.lock`; this moved test statistics by under 0.05% and changed the simulation-based DHARMa checks (see the diagnostics audit). |
| Genetics | `figures/43_final_snp_*`, `output/tables/43_final_snp_*` | Final host clusters from Rachael (2026-10-13). Symbiont read shares are rough. Since 2026-11-12 the RNA-seq covariate tables (`31_*`, `43_final_snp_rnaseq_covariates`) no longer carry the retired source-resilience scores. |
| Morphology scores | `data/raw/physio_morphology/data.csv`, `data/raw/microscope_physio/data.csv` | October 16 corrected exports. |
| Model checks | `docs/diagnostics/MODEL_DIAGNOSTICS_2026-11-11.md` | Model structure decisions and diagnostics for every reported model. |
| Figure catalog | `figures/FIGURE_INDEX.md` | |

## Historical: keep for provenance, don't cite

| What | Superseded by |
|---|---|
| Old scripts, figures, tables, the master results table and the daily yes/no models (deleted with `extras/` on 2026-11-13/14; tags `pre-prune-2026-11-13`, `pre-cleanup-2026-11-14`). Unused raw-data copies and templates (archived 2026-11-15; tag `pre-cleanup-2026-11-15`) | The 10 core scripts in `code/` (merged 2026-11-11) |
| Old stage-timing displays `39a`, `39c`, `40`, `41` (deleted 2026-11-13) | `39b` (Figure 1B) and `39d`. The older options can misrepresent the 12/12 vs 4/12 budding result (commit 0da201d). |
| `data/raw/*/data_before_2026-10-16.csv` | The October 16 corrected `data.csv` in the same folder |

## Why files show as changed after running the pipeline

`code/run_all.R` rewrites every table, figure and processed `.rds` file, and all of those are tracked in git. Two things follow.

1. **Real changes.** When code or raw data change, most outputs genuinely change. Between 2026-10-13 and 2026-10-27, three rounds of updates (final SNPs, October 16 scoring corrections, the fragment 119 label and the YSI refresh) were run but only partly committed. That is why about 170 files showed as changed. They were committed on 2026-11-10.
2. **Noise.** Some files change on every run even when nothing real did:
   - PDFs embed the time they were written
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

- **Legends back in the artwork (2026-11-17, Adrian):** every figure, including Figures 1B and 2, now shows a legend at the bottom explaining colours, point shapes and line types. Shared labels `LAB_TEMP`, `LAB_WOUND`, `LAB_GENO` in `code/00_setup.R`; Figure 2 collects its panel legends into one key and is 15 mm taller. The results-summary interactive plots also show legends, and a growth-plot grouping bug (levels still said clipped/unclipped after the wording change) was fixed. Replaces the 2026-11-12 no-legend rule.
- **Growth temperature test (2026-11-10, Adrian):** report the model test (temperature, F1,6 = 6.46, p = 0.044, from `03_physiology.R`). The tank permutation (p = 0.086, `05_buoyant_weight_tank_test.csv`) goes in the supplement as a check. (Numbers updated 2026-11-13 after the growth correction below; they were F = 6.33, p = 0.046 and p = 0.057.)
- **Figure 1B (2026-11-12):** follows Molly's working draft: four regeneration milestones (axial polyp, axial corallite forms, axial corallite extends, radial corallites), main experiment only. Healing and wound-surface remodelling go in the all-stage version (`39d`) for the supplement.
- **Old AI-written docs removed (2026-11-12):** RESULTS.docx, the local manuscript draft and its source extracts, session logs in `docs/provenance/`, email drafts, the roadmap, the diagnostic inventory and `notes/archive/`. Facts worth keeping went to `docs/DATA_NOTES.md`; the rest is in git history. Molly's Google Doc is the manuscript draft.

- **Oxygen and artwork check (2026-11-13):** raw YSI readings independently reproduce treatment oxygen concentration and saturation summaries. The summary no longer infers a small oxygen contribution; temperature and oxygen effects cannot be separated by these data. Color/shape keys and genetic-cluster sample sizes are in captions, with titles, subtitles, captions and legends suppressed by the figure theme. Full pipeline and HTML rebuild completed; numerical output tables were unchanged.
- **Growth calculation corrected (2026-11-13):** `01_clean_data.R` now follows Molly's buoyant-weight sheet: seawater density from the glass reference object at each weighing, and the aragonite plug subtracted from Day 0 and Day 15 dry mass. Earlier code used the freshwater density curve as seawater density and kept the plug (~60% of starting mass), which shrank % growth about 2.6-fold. Growth is now 16.0% (28 °C) vs 10.2% (31 °C), 36% lower under heat; the code checks every fragment against the metadata column `percent_growth_bw`. Only growth outputs changed.
- **Final cleanup (2026-11-14):** deleted the rest of `extras/`, the Nov 11 audit scripts (`docs/diagnostics/scripts/`; they no longer match the current models), `manuscript/references.bib` (unused; the manuscript is Molly's Google Doc) and `.Renviron`. The 1.8 MB Apex duplicate list was replaced by `08_apex_duplicate_summary.csv`. All recoverable from tag `pre-cleanup-2026-11-14`.
- **Accessibility pass (2026-11-14):** added `docs/GLOSSARY.md`, `docs/HOW_TO_CHANGE.md` and `CITATION.cff`; the README now shows Figures 1B and 2 with plain-language captions and screen-reader descriptions; the results summary is served as a web page through GitHub Pages (from `docs/`). Deleted 49 stale local plots in `figures/diagnostics/` (untracked; written by the removed audit scripts). No code or outputs changed.
- **Results summary revised after a Craig-style pre-review (2026-11-14):** written in "we", with provenance removed. Terms are now consistent throughout:
  - "clipped / unclipped" for the wounding treatment;
  - healing, then regeneration stages 1–5, using the stage names in figures 39b and 39d;
  - "host cluster" for the genetic groups;
  - symbiont reference groups renamed S1–S4 (Rachael's A–D), so they don't collide with the patch letters.

  Facts added from the data and Molly's methods notes:
  - the heat ramp, 28–31 May, from the tank sheet and the Apex logs;
  - fragment size and how much was clipped;
  - symbiont biopsy and surface-area methods;
  - details of the microscope experiment and its 8/8 vs 2/8 result;
  - a flatworm cross-tab: among heated clipped fragments with no flatworm record, 3 of 9 formed radial corallites;
  - the Day 0 weighing order;
  - tank 3 water conditions.

  A "Methods details still to add" list holds what only Molly can supply. No reported numbers changed.

  Follow-up: `code/07_supplement.R` calls tank 3 "slow to reach temperature", but the Apex record from 28 May shows it at 28.1 °C throughout. (The figure 36/36b stage-label mismatch was fixed on 2026-11-15.)
- **Draft Methods and Results text revised after a second Craig-style pre-review (2026-11-15):** the eight draft boxes in the summary's Statistics drop-downs now give:
  - the modelled days (Days 0–14 for Fv/Fm and color);
  - the reasons for the main model choices, and R and package versions;
  - the growth caveat (permutation p and tank 3) in the Results text itself;
  - the clipped-vs-unclipped symbiont ratios;
  - the flatworm check and a Fisher's test for the microscope repeat (p = 0.007);
  - the oxygen drop against what warmer water alone would cause;
  - the multiple-test adjustment and adjusted p for patch C.

  "Similar" and "did not differ" now read "did not detectably differ". Items only Molly can answer are marked [VERIFY] and added to the "Methods details still to add" list. No reported numbers changed, apart from rounding growth to one decimal and showing two decimals for the near-zero growth F values.
- **Summary: plain-language, biology-first pass (2026-11-15):**
  - removed the "How much confidence to place in the results" section;
  - the banner, sections 1–5 and the take-homes now lead with what the corals did, with the statistics after;
  - the oxygen text no longer says "mostly because warm water holds less oxygen", since solubility explains only 0.29 of the 0.63 mg/L drop;
  - stage 5 is described as "fewer, not later";
  - per-patch Fv/Fm, color and symbiont losses are now given as numbers in the main text.

  No reported numbers changed.
- **Interactive figures in the summary (2026-11-15, Adrian):** nine summary figures are now interactive (plotly; hover shows the tank or fragment, its patch and the exact value; drag to zoom):
  - tank temperature and water conditions;
  - growth and symbiont density by group;
  - the six stage curves and each fragment's first day per stage;
  - the SNP map, symbionts within host clusters and symbiont mix.

  They are built in the summary's `interactive-figures` chunk from the same data, colours and theme. The static files in `figures/` stay the manuscript versions. Figures 1B, 2 and the photo panels stay static. The page grew from 7.4 to 9.4 MB.
- **Wounded / unwounded throughout (2026-11-15, Adrian):** figures, the results summary, README, website, glossary and Figure 2 legend now call the treatment wounded / unwounded, matching the manuscript draft. This replaces the clipped / unclipped wording from earlier the same day. "Clipping" now only describes how the wound was made (about 1 cm clipped off the tip). Axis labels read "Day relative to wounding" and "Days after wounding". Only figure labels and the panel names in `47_other_daily_observations.csv` changed; no data or statistics changed.
- **README and website front page (2026-11-15, Adrian):** rewrote `README.md` around the question, findings and design, using clipped/unclipped throughout. Added `docs/index.html` as the website's front page, in the summary's style, linking to the summary, the code, the glossary and the how-to. The summary's "Methods details still to add" list is now in plainer language, split into questions for Molly and items for us to write.
- **Larger text on the summary's figures (2026-11-15, Adrian):**
  - `save_fig()` now also writes a web copy of every figure to `figures/web/`. It's the same plot drawn at 75% of the print size with the same point sizes, so the text is about a third larger once the page shrinks it to the text column. Figure 2 uses 85% and the patch panel 90%, so their labels and titles fit.
  - The summary uses these copies; the print files in `figures/` stay the manuscript versions.
  - The interactive plots' base text size went from 9–10 to 13.
  - Figure 2 panel A no longer labels Day −1, because its label ran into Day 0. The Day −1 points are still plotted.

  No data or statistics changed.
- **Summary page redesigned (2026-11-15, Adrian):**
  - the opening is now a lead paragraph plus a 28 °C vs 31 °C table of six key measures, all read from `output/paper/`;
  - type is Literata for text and Instrument Sans for headings, captions and navigation, with text kept to about 68 characters per line;
  - the only accent colours are the figures' two temperature colours;
  - figures lose their card frames, and the Statistics drop-downs are quieter;
  - the figure key is a list, and the layout works on phones.

  No text or numbers changed apart from the figure key's layout.
- **Figures restyled (2026-11-15, Adrian):** every figure now uses an "Information is Beautiful" look, set in `theme_pub()`, with black x and y axis lines on every panel. Labels are unified across figures:
  - 28 °C / 31 °C;
  - clipped / unclipped;
  - "Days after clipping" and "Sampling day";
  - stage names matching `39b` and `39d`;
  - "Cluster N" panel labels.

  In figure 06, the dots now sit over their own temperature's box. No data or statistics changed; only figures and the panel names in `47_other_daily_observations.csv` changed. Manuscript Figures 1B and 2 now say "clipped" rather than "wounded"; if the manuscript keeps "wounded", that is a one-line label change in `03_physiology.R` and `04_regeneration.R`.
- **Hochberg replaces Holm (2026-11-15, Adrian):** every multiple-test adjustment now uses the Hochberg method. This covers the six per-patch heat contrasts per measure (`03_physiology.R`, `12_genet_treatment_effects.csv`) and the per-cluster symbiont contrasts (`06_genetics.R`, `43_final_snp_symbiont_contrasts.csv`, `43_final_snp_leave_one_tank_out.csv`). The column `p_holm` is now `p_hochberg`. Hochberg is valid when the tests are independent or positively related, as contrasts from one model usually are, and is slightly less strict. Only adjusted p-values changed, and only slightly (the largest drop is 0.187 to 0.141 for growth, patch A and D unclipped); the patch C values quoted in the summary are unchanged, and no conclusion changed.
- **Unused files archived (2026-11-15, Adrian):** removed files nothing reads: the blank printable field sheets, the example chlorophyll workbook, the old buoyant-weight working copy, two CSV copies of buoyant-weight workbook tabs, two symbiont-count spreadsheet intermediates, the Nov 11 stats call scan, the never-filled daily health log, a duplicate wax standard curve and the out-of-date `data/DATA_DICTIONARY.docx` (June 12; the codebooks in `data/metadata/` and `docs/GLOSSARY.md` replace it). Their codebooks and `source_inventory.csv` rows were updated (inventory rows marked `archived`). The originals stay on Drive; everything is in tag `pre-cleanup-2026-11-15`. No code or outputs changed.
- **Literature PDFs removed from git history (2026-11-14):** the 101 PDFs (moved to Drive on 2026-11-13) were stripped from every past commit with `git filter-repo`, so the public history no longer carries copyrighted papers. All commit IDs and the recovery tags changed; files at every commit are otherwise identical. The repo's history shrank from 571 MB to 97 MB. A full backup of the old history is at `~/repo-backups/acropora-lth-before-history-rewrite-2026-11-14.bundle` on Adrian's laptop. Anyone with an older clone should re-clone.
