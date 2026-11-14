# Which outputs are current

Updated 2026-11-14. Use this page to tell current results from historical files before reading figures or tables.

## Current

| What | Where | Notes |
|---|---|---|
| Every reported number | `output/paper/*_numbers.csv` | Single source of truth, written by the core scripts. The team summary reads its numbers from here. |
| Plain-language results | `docs/team_summary/LTH_results_summary.html` | Rebuilt 2026-11-14. |
| Manuscript figure map | `manuscript/MAIN_FIGURE_PLAN.md` | Each figure in Molly's working draft and supplement, with the file and script that make it. |
| Figure 1B | `figures/39b_stage_timing_cumulative_option.*` | The four regeneration milestones, main experiment. All six stages: `39d_stage_cumulative_all_stages`. |
| Figure 2 | `figures/45_physiology_figure2.*` | Legend: `manuscript/FIGURE_2_LEGEND.md`. |
| Figure artwork | `figures/*` | Rebuilt 2026-11-12 with no embedded titles, subtitles, captions or legends; the colour/shape key is in `figures/FIGURE_INDEX.md`. All outputs regenerated the same day under the package versions pinned in `renv.lock`; this moved test statistics by under 0.05% and changed the simulation-based DHARMa checks (see the diagnostics audit). |
| Genetics | `figures/43_final_snp_*`, `output/tables/43_final_snp_*` | Final host clusters from Rachael (2026-10-13). Symbiont read shares are rough. Since 2026-11-12 the RNA-seq covariate tables (`31_*`, `43_final_snp_rnaseq_covariates`) no longer carry the retired source-resilience scores. |
| Morphology scores | `data/raw/physio_morphology/data.csv`, `data/raw/microscope_physio/data.csv` | October 16 corrected exports. |
| Model checks | `docs/diagnostics/MODEL_DIAGNOSTICS_2026-11-11.md` | Model structure decisions and diagnostics for every reported model. |
| Figure catalog | `figures/FIGURE_INDEX.md` | |

## Historical: keep for provenance, don't cite

| What | Superseded by |
|---|---|
| Old scripts, figures, tables, the master results table and the daily yes/no models (deleted with `extras/` on 2026-11-13/14; tags `pre-prune-2026-11-13`, `pre-cleanup-2026-11-14`) | The 10 core scripts in `code/` (merged 2026-11-11) |
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

- **Growth temperature test (2026-11-10, Adrian):** report the model test (temperature, F1,6 = 6.46, p = 0.044, from `03_physiology.R`). The tank permutation (p = 0.086, `05_buoyant_weight_tank_test.csv`) goes in the supplement as a check. (Numbers updated 2026-11-13 after the growth correction below; they were F = 6.33, p = 0.046 and p = 0.057.)
- **Figure 1B (2026-11-12):** follows Molly's working draft: four regeneration milestones (axial polyp, axial corallite forms, axial corallite extends, radial corallites), main experiment only. Healing and wound-surface remodelling go in the all-stage version (`39d`) for the supplement.
- **Old AI-written docs removed (2026-11-12):** RESULTS.docx, the local manuscript draft and its source extracts, session logs in `docs/provenance/`, email drafts, the roadmap, the diagnostic inventory and `notes/archive/`. Facts worth keeping went to `docs/DATA_NOTES.md`; the rest is in git history. Molly's Google Doc is the manuscript draft.

- **Oxygen and artwork check (2026-11-13):** raw YSI readings independently reproduce treatment oxygen concentration and saturation summaries. The summary no longer infers a small oxygen contribution; temperature and oxygen effects cannot be separated by these data. Color/shape keys and genetic-cluster sample sizes are in captions, with titles, subtitles, captions and legends suppressed by the figure theme. Full pipeline and HTML rebuild completed; numerical output tables were unchanged.
- **Growth calculation corrected (2026-11-13):** `01_clean_data.R` now follows Molly's buoyant-weight sheet: seawater density from the glass reference object at each weighing, and the aragonite plug subtracted from Day 0 and Day 15 dry mass. Earlier code used the freshwater density curve as seawater density and kept the plug (~60% of starting mass), which shrank % growth about 2.6-fold. Growth is now 16.0% (28 °C) vs 10.2% (31 °C), 36% lower under heat; the code checks every fragment against the metadata column `percent_growth_bw`. Only growth outputs changed.
- **Final cleanup (2026-11-14):** deleted the rest of `extras/`, the Nov 11 audit scripts (`docs/diagnostics/scripts/`; they no longer match the current models), `manuscript/references.bib` (unused; the manuscript is Molly's Google Doc) and `.Renviron`. The 1.8 MB Apex duplicate list was replaced by `08_apex_duplicate_summary.csv`. All recoverable from tag `pre-cleanup-2026-11-14`.
- **Accessibility pass (2026-11-14):** added `docs/GLOSSARY.md`, `docs/HOW_TO_CHANGE.md` and `CITATION.cff`; the README now shows Figures 1B and 2 with plain-language captions and screen-reader descriptions; the results summary is served as a web page through GitHub Pages (from `docs/`). Deleted 49 stale local plots in `figures/diagnostics/` (untracked; written by the removed audit scripts). No code or outputs changed.
- **Literature PDFs removed from git history (2026-11-14):** the 101 PDFs (moved to Drive on 2026-11-13) were stripped from every past commit with `git filter-repo`, so the public history no longer carries copyrighted papers. All commit IDs and the recovery tags changed; files at every commit are otherwise identical. The repo's history shrank from 571 MB to 97 MB. A full backup of the old history is at `~/repo-backups/acropora-lth-before-history-rewrite-2026-11-14.bundle` on Adrian's laptop. Anyone with an older clone should re-clone.
