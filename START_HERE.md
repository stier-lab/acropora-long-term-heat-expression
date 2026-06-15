# START HERE 👋

> 🗂️ **Front door / onboarding tour** · Updated 2026-06-12 · Full document list: [`DOCS_INDEX.md`](DOCS_INDEX.md) · next stops: `README.md`, `RESULTS.md`.

**New to this repo? Read this page first — it's a 5-minute guided tour.**
For the full scientific detail, `README.md` is the next stop.

---

## What this project is

A heat × wounding experiment on the branching coral *Acropora pulchra* (Mo'orea, French Polynesia,
2025). We clipped the growing tip off half the fragments, held corals at **28 °C (ambient)** or
**31 °C (heated, +3 °C)**, and tracked how they healed and regrew over 15 days — alongside
photochemistry, pigmentation, symbionts, and growth. Tissue samples were also taken for gene
expression (RNA-seq, still in progress).

**The paper's headline is the transcriptomic mechanism** (S. Banerjee, lead author; RNA-seq analysis
pending). The phenotype experiment supplies the organismal context for it. **The phenotype result:**
heat does **not** slow wound recovery uniformly — corals seal the wound (tissue healing) at the same
rate whether hot or not, but heated corals fail to rebuild skeleton at the tip (regeneration). Heat
hits one phase of recovery, not the other.

## Who's who

| Person | Owns |
|---|---|
| **Shreya Banerjee** (UC Davis, Bay lab) | **Leads the paper** — gene-expression analysis (RNA-seq DE + matching our 3 thickets to Cunning's genotyped genets) and the narrative (Introduction, Discussion, Abstract) → see **`docs/for_shreya/`** |
| **Molly Brzezinski** (UCSB) | Fieldwork, the phenotype data (PAM, color, morphology, growth, symbionts), the original exploratory scripts (`code/archive/molly_original/`) |
| **Stier lab** (UCSB) | Contributed the phenotype Methods + Results and the reproducible pipeline (`code/`) as a drop-in foundation |
| **Adrian Stier** (UCSB, PI) | Project PI; co-author |

## The one mental model you need

```
 Google Drive  ── source of truth ──►  THIS REPO  ── reproducible analysis ──►  results + manuscript
 (raw Sheets,                          (data/raw/ is exported FROM Drive;
  field notes,                          the repo NEVER writes back to Drive)
  manuscript Doc)
```

This repo is the **analysis layer**. Raw data is exported *from* Drive into `data/raw/`. When Drive
data changes, re-export and re-run the pipeline. Details + the Drive folder ID are in `README.md`
and `PROJECT_GUIDE.md`.

## Run it in 3 commands

```bash
cd ~/Stier-LTH-expression-by-temperature-2025
Rscript -e 'renv::restore()'     # one-time: install the exact package versions (renv.lock)
Rscript code/_run_all.R          # runs the full pipeline in order (~3 min); regenerates every figure + table
```

Everything is reproducible from `code/_run_all.R`. **Never hand-edit numbers** — every statistic
lives in `output/tables/20_master_results.csv` and regenerates with the pipeline.

## Read the docs in this order

*(For the complete catalog of every document — with timestamps and how each connects to the code and
data — see **[`DOCS_INDEX.md`](DOCS_INDEX.md)**. This list is just the main reading path.)*

1. **`START_HERE.md`** ← you are here
2. **`docs/ANALYSIS_SUMMARY.md`** — a visual, first-time-reader walkthrough of the analysis (what was measured, how, and what it found)
3. **`README.md`** — research question, design table, key findings, repo map
4. **`RESULTS.md`** — the full results narrative (all responses, genet effects, thermal context, caveats)
5. **`NEXT_STEPS.md`** — what's still pending (chl-a values, RNA-seq, genet matching)
6. **`literature/KNOWN_UNKNOWN_synthesis.md`** — what the field knows/doesn't about *A. pulchra*
   wounding & thermal tolerance, and where this study is the first to answer a question
7. **`manuscript/Manuscript_LTH.md`** — the working paper draft

## Read the code in this order

Each script starts with a header block stating its **Purpose / Input / Output** — skim those headers
to navigate. The numbering is the run order:

| Scripts | What they do |
|---|---|
| `00_setup.R` | Packages, file paths, `theme_pub()`, color palettes — sourced by every other script |
| `01` | Load + clean the master metadata (one row per coral) |
| `02`–`10` | One response variable each: PAM, color, morphology, buoyant weight, symbionts, wax SA, Apex temp, YSI, worms |
| `11`, `16` | Publication multi-panel figures |
| `12`–`15` | Core statistics: mixed models, genet interactions, survival (Kaplan–Meier/Cox), multivariate PCA |
| `12b`, `12c` | Robustness refits (ordinal color model; penalized morphology GLMMs) |
| `17`–`19` | Figure audit, data validation, genet resilience dashboard |
| `20` | **Master results table** — aggregates every statistic into one CSV |
| `21` | RNA-seq stub (placeholder — not run until count data lands) |
| `22`–`29` | Sensitivity & diagnostic suite (flagged samples, time-series, model coverage, thermal context, variance partitioning, multiple-testing, probability-scale contrasts) |
| `diagnostics/A–H` | Model-diagnostic suites (DHARMa residuals, design alignment, reproducibility) |

## Where everything lives

| Folder | Contents |
|---|---|
| `code/` | All analysis scripts (run order = file number) + `archive/molly_original/` (preserved) |
| `data/raw/` | Exported from Drive — never edited by hand |
| `data/processed/` | Cleaned `.rds` files the pipeline produces |
| `data/external/` | Cunning et al. 2024 CBASS ED50 reference data (+ provenance README) |
| `output/tables/` | Every result as CSV (`20_master_results.csv` is the single source of truth) |
| `output/models/` | Saved fitted models (`.rds`) |
| `output/diagnostics/` | Model-diagnostic reports |
| `figures/` | All figures (`.pdf` + `.png`) |
| `literature/` | PDF library + `LIBRARY_MAP.md`, `LITERATURE.md`, `KNOWN_UNKNOWN_synthesis.md`, `CITATION_AUDIT.md` |
| `manuscript/` | The working paper |
| `notes/` | Field notes, experimental plans, QA/QC flags, progress notes (exported from Drive) |
| `docs/for_shreya/` | RNA-seq + genet-matching brief |

## Status at a glance

- ✅ **Phenotype analysis: complete and reproducible.** Physiology, morphology, growth, genet
  variation, thermal context — all in `RESULTS.md` and `output/tables/20_master_results.csv`.
- ✅ **Manuscript: phenotype Methods + Results drafted** (Stier lab). Introduction, Discussion, and
  Abstract are Shreya's to write — the manuscript has placeholders, with *optional* draft notes in
  `docs/for_shreya/optional_intro_discussion_notes.md` (take or leave).
- ⏳ **Pending:** chlorophyll-a values (assay), and the **RNA-seq** (144 libraries at the Bay lab) →
  see `NEXT_STEPS.md` and `docs/for_shreya/`.

## Glossary — the recovery terminology (used consistently throughout)

The lab uses a **biphasic** framework (see Methods **Table 1** in the manuscript):

- **Tissue-healing phase** = re-epithelialization / **coenosarc coverage** that seals the wound
  (traits: `hole_in_center`, `polyp_in_hole`, `wound_smoothed`).
- **Regeneration phase** = reappearance of polyps + **skeletal/calyx/corallite** regrowth
  (traits: `tip_exist`, `tip_extension`, `new_corallites_on_tip`).
- **Wound bed** = the active recovery region (here, the excised apical tip).
- **Healing-to-regeneration lag** = the gap between sealing the wound and rebuilding skeleton
  (`code/14`). At 31 °C, two-thirds of corals close the wound but never cross into regeneration.
- **Genet / thicket (A, C, D)** = the three parent colonies; modeled as a fixed effect (only 3).

## Getting help

- Questions on the analysis pipeline or a specific script → read its header, then `RESULTS.md`.
- Questions on a citation → `literature/LIBRARY_MAP.md` (every PDF mapped to its relevance).
- Anything about Drive ↔ repo, permits, or terminology conventions → `PROJECT_GUIDE.md`.
