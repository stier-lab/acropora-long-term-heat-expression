# Document Index — LTH (Project #17)

> **What this is** — the master list of every document in this repository: what it is, when it was
> last updated, and how it connects to the README, the code, and the data. If you're new, start with
> [`START_HERE.md`](START_HERE.md); this index is the map of everything else.
>
> **Updated** 2026-06-12 · **Maintainer note:** when you add or rename a document, add a row here.

## Legend

- **Updated** = date of the last change in git (`git log -1 --format=%cs -- <file>`).
- **Type:**
  - 📘 *manual* — hand-maintained prose; edit directly.
  - ⚙️ *generated* — written by the pipeline; **do not hand-edit** (it is overwritten on the next
    `Rscript code/_run_all.R`). Its "timestamp" is "regenerated every run."
  - 🗄️ *provenance* — raw export from the Google Drive project folder; kept read-only for the record.

---

## 1. Start here — orientation (repo root)

| Document | Type | Purpose | Updated | Connects to |
|---|---|---|---|---|
| [`START_HERE.md`](START_HERE.md) | 📘 | 5-minute guided tour; the front door | 2026-06-12 | → `README.md`, `RESULTS.md`; the reading path |
| [`README.md`](README.md) | 📘 | Full reference: question, design, findings, repo map, how to reproduce | 2026-06-12 | → all code/data/docs; links this index |
| [`RESULTS.md`](RESULTS.md) | 📘 | Authoritative phenotype results narrative (all responses, genet effects, thermal context, caveats) | 2026-06-12 | ← `code/01–29`; every number → `output/tables/20_master_results.csv` |
| [`NEXT_STEPS.md`](NEXT_STEPS.md) | 📘 | Pending work: chl-a, RNA-seq, genet matching, validations | 2026-06-12 | → `code/06`, `docs/for_shreya/`, `data/raw/` |
| [`CLAUDE.md`](CLAUDE.md) | 📘 | Repo guide + **Authorship & scope** (who leads, do-not-overreach rules), Drive↔repo mapping, conventions | 2026-06-12 | → README, `docs/for_shreya/`, Drive |
| [`SESSION_SUMMARY.md`](SESSION_SUMMARY.md) | 🗄️ | Internal/historical build log (2026-05-23) — provenance only, not a usage guide | 2026-06-12 | → `START_HERE.md` to actually use the repo |

## 2. Analysis walkthrough (`docs/`)

| Document | Type | Purpose | Updated | Connects to |
|---|---|---|---|---|
| [`docs/ANALYSIS_SUMMARY.md`](docs/ANALYSIS_SUMMARY.md) | 📘 | Plain-language, figure-by-figure walkthrough of the phenotype analysis | 2026-06-12 | ← figures in `figures/`; numbers → `output/tables/20_master_results.csv` |

## 3. Lead-author resources (`docs/for_shreya/`)

The hand-off landing zone for the lead author, **S. Banerjee** (gene-expression analysis + the
Introduction/Discussion/Abstract). These are *suggestions/resources*, not a prescribed pipeline.

| Document | Type | Purpose | Updated | Connects to |
|---|---|---|---|---|
| [`docs/for_shreya/README.md`](docs/for_shreya/README.md) | 📘 | Collaboration brief — who owns what, what's reproducible, where things live | 2026-06-12 | → all of `docs/for_shreya/`, `RESULTS.md` |
| [`docs/for_shreya/analysis_proposal.md`](docs/for_shreya/analysis_proposal.md) | 📘 | RNA-seq *questions/goals* + the genet-matching goal (A/C/D ↔ Cunning genets) | 2026-06-12 | → `data/external/`, `output/tables/19_genet_resilience_summary.csv`, `notes/sequencing-plan-keck-LTH.md` |
| [`docs/for_shreya/gene_expression_integration_map.md`](docs/for_shreya/gene_expression_integration_map.md) | 📘 | Phenotype → expression hypotheses (bidirectional), mapped to the Day 1/3/10 design | 2026-06-12 | → `output/tables/31_rnaseq_phenotype_covariates.csv`, `RESULTS.md` |
| [`docs/for_shreya/optional_intro_discussion_notes.md`](docs/for_shreya/optional_intro_discussion_notes.md) | 📘 | Optional phenotype-side framing angles (in tension) + verified citation bank | 2026-06-12 | → `manuscript/references.bib`, `literature/CITATIONS_INDEX.md` |

## 4. Data documentation (`data/`)

| Document | Type | Purpose | Updated | Connects to |
|---|---|---|---|---|
| [`data/DATA_DICTIONARY.md`](data/DATA_DICTIONARY.md) | 📘 | Variable definitions across all raw streams (one place to decode columns) | 2026-06-12 | ← `data/raw/`, `data/metadata/*_codebook.csv`; used by `code/01` |
| [`data/external/README.md`](data/external/README.md) | 📘 | Provenance for the Cunning 2024 CBASS ED50 reference data | 2026-06-07 | ← `data/external/cunning2024_apulchra_ed50.csv`; used by `code/26` |
| `data/metadata/*_codebook.csv` (6) | 📘 | Per-stream column codebooks (metadata, PAM, color, morphology, buoyant weight, wax) | 2026-06-12 | ← `data/raw/`; first draft for BCO-DMO |

## 5. Figures (`figures/`)

| Document | Type | Purpose | Updated | Connects to |
|---|---|---|---|---|
| [`figures/FIGURE_INDEX.md`](figures/FIGURE_INDEX.md) | 📘 | Catalog of every figure (PDF+PNG) → which script makes it, what it shows | 2026-06-12 | ← `code/11,13–16,19,26`; figures in `figures/` |

## 6. Literature (`literature/`)

| Document | Type | Purpose | Updated | Connects to |
|---|---|---|---|---|
| [`literature/LITERATURE.md`](literature/LITERATURE.md) | 📘 | Overview of the *A. pulchra* literature relevant to LTH | 2026-06-12 | → `literature/pdfs/` (98 PDFs), `manuscript/references.bib` |
| [`literature/LIBRARY_MAP.md`](literature/LIBRARY_MAP.md) | 📘 | Each paper → its relevance → which model/project it connects to | 2026-06-12 | → `literature/pdfs/`, `library_map.csv` |
| [`literature/CITATIONS_INDEX.md`](literature/CITATIONS_INDEX.md) | 📘 | Citation → PDF index for the manuscript | 2026-06-12 | → `manuscript/references.bib`, `manuscript/Manuscript_LTH.md` |
| [`literature/CITATION_AUDIT.md`](literature/CITATION_AUDIT.md) | 📘 | Backward-snowball citation audit of the library | 2026-06-12 | → `literature/pdfs/` |
| [`literature/KNOWN_UNKNOWN_synthesis.md`](literature/KNOWN_UNKNOWN_synthesis.md) | 📘 | What the field knows/doesn't about *A. pulchra* wounding & thermal tolerance | 2026-06-12 | → `literature/pdfs/`; framing for the Intro |
| `literature/library_map.csv` | 📘 | Machine-readable version of `LIBRARY_MAP.md` | 2026-06-12 | ← `literature/pdfs/` |

## 7. Manuscript (`manuscript/`)

| Document | Type | Purpose | Updated | Connects to |
|---|---|---|---|---|
| [`manuscript/Manuscript_LTH.md`](manuscript/Manuscript_LTH.md) | 📘 | Working draft. Phenotype Methods + Results drafted; Intro/Discussion/Abstract are the lead author's | 2026-06-12 | ← `RESULTS.md`; checked by `code/30`; refs in `references.bib` |
| `manuscript/references.bib` | 📘 | Verified BibTeX references (DOIs) | 2026-06-12 | ← `literature/CITATIONS_INDEX.md` |

## 8. Field notes & plans (`notes/`) — raw Drive exports (provenance)

Read-only exports from the Google Drive project folder, kept for the record. Listed here so they're
discoverable; not part of the active hand-off doc set.

| Document | Type | Purpose | Updated |
|---|---|---|---|
| [`notes/project_README_from_drive.md`](notes/project_README_from_drive.md) | 🗄️ | Original project README from Drive | 2026-05-23 |
| [`notes/Progress_Notes.md`](notes/Progress_Notes.md) | 🗄️ | Running progress notes (lab meetings) | 2026-05-23 |
| [`notes/Experimental_Plan_Gene_Expression.md`](notes/Experimental_Plan_Gene_Expression.md) | 🗄️ | Gene-expression sampling plan | 2026-05-23 |
| [`notes/Experimental_Plan_Microscope_photographs.md`](notes/Experimental_Plan_Microscope_photographs.md) | 🗄️ | Microscope photo protocol | 2026-05-23 |
| [`notes/Field_Notes_A._pulchra_healing_notes.md`](notes/Field_Notes_A._pulchra_healing_notes.md) | 🗄️ | Daily healing field observations | 2026-05-23 |
| [`notes/sequencing-plan-keck-LTH.md`](notes/sequencing-plan-keck-LTH.md) | 🗄️ | Sequencing / plate-layout plan | 2026-05-23 |
| [`notes/LTH_Color_Card_Photos.md`](notes/LTH_Color_Card_Photos.md) | 🗄️ | Color-card photo index (NAS) | 2026-05-23 |
| [`notes/LTH_Microscope_Characterization_Photos.md`](notes/LTH_Microscope_Characterization_Photos.md) | 🗄️ | Microscope photo index (NAS) | 2026-05-23 |
| [`notes/growth_allometry.md`](notes/growth_allometry.md) | 📘 | Why areal calcification is the growth metric (justifies `code/05`) | 2026-06-04 |
| [`notes/QAQC_flagged_samples.md`](notes/QAQC_flagged_samples.md) | 📘 | Flagged samples/tanks (feeds the `code/22` sensitivity check) | 2026-06-04 |

## 9. Auto-generated diagnostic reports (`output/diagnostics/`) — ⚙️ do not hand-edit

Regenerated on every `Rscript code/_run_all.R`. Each is produced by a diagnostic script and reports
model assumptions/coverage. "Timestamp" = whenever the pipeline last ran.

| Document | Produced by | Purpose |
|---|---|---|
| `A_continuous_report.md` | `code/diagnostics/A_continuous_lmm.R` | Continuous-response LMM diagnostics (DHARMa) |
| `B_morphology_report.md` | `code/diagnostics/B_morphology_glmm.R` | Morphological GLMM diagnostics |
| `C_cox_report.md` | `code/diagnostics/C_cox_diagnostics.R` | Cox PH (Schoenfeld) diagnostics |
| `D_pca_lrt_report.md` | `code/diagnostics/D_pca_lrt.R` | PCA + genet LRT diagnostics |
| `E_design_alignment_report.md` | `code/diagnostics/E_design_alignment.R` | Design/coding alignment audit |
| `F_model_reproducibility_report.md` | `code/diagnostics/F_model_reproducibility.R` | Re-fit reproducibility check |
| `G_inventory_report.md` | `code/diagnostics/G_diagnostic_plots.R` | Diagnostic-plot inventory |
| `H_coverage_report.md` | `code/diagnostics/H_spreadsheet_coverage.R` | Master-table coverage check |
| `I_timeseries_report.md` | `code/23_timeseries_diagnostics.R` | Repeated-measures/time-series diagnostics |
| `J_headline_model_comparison.md` | `code/24_headline_model_comparison.R` | Is a richer primary model worth it? (no) |
| `K_model_coverage_report.md` | `code/25_model_diagnostic_coverage.R` | Every model has a diagnostic (34/34) |

## 10. Code & archive documentation

| Document | Type | Purpose | Updated | Connects to |
|---|---|---|---|---|
| [`code/archive/molly_original/README.md`](code/archive/molly_original/README.md) | 📘 | Provenance note for Molly's original exploratory scripts | 2026-06-04 | ← `code/archive/molly_original/*.R` |

> **Script headers:** every file in `code/` opens with a `Purpose / Input / Output` header block —
> that is the per-script documentation; `README.md` (Pipeline stages) and `START_HERE.md` (Read the
> code in this order) give the run order. The single source of truth for every statistic is
> `output/tables/20_master_results.csv`.
