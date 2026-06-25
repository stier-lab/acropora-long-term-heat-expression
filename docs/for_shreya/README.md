# For Shreya — collaboration brief

> 🗂️ **Lead-author landing page** · Updated 2026-06-12 · Index: [`README.md`](../../README.md) · front door to everything in `docs/for_shreya/`.

Welcome, and thanks for leading this one. This folder is your front door to everything the Stier lab
has built on the LTH project (#17: heat × wound × *A. pulchra*, Mo'orea 2025). The short version:

> **You lead the paper and the gene-expression analysis. We've finished and handed over the
> phenotype half — the experiment's Methods and Results — as a clean, reproducible, drop-in
> foundation, with every number traceable. The narrative (Introduction, Discussion, Abstract) is
> yours.** Nothing here is meant to box you in.

## Who owns what

| | Lead | Notes |
|---|---|---|
| **Gene-expression analysis** (DE, WGCNA, SNP/genet-matching) | **Shreya** | the paper's lead result; open *questions/goals* (not a prescribed pipeline) in `analysis_proposal.md` |
| **Introduction / Discussion / Abstract** (the narrative) | **Shreya** | manuscript has placeholders; *optional* phenotype-side framing angles (in tension) + a citation bank in `optional_intro_discussion_notes.md` — material to mine, not a recommended spine |
| **Phenotype Methods + Results** (physiology, morphology, growth, genet variation, thermal context) | Stier lab (done) | the organismal context for your mechanism — written, number-rich, reproducible; drop in and trust the numbers |
| **Data / code / figures / stats** | Stier lab (done) | one-command reproducible; see below |
| **Author order** | **settled** | you lead and are corresponding author; R. Bay + A. Stier are co-senior. (Some of your correspondence fields are still `[to confirm]` in the masthead.) |
| Target journal, review process | you + Adrian | still open — your call to drive |

## What's done and reproducible (so you don't have to reconstruct it)

- **The whole phenotype pipeline runs with one command** and regenerates every figure and table:
  `Rscript code/_run_all.R` (~4 min). Start at **`README.md`** (repo root) for the tour.
- **Every statistic lives in one table:** `output/tables/20_master_results.csv` (description +
  effect size + test stat + df + p + CI per row; formatted version `_paper_ready.csv`). Never
  hand-copy a number — cite the table.
- **The phenotype numbers won't silently drift out of sync:** `code/30_manuscript_audit.R` (the last
  pipeline step) recomputes every phenotype number and **warns** (advisory only — it does *not* fail
  the run) if the manuscript no longer matches. It is scoped to the Stier-lab **phenotype
  Methods/Results** — it does not police your Introduction, Discussion, Abstract, or the
  transcriptomics. So when an analysis is edited, re-running tells you what to update, but nothing
  blocks you. (Currently 15/15 phenotype checks pass.)
- The Results section of `manuscript/Manuscript_LTH.md` is **already written** for the phenotype side
  — photochemistry, pigmentation, calcification, symbionts, the healing-vs-regeneration morphology,
  genet variation, and the Cunning thermal-context placement.

## Where everything lives

| Need | Path |
|---|---|
| Orientation / how to run | `README.md` (repo root) |
| Visual first-read walkthrough of the analysis | `README.md` (Findings section) |
| Full phenotype results narrative | `RESULTS.md` (incl. §10 limitations/caveats) |
| Every statistic (source of truth) | `output/tables/20_master_results.csv` |
| Your RNA-seq questions/goals + the genet-matching goal | `docs/for_shreya/analysis_proposal.md` |
| **Phenotype → expression hypotheses (suggested)** | `docs/for_shreya/gene_expression_integration_map.md` |
| **Optional Intro/Discussion draft (take or leave)** | `docs/for_shreya/optional_intro_discussion_notes.md` |
| Verified references + cite→PDF index | `manuscript/references.bib`, `literature/CITATIONS_INDEX.md` |
| Literature library (98 PDFs, mapped) | `literature/LIBRARY_MAP.md`, `literature/LITERATURE.md` |
| What the field knows/doesn't (A. pulchra) | `literature/KNOWN_UNKNOWN_synthesis.md` |
| RNA-seq plate layout / sequencing plan | `notes/sequencing-plan-keck-LTH.md`, `code/plate_fig.R` |

## How the two halves connect

The phenotype results set up specific, testable questions for your RNA-seq — the healing→regeneration
transition, the genet-c resilience axis, the local wound response, the chronic-vs-acute signature,
and the SNP-based genet-matching to Cunning's CBASS genets. These are laid out (as suggestions, your
call) in **`gene_expression_integration_map.md`**, mapped onto your actual Day 1 / 3 / 10 margin
design. Each link runs both ways — the expression data can confirm, extend, **or revise** the
phenotype-side hypothesis. The manuscript Results marks the integration points with ⟶ tags.

## A few honest things worth knowing

- **The regeneration result has a soft spot:** it is strongest for interval-censored new-corallite
  onset (p = 1.4e-7); tip-exist is also delayed, while tip-extension remains n.s. The
  censored fraction (67% closed-but-never-regenerated vs 0%) and the per-coral lag are the strongest
  framing. See `RESULTS.md` §10.
- **Chlorophyll-a was not run**; the metadata slot remains for provenance, but the analysis uses
  PAM, color-card scores, and symbiont counts instead.
- **Three genets** means we can detect genotype variation but not its genetic architecture — your
  SNP work is what would resolve that and link it to Cunning's ED50s.

## What we can still do for you

A per-library phenotype covariate table is **already done** —
`output/tables/31_rnaseq_phenotype_covariates.csv` (one row per RNA-seq library, joined by
`Fragment_ID`; design factors + symbiont density + per-genet resilience covariates). Built by
`code/31_rnaseq_covariate_table.R`. Beyond that, anything that lowers friction — just ask.
