# PROJECT_GUIDE.md — LTH (project #17) repo guide

> 🗂️ **Repo guide + authorship/scope rules** · Updated 2026-06-12 · Full document list: [`DOCS_INDEX.md`](DOCS_INDEX.md).

Heat × wound × gene expression in *Acropora pulchra*, Mahana/Tiahura, Mo'orea (2025).
Repo: `stier-lab/17.-LTH_expression_by_temperature_2025`. Local: `~/Stier-LTH-expression-by-temperature-2025/`.

## Authorship & scope — read this before editing prose or analysis docs

**Shreya Banerjee (UC Davis, Bay lab) leads this paper.** She owns the gene-expression analysis
(RNA-seq DE / WGCNA / SNP genet-matching) and the **Introduction, Discussion, and Abstract**, and is
the **corresponding author**. **R. A. Bay and A. C. Stier are co-senior authors.** Byline order:
Banerjee, Brzezinski, Diminuco, Seifert, Osenberg, Bay, Stier.

The Stier-lab contribution is the **phenotype experiment's Methods + Results and a reproducible
pipeline** — this is the *organismal context* for the paper, **not** its headline. The paper's headline
is the transcriptomic mechanism (hers, pending).

When working in this repo, **do not over-reach into the lead author's domain**:
- Don't write or pre-empt her Intro/Discussion/Abstract, and don't bake Discussion-level interpretation
  into the Results (effect sizes + direction only there).
- Don't prescribe her RNA-seq pipeline. The `docs/for_shreya/` files pose *questions/goals* and offer
  *suggestions* (no fixed DE model, no candidate-gene symbols, no step-by-step SNP recipe); phrase every
  phenotype↔expression link bidirectionally ("expression can test, extend, **or revise** this").
- "Energetic triage" and similar framings are candidate labels for her to choose, **not** house
  vocabulary.
- The manuscript reproducibility check (`code/30_manuscript_audit.R`) is **advisory** (warns, never
  fails the build) and scoped to the Stier-lab phenotype Methods/Results only — it does not police her
  narrative.

Still open for PI/institutional decision: target journal & review process;
Banerjee's correspondence department/address/email (still `[to confirm]` in the masthead); and
rehoming/renaming the repo to a neutral or co-owned home.

## Google Drive project folder (source of truth for raw data & docs)

This repo is the **analysis** layer; the **project of record** lives in Google Drive:

- **Folder:** `17. LTH_expression_by_temperature_2025`
- **Drive ID:** `1sXfnHN-vmSBuwMEfERiYOWDeKRmjFWJP`
- **Path:** `…/My Drive/Stier Lab/People/Adrian Stier/Projects/In Progress/Coral-Regeneration/Projects/17. LTH_expression_by_temperature_2025/`
- Access via the `gws` CLI (authenticated as astier@ucsb.edu), e.g.
  `gws drive files list --params '{"q":"'\''1sXfnHN-vmSBuwMEfERiYOWDeKRmjFWJP'\'' in parents and trashed=false","fields":"files(id,name,mimeType)"}'`

**Repo ↔ Drive mapping:**
- Drive `data/` (Google Sheets: metadata, PAM, color, morphology, buoyant weight, wax, symbionts, worms, APEX/YSI) → exported to `data/raw/` here.
- Drive `notes/` (Progress Notes, Field Notes, Experimental Plans, sequencing plan, plate layout) → exported to `notes/` here (markdown).
- Drive `Manuscript_LTH` (Google Doc) ↔ `manuscript/Manuscript_LTH.md`.
- Raw data is exported FROM Drive; this repo never writes back to Drive. When Drive data changes, re-export into `data/raw/` and re-run `code/_run_all.R`.

## How the analysis works

- Run everything: `Rscript code/_run_all.R` (the full pipeline in dependency order, ~3 min; `21_rnaseq_stub.R` is excluded — it waits on sequencing).
- **Single source of truth for every statistic:** `output/tables/20_master_results.csv` (+ `_paper_ready.csv`). Never hardcode numbers in prose — cite the table; it regenerates with the pipeline.
- Full narrative: `RESULTS.md`. Pending work: `NEXT_STEPS.md`. RNA-seq + genet-matching brief: `docs/for_shreya/`.
- Genet (thicket A/C/D) is a **fixed** effect; continuous LMMs use `lmerTest::lmer` (type-III Satterthwaite). Growth = **areal calcification** (mg CaCO₃ cm⁻² d⁻¹). Confirmatory a-priori tests (literature-grounded) are reported unadjusted; exploratory tests are BH-corrected (`code/28`).

## Regeneration-staging terminology (the lab's biphasic vocabulary — a starting convention, not a mandate)

The **biphasic** framework below is the lab's working vocabulary, broadly aligned with the wound-type
manuscript so the sibling projects are easy to read together. Use it as a sensible default for the
*staging* descriptors, but it is **not** binding on this paper: Shreya leads the narrative, and the
framing vocabulary (including interpretive labels) is hers to set. In particular, **"energetic triage"
is a candidate interpretive label, not house vocabulary** — it belongs in the Discussion only if the
lead author chooses it. In prose (not data/column names, which stay fixed):

- **Tissue-healing phase** = re-epithelialization / **coenosarc coverage** that seals the wound (LTH traits: `hole_in_center`, `polyp_in_hole`, `wound_smoothed`). Say "tissue healing / coenosarc coverage", not just "wound closure".
- **Regeneration phase** = reappearance of polyps + **skeletal/calyx/corallite** regrowth (LTH traits: `tip_exist`, `tip_extension`, `new_corallites_on_tip`).
- Use **"wound bed"** for the active region; **"algal colonization of the wound bed"** (or "algal plug" when dense) rather than "algae on wound".
- The healing→regeneration interval is the **"healing-to-regeneration lag"** (`code/14`).
- Geometric caveat to state explicitly: LTH wounds are apical-**tip excisions** (regeneration scored as new skeleton at the branch tip), whereas the wound-type paper scores polyp reappearance in a **surface wound bed**. Don't imply identical assays.

## Conventions

- R ≥ 4.3 (dev 4.5.2); packages pinned in `renv.lock`; `code/00_setup.R` defines `theme_pub()`, Okabe-Ito palettes, paths.
- Citation grounding: verify any new citation via NotebookLM/Zotero before it enters the manuscript.
- Related: wound-type repos (`stier-wound-type-*`) and the wound-healing model (`coral-wound-healing-model`) are conceptual siblings, not formal comparators in this analysis.
