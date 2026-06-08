# CLAUDE.md — LTH (project #17) repo guide

Heat × wound × gene expression in *Acropora pulchra*, Mahana/Tiahura, Mo'orea (2025).
Repo: `stier-lab/17.-LTH_expression_by_temperature_2025`. Local: `~/Stier-LTH-expression-by-temperature-2025/`.

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

- Run everything: `source("code/_run_all.R")` (scripts 01–29, ~2.5 min).
- **Single source of truth for every statistic:** `output/tables/20_master_results.csv` (+ `_paper_ready.csv`). Never hardcode numbers in prose — cite the table; it regenerates with the pipeline.
- Full narrative: `RESULTS.md`. Pending work: `NEXT_STEPS.md`. RNA-seq + genet-matching brief: `docs/for_shreya/`.
- Genet (thicket A/C/D) is a **fixed** effect; continuous LMMs use `lmerTest::lmer` (type-III Satterthwaite). Growth = **areal calcification** (mg CaCO₃ cm⁻² d⁻¹). Confirmatory a-priori tests (literature-grounded) are reported unadjusted; exploratory tests are BH-corrected (`code/28`).

## Regeneration-staging terminology (keep consistent with the wound-type manuscript)

Use the **biphasic** framework standard in the lab's wound-type paper. In prose (not data/column names, which stay fixed):

- **Tissue-healing phase** = re-epithelialization / **coenosarc coverage** that seals the wound (LTH traits: `hole_in_center`, `polyp_in_hole`, `wound_smoothed`). Say "tissue healing / coenosarc coverage", not just "wound closure".
- **Regeneration phase** = reappearance of polyps + **skeletal/calyx/corallite** regrowth (LTH traits: `tip_exist`, `tip_extension`, `new_corallites_on_tip`).
- Use **"wound bed"** for the active region; **"algal colonization of the wound bed"** (or "algal plug" when dense) rather than "algae on wound".
- The healing→regeneration interval is the **"healing-to-regeneration lag"** (`code/14`).
- Geometric caveat to state explicitly: LTH wounds are apical-**tip excisions** (regeneration scored as new skeleton at the branch tip), whereas the wound-type paper scores polyp reappearance in a **surface wound bed**. Don't imply identical assays.

## Conventions

- R ≥ 4.3 (dev 4.5.2); packages pinned in `renv.lock`; `code/00_setup.R` defines `theme_pub()`, Okabe-Ito palettes, paths.
- Cite the lab figure/writing standards in `~/.claude/CLAUDE.md`. Citation grounding: verify any new citation via NotebookLM/Zotero before it enters the manuscript.
- Related: wound-type repos (`stier-wound-type-*`) and the wound-healing model (`coral-wound-healing-model`) are conceptual siblings, not formal comparators in this analysis.
