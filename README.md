# LTH: Long-Term Heating × Wounding × Gene Expression in *Acropora pulchra*

> 👋 **First time here? Read [`START_HERE.md`](START_HERE.md) first** — a 5-minute guided tour for
> new readers (Shreya, Molly, and anyone picking up the repo). This README is the fuller reference.

Project #17 of the Coral Regeneration program. May–July 2025, Gump Station, Mo'orea, French Polynesia.

## Research question

Branching corals like *Acropora pulchra* regenerate readily after physical damage to the apical (growing) tip, and elevated temperatures are known to slow growth and bleach symbionts. **How do those two stressors interact?** Specifically: does mild long-term heating (+3 °C above ambient) compromise wound regeneration, and what gene-expression machinery underlies any response?

## Key phenotype findings (organismal context; the transcriptomic mechanism is the lead result, pending)

> **Shreya Banerjee leads this paper.** Its headline is the **transcriptomic mechanism** (S. Banerjee,
> Bay lab; RNA-seq analysis pending). The phenotype findings below — the Stier-lab contribution — are
> the **organismal context** for that mechanism, not the paper's lead result.

1. **Heat impairs *regeneration*, not wound *closure*.** Wounds seal at the same rate in both temperatures, but new-corallite formation at the regenerating tip collapses under heat (Cox HR = 0.22, p = 0.010).
2. **Sustained 31 °C broadly compromises physiology** — photochemistry (Fv/Fm treatment×day F = 106.6), pigmentation (58–67 % paled vs 0–8 % ambient by D14), symbiont density, and **areal calcification (38 % reduction)**.
3. **Genotype-level thermal tolerance: C > D > A.** Genet C is 2.7–3.5× less heat-sensitive across photochemistry, pigmentation, and symbiont retention.
4. **The 31 °C stress is chronic-sublethal** — ~4.4 °C *below* this population's acute Fv/Fm ED50 (35.4 °C; Cunning et al. 2024, same Mahana site), so the responses are accumulated sublethal stress, not acute photoinhibition.

Full numbers and prose: **`RESULTS.md`**. Every statistic (p-values, effect sizes, test stats, df, CIs) lives in one regenerable table: **`output/tables/20_master_results.csv`** (+ `_paper_ready.csv`).

## Design

| Factor | Levels |
|---|---|
| **Temperature** | 28 °C (ambient) vs 31 °C (heated). Ramped 1 °C/day from ambient. |
| **Wound** | Wounded (~1 cm clipped off growing tip with band-saw + caliper) vs unwounded sham. Applied 7 days after target temperature was reached. |
| **Genotype (thicket)** | A, C, D — modeled as a **fixed effect** (only 3 field-collected genets; too few for a variance component). |
| **Tank** | 8 total: 4 per temperature (28 °C: 3, 6, 9, 12; 31 °C: 4, 5, 10, 11). Random effect. |
| **Time** | Daily non-destructive obs. Destructive biopsies at D0, D1, D3, D10, D15. |

**Coral fragments (n = 208 total):**
- **192** for gene expression + destructive physiology (chl-a, symbionts) — 3 genotypes × 4 tanks × 2 temperatures × 2 wound × ~4 timepoints
- **48** of these 192 also tracked for non-destructive metrics (PAM, color card, buoyant weight)
- **16** dedicated to daily microscope photography (genotypes A and C only)

**Collection site:** Mahana / Tiahura, NW Mo'orea — three parent thickets at A 17.49735 °S / 149.91557 °W, C 17.49808 °S / 149.91595 °W, D 17.49726 °S / 149.91581 °W (from `metadata.csv` `coord_lat`/`coord_long`). Same named site as Cunning et al. 2024's CBASS genets.

The full README from Drive lives at `notes/project_README_from_drive.md`. Detailed protocols: `notes/Experimental_Plan_Gene_Expression.md`, `notes/Experimental_Plan_Microscope_photographs.md`. Plate layout / sequencing plan: `notes/sequencing-plan-keck-LTH.md`.

## Google Drive project folder

This repo is the **analysis layer**; the **project of record** (raw Sheets, field notes, manuscript Doc, slides) lives in Google Drive:

- **Folder:** `17. LTH_expression_by_temperature_2025` · **Drive ID:** `1sXfnHN-vmSBuwMEfERiYOWDeKRmjFWJP`
- **Path:** `…/My Drive/Stier Lab/People/Adrian Stier/Projects/In Progress/Coral-Regeneration/Projects/17. LTH_expression_by_temperature_2025/`
- **Mapping:** Drive `data/` Sheets → `data/raw/` here; Drive `notes/` Docs → `notes/` here; Drive `Manuscript_LTH` Doc ↔ `manuscript/Manuscript_LTH.md`. Raw data is exported *from* Drive; the repo never writes back. When Drive data changes, re-export into `data/raw/` and re-run `code/_run_all.R`.
- Access programmatically with the `gws` CLI (authenticated as astier@ucsb.edu). See `CLAUDE.md` for the exact query and the repo↔Drive details.

Regeneration-staging terminology uses the lab's biphasic vocabulary (tissue healing / coenosarc coverage → regeneration / skeletal regrowth) as a working default — a convention for the phenotype descriptors, not a mandate on the paper's framing, which the lead author sets; see `CLAUDE.md`.

## Repository documents — where to find things

| Document | What it is |
|---|---|
| **`START_HERE.md`** | 5-minute onboarding tour — read this first |
| **`RESULTS.md`** | Full results narrative (all responses, genet effects, thermal context, limitations) |
| **`output/tables/20_master_results.csv`** | **Single source of truth** — every p-value, effect size, test stat, df, CI, with source script per row. Regenerated by the pipeline. `_paper_ready.csv` is the formatted version. |
| **`NEXT_STEPS.md`** | Pending work: chl-a values, RNA-seq integration, genet matching, validations for Adrian |
| **`SESSION_SUMMARY.md`** | Handoff log of the build-out |
| **`docs/for_shreya/`** | Lead-author resources for Shreya (Bay lab): RNA-seq *questions/goals* (not a prescribed pipeline) and the goal of matching thickets A/C/D to Cunning's genotyped genets — suggestions to adopt, change, or discard |
| **`literature/KNOWN_UNKNOWN_synthesis.md`** | What the field knows/doesn't about *A. pulchra* wounding, healing, regeneration & thermal tolerance — and where LTH is first to answer |
| **`literature/LIBRARY_MAP.md`** · **`LITERATURE.md`** | Every PDF mapped to its *A. pulchra* relevance; discovery index |
| **`manuscript/Manuscript_LTH.md`** | Working manuscript draft — Stier-lab phenotype Methods + Results drafted; Introduction, Discussion, and Abstract are the lead author's (S. Banerjee) to write |
| **`output/diagnostics/`** | Model-diagnostic reports (continuous, morphology, Cox, PCA, time-series, model-coverage) — all regenerated by the pipeline |
| **`code/archive/molly_original/`** | Molly's original exploratory R scripts, preserved unmodified for provenance |
| **`data/external/`** | Cunning et al. 2024 *A. pulchra* CBASS ED50 reference data (+ provenance README) |
| **`notes/`** | Drive exports: experimental plans, field notes, progress notes, sequencing plan, photo indices |

## Data inventory

| Stream | Source file | Rows | What it captures |
|---|---|---|---|
| Coral metadata (one per fragment) | `data/raw/metadata/metadata.csv` | 208 | thicket, id, tank, treatment, wound, biopsy day/date, **coord_lat/long**, calculated SA, chl-a, zoox |
| PAM (Fv/Fm) | `data/raw/pam/PAM_data.csv` | 672 obs | F, M, Y, E, Fv/Fm by date × tank × sample × location (top/bottom) |
| Color card (Siebeck D-scale) | `data/raw/color_card/data.csv` | 336 obs | health_status, color (D-scale, split scores averaged), paling, hole_at_center |
| Morphology (9 binary traits) | `data/raw/physio_morphology/data.csv` | 768 obs | polyps_out, hole_in_center, polyp_in_hole, wound_smoothed, pigment_over_wound, tip_exist, tip_extension, new_corallites_on_tip, algae_on_wound |
| Buoyant weight (growth) | `data/raw/buoyant_weight/data.csv` | 48 corals | initial + final coral & plug weights, water-density correction → areal calcification |
| Wax dipping (surface area) | `data/raw/wax_dipping/data.csv` | 192 corals | diameter, height, dry weight, wax weights, SA from standard curve |
| Standard curve for wax | `data/raw/wax_dipping/Standard_curve.csv` | 19 | wax mass vs SA from cylinders |
| Symbiont counts | `data/raw/symbiont_counts/Raw_counts.csv` | 768 (4 reps × 192) | hemocytometer counts → cells/cm² via wax SA |
| Worm presence | `data/raw/worm_presence/Sheet1.csv` | 192 corals × 3 dates | AEFW surveillance 06/07, 06/08, 06/12 |
| APEX temperature/pH | `data/raw/apex/datalog*.xml` | continuous | per-tank temperature and pH, May–June |
| YSI daily spot checks | `data/raw/ysi/Sheet1.csv` | 72 obs | TEMP, DO%, DO mg/L, SAL, pH daily |
| Sequencing plate layout | `data/raw/plate_layout/Plate_{1,2}.csv` | 96 wells × 2 | wells, sample ID, day, temp, wound, tank, genotype |
| Sample shipment manifest | `data/raw/shipping/*.csv` | varies | shipped to UC Davis Bay lab |
| **External: Cunning 2024 ED50** | `data/external/cunning2024_apulchra_ed50.csv` | 20 genets | acute CBASS Fv/Fm ED50 for *A. pulchra*, Mahana (thermal-tolerance benchmark) |

**Photographs** (color-card, microscope) live on the Stier Lab NAS (`smb://stier-nas1.eemb.ucsb.edu`); see `notes/LTH_*_Photos.md`. **RNA-seq reads** are processed at UC Davis (Bay lab); NCBI BioProject TBD.

## How to reproduce

### Prerequisites
- R ≥ 4.3.0 (developed on 4.5.2)
- Packages (see `code/00_setup.R`): tidyverse, lubridate, janitor, readxl, patchwork, scales, broom, broom.mixed, lme4, lmerTest, emmeans, DHARMa, MuMIn, car, survival, ordinal, blme, nlme, mgcv, ggrepel, here. `renv.lock` pins versions.

### Run the whole pipeline
```r
renv::restore()              # restore pinned packages
source("code/_run_all.R")    # runs the full pipeline in dependency order (~3 min)
```
Figures land in `figures/`, tables in `output/tables/`, diagnostics in `output/diagnostics/`, processed data in `data/processed/`, fitted models in `output/models/`. The final step rebuilds the master results table.

### Pipeline stages
`01` metadata → `02` PAM → `03` color → `04` morphology → `05` calcification → `06` symbionts/chl → `07` wax SA → `08` APEX temp → `09` YSI → `10` worms → `11` combined figure → `12` primary mixed models (+ `12b` color CLMM robustness, `12c` blme morphology) → `13` genet interaction → `14` Kaplan–Meier / Cox → `15` PCA → `16` manuscript figure → `17` figure audit → `18` data validation → `19` genet dashboard → `22` flagged-sample sensitivity → `23` time-series diagnostics → `24` headline-model comparison → `26` thermal context → `20` master results table → `25` model-diagnostic coverage. (`21_rnaseq_stub.R` waits on sequencing; not in the run-all.)

## Directory structure

```
.
├── code/
│   ├── _run_all.R               # runs the full pipeline in order
│   ├── 00_setup.R               # packages, theme_pub(), palettes, paths
│   ├── 01–11_*.R                # per-stream cleaning + first-pass analysis/figures
│   ├── 12_extended_stats.R      # primary mixed models (genet fixed; full interactions)
│   ├── 12b/12c_*.R              # color CLMM + penalized (blme) morphology robustness
│   ├── 13_genet_interaction.R   # genet × treatment LRTs + reaction norms
│   ├── 14_morphology_kaplan.R   # KM curves + Cox PH (+ PH diagnostics)
│   ├── 15_multivariate.R        # physiology PCA
│   ├── 16_main_figure.R         # publication Figure 1
│   ├── 17–18_*.R                # figure audit, data validation
│   ├── 19_genet_dashboard.R     # cross-response resilience ranking
│   ├── 20_master_results_table.R# aggregate every stat → one table
│   ├── 21_rnaseq_stub.R         # counts+metadata import example (awaiting data; DE design is the lead author's)
│   ├── 22_sensitivity_flagged.R # drop Molly's flagged samples → robustness
│   ├── 23_timeseries_diagnostics.R # autocorrelation, nonlinearity, random slopes
│   ├── 24_headline_model_comparison.R # is a richer model worth it? (no)
│   ├── 25_model_diagnostic_coverage.R # audit: every model has a diagnostic
│   ├── 26_thermal_context.R     # LTH treatments vs Cunning 2024 acute ED50
│   ├── archive/molly_original/  # Molly's original scripts (provenance)
│   ├── diagnostics/             # diagnostic-suite scripts (A–H)
│   └── functions/               # shared helpers
├── data/
│   ├── raw/                     # immutable — never edit
│   ├── processed/               # regenerable
│   ├── metadata/                # per-stream codebooks
│   └── external/                # Cunning 2024 ED50 reference (+ README)
├── docs/for_shreya/             # RNA-seq + genet-matching collaborator brief
├── figures/                     # publication figures (PDF + PNG) + diagnostics/
├── output/
│   ├── tables/                  # CSV summaries incl. 20_master_results.csv
│   ├── models/                  # saved model objects (.rds)
│   └── diagnostics/             # diagnostic reports (A–K)
├── manuscript/Manuscript_LTH.md # working draft
├── notes/                       # Drive exports + QA/QC + allometry notes
├── RESULTS.md  NEXT_STEPS.md  SESSION_SUMMARY.md
└── renv.lock  *.Rproj
```

## Statistical approach

Genet (thicket) is treated as a **fixed effect** throughout (only 3 field-collected genets — too few for a reliable variance component; Bolker 2008, Gelman 2005), which also surfaces per-genet effects directly.

- **Continuous responses** (PAM Fv/Fm, color D-scale, log symbiont density, areal calcification): linear mixed models `response ~ treatment * wound * day * thicket + (1|tank) + (1|id)` (`lme4`/`lmerTest`; `(1|id)` dropped for single-observation responses).
- **Growth = areal calcification** (mg CaCO₃ cm⁻² d⁻¹, surface-area-normalized via wax SA) is the primary metric; % mass change and SGR are retained as robustness (see `notes/growth_allometry.md`).
- **Morphology** (9 binary wound-healing traits): binomial GLMMs; 7 traits refit with Cauchy(0,2.5) priors (`blme`) for separation.
- **Color** D-scale also refit as an ordinal CLMM (`ordinal`) robustness check.
- **Healing milestones:** Kaplan–Meier + Cox PH (`survival`), thicket-stratified, with full PH (Schoenfeld) diagnostics.
- **Multivariate:** centered/scaled PCA on the four endpoint responses.
- **Diagnostics:** DHARMa for every model + dedicated suites for time-series (AR(1), nonlinearity, random slopes) and Cox PH; coverage audited at 34/34 models. All contrasts Tukey-adjusted (`emmeans`).
- **Reference levels:** treatment = `28C`, wound = `no`.

Every numerical claim traces to `output/tables/20_master_results.csv`; reproduce with `code/_run_all.R`.

## Authors

**S. Banerjee leads the paper** (RNA-seq analysis + the Introduction, Discussion, and Abstract) and is
the **corresponding author**. **R. A. Bay and A. C. Stier are co-senior authors.** Byline order:
Banerjee, Brzezinski, Diminuco, Seifert, Osenberg, Bay, Stier. (Banerjee's correspondence
department/address/email are still `[to confirm]` in the manuscript masthead.)

- **Shreya Banerjee** (lead, corresponding) — UC Davis Bay Lab; RNA-seq differential expression + matching thickets A/C/D to Cunning's genotyped genets; manuscript narrative (Intro/Discussion/Abstract)
- **Molly Brzezinski** — UCSB Stier Lab; experimental execution, data management, microscope analysis. ORCID 0000-0002-0417-3406
- **Michelle Diminuco** — UGA Osenberg Lab; field collection
- **Ashley W. Seifert** — University of Kentucky; regeneration biology PI. ORCID 0000-0001-6576-3664
- **Craig W. Osenberg** — UGA; ORCID 0000-0003-1918-7904
- **Rachael A. Bay** (co-senior) — UC Davis Bay Lab; gene-expression PI
- **Adrian C. Stier** (co-senior) — UCSB Stier Lab; astier@ucsb.edu, ORCID 0000-0002-4704-4145

## License

Code: MIT. Data: CC-BY 4.0 (until publication, then DOI). Manuscript: All rights reserved.

## Funding

NSF support to A.C. Stier and collaborators. Field work at the UC Berkeley Gump Research Station, Mo'orea, French Polynesia.
