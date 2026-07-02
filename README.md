# LTH: Long-Term Heating × Wounding × Gene Expression in *Acropora pulchra*

> Project #17 of the Coral Regeneration program · Gump Station, Mo'orea, French Polynesia · May–July 2025
> · Repo: `stier-lab/17.-LTH_expression_by_temperature_2025` · Last updated 2026-07-02

<p align="center">
  <img src="docs/team_summary/img/healing_ambient_vs_heated.gif" width="720"
       alt="Wound-healing time-lapse, ambient 28 °C vs heated 31 °C, Day 0 to Day 15"><br>
  <sub><b>Wound healing over 15 days.</b> <b><span>Ambient&nbsp;28&nbsp;°C</span></b> (left, blue) closes the wound <b>and</b> regrows the tip;
  <b>Heated&nbsp;31&nbsp;°C</b> (right, orange) closes the wound but never regenerates.</sub>
</p>

A heat × wounding experiment on the branching coral *Acropora pulchra*. We clipped the growing tip off
half the fragments, held corals at **28 °C (ambient)** or **31 °C (heated, +3 °C)**, and tracked
healing and regrowth over 15 days — with photochemistry, pigmentation, symbionts, and growth. We also
took tissue for gene expression (RNA-seq, in progress).

**The paper's lead result is the transcriptomic mechanism** (RNA-seq pending);
the phenotype experiment here supplies its **organismal context**. **The phenotype result:** heat does
*not* slow wound recovery uniformly — corals seal the wound (tissue healing) at the same rate hot or
not, but heated corals fail to rebuild skeleton at the tip (regeneration). Heat affects one phase of
recovery, not the other.

For the full results narrative see **`RESULTS.md`**; for the RNA-seq / lead-author brief see
**`docs/rnaseq/`**.

## Contents

- [Quick start](#quick-start)
- [Working in this repo](#working-in-this-repo)
- [How the repo relates to Google Drive](#how-the-repo-relates-to-google-drive)
- [Research question & design](#research-question--design)
- [Findings (the phenotype half)](#findings-the-phenotype-half)
- [Statistical approach](#statistical-approach)
- [Reproducing the analysis](#reproducing-the-analysis)
- [Repository map](#repository-map)
- [Data inventory](#data-inventory)
- [Glossary](#glossary)
- [Regeneration-staging terminology](#regeneration-staging-terminology)
- [Pending work](#pending-work)
- [License and funding](#license-and-funding)

---

## Quick start

```bash
cd ~/Stier-LTH-expression-by-temperature-2025
Rscript -e 'install.packages(c("tidyverse","here","readxl","janitor","scales","patchwork","lme4","lmerTest","blme","nlme","MuMIn","car","emmeans","DHARMa","broom","broom.mixed","influence.ME","lmtest","survival","survminer","xml2"))'   # one-time: install packages (current CRAN versions)
Rscript code/_run_all.R          # full pipeline in dependency order (~4 min); regenerates every figure + table
```

Everything is reproducible from `code/_run_all.R`. **Never hand-edit numbers** — every statistic lives
in `output/tables/20_master_results.csv` (`+ _paper_ready.csv`) and regenerates with the pipeline. The
last step (`code/30_manuscript_audit.R`) recomputes every phenotype number and *warns* if the
manuscript has drifted out of sync (advisory only — never fails the run). Requires R ≥ 4.3 (developed
on 4.5.2).

## Working in this repo

The gene-expression analysis and the manuscript narrative are out of scope for this repo — keep to these scope rules:
- Don't write or pre-empt the Introduction / Discussion / Abstract, and don't bake Discussion-level
  interpretation into the Results (effect sizes + direction only there).
- Don't prescribe the RNA-seq pipeline. The `docs/rnaseq/` files pose *questions/goals* and offer
  *suggestions*; phrase every phenotype↔expression link bidirectionally ("expression can test, extend,
  **or revise** this").
- Interpretive labels (e.g. "energetic triage") are out of scope here, not house vocabulary.
- The manuscript audit (`code/30_manuscript_audit.R`) is advisory and scoped to phenotype
  Methods/Results only — it does not police the narrative.

Still open for PI/institutional decision: target journal & review process; and rehoming/renaming the repo.

## How the repo relates to Google Drive

```
 Google Drive  ── source of truth ──►  THIS REPO  ── reproducible analysis ──►  results + manuscript
 (raw Sheets,                          (data/raw/ is exported FROM Drive;
  field notes,                          the repo NEVER writes back to Drive)
  manuscript Doc)
```

This repo is the **analysis layer**; the **project of record** (raw Sheets, field notes, manuscript
Doc, slides) lives in Drive:

- **Folder:** `17. LTH_expression_by_temperature_2025` · **Drive ID:** `1sXfnHN-vmSBuwMEfERiYOWDeKRmjFWJP`
- **Path:** `…/My Drive/Stier Lab/People/Adrian Stier/Projects/In Progress/Coral-Regeneration/Projects/17. LTH_expression_by_temperature_2025/`
- **Access** via the `gws` CLI (authenticated as astier@ucsb.edu), e.g.
  `gws drive files list --params '{"q":"'\''1sXfnHN-vmSBuwMEfERiYOWDeKRmjFWJP'\'' in parents and trashed=false","fields":"files(id,name,mimeType)"}'`
- **Mapping:** Drive `data/` Sheets → `data/raw/` here; Drive `notes/` Docs → `notes/` here; Drive
  `Manuscript_LTH` Doc ↔ `manuscript/Manuscript_LTH.md`. When Drive data changes, re-export into
  `data/raw/` and re-run `code/_run_all.R`.

## Research question & design

Branching corals like *A. pulchra* regenerate after damage to the apical (growing) tip. Elevated
temperatures slow growth and bleach symbionts. **How do those two stressors interact?** Does mild
long-term heating (+3 °C) weaken wound regeneration, and what gene-expression machinery underlies the
response?

| Factor | Levels |
|---|---|
| **Temperature** | 28 °C (ambient) vs 31 °C (heated). Ramped 1 °C/day from ambient. |
| **Wound** | Wounded (~1 cm clipped off growing tip with band-saw + caliper) vs unwounded sham. Applied 7 days after corals reached the target temperature. |
| **Genotype (thicket)** | A, C, D — modeled as a **fixed effect** (only 3 field-collected genets; too few for a variance component). |
| **Tank** | 8 total: 4 per temperature (28 °C: 3, 6, 9, 12; 31 °C: 4, 5, 10, 11). Random effect. |
| **Time** | Daily non-destructive obs. Destructive biopsies at D0, D1, D3, D10, D15. |

**Coral fragments (n = 208 total):** 192 for gene expression + destructive physiology (3 genotypes ×
4 tanks × 2 temperatures × 2 wound × ~4 timepoints); 48 of those also tracked for non-destructive
metrics (PAM, color, buoyant weight); 16 dedicated to daily microscope photography (genets A and C only).

**Collection site:** Mahana / Tiahura, NW Mo'orea — three parent thickets at A 17.49735 °S /
149.91557 °W, C 17.49808 °S / 149.91595 °W, D 17.49726 °S / 149.91581 °W. Same named site as Cunning
et al. 2024's CBASS genets. Protocols: sequencing / gene-expression design in `docs/rnaseq/`;
earlier field and microscope-photo plans in `notes/archive/`.

## Findings (the phenotype half)

The four take-homes below are the phenotype results — organismal context, not the
paper's lead result. Every number traces to `output/tables/20_master_results.csv`; the full narrative
(with caveats) is in `RESULTS.md`. Summary figure: `figures/16_manuscript_fig1.pdf`.

**1 — Sustained heat broadly compromises physiology.** At 31 °C, photochemistry, pigmentation,
symbiont density, and growth all declined progressively while 28 °C held steady. By Day 14
heated corals had paled (58–67 % vs 0–8 % ambient) and grew **34 % less** (6.10 → 4.03 % skeletal
mass change). *How asked:* linear mixed model per response (`response ~ treatment × wound × day × genet`
+ random tank & coral); the heat signal is the **treatment × day interaction** (the rate of
divergence). All significant — Fv/Fm F = 111.7, color F = 287.6, symbionts F = 94.0; all p < 0.001
(type-III ANOVA).

**2 — Heat blocks *regeneration*, not *healing* (the phenotype result).** Wounds *sealed* at the
same rate in both temperatures, but new corallites formed in 100 % of ambient wounded corals vs only
33 % of heated ones — and at 28 °C every coral that healed went on to regenerate, while at 31 °C 67 %
healed but never rebuilt skeleton. *How asked:* each binary recovery trait → a time-to-event interval;
primary test is an **interval-censored Weibull AFT** (new-corallite onset time ratio = 1.32, 95 % CI
1.19–1.47, p = 1.4e-7), with Kaplan–Meier / Cox as first-observed-day summaries (Cox HR = 0.22,
p = 0.010). Survival analysis fits because many corals never reach the milestone in 15 days ("censored"). See
`figures/14_morphology_KM.pdf`.

**3 — Genotype matters: a resilience gradient C > D > A.** Genet C consistently defended its
physiology and was likeliest to regenerate under heat; A and D were sensitive. In multivariate
physiology space, genet C's state shifted **3.6× less** under heat than genet A's (PCA centroid
displacement 1.02 vs 3.71). *How asked:* genet as a fixed effect (only 3 levels); likelihood-ratio
tests on the genet × treatment interaction (Fv/Fm, color, symbionts χ² = 78.1, 149.8, 64.3; all
p < 0.001). A composite standardized sensitivity across 11 responses (A = 0.56, D = 0.43, C = 0.01)
and the PCA displacement agree. See `figures/19_genet_dashboard.pdf`.

**4 — This is chronic-sublethal stress, not acute bleaching.** Benchmarked against an independent
acute heat-tolerance assay for the same species and reef (Cunning et al. 2024 CBASS Fv/Fm ED50), the
mean ED50 is **35.4 °C**, so 31 °C sits **~4.4 °C below** the acute threshold. The declines reflect
sub-bleaching stress accumulated over weeks — the realistic regime for recurrent moderate warming —
not acute photoinhibition. See `figures/26_thermal_context.pdf` (or `figures/08_apex_temperature.pdf`
for treatment validation: tanks held within ~0.3 °C of setpoints).

## Statistical approach

Genet (thicket) is a **fixed effect** throughout (only 3 field-collected genets — too few for a
reliable variance component; Bolker 2008, Gelman 2005), which also exposes per-genet effects directly.
Reference levels: treatment = `28C`, wound = `no`.

- **Continuous responses** (PAM Fv/Fm, color D-scale, log symbiont density, growth):
  LMMs `response ~ treatment * wound * day * thicket + (1|tank) + (1|id)` (`lme4`/`lmerTest`,
  type-III Satterthwaite; `(1|id)` dropped for single-observation responses).
- **Growth = % skeletal mass change** is the primary metric; specific growth rate (SGR) is the
  robustness check. No areal calcification rate — we did not measure whole-fragment surface area (only
  the symbiont sub-fragment was wax-dipped), so we report SA-free metrics.
- **Morphology** (8 binary wound-healing traits): binomial GLMMs; traits with separation refit with
  Cauchy(0,2.5) priors (`blme`). **Color** D-scale also refit as an ordinal CLMM (`ordinal`).
- **Healing milestones:** interval-censored Weibull AFT for inference; Kaplan–Meier + Cox PH
  first-observed-day summaries with full Schoenfeld PH diagnostics.
- **Multivariate:** centered/scaled PCA on the four endpoint responses.
- **Multiplicity:** confirmatory a-priori (literature-grounded) tests reported unadjusted; exploratory
  tests BH-corrected (`code/sensitivity/28_multiple_testing.R`).
- **Diagnostics:** DHARMa for every model + dedicated time-series (AR(1), nonlinearity, random slopes)
  and Cox PH suites; coverage audited across all models. Contrasts Tukey-adjusted (`emmeans`).

## Reproducing the analysis

```r
# one-time: install the packages listed in Quick start (current CRAN versions), then:
source("code/_run_all.R")    # full pipeline in dependency order (~4 min)
```

Outputs: figures in `figures/`, tables in `output/tables/`, diagnostics in `output/diagnostics/`,
processed data in `data/processed/`, fitted models in `output/models/`. The final step rebuilds the
master results table and runs the manuscript audit.

**Pipeline order:** `01` metadata → `02` PAM → `03` color → `04` morphology → `07` wax SA → `05`
growth → `06` symbionts → `08` APEX temp → `09` YSI → `10` worms → `11` combined figure → `12`
models (primary + color-CLMM & morphology-blme robustness) → `13` genet interaction → `14`
interval/KM/Cox timing → `15` PCA → `16` manuscript figure → `17` figure audit → `18` data validation
→ `19` genet dashboard → `sensitivity/22–29` (flagged-sample, time-series, headline comparison,
thermal context, variance partitioning, multiple testing, prob contrasts) → `diagnostics/A–H` →
`20` master results table → `31` RNA-seq covariate table → `30` manuscript audit. (`21_rnaseq_stub.R`
waits on sequencing; not in the run-all. Exact order in `code/_run_all.R`.)

## Repository map

| Path | Contents |
|---|---|
| `code/` | All analysis scripts; run order = file number (see [Reproducing](#reproducing-the-analysis)). Each script opens with a Purpose / Input / Output header. |
| `code/12_models.R` | Primary mixed models + color-CLMM and penalized-morphology robustness (merged from the former 12/12b/12c). |
| `code/sensitivity/` | Robustness & sensitivity analyses (22–29): flagged samples, time-series, headline-model comparison, model coverage, thermal context, variance partitioning, multiple testing, probability-scale contrasts. |
| `code/diagnostics/` | Model-diagnostic suites (A–H): DHARMa residuals, design alignment, reproducibility, spreadsheet coverage. |
| `code/archive/molly_original/` | Molly's original exploratory scripts — historical, **not part of the pipeline**; preserved for provenance. |
| `data/raw/` | Exported from Drive — never edited by hand. |
| `data/processed/` | Cleaned `.rds` files the pipeline produces (regenerable). |
| `data/external/` | Cunning et al. 2024 CBASS ED50 reference data (+ provenance README). |
| `output/tables/` | Every result as CSV. `20_master_results.csv` is the single source of truth. |
| `output/models/` | Saved fitted models (`.rds`). |
| `output/diagnostics/` | Model-diagnostic reports (A–K), regenerated by the pipeline. |
| `figures/` | All figures (`.pdf` + `.png`) + `diagnostics/`; catalogued in `figures/FIGURE_INDEX.md`. |
| `literature/` | PDF library (101 PDFs) + `LITERATURE.md` (consolidated bibliography + know/don't-know synthesis). |
| `manuscript/Manuscript_LTH.md` | Working draft — phenotype Methods + Results drafted; Intro/Discussion/Abstract are out of scope here. |
| `notes/` | Sequencing / gene-expression plan. `notes/archive/` holds field notes, QA/QC flags, photo indices, and superseded plans/`growth_allometry.md`. |
| `docs/rnaseq/README.md` | Single lead-author brief: RNA-seq design, analysis proposal + genet-matching, phenotype↔expression map, candidate genes (suggestions, not a prescribed pipeline). |
| `RESULTS.md` | Full results narrative (all responses, genet effects, thermal context, §10 limitations). |

## Data inventory

| Stream | Source file | Rows | What it captures |
|---|---|---|---|
| Coral metadata (one per fragment) | `data/raw/metadata/metadata.csv` | 208 | thicket, id, tank, treatment, wound, biopsy day/date, **coord_lat/long**, calculated SA, planned chl-a slot, zoox |
| PAM (Fv/Fm) | `data/raw/pam/PAM_data.csv` | 672 obs | F, M, Y, E, Fv/Fm by date × tank × sample × location (top/bottom) |
| Color card (Siebeck D-scale) | `data/raw/color_card/data.csv` | 336 obs | health_status, color (split scores averaged), paling, hole_at_center |
| Morphology (9 raw columns → 8 analyzed traits) | `data/raw/physio_morphology/data.csv` | 768 obs | polyps_out, hole_in_center, polyp_in_hole, wound_smoothed, pigment_over_wound, tip_exist, tip_extension, new_corallites_on_tip, algae_on_wound (the byte-identical hole_in_center/polyp_in_hole are combined to `axial_polyp_formation` in analysis — see `data/raw/physio_morphology/SCORING_NOTES.md`) |
| Buoyant weight (growth) | `data/raw/buoyant_weight/data.csv` | 48 corals | initial + final coral & plug weights, water-density correction → % skeletal mass change |
| Wax dipping (surface area) | `data/raw/wax_dipping/data.csv` | 192 corals | diameter, height, dry weight, wax weights, SA from standard curve |
| Standard curve for wax | `data/raw/wax_dipping/Standard_curve.csv` | 19 | wax mass vs SA from cylinders |
| Symbiont counts | `data/raw/symbiont_counts/Raw_counts.csv` | 768 (4 reps × 192) | hemocytometer counts → cells/cm² via wax SA |
| Worm presence | `data/raw/worm_presence/Sheet1.csv` | 192 × 3 dates | AEFW surveillance 06/07, 06/08, 06/12 |
| APEX temperature/pH | `data/raw/apex/datalog*.xml` | continuous | per-tank temperature and pH, May–June |
| YSI daily spot checks | `data/raw/ysi/Sheet1.csv` | 72 obs | TEMP, DO%, DO mg/L, SAL, pH daily |
| Sequencing plate layout | `data/raw/plate_layout/Plate_{1,2}.csv` | 96 wells × 2 | wells, sample ID, day, temp, wound, tank, genotype |
| Sample shipment manifest | `data/raw/shipping/*.csv` | varies | shipped to UC Davis Bay lab |
| **External: Cunning 2024 ED50** | `data/external/cunning2024_apulchra_ed50.csv` | 20 genets | acute CBASS Fv/Fm ED50, Mahana (thermal-tolerance benchmark) |

**Photographs** (color-card, microscope) live on the Stier Lab NAS (`smb://stier-nas1.eemb.ucsb.edu`);
see `notes/archive/LTH_*_Photos.md`. **RNA-seq reads** are processed at UC Davis (Bay lab); NCBI BioProject TBD.

## Glossary

| Term | What it means here |
|---|---|
| **Linear mixed model (LMM)** | Regression for the continuous responses that accounts for repeated measures (same coral over time) and grouping (corals in tanks) via random effects. |
| **treatment × day interaction** | The *rate* at which heated and ambient corals diverge — the heat signal, not the day-0 main effect. |
| **Interval survival / Kaplan–Meier / Cox** | Survival analysis for "time until a recovery milestone." Interval models handle discrete scoring days; KM/Cox summarize first observed onset. |
| **genet as a fixed effect** | With only 3 genets, each one's effect is estimated directly rather than as a variance component (more reliable with so few groups). |
| **type-III ANOVA** | Tests each term adjusting for all others (the correct default for interaction models). |
| **confirmatory vs exploratory** | A-priori predictions (literature-grounded) reported unadjusted; exploratory tests BH-corrected (`code/sensitivity/28`). |
| **Tissue-healing phase** | Re-epithelialization / coenosarc coverage that seals the wound (traits `axial_polyp_formation`, `wound_smoothed`). |
| **Regeneration phase** | Reappearance of polyps + skeletal/calyx/corallite regrowth (traits `tip_exist`, `tip_extension`, `new_corallites_on_tip`). |
| **Healing-to-regeneration lag** | The gap between sealing the wound and rebuilding skeleton (`code/14`). At 31 °C, two-thirds of corals close the wound but never cross into regeneration. |

## Regeneration-staging terminology

The **biphasic** framework above is the lab's working vocabulary, aligned with the sibling wound-type
project so they read together. Use it as a default for the *staging* descriptors, but it is **not**
binding on this paper — the framing and interpretive labels (e.g. "energetic
triage") are out of scope here. In prose (not data/column names, which stay fixed):

- Say "tissue healing / coenosarc coverage", not just "wound closure".
- Use **"wound bed"** for the active region; **"algal colonization of the wound bed"** (or "algal
  plug" when dense) rather than "algae on wound".
- State the geometric caveat: LTH wounds are apical-**tip excisions** (regeneration =
  new skeleton at the branch tip), unlike a surface wound bed — don't imply identical assays.

## Pending work

See `docs/rnaseq/` for the RNA-seq resources. Outstanding items:

- **RNA-seq (UC Davis Bay lab):** 144 libraries shipped (`data/raw/shipping/Box_list.csv`). When
  counts land, they go in `data/raw/sequencing/` (`counts.csv` + `sample_metadata.csv`, matched to
  `data/raw/plate_layout/Plate_{1,2}.csv`); gitignore the FASTQ stage, commit the count matrix. The DE
  / WGCNA / SNP design is out of scope here. A per-library phenotype covariate table is already built:
  `output/tables/31_rnaseq_phenotype_covariates.csv` (`code/31_rnaseq_covariate_table.R`).
- **Genet matching:** call SNPs from host RNA-seq → per-thicket genotypes → match A/C/D to Cunning's
  CBASS genets, enabling an acute-ED50-vs-chronic-resilience test (`docs/rnaseq/README.md`).
- **Chlorophyll-a:** tissue slurries frozen but assay not run; the metadata slot remains for provenance.
- **Validations before submission** (PI judgment, not code): Day-1 symbiont gap; the
  non-significant wound × treatment color interaction (type-III F₁,₅₉ = 1.04, p = 0.31); genet-C
  resilience with only 3 genets; the singular `(1|tank)` fit in the growth model.
- **Manuscript narrative:** Introduction, Discussion, Abstract (out of scope here). Phenotype Methods +
  Results are drafted in `manuscript/Manuscript_LTH.md`.

## License and funding

**License:** Code MIT · Data CC-BY 4.0 (until publication, then DOI) · Manuscript all rights reserved.
**Funding:** NSF support to the project investigators. Field work at the UC Berkeley Gump Research
Station, Mo'orea, French Polynesia.
