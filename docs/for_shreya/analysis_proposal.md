# LTH RNA-seq: analysis proposal & genet-matching plan

**For:** Shreya Banerjee (UC Davis, Bay lab)
**From:** Adrian Stier lab (LTH project #17 — *Expression by Temperature 2025*)
**Re:** what the 144 RNA-seq libraries should test, and how to link the three
LTH thickets (A, C, D) to Ross Cunning's genotyped *A. pulchra* CBASS genets.
**Status:** phenotype analysis complete; this is the open genomics piece.
**Last updated:** 2026-06-12 · **Document index:** [`README.md`](../../README.md)

> **You own the gene-expression analysis — this is a resource, not a plan you have to follow.**
> Everything below is the Stier lab's *suggestion* given the phenotype results. The model, the
> contrasts, and the genet-matching strategy are yours to adopt, change, or discard; we wrote it down
> so the two halves of the paper connect and so you have a starting point if it's useful, not to
> prescribe your pipeline. Section 5's "next steps" is one possible division of labor — please
> rewrite it however fits how you want to work.

---

## 1. One-paragraph context

The LTH experiment crossed **wounding** (clip the growing tip, yes/no) with
**chronic heating** (28 °C vs 31 °C, weeks) in *Acropora pulchra* from Mahana,
Mo'orea, using **three field-collected parent thickets (A, C, D)**. The
phenotypic analysis is done and points to three results the transcriptomics
should explain mechanistically:

1. **Heat impairs *regeneration*, not wound *closure*.** Wounds seal at the same
   rate in both temperatures, but new-corallite formation at the regenerating
   tip is delayed/suppressed under heat (interval-censored time ratio = 1.32,
   p = 1.4e-7; first-observed Cox HR = 0.22). The lesion is in the
   skeletal-regrowth program, not re-epithelialization.
2. **Genotype-level thermal tolerance: C > D > A.** Genet C is the most resilient
   across photochemistry, pigmentation, symbiont retention, and calcification;
   the heat effect on physiology is 2.7–3.5× weaker in C than in A or D.
3. **The 31 °C stress is chronic-sublethal.** It sits ~4.4 °C *below* the acute
   Fv/Fm ED50 for this population (35.4 °C; Cunning et al. 2024), so the
   responses are accumulated sublethal stress, not acute photoinhibition.

Full numbers: `RESULTS.md` and `output/tables/20_master_results.csv` (single
source of truth). Sampling design: `notes/sequencing-plan-keck-LTH.md`.

---

## 2. RNA-seq design (as built)

144 libraries (`notes/LTH_PlateLayout_with_IDs`, plate maps in `data/raw/plate_layout/`):

- **Margin biopsies** (wound margin vs matched unwounded margin) at **Days 1, 3, 10**
- **Factors:** temperature (28/31) × wound (W/U) × genotype (A/C/D) × day, with
  4 tanks/temperature; fully interleaved across 2 plates / 1 NovaSeq run
  (temperature, wound, day, tank, genotype all orthogonal to plate/lane).

This balance means the design cleanly supports the contrasts below without
batch confounding.

## 3. Questions the differential-expression analysis could address (phenotype-anchored)

These are **goals/questions**, not a prescribed model or pipeline — the design (factors, normalization,
fixed/random structure, tool) is yours to specify. The phenotype results raise the questions below; the
expression data can **test, extend, or revise** each one. Predicted *processes* are listed as
hypotheses to confirm or overturn — we deliberately do **not** name candidate gene symbols, so the gene
sets you interrogate are your call.

- **Does heat suppress the regeneration program while sparing the healing program?** Phenotype: wounds
  seal at the same rate in both temperatures, but new-corallite formation is delayed/suppressed under heat
  (interval-censored time ratio = 1.32, p = 1.4e-7; first-observed Cox HR = 0.22). The expression data can ask whether, at the wound margin, the
  early-healing processes (re-epithelialization, immune, ECM remodeling, proliferation) are induced
  similarly in both temperatures while skeletal/biomineralization and corallite-patterning processes
  are specifically suppressed under heat — **or** whether heat suppresses healing too, which would
  *revise* the phase-decoupling story. This is the crux and could go either way.
- **What distinguishes the resilient genet (C) from the sensitive ones (A, D)?** Phenotype: genet C
  defends photochemistry, pigmentation, and symbiont retention far better than A/D, with the spread
  sharpest in unwounded corals. The expression data can ask whether genet C shows a smaller heat-induced
  transcriptional shift and/or a constitutively frontloaded stress-tolerance signature (proteostasis,
  antioxidant/ROS handling, symbiosis) — and may identify, extend, or contradict candidate
  thermal-tolerance processes underlying that heritable variation.
- **Is the chronic 31 °C signature sustained/constitutive rather than an acute heat-shock spike?**
  Phenotype/design: 31 °C sits ~4.4 °C below the acute ED50 and wounding came after 7 days at
  temperature, so a transient heat-shock response has likely subsided. The expression data can test
  whether the signature is accumulating sublethal stress vs acute photoinhibition.
- **Does the wound response narrow the genet differences (matching the physiology)?** Phenotype: the
  genet spread in heat sensitivity is large when unwounded and compresses when wounded. The expression
  data can test whether the genet (and genet × temperature) effect is correspondingly smaller in
  wounded than unwounded margins.

Tie-in: the **Day 10** tip biopsies are where new-corallite divergence is sharpest in the phenotype
data, so they are likely the highest-yield timepoint for the first question — but that is a phenotype
expectation the expression data can confirm or move.

Symbiont (*Symbiodiniaceae*) reads, if retained, could corroborate (or qualify) the measured
symbiont-density loss and genet C's retention.

---

## 4. The genet-matching question (A, C, D ↔ Cunning genets)

**Why it matters.** Cunning et al. 2024 (*Coral Reefs*, doi:10.1007/s00338-024-02577-7)
measured **acute** CBASS Fv/Fm **ED50** for **20 genotyped *A. pulchra* genets
from Mahana, Mo'orea** (range 34.4–36.6 °C; ED50 predicts classic bleaching,
R = 0.74). Our experiment measured **chronic** resilience for 3 thickets from the
same site. If we can match thickets A/C/D to Cunning genets, we can ask a novel
question: **does an 18-hour acute assay predict chronic, wound-context
resilience?** That would be a strong, publishable cross-method validation.

**The obstacle.** The LTH thicket labels A, C, D are **arbitrary** — there is no
a priori link to Cunning's numbered genets.

### What we already have to work with

- **Exact GPS for each LTH parent thicket** (decimal degrees, from
  `data/raw/metadata/metadata.csv`, columns `coord_lat`/`coord_long`):

  | Thicket | Latitude (°S) | Longitude (°W) | n fragments |
  |---|---|---|---|
  | **A** | 17.49735 | 149.91557 | 72 |
  | **C** | 17.49808 | 149.91595 | 72 |
  | **D** | 17.49726 | 149.91581 | 64 |

  The three thickets are ~40–90 m apart, all within the Mahana / Tiahura
  *A. pulchra* stand on NW Mo'orea — the **same named site** Cunning sampled
  ("mahana", collected Dec 2022; see his repo `data/classic_cbass/collection_metadata.csv`).

- **Host RNA-seq reads from all three thickets** (forthcoming) — enough to call
  genotype-distinguishing SNPs.

- **Cunning's genets are genotyped** — `data/reproducibility/genet_map.xlsx` and
  the CBASS data in `github.com/jrcunning/CBASS_methods`; per-genet ED50 is
  cached for us at `data/external/cunning2024_apulchra_ed50.csv`.

### The goal, and the one external ask that unlocks it

**Goal.** Link thickets A/C/D to Cunning's genotyped genets well enough to ask whether **acute CBASS
ED50 predicts our chronic, wound-context resilience ranking (C > D > A)** — and, if matched, to relate
that to the genet × temperature expression signal. The *how* (reference genome, variant caller, the
relatedness/identity metric, how much to lean on the spatial check) is yours to decide; we have no
prescription here.

**The one thing that actually requires us / an external request.** The pivotal, repo-external dependency
is **Cunning's per-genet host SNP genotypes** (not just ED50 + genet number). Likely held by
Cunning / Putnam — the CBASS_methods repo / associated genomic archive, or a direct request (co-authors
Detmer & Moeller are in the UCSB/Mo'orea network — an easy intro). We can chase this for you. Everything
downstream of obtaining comparable genotype data is analysis you own.

**Supporting data already in hand.** The GPS table above places all three thickets in the same named
Mahana/Tiahura stand Cunning sampled; spatial proximity is **suggestive, not conclusive** (*A. pulchra*
forms clonal thickets, so adjacent colonies can be the same genet and distant ones clonemates). The
cached per-genet ED50 is at `data/external/cunning2024_apulchra_ed50.csv`; our per-thicket chronic
resilience scores are in `output/tables/19_genet_resilience_summary.csv`.

### Honest fallback

If genotype data can't be obtained or the panels aren't comparable, we fall back
to a **population-level** statement (already in the manuscript): both acute
(CBASS) and chronic (LTH) methods independently detect substantial heritable
thermal-tolerance variation in Mahana *A. pulchra*. The genet-level correlation
is the upside if matching succeeds.

---

## 5. A possible division of labor (yours to redraw)

This is just a sketch of how the pieces *could* split — please adjust to match how you want to run
the analysis. The only items that are genuinely ours to chase are the data requests and the site
label.

- The **gene-expression analysis** (DE design, SNP-calling, matching) is yours to lead and structure
  however you see fit — Sections 3–4 are suggestions, not a checklist.
- **We can help with:** requesting Cunning/Putnam host genotype data + exact genet collection coords
  and confirming the reference genome (the pivotal external ask); and confirming the LTH collection
  site label with the Gump/field team (coords say Mahana/Tiahura).
- **Anything you need from us** to make your side easier — a per-library phenotype covariate table,
  the cached Cunning ED50 data, clarification on the design — just ask.

Questions or data requests → Adrian. Phenotype data, code, and the cached Cunning ED50 table are
all in this repo.
