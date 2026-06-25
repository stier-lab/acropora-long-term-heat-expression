# Stier-lab phenotype framing angles — Intro/Discussion (optional, not the spine)

> 🗂️ **Optional framing angles + citation bank** · Updated 2026-06-12 · Index: [`README.md`](../../README.md) · references → `manuscript/references.bib`, `literature/CITATIONS_INDEX.md`.

**Shreya — you lead this paper. The recommended spine is yours and should lead with the
transcriptomics.** This file is *not* the manuscript text and *not* a recommendation about how to
frame the paper. It is the Stier lab's phenotype-side thinking: candidate framings (deliberately left
in tension, not resolved into one) and a verified citation bank you can draw on if any of them is
useful. Take a framing, combine them, or set them all aside once the gene-expression results are in.

These are **phenotype-first** angles, written before your gene-expression results exist. The real
paper is led by the transcriptomics, so none of these is "the" framing — they're raw material.

---

## Candidate framings (in tension — pick, blend, or discard)

**A. Energetic-triage / phase-decoupling (phenotype-first).**
- Spine: heat does not slow recovery uniformly; it spares the tissue-healing phase (coenosarc
  closure) and arrests the regeneration phase (skeletal regrowth at the tip).
- Hook: survival of the tissue is a poor proxy for recovery of structure — a *restoration* angle
  (selecting genets on bleaching survival may propagate corals that persist without regenerating).
- Tension: this is an *inference* about energy allocation that the phenotype alone cannot prove; it
  pre-commits the paper to an allocation story your expression data may support, extend, or overturn.

**B. Molecular basis of the healing→regeneration transition (transcriptomics-led).**
- Spine: the phenotype locates *where* recovery fails (the healing→regeneration seam under heat);
  the expression data ask *why* — which programs are induced vs suppressed at that seam.
- Tension with A: leads with mechanism, not allocation; lets the molecular result, not the
  energetic narrative, set the paper's claim. Likely the stronger spine, but it's yours to judge.

**C. Genotype-dependent thermal tolerance of regeneration (cross-method).**
- Spine: regenerative capacity under heat is heritable (genet C resilient vs A/D); pairs with an
  independent acute CBASS assay on the same population (Cunning et al. 2024) and with candidate
  thermal-tolerance genes from the expression contrast.
- Tension with A/B: foregrounds heritable variation and the acute↔chronic validation rather than the
  within-coral phase story; could be the lead, a major thread, or a Discussion section.

**Note:** A and B make different *claims* about the same data; C can sit under either. They are listed
in tension on purpose — resolving them is a framing decision that is yours to make with the expression
results in hand.

---

## Phenotype results these framings rest on (numbers → `output/tables/20_master_results.csv`)

- Tissue healing (coenosarc closure) statistically indistinguishable between 28/31 °C; new-corallite
  regeneration is delayed/suppressed under heat in the interval-censored model (time ratio = 1.32,
  95% CI 1.19–1.47, p = 1.4e-7; first-observed Cox HR = 0.22, 95% CI 0.07–0.69).
- Closed-but-never-regenerated: 67% at 31 °C vs 0% at 28 °C (median healing-to-regeneration lag 10 vs 8 d).
- Whole-colony physiology responds to temperature, barely to wounding; areal calcification 38% lower at 31 °C.
- Genet C multivariate displacement 1.06 vs 3.72 (A) / 3.34 (D); most likely to regenerate at 31 °C.
- 31 °C sits ~4.4 °C below the population acute ED50 (35.4 °C; Cunning et al. 2024) — chronic/sublethal.

**Fragility to carry forward (keep honest):** the regeneration result is strongest for interval-censored
new-corallite onset (p = 1.4e-7); tip-exist is also delayed, while tip-extension points the same way but
is n.s. The first-observed Cox model for new corallites is close to the PH diagnostic boundary, so the
interval model and censored-fraction result are the cleaner anchors. Three genets resolve variation but not its architecture.
Apical-tip excision is not identical to a surface wound bed (Munk 2024). See `RESULTS.md` §10.

---

## Verified citation bank

Every entry below has a DOI in `manuscript/references.bib` (index: `literature/CITATIONS_INDEX.md`),
so any sentence you build on one is on solid ground. Grouped by what they support.

**The biphasic healing↔regeneration distinction, growth & injury ecology of branching *Acropora***
- Henry & Hart 2005 — wound closure precedes regeneration across animals (the two-phase logic).
- Yap & Gomez 1984 — *A. pulchra* branch extension 13–16 cm yr⁻¹.
- Highsmith 1982; Madin et al. 2014 — fragmentation, breakage, and the injury-prone arborescent skeleton.

**Thermal sensitivity of *Acropora* and *A. pulchra***
- Hoegh-Guldberg 1999; Hughes et al. 2017 — *Acropora* among first to bleach/die under warming.
- Berg et al. 2020 — *A. pulchra* persistent photosystem damage and mortality under sustained heat.

**Heat drains the coral energy budget (photochemistry, symbionts, calcification)**
- Warner et al. 1999 — PSII damage lowers photochemical efficiency under heat.
- Hoegh-Guldberg 1999; Jokiel & Coles 1990 — pigment/symbiont loss, fixed-carbon supply.
- Jokiel & Coles 1977; Comeau et al. 2014 — warming suppresses calcification.

**Prior heat × injury work (reports recovery as a single rate — the gap this study addresses)**
- Meesters & Bak 1993 — bleaching reduces lesion regeneration.
- Bonesso et al. 2017 — sub-bleaching warming slows post-injury recovery.
- Traylor-Knowles et al. 2016 — temperature reshapes branching-coral wound-healing response.

**Heritable, genotype-level thermal tolerance & the cross-method link**
- Dixon et al. 2015; Shaw et al. 2016 — thermal tolerance is heritable and varies among genotypes.
- Cunning et al. 2024 — CBASS *Fv/Fm* ED50 for Mahana *A. pulchra* (the acute assay to match against).

**Wound-geometry caveat**
- Munk 2024 — surface wound-bed polyp reappearance (not identical to apical-tip excision).
