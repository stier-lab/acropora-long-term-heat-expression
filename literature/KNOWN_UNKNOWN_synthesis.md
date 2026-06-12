# What we know and don't know — *A. pulchra* wounding, healing, regeneration & thermal tolerance

A reading synthesis across **this LTH library** (`literature/`, 98 PDFs), the **wound-healing-model
library** (`~/coral-wound-healing-model/literature/`), and the curated NotebookLM corpora
("Acropora pulchra" notebook `f9a97210`; "Coral regeneration all sources" notebook `bb37eb1a`).
Built 2026-06-08 to orient the LTH manuscript Introduction/Discussion and flag where the LTH study
is the first to answer a standing question.

**The one-paragraph version.** The literature *agrees* that heat damages *A. pulchra* physiology,
*disagrees* on whether heat impairs its wound recovery, and has *almost never* (i) separated tissue
healing from skeletal regeneration within the same coral or (ii) linked recovery capacity to a
genotype's measured thermal tolerance. That unresolved core is exactly what LTH targets. The
apparent contradiction in the field is best read as **dose-dependence**: at +1.65 °C *A. pulchra*
recovers fine (Munk 2024), at +3 °C healing proceeds but regeneration stalls (LTH, this study), and
at sub-bleaching +5–6 °C skeletal regrowth halts in a congener (Bonesso 2017, *A. aspera*).

> ⚠ **Library species-attribution corrections found during this synthesis** (now fixed in
> `LITERATURE.md` and `LIBRARY_MAP.md`): four papers were over-attributed to *A. pulchra*.
> **Zhang 2025** is the macroalga *Halimeda macroloba* (NOT a coral); **Almeida 2024** is
> *Pocillopora* cf. *damicornis* (Kenya); **Denis 2024** is *A. spathulata* (GBR); **Raymundo 2025**
> is *A. aspera* outplants (Guam, with *A. pulchra* only as the background congener). Net effect:
> *A. pulchra* proteomic/metabolomic thermal acclimation is essentially **unstudied** (the Zhang
> slot was the only claimed source).

---

## 1. Wounding, healing & regeneration of *A. pulchra*

### Known (direct *A. pulchra* data)
- **Growth & breakage.** Branches extend **13.1–15.8 cm yr⁻¹**; broken branches regenerate
  ~**12.6 cm yr⁻¹** *if they survive*, but **36.4 % of broken branches die** in the wild, with
  zero-growth concentrated in the warm season (Yap & Gomez 1984). Transplantation, prone
  positioning, and sedimentation all retard growth, worst in warm months (Yap 1985).
- **Biphasic recovery is documented** (Munk 2024, Mo'orea MSc thesis — the single closest
  comparator to LTH): healing scored as polyp, corallite, and tissue regeneration; **92 % (13/14)**
  of fragments showed ≥2 of 3 signs by day 19; mean calcification 1.62 mg cm⁻² d⁻¹. **Tissue
  regeneration was faster in abrasion than fragmentation wounds** — a "phoenix effect," where
  residual tissue in the wound subsidizes repair.
- **Restoration workhorse.** Among the top-20 coral species used in restoration globally; forms
  extensive back-reef thickets in Mo'orea (Conetta 2021; Soong & Chen 2003).

### Known (congener, same apical-tip assay as LTH)
- **Bonesso 2017** (*A. aspera*) — sub-bleaching **32 °C halted skeletal regrowth at the tip**
  (~0.5 mm vs ~2.2 mm at ambient over 12 d) and suppressed GFP at the wound. The closest precedent
  to the LTH phase-specific result.

### Not known
- No high-resolution quantitative tissue-closure rates (mm² d⁻¹) for *A. pulchra*.
- Tissue architecture (perforate/imperforate), reserve thickness, and **internal resource
  translocation distance** are unmeasured — directly the model's weakly-grounded `D_E`/`resid`.
- Heterotrophic-compensation capacity during repair is unquantified.
- Whether a colony abandons badly damaged branches; algal/turf colonization of *A. pulchra* wound
  beds (scored but never observed in LTH).
- The **molecular basis** of healing vs regeneration (the LTH RNA-seq, pending).

---

## 2. Thermal tolerance of *A. pulchra*

### Known (direct)
- **Acute limits, your exact population.** CBASS *Fv/Fm* **ED50 = 34.3–36.6 °C** across 20
  Mahana/Mo'orea colonies; ED50 predicts performance in a longer heatwave-like assay (R = 0.74)
  (Cunning 2024).
- **Photophysiology.** Heat causes immediate *Fv/Fm* collapse and photodamage that takes **3× the
  stress duration to repair**, with the worst effects *during* recovery; shading buffers it
  (Berg 2020, *A.* cf. *pulchra*, Guam).
- **Calcification.** Thermal optimum **~28 °C** (asymmetric parabola); at 29.8 °C net calcification
  falls **18 % (high light) to 50 % (low light)** at ambient pCO₂ (Comeau 2014, 2016).
- **Heritable variation.** Intraspecific differences in calcification under temperature × pCO₂, with
  a **growth–tolerance trade-off** — the fastest growers decline most under stress (Shaw 2016, Mo'orea).
- **Symbiont basis.** Fidelity for thermotolerant *Cladocopium* C40 + *Durusdinium* D1 lets
  *A. pulchra* persist through repeated bleaching, with temporally stable communities (Anthony 2023;
  Rouzé 2017/2019); cell density drops sharply after transplantation but communities don't shift
  (Lock 2025).

### Not known
- **Proteomic/metabolomic thermal acclimation** of *A. pulchra* — essentially no data (the prior
  "Zhang 2025" source is a macroalga; see correction above).
- Landscape-scale ED50 mapping (exists for *A. spathulata*, not *A. pulchra*); cold-bleaching
  thresholds for Pacific populations.
- Whether ED50 predicts anything **beyond** *Fv/Fm* — e.g. calcification or regeneration under
  chronic heat (Cunning's genets are not trait-mapped to those outcomes).

---

## 3. The interaction (heat × wounding) — the genuinely unresolved core

### Heat *impairs* recovery
- Meesters & Bak 1993 — bleached colonies regenerate slower, sometimes cease (symbiont loss blocks
  photosynthate/translocation; Fine et al. 2002).
- Bonesso 2017 — skeletal regrowth halted at 32 °C (*A. aspera*, apical tip).
- Rice 2019 — *Pocillopora* (Moorea): wound healing **−66 % at 29 °C**, while growth is maintained →
  a growth-over-repair triage.
- Paradis 2019 — *A. cervicornis*: abrasion + heat drives **P:R < 1** (negative energy balance).
- Kaufman 2021 — *A. cervicornis*: only **35 % healed at 31.5 °C vs 99 % at 28 °C**, and donor-reef
  thermal history predicts healing under heat → a genotype effect.

### Heat *doesn't* impair / even *enhances* recovery
- Dias 2018 — regeneration rates generally *increased* with temperature up to 32 °C in several
  Indo-Pacific species.
- Burmester 2017 — temperate *Astrangia*: recovery rose from 9 → 24 °C (cold is the constraint).
- **Munk 2024 — *A. pulchra* at +1.65 °C (29.5 °C): regeneration and calcification unaffected,
  "robust," warming even raised daily P:R productivity.**
- Traylor-Knowles 2016 — no healing-rate difference between thermal-history regimes (*A. hyacinthus*).

### The reconciliation the field hasn't stated — dose-dependence
The contradiction largely dissolves on temperature dose: **+1.65 °C (Munk) → fine; +3 °C / 31 °C
(LTH) → healing proceeds, regeneration blocked; sub-bleaching +5–6 °C (Bonesso) → skeleton halts.**
The threshold between "tolerable warming" and "regeneration-blocking heat" for *A. pulchra* is
**unmapped**, and the LTH 31 °C treatment sits on that seam.

### Mechanistic candidates (why regeneration is the fragile phase)
- Symbiont-loss-driven **photosynthate/translocation blockage** (Fine 2002; Meesters & Bak 1993).
- **Negative energy balance** (P:R < 1) under combined heat + injury (Paradis 2019).
- **Immune/antioxidant/HSP costs** diverting resources from repair (Madeira 2022; Lock 2022).
- Explicit **growth-vs-regeneration triage** (Rice 2019, citing *Oculina*; Bonesso skeleton-specific;
  Denis 2013, *A. muricata* — fast growth trades off against regeneration).

### What LTH adds (the gaps it closes)
1. **Phase-specificity** — whether heat acts uniformly or **decouples healing from regeneration
   within the same coral**. Almost every prior study measures a single recovery rate; LTH resolves
   it per coral (healing on schedule; 67 % closed-but-never-regenerated at 31 °C vs 0 % at 28 °C;
   new-corallite Cox HR 0.22).
2. **The ED50 ↔ regeneration link** — whether thermally tolerant genets retain regenerative
   capacity. Cunning gives ED50s, Kaufman gives donor-history effects, but no one has correlated a
   genet's thermal threshold with its regeneration under heat. LTH approaches this; the SNP-matching
   of thickets A/C/D to Cunning's genotyped genets is the pending step (see `NEXT_STEPS.md`).
3. **Mechanism** (RNA-seq, pending) and **field realism** (healing under corallivory/algal
   competition + heat) remain open for *A. pulchra*.

---

## Source pointers
Direct *A. pulchra*: Yap & Gomez 1984; Yap 1985; Munk 2024; Berg 2020; Comeau 2014/2016; Shaw 2016;
Cunning 2024; Anthony 2023; Lock 2025; Rouzé 2017/2019; Conetta 2021; Conn 2025 (genome).
Congener / cross-species: Bonesso 2017; Traylor-Knowles 2016; Burmester 2017/2018; Dias 2018;
Meesters & Bak 1993; Meesters 1994/1997; Rice 2019; Paradis 2019; Kaufman 2021; Denis 2013;
Rinkevich 1996; Henry & Hart 2005; Lock 2022; Madeira 2022; Fine 2002.
Full per-paper relevance and file paths: [`LIBRARY_MAP.md`](LIBRARY_MAP.md) ·
discovery index: [`LITERATURE.md`](LITERATURE.md) · grounding audit: [`CITATION_AUDIT.md`](CITATION_AUDIT.md).
