# Gene-expression ↔ phenotype integration map (suggestions)

> 🗂️ **Phenotype → expression hypotheses (bidirectional)** · Updated 2026-06-12 · Index: [`README.md`](../../README.md) · covariates → `output/tables/31_rnaseq_phenotype_covariates.csv`.

**Shreya — these are suggestions, not a plan you have to follow.** You own the gene-expression
analysis and the paper's framing. This file just lays out where the phenotype results we've finished
*set up* questions your RNA-seq can answer, so the two halves of the paper reinforce each other. Pick
up what's useful, drop what isn't. Every link below runs **both ways**: the expression data can confirm,
extend, **or revise** the phenotype-side hypothesis — none of these is a result the phenotype has already
decided for you.

Every phenotype number below traces to `output/tables/20_master_results.csv` (single source of truth);
the narrative is in `RESULTS.md`; the DE pipeline + genet-matching plan you already drafted is in
`analysis_proposal.md` (this file complements it — it's the "what to test, and why, given the
phenotype" layer).

## Your RNA-seq design (for reference)

144 libraries = **3 genets (A, C, D) × 2 wound (U, W) × 2 temp (28, 31 °C) × 4 tanks × 3 days**,
**wound-margin tissue only**, **Days 1, 3, 10** post-wound (Day 1 = 24 h). Fully balanced, plates
interleaved so no factor is confounded with plate/lane. That timecourse lines up almost perfectly
with the phenotype trajectory:

| Day | Phenotype state (from our data) | What the expression snapshot likely captures |
|---|---|---|
| **1** | wound fresh; tissue-healing just starting | early wound response (ECM remodeling, immune, proliferation) |
| **3** | coenosarc coverage advancing (≥90% by D5–7) | healing program at peak; regeneration not yet engaged |
| **10** | **new-corallite formation diverges by temp** | the healing→regeneration transition — where heat bites |

---

## Phenotype anchors → suggested expression contrasts

### 1. The central phenotype anchor: healing proceeds, regeneration stalls under heat
**Phenotype.** Tissue healing (coenosarc closure) is identical in both temperatures, but new-corallite
regeneration is delayed/suppressed under heat (interval-censored time ratio = 1.32, 95% CI 1.19–1.47;
first-observed Cox HR = 0.22, 95% CI 0.07–0.69; 67% of heated corals closed but never rebuilt skeleton
vs 0% at ambient). The divergence appears by **Day 10**.
**Hypothesis (expression can confirm, extend, or revise).** At the margin, the *healing* program
(Day 1–3) is induced similarly in both temps, while by **Day 10** heat specifically suppresses the
*biomineralization / skeletal-organic-matrix* program that drives regeneration.
**Question.** Temperature × timepoint (wounded margins): are there genes/modules induced at Day 10 in
28 °C but not 31 °C? The predicted *processes* (skeletal-organic-matrix/biomineralization and
corallite-patterning developmental signaling on the regeneration side; ECM/immune/cytoskeletal/
proliferation on the early-healing side) are hypotheses to confirm or overturn — the specific gene sets
you interrogate are your call.
**Refutes the decoupling story if:** heat also suppresses the Day 1–3 healing program — then the
impairment is uniform, not phase-specific. (Worth checking honestly; it's the crux, and the expression
data, not the phenotype, settles it.)

### 2. Genet-c resilience → candidate thermal-tolerance genes
**Phenotype.** Genet c defends photochemistry, pigmentation, and symbiont density far better than
a and d (multivariate displacement 1.06 vs 3.72/3.34) and is the most likely to regenerate at 31 °C.
The genet spread is **largest in unwounded corals**.
**Hypothesis (expression can confirm, extend, or revise).** Genet c has a distinct heat response — plausibly a smaller heat-induced
transcriptional shift and/or *constitutively frontloaded* stress-tolerance transcripts (the Barshis
et al. 2013 frontloading pattern).
**Contrast.** Genet (c vs a, d) × temperature, **prioritizing the unwounded margins** (where the
phenotype genet difference is sharpest) at each day. This is the most direct route to candidate
genes underlying heritable thermal tolerance in *A. pulchra*. A WGCNA module that tracks the
per-genet resilience score (`19_genet_resilience_summary.csv`) would be a clean way in.

### 3. Wounding homogenizes the genet response
**Phenotype.** The genet spread in heat sensitivity is large when unwounded (a = 0.99, d = 0.87,
c = 0.44) but compresses to near-uniform when wounded (a = 0.36, d = 0.26, c = 0.28).
**Hypothesis (expression can confirm, extend, or revise).** Wounding induces a shared metabolic/wound program that narrows genet differences.
**Contrast.** Genet × wound: expect the magnitude of the genet effect (and the genet × temp effect)
to be *smaller* in wounded than unwounded margins. If so, the molecular story matches the physiology.

### 4. The wound response is local (validates margin sampling)
**Phenotype.** Whole-colony physiology responds strongly to temperature but barely to wounding —
the wound response is concentrated locally at the site.
**Hypothesis (expression can confirm, extend, or revise).** At the *margin*, a strong wound main effect should be present at every day even
though bulk physiology shows almost none — confirming that the margin captures the local program
the colony-level metrics average away.
**Contrast.** Wound main effect at the margin (per day). A robust wound signature here is a good
positive control that the sampling worked.

### 5. Expect a chronic/constitutive signature, not an acute heat-shock spike
**Phenotype/design note.** 31 °C sits ~4.4 °C *below* the population's acute ED50 (Cunning 2024),
and wounding was applied **after 7 days at temperature** — so the acute heat-shock response has
likely subsided before Day 1. **Interpretation cue:** expect *sustained/constitutive* heat
signatures (altered baseline, frontloaded transcripts) rather than a classic transient HSP70 burst.
Distinguishing acute transcriptional plasticity from altered constitutive expression is the right
lens for this chronic-sublethal design.

### 6. Close the loop to Cunning's CBASS genets (high-value, near-publishable)
**Phenotype.** Our thicket labels A/C/D are **not yet matched** to Cunning et al. (2024)'s
genotyped CBASS genets (same Mahana population; ED50 34.3–36.6 °C, range 2.2 °C).
**Action.** Call SNPs from these RNA-seq libraries, match thickets A/C/D to Cunning's genotype
reference, then test whether **acute CBASS ED50 predicts our chronic, wound-context resilience
ranking (c > d > a)**. A positive result links a rapid acute assay to chronic regenerative capacity —
a clean, citable cross-method validation. (Details in `analysis_proposal.md`.)

---

## Possible entry points (no required order — sequence is yours)

These are angles the phenotype sets up, not a ranked work plan:

- **QC + the wound main effect at the margin** (Anchor 4) — a fast positive control that the margin
  sampling captured the local program.
- **Temperature × timepoint, wounded** (Anchor 1) — the D3→D10 shift; the question most central to the
  phenotype decoupling story.
- **Genet-c vs a/d × temperature, unwounded** (Anchor 2) — toward candidate thermal-tolerance processes.
- **Module-level analyses** (e.g. WGCNA) vs the phenotype axes (PC1 heat-stress axis; per-coral
  regeneration outcome; genet resilience score) — to tie modules to organismal outcomes.
- **SNP calling + Cunning genet-matching** (Anchor 6) — can run in parallel; high payoff.

Phenotype covariates you can pull per coral to correlate against expression (all in
`data/processed/` / `output/tables/`). **A tidy per-library covariate table is ready for you:**
`output/tables/31_rnaseq_phenotype_covariates.csv` — one row per RNA-seq library, joined by
`Fragment_ID` (= `id`), with the design factors (treatment/wound/genet/tank/day/plate, harmonized to
the phenotype coding), the per-fragment symbiont density, and the per-genet resilience covariates
(mean heat sensitivity, PCA displacement, rank). Built by `code/31_rnaseq_covariate_table.R`; drop in
more per-coral covariates there if you need them. A **raw, un-recoded** lookup
(`output/tables/31_rnaseq_library_lookup_raw.csv`) ships alongside it — library_id ↔ fragment_id ↔ the
original design fields verbatim from the plate layout (no 28C/31C, a/c/d, or no/yes recoding) — so you
set your own factor levels and reference categories.
