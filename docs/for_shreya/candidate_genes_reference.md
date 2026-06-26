# Candidate genes — literature + lab-data reference (optional)

> 🗂️ **Optional candidate-gene reference** · Updated 2026-06-26 · Index: [`README.md`](../../README.md) · companion to `analysis_proposal.md` and `gene_expression_integration_map.md` · machine-readable twin: [`candidate_genes_reference.csv`](candidate_genes_reference.csv).

**Read this as a resource to mine, not a prescribed panel.** The rest of `docs/for_shreya/`
deliberately frames the RNA-seq as *hypothesis-led, not candidate-gene-led* (see
`analysis_proposal.md`: "we deliberately do not name candidate gene symbols, so the gene discovery
stays unbiased"). This table does name genes — so treat it as **prior context for interpretation**
(what the wound and heat literature, plus the lab's own 2023 data, flag as likely players), **not** a
target list that should bias the differential-expression analysis. Keep discovery genome-wide; use
this to *interpret* and *cross-check* the output, and to seed a confirmatory qRT-PCR panel if wanted.

**Relevance to LTH.** Shreya's LTH design is heat (28 vs 31 °C) × wound × genet, biopsies at
D0/D1/D3/D10/D15. The headline phenotype is that **heat blocks regeneration, not wound closure**
(`RESULTS.md`), with a **genet resilience gradient C > D > A**. So the two tables below split the
candidates by axis — *wound healing / regeneration* (do these fire on schedule under heat?) and
*thermal stress* (are they diverting resources, and are they "frontloaded" in the resilient genet?).
The **heat × wound interaction** is the molecular version of the LTH phenotype; the **per-genet
baseline ("frontloading," Barshis 2013)** is the molecular version of the C > D > A gradient.

**Provenance & verification.** Consolidated from the lab's planning docs (Mar-15-2026 Keck signaling
notes; Porites Common-Garden design H10–H12; NSF BIO/OCE Aim-2 panel) and grounded in the cited
papers. Thermal-gene rows were extracted directly from the repo's thermal-stress PDFs. All citations
below were checked (local PDFs, repo DOIs, or NotebookLM); see **References**. Two citations in the
planning docs were wrong and have been corrected here (see References footnotes).

---

## Table 1 — Wound-healing / regeneration candidates

| Gene / pathway | Functional role | Axis | Why it matters for LTH | Source |
|---|---|---|---|---|
| **c-Fos** | Immediate-early transcription factor; peaks ~6 h post-wound | Wound (acute) | Cleanest acute wound marker (no known positional confound in *Acropora*); tests whether the acute response itself is intact under heat | Xu 2023 |
| **JNK** | Stress-activated MAP-kinase; early wound signaling | Wound (acute) | Part of the immediate-early injury cascade; pairs with c-Fos | Xu 2023 |
| **NF-κB** | Innate-immune transcription factor; up early | Wound / immunity | Immune phase of healing; may be heat-sensitive | Signaling notes |
| **Wnt / Wnt3 / β-catenin** | Axis patterning, tissue regeneration | Regeneration | Core **regeneration** pathway — the phase heat impairs in LTH; proliferative/Phoenix marker | Xu 2023; NSF Aim-2 |
| **FGF ligands / FGF-receptor / Sprouty** | Growth-factor signaling; wound-specific | Regeneration | Proliferative regrowth signal; expected suppressed under heat if regeneration is blocked | Xu 2023; NSF Aim-2 |
| **TGF-β** | Wound repair, cell migration | Wound → regeneration | Bridges closure and tissue remodeling | Xu 2023; NSF Aim-2 |
| **MMPs (matrix metalloproteinases)** | Tissue/ECM remodeling | Wound | Remodeling during closure and regrowth | Common-Garden H10 |
| **ADAMTS metalloproteases** | ECM remodeling; wound-specific | Wound | ECM turnover at the wound bed | Signaling notes |
| **PCNA** | Proliferating-cell nuclear antigen | Regeneration | Reads out cell proliferation; key for the "regeneration from within the skeleton" (Phoenix) test | NSF Aim-2 |
| **Carbonic anhydrase, galaxin, SLC4γ** | Calcification / biomineralization | Skeletal repair | New-skeleton/corallite formation — the exact LTH regeneration milestone (`new_corallites_on_tip`) | Signaling notes |
| **Immune cells / wound-healing machinery** | Granular amoebocytes, melanin-synthesis (PO), phagocytosis | Wound / immunity | Corals use conserved immune + wound-healing processes; baseline for the healing phase | Palmer 2011 |

## Table 2 — Thermal-stress candidates (grounded in repo PDFs)

| Gene / family | Functional role | Direction under heat | Why it matters for LTH | Source |
|---|---|---|---|---|
| **HSP70 / HSP75 / HSP90** | Molecular chaperones (core heat-shock response) | Up acutely; **higher constitutive baseline ("frontloaded") in heat-tolerant corals** | The thermal axis, and the genet test: frontloading predicts the resilient genet (C) carries higher baseline HSPs | Barshis 2013; Bay & Palumbi 2015 |
| **Small HSPs (HSPB1/Hsp23, Hsp16.2)** | Small heat-shock chaperones | Up under acute heat; lower in tolerant genotypes | Acute-stress reporters; contrast with frontloaded HSP70 | Barshis 2013 |
| **HSF (heat-shock transcription factor)** | Master regulator of the HSP response | Activates HSPs under heat | Upstream switch for the whole HSP program; CRISPR knockout reduces thermal tolerance | Cleves & Tinoco 2020 |
| **Antioxidants: Cu/Zn-SOD, catalase, peroxiredoxin/peroxidasin** | Reactive-oxygen-species detoxification | Catalase up in bleaching; peroxidasin frontloaded in tolerant | Sub-bleaching oxidative load — the likely cost diverting resources from regeneration at 31 °C | Barshis 2013; Seneca & Palumbi 2010 |
| **Apoptosis / TNF pathway: TNFR, TRAF3** | Programmed cell death + immune signaling | TRAF3 among top up-regulated under heat; frontloaded in tolerant | Cell-fate decisions under heat; ties stress to tissue loss | Barshis 2013 |
| **Ubiquitin–proteasome (E3/RING-finger ligases, proteasome)** | Protein quality control / turnover of damaged protein | Up under acute heat; elevated pre-stress in tolerant | Proteostasis under thermal load | Bay & Palumbi 2015; Dixon 2015 |
| **Ion / Ca²⁺ transporters + oxidoreductases** | Ion homeostasis, redox balance | Higher pre-stress expression in heat-tolerant larvae (heritable) | Heritable heat-tolerance correlate; Ca²⁺ handling also links to calcification | Dixon 2015 |
| **Chromoprotein / GFP-like pigments** | Photoprotection (screening) | Up in bleaching | Pigmentation response — the molecular side of the LTH color-card paling | Seneca & Palumbi 2010 |
| **C-type lectins** | Symbiont recognition / innate immunity | Up in bleaching | Host–symbiont signaling under heat (LTH symbiont-density loss) | Seneca & Palumbi 2010; Barshis 2013 |
| **Signaling / TF network modules** (GPCR signaling, sequence-specific TFs up; ECM modules down) | Co-expressed regulatory modules | Rapidly modulated under heat; faster recovery in resilient species | Module-level (WGCNA) view; resilience = rapid signaling modulation | Thomas 2019 |

**Two interpretation handles for LTH:**
- **Frontloading (Barshis 2013):** thermal tolerance shows up as *higher constitutive baseline* of
  stress genes, not bigger acute induction. → Compare **per-genet baselines** (C vs A) at D0, not
  just heat-induced fold-changes — this is the molecular test of the C > D > A gradient.
- **Heat-tolerance ↔ growth tradeoff (Cornwell 2021):** thermally resilient colonies grow less /
  carry fewer symbionts. → A candidate mechanism for the LTH 38 % calcification reduction under heat.

## Caveats

- **The list is prior context, not a target panel** (see the framing note up top). Don't filter the
  DE analysis to these genes.
- **Mostly congeners.** Wound genes are from *Acropora millepora* (Xu 2023) and the lab's *A. pulchra*
  near/far pilot; thermal genes are from *Acropora* spp. and *Porites astreoides* (Dixon 2015).
  Orthology to the *A. pulchra* Conn 2025 genome must be established before these map to specific loci
  (a planned next step — gene-ID mapping was deferred).
- **Pathways, not single loci.** Many rows are families/pathways; the actual count matrix will have
  multiple paralogs per family.
- **Positional confound** (only if any spatial contrast is attempted): in *Acropora*, tip vs. base
  tissue differs by ~2,200 transcripts constitutively (Hemond et al. 2014) — minor for LTH's
  whole-fragment biopsies.

## References (verified)

| Citation | Title (short) | Venue | DOI | Verification |
|---|---|---|---|---|
| Xu et al. 2023 | Wound healing & regeneration in *Acropora millepora* | Front. Ecol. Evol. | 10.3389/fevo.2022.979278 | ✅ local PDF |
| Palmer et al. 2011 | Corals use similar immune cells & wound-healing processes | PLoS ONE | 10.1371/journal.pone.0023992 | ✅ local PDF |
| Barshis et al. 2013 | Genomic basis for coral resilience (frontloading) | PNAS | 10.1073/pnas.1210224110 | ✅ repo |
| Seneca & Palumbi 2010 | Gene expression in a coral undergoing natural bleaching | Mar. Biotechnol. | 10.1007/s10126-009-9247-5 | ✅ repo |
| Bay & Palumbi 2015 | Rapid acclimation via transcriptome change | Genome Biol. Evol. | 10.1093/gbe/evv085 | ✅ repo |
| Dixon et al. 2015 | Genomic determinants of heat tolerance across latitudes | Science | 10.1126/science.1261224 | ✅ repo |
| Thomas et al. 2019 | Transcriptomic resilience + symbiont shuffling | Mol. Ecol. | 10.1111/mec.15143 | ✅ repo |
| Cornwell et al. 2021 | Heat tolerance × symbiont load → growth tradeoff | eLife | 10.7554/eLife.64790 | ✅ repo |
| Cleves & Tinoco et al. 2020 | CRISPR knockout of a heat-shock TF reduces thermal tolerance | PNAS | 10.1073/pnas.1920779117 | ✅ NotebookLM |
| Hemond, Kaluziak & Vollmer 2014 | Genetics of colony form & function in Caribbean *Acropora* | BMC Genomics | 10.1186/1471-2164-15-1133 | ✅ NotebookLM |

*Lab-internal sources (not citations): the 2023 *A. pulchra* near/far DEG pilot (76 near / 165 far
DEGs at Day 7; 36 TRIzol samples at −80 °C, unprocessed); the NSF BIO/OCE Aim-2 candidate panel; and
the Mar-15-2026 Keck signaling meeting notes.*

> **Corrections made during verification (2026-06-26):** the planning docs cited "Tinoco et al. 2023"
> (does not exist — the only Tinoco paper is Cleves & Tinoco 2020, on a heat-shock factor) and
> "Hemond & Vollmer 2015" (correct paper is Hemond, Kaluziak & Vollmer **2014**). Both fixed above.
