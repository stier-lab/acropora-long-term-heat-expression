# Candidate genes — literature + lab-data reference (optional)

> 🗂️ **Optional candidate-gene reference** · Updated 2026-06-26 · Index: [`README.md`](../../README.md) · companion to `analysis_proposal.md` and `gene_expression_integration_map.md` · machine-readable twin: [`candidate_genes_reference.csv`](candidate_genes_reference.csv).

**Read this as a resource to mine, not a prescribed panel.** The rest of `docs/for_shreya/`
deliberately frames the RNA-seq as *hypothesis-led, not candidate-gene-led* (see
`analysis_proposal.md`). This table names genes — so treat it as **prior context for
interpretation**, not a target list that should bias the differential-expression analysis. Keep
discovery genome-wide; use this to *interpret* and *cross-check* the output, and to seed a
confirmatory qRT-PCR panel if wanted.

**Every gene–citation pair below was verified by reading the actual paper** (close reads of the
PDFs; Cleves & Tinoco 2020 via NotebookLM). The verification corrected several errors carried over
from the planning docs — see the **Verification ledger** at the bottom. Source tags: a paper name =
the gene is named *in that paper*; **(lab)** = from the lab's own planning docs/panels (no external
primary confirmed); ⚠ = a flagged discrepancy.

**Relevance to LTH.** Heat (28 vs 31 °C) × wound × genet, biopsies D0/D1/D3/D10/D15. Headline:
**heat blocks regeneration, not closure**, with a genet gradient **C > D > A**. Tables split by axis;
the **heat × wound interaction** is the molecular LTH phenotype, and **per-genet baseline
("frontloading," Barshis 2013)** is the molecular version of C > D > A.

---

## Table 1 — Wound-healing / regeneration candidates

| Gene / pathway | Functional role | Axis | Why it matters for LTH | Source (verified) |
|---|---|---|---|---|
| **c-Fos** | Immediate-early transcription factor | Wound (acute) | Up in the *earliest* regeneration stage; cleanest acute marker — is the acute response intact under heat? | **Xu 2023** ✓ (up, earliest stages) |
| **JNK cascade** | Stress-activated MAP-kinase signalling | Wound (acute) | Implicated via GO enrichment; regulates c-Fos | **Xu 2023** ✓ (GO) |
| **NF-κB** | Innate-immune signalling | Wound / immunity | "NF-κB-inducing kinase" in GO enrichment | **Xu 2023** ✓ (GO) |
| **Wnt ligands + β-catenin** (APC, Wntless) | Axis patterning, regeneration | Regeneration | Core **regeneration** pathway — the phase heat impairs | **Xu 2023** ✓ (Wnt-1/2b/4/7b/8/8a/16/10; β-catenin). ⚠ *Wnt3* specifically is **not** in Xu — it is the NSF panel's pick (lab) |
| **FGF ligands / FGF-receptors / Sprouty** | Growth-factor signalling; wound-specific | Regeneration | Proliferative regrowth signal; expected suppressed if regeneration is blocked under heat | **Xu 2023** ✓ (FGF1/2/7/10/15/18; 3 FGF-Rs; sprouty) |
| **ADAMTS metalloproteases** | ECM remodelling | Wound | ECM turnover at the wound bed (ADAMTS18-like) | **Xu 2023** ✓ |
| **Galaxin** | Skeletal organic-matrix / biomineralization | Skeletal repair | New-skeleton/corallite formation — the LTH regeneration milestone (`new_corallites_on_tip`) | **Xu 2023** ✓ |
| **Ca²⁺-transport genes** (carbonic anhydrase, SLC4γ) | Calcification ion transport | Skeletal repair | Calcification under heat (LTH 38 % drop) | Xu 2023 ✓ for "calcium-ion transmembrane transport" (GO); *carbonic anhydrase / SLC4γ* named only in **(lab)** signaling notes |
| **MMPs (matrix metalloproteinases)** | Tissue/ECM remodelling | Wound | Remodelling during closure/regrowth | **(lab)** Common-Garden H10. ⚠ Xu reports **ADAMTS, not MMPs** |
| **TGF-β** | Wound repair / cell migration (vertebrate paradigm) | (see note) | Listed in the NSF panel — **but test, don't assume** | **(lab)** NSF Aim-2. ⚠ **Xu 2023 explicitly found TGF-β *not* involved in coral regeneration** |
| **PCNA** | Proliferation marker (Phoenix-effect test) | Regeneration | Reads out proliferation | **(lab)** NSF Aim-2 — ⚠ **not** named in Xu 2023; the planning-doc co-cite ("Tinoco 2023") does not exist. Substitute a verified proliferation marker before use |
| **Granular amoebocytes + melanin-bearing immune cells** *(cellular — no gene symbols)* | Four-phase wound response: degranulation → immune infiltration (6 h) → proliferation (24 h) → maturation/apoptosis (48 h) | Wound / immunity | Cellular framework + timeline for the healing phase | **Palmer 2011** ✓ (cells/timeline). ⚠ names **no genes**; phenoloxidase & phagocytosis are discussed only re *other* organisms, not shown in coral |

## Table 2 — Thermal-stress candidates

| Gene / family | Functional role | Direction under heat | Why it matters for LTH | Source (verified) |
|---|---|---|---|---|
| **HSP70 (+ HSP75)** | Core heat-shock chaperones | Up acutely; **higher constitutive baseline ("frontloaded") in tolerant** | Thermal axis + genet test: frontloading predicts resilient genet C carries higher baseline | **Barshis 2013** ✓ (HSP70/HSPA5); **Bay & Palumbi 2015** ✓ (HSP75). (HSP90 not named in either — omitted) |
| **Small HSPs** (HSPB1/Hsp23, Hsp16.2) | Small heat-shock chaperones | Up acute; lower in tolerant genotypes | Acute-stress reporters; contrast with frontloaded HSP70 | **Barshis 2013** ✓ |
| **HSF (heat-shock transcription factor)** | Master regulator of the HSP response | Activates HSPs under heat | Upstream switch; CRISPR knockout reduces thermal tolerance | **Cleves & Tinoco 2020** ✓ (NotebookLM-verified; no local PDF) |
| **Cu/Zn-SOD; catalase; peroxidasin** | Reactive-oxygen-species detox | Catalase up in bleaching; SOD & peroxidasin in tolerant | Sub-bleaching oxidative load — likely cost diverting resources from regeneration at 31 °C | **Barshis 2013** ✓ (Cu/Zn-SOD, peroxidasin); **Seneca & Palumbi 2010** ✓ (catalase, AmCat) |
| **TNFR; TRAF3 (TNF/apoptosis)** | Programmed cell death + immune signalling | TRAF3 among top up-regulated; frontloaded in tolerant | Cell-fate decisions under heat; links stress to tissue loss | **Barshis 2013** ✓ |
| **Ubiquitin–proteasome** (E3/RING-finger ligases; proteasome) | Protein quality control | Up acute; elevated pre-stress in tolerant | Proteostasis under thermal load | **Bay & Palumbi 2015** ✓ (ubiquitin/RING-finger); **Dixon 2015** ✓ (proteasome) |
| **Ion / Ca²⁺ transporters + oxidoreductases** | Ion homeostasis, redox balance | Higher pre-stress in heat-tolerant larvae (heritable) | Heritable heat-tolerance correlate; Ca²⁺ handling also links to calcification | **Dixon 2015** ✓ |
| **Chromoprotein / GFP-like pigments** | Photoprotection (screening pigments) | Up in bleaching | Molecular side of the LTH color-card paling | **Seneca & Palumbi 2010** ✓ (AmCh) |
| **C-type lectin** | Symbiont recognition / innate immunity | Up in bleaching | Host–symbiont signalling under heat (LTH symbiont loss) | **Seneca & Palumbi 2010** ✓ (AmCTL). (Barshis named *mannose-binding* lectin, down — different lectin) |
| **Signalling / TF network modules** (GPCR signalling, sequence-specific TFs up; ECM modules down) | Co-expressed regulatory modules (WGCNA) | Rapidly modulated; faster recovery in resilient species | Module-level resilience view | **Thomas 2019** ✓ (module-level; no individual gene names) |

**Two interpretation handles for LTH:**
- **Frontloading (Barshis 2013):** tolerance = *higher constitutive baseline* of stress genes, not
  bigger acute induction. → Compare **per-genet baselines** (C vs A) at D0 — the molecular test of C > D > A.
- **Heat-tolerance ↔ growth tradeoff (Cornwell 2021, phenotypic — no genes):** resilient colonies
  grow less / carry fewer symbionts. → Candidate mechanism for the LTH 38 % calcification reduction.

## Caveats

- **Prior context, not a target panel** — don't filter the DE analysis to these genes.
- **Mostly congeners.** Wound genes from *A. millepora* (Xu 2023) + the lab's *A. pulchra* near/far
  pilot; thermal genes from *Acropora* spp. and *Porites astreoides* (Dixon 2015). Orthology to the
  *A. pulchra* Conn 2025 genome must be established before any of these map to specific loci
  (deferred next step).
- **Pathways, not single loci** — many rows are families; the count matrix will have multiple
  paralogs per family.

## Verification ledger (2026-06-26)

Every cited paper was read and each gene checked against its text. Corrections applied:

| Item | Planning-doc claim | What the paper actually says | Action |
|---|---|---|---|
| **TGF-β** | regeneration candidate, cited to Xu 2023 | Xu 2023 states Wnt and FGF "**but not TGF-β**" pathways are involved | Moved to **(lab)** NSF-panel only; flagged ⚠ as contradicted by Xu |
| **Wnt3** | specific candidate (NSF panel, "Xu/Tinoco") | Xu names many Wnts (1/2b/4/7b/8/8a/16/10) + β-catenin, **not Wnt3** | Relabeled "Wnt ligands"; Wnt3 marked as the lab panel's pick |
| **PCNA** | proliferation marker (NSF panel, "Xu/Tinoco") | **not** in Xu 2023; co-cite "Tinoco 2023" does not exist | Marked **(lab)**, unverified; advise substituting a verified marker |
| **"Tinoco et al. 2023"** | cited for Wnt/FGF/PCNA/TGF-β | does not exist; only Cleves & **Tinoco 2020** (HSF, PNAS) | Removed; added Cleves & Tinoco 2020 → HSF |
| **"Hemond & Vollmer 2015"** | positional-confound caveat | correct paper is **Hemond, Kaluziak & Vollmer 2014** (BMC Genomics) | Citation corrected |
| **Palmer 2011 (PO, phagocytosis)** | wound machinery incl. phenoloxidase, phagocytosis | both discussed only re *other* organisms; coral findings are amoebocytes + melanin cells; **no gene symbols** | Row reframed as cellular/timeline; PO/phagocytosis removed |
| **MMPs** | wound remodelling (would-be Xu) | Xu reports **ADAMTS**, not MMPs | MMPs kept as **(lab)** Common-Garden; ADAMTS credited to Xu |
| **peroxiredoxin** | antioxidant (Barshis) | Barshis names **peroxidasin** (a distinct ECM peroxidase) | Corrected to peroxidasin |
| **HSP90** | implied HSP candidate | not named in Barshis 2013 or Bay & Palumbi 2015 | Omitted (kept HSP70/HSP75 only) |
| **C-type lectin → Barshis** | lectin (Barshis + Seneca) | Barshis named *mannose-binding* lectin (down); C-type lectin is Seneca's (AmCTL, up) | C-type lectin credited to Seneca only |

## References (verified)

| Citation | Title (short) | Venue | DOI | Verification |
|---|---|---|---|---|
| Xu et al. 2023 | Wound healing & regeneration in *Acropora millepora* | Front. Ecol. Evol. | 10.3389/fevo.2022.979278 | ✅ PDF read |
| Palmer et al. 2011 | Corals use similar immune cells & wound-healing processes | PLoS ONE | 10.1371/journal.pone.0023992 | ✅ PDF read |
| Barshis et al. 2013 | Genomic basis for coral resilience (frontloading) | PNAS | 10.1073/pnas.1210224110 | ✅ PDF read |
| Seneca & Palumbi 2010 | Gene expression in a coral undergoing natural bleaching | Mar. Biotechnol. | 10.1007/s10126-009-9247-5 | ✅ PDF read |
| Bay & Palumbi 2015 | Rapid acclimation via transcriptome change | Genome Biol. Evol. | 10.1093/gbe/evv085 | ✅ PDF read |
| Dixon et al. 2015 | Genomic determinants of heat tolerance across latitudes | Science | 10.1126/science.1261224 | ✅ PDF read |
| Thomas et al. 2019 | Transcriptomic resilience + symbiont shuffling | Mol. Ecol. | 10.1111/mec.15143 | ✅ PDF read |
| Cornwell et al. 2021 | Heat tolerance × symbiont load → growth tradeoff (phenotypic) | eLife | 10.7554/eLife.64790 | ✅ PDF read (no genes) |
| Cleves & Tinoco et al. 2020 | CRISPR knockout of a heat-shock TF reduces thermal tolerance | PNAS | 10.1073/pnas.1920779117 | ✅ NotebookLM (no local PDF) |
| Hemond, Kaluziak & Vollmer 2014 | Genetics of colony form & function in Caribbean *Acropora* | BMC Genomics | 10.1186/1471-2164-15-1133 | ✅ NotebookLM |

*Lab-internal sources (not external citations): the 2023 *A. pulchra* near/far DEG pilot (76 near /
165 far DEGs at Day 7; 36 TRIzol samples at −80 °C, unprocessed); the NSF BIO/OCE Aim-2 candidate
panel; the Mar-15-2026 Keck signaling notes.*
