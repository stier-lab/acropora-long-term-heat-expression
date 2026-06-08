# Citation audit — *Acropora pulchra* library (backward snowball)

> **ACQUISITION STATUS (2026-06-07):** Acted on this audit — **37 of the candidate papers acquired**
> (library 37 → **74 PDFs**). All 16 Tier-A (lab-library copies) + 21 of 44 Tier-B/C (PMC/EuropePMC,
> Springer `content/pdf`, green-OA repositories: JCU, ANU). **23 still missing** — genuinely paywalled
> Wiley/Elsevier/Science/Nature/Annual-Reviews/Cell/JSTOR with no open-access route and an
> unauthenticated browser. Of those, **6 are open-access but Cloudflare-blocks headless curl**
> (Anthony 2007 & 2009, Omori 2001, Voolstra 2020, Evensen 2023, Thomas 2019) — retrievable with a
> normal interactive browser download or an authenticated UCSB session. Highest-value still-missing:
> Grottoli 2006 (`h`, Nature-closed), Voolstra 2020 (CBASS, OA-blocked), LaJeunesse 2018 (Cell),
> Anthony 2009 (energetics, OA-blocked), Dixon 2015 & Palumbi 2014 (Science-closed). Yap 1981
> (*A. pulchra* thesis) has no DOI/online copy — would need the U. Philippines library.



**Date:** 2026-06-07. **Method:** For every paper in the library, we mined its reference list
two ways: (1) extracted embedded DOIs from all 37 PDFs and counted cross-paper citation
frequency; (2) queried the "Acropora pulchra" NotebookLM notebook (26 sources) thematically
(wound healing, thermal tolerance/expression, symbiont/microbiome, energetics, reproduction)
and pulled the exact reference-list entries. Candidates were de-duplicated against the existing
library + LITERATURE.md, and every recommendation's DOI was verified via Crossref/OpenAlex
(none guessed; the few stub hits were corrected by hand and noted).

**Scope filter:** the corpus cites a large neighborhood (~250 distinct works). We kept only papers
relevant to *this* library's purpose — **A. pulchra biology, the wound-healing model's mechanisms
and parameters, the LTH expression-by-temperature project, and the Mo'orea field site.** We
deliberately **excluded** the ~100-paper coral-microbiome taxonomy/disease/probiotics neighborhood
(Endozoicomonas genomics, Red Sea 16S surveys, white-band/black-band disease, marine probiotics,
sponge/diatom molecular papers) — that belongs in a microbiome-specific notebook, not here. Lock
2025, Miller 2024, and de la Vega 2023 already anchor the microbiome angle for *A. pulchra*.

---

## Tier A — already in the `coral-regen-review` lab library → copy in now (zero cost, highest relevance)

These are the **wound-healing × temperature core of the model** and are already PDFs in
`~/coral-regen-review/literature/pdfs/`. They are the single most valuable gap because the model
*is* a wound-healing model and these directly inform `kc_mult`, the temperature×healing coupling,
and the energy-allocation tier.

| Paper | DOI | Why it matters → model link |
|---|---|---|
| **Traylor-Knowles 2016** — distinctive wound-healing in *Pocillopora damicornis* & *Acropora hyacinthus* in two temperature regimes | 10.1186/s12862-016-0791-0 | ⭐ wound healing × temperature × the model's two *existing* species |
| **Bonesso 2017** — elevated SST below bleaching threshold impairs recovery/regeneration after injury | 10.7717/peerj.3719 | sub-bleaching heat slows healing → temp-coupling of `kc_mult` |
| **Burmester 2017** — temperature & symbiosis affect lesion recovery | 10.3354/meps12114 | symbiotic state × temp × healing |
| **Dias 2018** — mortality, growth, regeneration after fragmentation under thermal stress | 10.1016/j.seares.2018.08.008 | fragmentation × heat → `m_d`, healing |
| **Edmunds & Lenihan 2010** — sub-lethal damage to juvenile *Porites* under temperature/flow | 10.1007/s00227-009-1372-1 | *Porites* (model species) damage × temp |
| **Rinkevich 1996** — do reproduction and regeneration compete for energy allocation? | 10.3354/meps143297 | ⭐ energy trade-off → model Tier-3 energy budget |
| **Henry & Hart 2005** — regeneration from injury & resource allocation (review) | 10.1002/iroh.200410759 | energy-allocation framing for healing |
| **Meesters 1993** — coral bleaching reduces tissue regeneration potential | 10.3354/meps096189 | ⭐ bleaching × regeneration (we only had Meesters 1997) |
| **Meesters 1994** — damage & regeneration link to growth (*Montastraea*) | 10.3354/meps112119 | growth↔regeneration trade-off |
| **Lirman 2000b** — lesion regeneration in branching *Acropora palmata* | 10.3354/meps197209 | lesion-area recovery (companion to Lirman 2000a we have) |
| **DeFilippo 2016** — surface lesion recovery in *Astrangia poculata* | 10.1016/j.jembe.2016.03.016 | lesion-recovery kinetics |
| **Counsell 2019** — colony size & depth affect wound repair in a branching coral | 10.1007/s00338-019-01807-7 | size-dependent healing |
| **Rice 2019** — nitrogen sources speed recovery from corallivory + microbiome | 10.7717/peerj.8056 | corallivory recovery (browsing axis) |
| **Bak 1983** — neoplasia, regeneration & growth in *Acropora palmata* | 10.1007/BF00395810 | foundational *Acropora* regeneration |
| **Bak & Steward-Van Es 1980** — regeneration of superficial damage (*Agaricia*, *Porites*) | (Bull Mar Sci, no DOI) | foundational regeneration |

*(Bonus also in lab: Burmester 2018 — autotrophy vs heterotrophy in lesion recovery.)*

---

## Tier B — must-acquire (directly on model mechanism / A. pulchra / Mo'orea / LTH expression)

### B1. Wound / mechanical damage / microfragmentation (not in lab library)
| Paper | DOI | Model link |
|---|---|---|
| Van de Water 2015 — elevated temps & coral immune response after physical damage | 10.1007/s10750-015-2243-z | immune/healing × temp |
| Edmunds & Yarid 2017 — ocean acidification & wound repair in *Porites* | 10.1016/j.jembe.2016.10.001 | wound repair (model species) |
| Marshall 2000 — skeletal damage vs colony morphology | 10.3354/meps200177 | `m_d` (breakage resistance) |
| Tunnicliffe 1981 — breakage & propagation of *Acropora cervicornis* | 10.1073/pnas.78.4.2427 | fragmentation → `m_d` |
| Madin & Connolly 2006 — ecological consequences of hydrodynamic disturbance | 10.1038/nature05328 | mechanical vulnerability → `m_d` |
| Lock 2022 — Ca-homeostasis disruption drives rapid growth after microfragmentation (*Porites lobata*) | 10.1002/ece3.9345 | ⭐ microfragmentation → model Tier-3 microfrag |
| Page 2018 — microfragmenting for restoration of slow-growing massive corals | 10.1016/j.ecoleng.2018.08.017 | microfragmentation |
| Soong & Chen 2003 — regeneration & growth of *Acropora* fragments in a nursery | 10.1046/j.1526-100x.2003.00100.x | fragment regeneration/growth |

### B2. Energetics / reserves / trade-offs (model Tier-3 energy budget `B,R,L`)
| Paper | DOI | Model link |
|---|---|---|
| Anthony 2009 — energetics approach to predicting mortality risk from stress | 10.1111/j.1365-2435.2008.01531.x | ⭐ bleaching energetics → energy-budget mortality |
| Anthony 2007 — bleaching, energetics & mortality risk (temp/light/sediment) | 10.4319/lo.2007.52.2.0716 | energy reserves under stress |
| Grottoli 2006 — heterotrophic plasticity & resilience in bleached corals | 10.1038/nature04565 | ⭐ `h` (heterotrophic buffer), foundational |
| Grottoli 2014 — cumulative annual bleaching turns winners into losers | 10.1111/gcb.12658 | ⭐ repeated-stress reserve depletion |
| Rodrigues & Grottoli 2007 — energy reserves/metabolism as recovery indicators | 10.4319/lo.2007.52.5.1874 | ⭐ `E_max` (reserve dynamics) |
| Levas 2013 — physiological traits of bleaching & recovery in *Porites lobata* | 10.1371/journal.pone.0063267 | *Porites* recovery physiology |

### B3. A. pulchra-specific / Mo'orea–French Polynesia field site
| Paper | DOI | Why |
|---|---|---|
| **Yap 1981** — growth, regeneration & transplantation of *A. pulchra* (Bolinao) | (M.Sc. thesis, UP Quezon City; no DOI) | ⭐⭐ rare *A. pulchra* regeneration data — direct param source; thesis, hard to get |
| Carroll 2006 — sexual reproduction of *Acropora* at Moorea, FP | 10.1007/s00338-005-0057-6 | ⭐⭐ *Acropora* reproduction AT Mo'orea → `B_mat` + field site |
| Rouzé 2017 — Symbiodiniaceae of Mo'orea scleractinians | 10.7717/peerj.2856 | Mo'orea symbiont identity (A. pulchra context) |
| Rouzé 2019 — quantitative Symbiodiniaceae signature of Mo'orea colonies | 10.1038/s41598-019-44017-5 | Mo'orea symbiont dynamics |
| Roff 2014 — *Porites* "Phoenix effect" recovery, Rangiroa (FP) | 10.1007/s00227-014-2426-6 | dramatic FP recovery → recovery dynamics |

### B4. Reproduction / fecundity / size-at-maturity (`B_mat`)
| Paper | DOI | Why |
|---|---|---|
| Wallace 1985 — reproduction/recruitment/fragmentation in nine *Acropora* | 10.1007/BF00392585 | ⭐ *Acropora* fecundity × size → `B_mat` |
| Baird 2009 — reproductive biology of scleractinian corals (review) | 10.1146/annurev.ecolsys.110308.120220 | reproduction framing |
| Omori 1999/2001 — drop in *Acropora* fertilization after bleaching | 10.4319/lo.2001.46.3.0704 | bleaching × reproduction trade-off |

### B5. Symbiont identity / shuffling / recovery threshold (`iZ`)
| Paper | DOI | Why |
|---|---|---|
| Baker 2003 — flexibility & specificity in coral–algal symbiosis (review) | 10.1146/annurev.ecolsys.34.011802.132417 | foundational symbiosis framing |
| Berkelmans & van Oppen 2006 — zooxanthellae role in thermal tolerance ("nugget of hope") | 10.1098/rspb.2006.3567 | ⭐ symbiont-mediated thermal tolerance (co-cited 3×) |
| Jones 2008 — endosymbiont community change after bleaching | 10.1098/rspb.2008.0069 | ⭐ symbiont shuffling → `iZ` |
| Silverstein 2017 — clade D *Symbiodinium* persist at temp extremes | 10.1242/jeb.148239 | thermally tolerant symbiont retention |
| LaJeunesse 2018 — systematic revision of Symbiodiniaceae | 10.1016/j.cub.2018.07.008 | essential symbiont taxonomy (co-cited 2×) |

### B6. Thermal-tolerance & gene-expression backbone (LTH "expression by temperature")
| Paper | DOI | Why |
|---|---|---|
| Voolstra 2020 — standardized short-term acute heat-stress assays (CBASS) | 10.1111/gcb.15148 | ⭐⭐ CBASS standard (co-cited 3×); LTH assay backbone |
| Evensen 2023 — CBASS low-cost portable system | 10.1002/lom3.10555 | CBASS methods |
| Cunning 2021 — heat-tolerance census, Florida staghorn *Acropora* | 10.1098/rspb.2021.1613 | acute-tolerance variation (companion to Cunning 2024 we have) |
| Barshis 2013 — genomic basis for coral resilience to climate change | 10.1073/pnas.1210224110 | ⭐⭐ expression basis of resilience (co-cited 2×) |
| Bay & Palumbi 2015 — rapid acclimation via transcriptome change | 10.1093/gbe/evv085 | expression plasticity |
| Dixon 2015 — genomic determinants of heat tolerance across latitudes | 10.1126/science.1261224 | heritable heat tolerance |
| Palumbi 2014 — mechanisms of reef-coral resistance to climate change | 10.1126/science.1251336 | ⭐ *A. hyacinthus* expression/acclimation |
| Seneca & Palumbi 2010 — gene expression in a coral undergoing natural bleaching | 10.1007/s10126-009-9247-5 | expression × bleaching |
| Parkinson 2018 — transcriptional variation challenges thermal biomarkers | 10.1111/mec.14517 | expression variation (verify DOI) |
| Thomas 2019 — transcriptomic resilience, symbiont shuffling, recurrent bleaching | 10.1111/mec.15143 | expression + shuffling |
| Cornwell 2021 — heat tolerance & symbiont load → growth trade-offs (*A. hyacinthus*) | 10.7554/eLife.64790 | ⭐ growth–tolerance trade-off |
| Fitt 2001 — interpretation of thermal tolerance limits & thresholds | 10.1007/s003380100146 | foundational bleaching thresholds |
| Warner, Fitt & Schmidt 1996 — elevated temp & photosynthetic efficiency of zooxanthellae | 10.1046/j.1365-3040.1996.d01-12.x | foundational photoinhibition (verify DOI) |

---

## Tier C — broad climate/bleaching context (optional; heavily co-cited, intro/discussion framing)

| Paper | DOI | Note |
|---|---|---|
| Hoegh-Guldberg 2007 — coral reefs under rapid climate change & OA | 10.1126/science.1152509 | co-cited 3×; classic intro cite |
| Hughes 2017 — global warming & recurrent mass bleaching of corals | 10.1038/nature21707 | co-cited 2× |
| Hughes 2017 — coral reefs in the Anthropocene | 10.1038/nature22901 | co-cited 2× |
| De'ath 2012 — 27-year decline of GBR coral cover | 10.1073/pnas.1208909109 | co-cited 2× |
| Buddemeier & Fautin 1993 — coral bleaching as an adaptive mechanism | 10.2307/1312064 | Adaptive Bleaching Hypothesis; co-cited 2× |

---

## Already tracked (not "missing") — flagged for context
These surfaced in the citation mining but are **already in LITERATURE.md** (acquired or known-paywalled):
Tanaka 2006 (organic-N translocation, *A. pulchra*), Houlbrèque & Ferrier-Pagès 2009, Schoepf 2013,
Bay 2016, Denis 2011/2013, Madin 2014, Conn 2025 (the bioRxiv preprint 10.1101/2025.03.27.645822
duplicates our Gigabyte version). The MEPS `meps09060` "Denis 2011" co-citation hit = our Denis 2011.

## Excluded from scope (intentionally)
~100 coral-microbiome papers (Endozoicomonas genomics, Red Sea/Thailand 16S surveys, white-/black-band
disease, probiotics/BMC, sponge & diatom molecular work) and broad symbiont-ecology surveys not tied
to *A. pulchra*, the model, or Mo'orea. Anchored already by Lock 2025 / Miller 2024 / de la Vega 2023.
