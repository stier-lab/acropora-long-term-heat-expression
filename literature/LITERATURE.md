# Literature — *Acropora pulchra* (LTH "expression by temperature" + wound-healing model)

> 🗂️ **Literature overview** · Updated 2026-06-12 · Index: [`README.md`](../README.md) · PDFs in `literature/pdfs/`; refs in `manuscript/references.bib`.

Curated library for the LTH heat × wounding project and for parameterizing a third
species (*A. pulchra*) in the coral wound-healing model (`~/coral-wound-healing-model`).
**Every PDF is mapped to its *A. pulchra* relevance + model/project connection in
[`LIBRARY_MAP.md`](LIBRARY_MAP.md)** (98 PDFs total — see `LIBRARY_MAP.md` for the Core /
includes-pulchra / grounding breakdown; sortable companion [`library_map.csv`](library_map.csv)).
Companion to the **"Acropora pulchra" NotebookLM notebook** (`f9a97210-…`, query for full text)
and the lab regeneration library `~/coral-regen-review/literature/` (191 PDFs).

> **Manuscript citations** are tracked separately: every work cited in `manuscript/Manuscript_LTH.md`
> has a Crossref-verified entry in **`manuscript/references.bib`** and a cite→DOI→PDF mapping in
> **[`CITATIONS_INDEX.md`](CITATIONS_INDEX.md)** (2026-06-09: +5 PDFs acquired — Jokiel & Coles
> 1977/1990, Cunning 2024, Davies 1989, Siebeck 2006; now **98 PDFs**). Six cited papers remain
> paywalled and are listed there for acquisition via authenticated Chrome.

**Status (2026-06-09):** discovery complete + **PDF acquisition done (now 98 PDFs in `pdfs/`)**
+ **NotebookLM ingestion done (notebook 26 sources; new round not yet ingested)**. A backward
citation-chain audit (`CITATION_AUDIT.md`) added 37 more papers — the wound-healing×temperature
core (Traylor-Knowles, Bonesso, Burmester, Dias, Meesters 1993/94, Rinkevich…), energetics, the
expression-by-temperature backbone (Barshis, Bay&Palumbi, Cornwell, Cunning 2021…), Mo'orea
symbionts/reproduction (Carroll, Rouzé), and symbiont-shuffling (Berkelmans, Jones). See
`CITATION_AUDIT.md` for the full audit. After a VPN-authenticated CDP round (in-page `fetch()`),
only **4 audit papers remain unobtained** (Anthony 2009, Edmunds & Yarid 2017, Page 2018,
Buddemeier 1993) — all with in-library substitutes; plus the original handful of A.-pulchra-index
papers that are paywalled with no OA route (mostly [NB] in the notebook or congener-substitutable). DOIs were tool-confirmed
(OpenAlex/Crossref); none guessed. PDFs are named `AuthorYear_ShortTitle.pdf`.

**Acquisition routes that worked (2026-06-07):** Springer `link.springer.com/content/pdf/{DOI}.pdf`
via plain curl (even for "paywalled" Mar Biol / Coral Reefs); Nature/PLOS/Frontiers/RSOS direct
PDF URLs; **EuropePMC `?pdf=render`** for OA papers that Cloudflare-blocked bare curl (PeerJ, ASM,
iScience, Wiley-with-PMC-deposit). int-res.com (MEPS), bepress theses, and `download-papers` CDP
all failed (JS anti-bot / unauthenticated session) — EuropePMC recovered every OA one of those.

Acquisition legend: **[NB]** already in the A. pulchra notebook · **[LAB]** already in
`~/coral-regen-review/literature/pdfs` (copy in) · **[OA]** open-access (curl) ·
**[PW]** paywalled (Chrome CDP) · relevance 1–5 (5 = core).

> ⚠ **Species-attribution corrections (2026-06-08).** A read-through against the source PDFs found
> four papers over-attributed to *A. pulchra*; rows below are corrected and flagged `⚠`:
> **Zhang 2025** = the macroalga *Halimeda macroloba* (NOT a coral) — so *A. pulchra* proteomic/
> metabolomic thermal acclimation is actually unstudied; **Almeida 2024** = *Pocillopora* cf.
> *damicornis* (Kenya); **Denis 2024** = *A. spathulata* (GBR, *A. pulchra* not in the study);
> **Raymundo 2025** = *A. aspera* outplants (Guam, *A. pulchra* only the background congener).
> These remain useful *grounding/framing* sources but are not *A. pulchra*-subject papers.

---

## 1. *A. pulchra*-specific papers

### Genomics / gene expression / proteomics (core for "expression by temperature")
| Paper | DOI | Acq | Rel | Note |
|---|---|---|---|---|
| Conn et al. 2025 — genome assembly & annotation (Mo'orea) | 10.46471/gigabyte.153 | [NB][OA] | 5 | reference genome |
| ⚠ Zhang et al. 2025 — protein & metabolite acclimation to temp variability | 10.3389/fmars.2025.1543591 | [OA] | 2 | **NOT A. pulchra — macroalga *Halimeda macroloba***; calcifier thermal-acclimation framing only |
| Anthony et al. 2023 — Symbiodiniaceae cellular plasticity / phenotypic change | 10.3389/fevo.2023.1288596 | [OA] | 5 | symbiont expression plasticity |
| de la Vega et al. 2023 — *Endozoicomonas* GU-1 genome from *A. pulchra* | 10.1128/mra.01355-22 | [OA] | 4 | microbiome resource |

### Thermal tolerance / bleaching / CBASS
| Paper | DOI | Acq | Rel | Note |
|---|---|---|---|---|
| Cunning et al. 2024 — rapid acute heat-tolerance (CBASS), *A. pulchra* | (DOI unconfirmed) | [NB] | 5 | Fv/Fm ED50 34.3–36.6 °C; Mo'orea; used in Pfab 2024 |
| Berg et al. 2020 — persistent photosystem damage fluorescence (A. cf. pulchra) | 10.1080/17451000.2021.1875245 | [NB][PW] | 5 | 58% mortality d34; recovery 3× stress |
| Shaw et al. 2016 — intraspecific variability, warming × acidification | 10.1007/s00227-016-2986-8 | [NB][PW] | 5 | genotype variation |
| Grottoli et al. 2021 — comparability among bleaching experiments | 10.1002/eap.2262 | [NB] | 4 | methods (Ecol. Appl. 31:e02262, 2021) |
| ⚠ Denis et al. 2024 — individual thermal-tolerance traits across the GBR | 10.1098/rspb.2024.0587 (preprint 10.1101/2024.01.28.576773) | [OA preprint] | 4 | **A. spathulata, GBR (not A. pulchra)**; CBASS intraspecific-variation framing |
| El-Khaled et al. 2025 — heat & cold bleaching vulnerability | 10.1038/s42003-025-09329-5 | [OA] | 5 | incl. *A. pulchra* |
| Hoadley et al. 2015 — temp × pCO2 physiology across genera | 10.1038/srep18371 | [OA] | 3 | |
| Puisay et al. 2018 — larval thermal resistance/acclimation | 10.1016/j.marenvres.2018.01.005 | [PW] | 4 | larvae |
| Babka 2017 — successive bleaching resilience | 10.60950/bgtl.2017.2.babka.sk | [OA] | 4 | |

### Symbiont community / microbiome
| Paper | DOI | Acq | Rel | Note |
|---|---|---|---|---|
| Lock et al. 2025 — Symbiodiniaceae + bacterial microbiome impacts on host | 10.1111/1462-2920.70175 | [OA] | 5 | |
| Miller et al. 2024 — seasonal tissue/mucus microbiome baseline | 10.7717/peerj.17421 | [OA] | 5 | symbiont/microbiome baseline |
| Alina et al. 2023 — low-tide exposure × microbial abundance (intertidal) | 10.1080/17451000.2023.2169464 | [PW→CDP] | 4 | intertidal microbiome; first author Dining Nika Alina (was misfiled as Nuñez Lendo); from completeness census |

### Growth / calcification / skeleton
| Paper | DOI | Acq | Rel | Note |
|---|---|---|---|---|
| Yap & Gomez 1984 — growth of *A. pulchra* | 10.1007/bf00393119 | [NB][PW] | 5 | 13.1–15.8 cm/yr; 36.4% branch mortality → `v_edge`,`m_d` |
| Yap 1985 — growth II (temp × daylength) | 10.1007/bf00539430 | [PW] | 5 | |
| Strømgren 1987 — light × growth, intertidal | 10.1007/bf00302211 | [PW] | 4 | |
| Comeau et al. 2014 — irradiance × calcification | 10.1016/j.jembe.2013.12.013 | [NB][PW] | 5 | tissue biomass 1.87–3.15 mg/cm² → `E_max` |
| Comeau et al. 2016 — calcification vs temp & pCO2 (parameterization) | 10.1007/s00338-016-1425-0 | [NB][PW] | 5 | |
| Roche et al. 2010a — skeletal porosity (X-ray microCT) | 10.1016/j.jembe.2010.10.006 | [PW] | 4 | perforate structure → `D_E`,`resid` |
| Roche et al. 2010b — spatial porosity variation | 10.1007/s00338-010-0679-1 | [PW] | 4 | |
| Tian et al. 2025 — internal skeletal hydrodynamics | 10.1016/j.isci.2025.111742 | [OA] | 3 | |

### Physiology / OA / nutrients / reproduction / ecology
| Paper | DOI | Acq | Rel | Note |
|---|---|---|---|---|
| Conetta 2021 — phenotypic plasticity (URI thesis) | 10.23860/thesis-conetta-dennis-2021 | [NB][OA] | 5 | |
| Tanaka et al. 2006 — organic-N translocation/conservation | 10.1016/j.jembe.2006.04.011 | [PW] | 4 | translocation → `D_E` |
| Tanaka et al. 2009 — net DOM release | 10.1016/j.jembe.2009.06.023 | [PW] | 3 | |
| Buckingham et al. 2022 — N+P enrichment, skewed N:P | 10.1007/s00338-022-02223-0 | [OA] | 3 | |
| Huang 2009 — oocyte development histology (Sanya) | (no DOI) | [PW] | 4 | only A. pulchra reproduction source found |
| ⚠ Almeida et al. 2024 — intertidal vs subtidal recovery/growth/survival | 10.1007/s00227-024-04546-8 | [OA] | 3 | **Pocillopora cf. damicornis, Kenya (not A. pulchra)**; recovery-dynamics framing |
| Ladd et al. 2025 — growth–predation tradeoffs, distribution | 10.1038/s41598-025-21028-z | [OA] | 4 | |
| Raymundo et al. 2025 — restoration in stressful environment (Guam) | 10.1016/j.isci.2025.112244 | [OA] | 4 | |
| Munk 2024 — host & symbiont physiology during wound regeneration (thesis) | (no DOI) | [NB] | 5 | **KEY**: abrasion vs fragmentation × 27.9/29.5 °C, Mo'orea → maps to model excav × bleaching |

## 2. *Acropora*-genus trait grounding (for model parameters lacking *A. pulchra* data)

| Trait → param | Paper | DOI | Acq | Rel |
|---|---|---|---|---|
| Heterotrophy `h` (low plasticity) | Conti-Jerpe & Baker 2020 (Sci Adv) | 10.1126/sciadv.aaz5443 | [LAB][OA] | 5 |
| Heterotrophy `h` | Houlbrèque & Ferrier-Pagès 2009 | 10.1111/j.1469-185x.2008.00058.x | [PW] | 5 |
| Heterotrophy `h` | Radice et al. 2019 (FA+SIA) | 10.1371/journal.pone.0222327 | [OA] | 4 |
| Reserves `E_max` / fast-fragile | Denis et al. 2013 (*A. muricata* growth↔regen) | 10.1371/journal.pone.0072618 | [LAB][OA] | 5 |
| Reserves `E_max` | Schoepf et al. 2013 (Acropora reserves) | 10.1371/journal.pone.0075049 | [OA] | 4 |
| Reserves `E_max` | Leinbach et al. 2021 (recovery costs) | 10.1038/s41598-021-02807-w | [OA] | 4 |
| Translocation `D_E` (perforate) | Oren et al. 1997 (oriented 14C transport) | 10.3354/meps161117 | [OA] | 5 |
| Translocation `D_E` | Fine & Loya 2002 (bleaching cuts translocation) | 10.3354/meps234119 | [OA] | 4 |
| Translocation `D_E` | Swain et al. 2018 (colony integration) | 10.3354/meps12445 | [OA] | 4 |
| Healing rate `kc_mult` | Meesters et al. 1997 (lesion-area model) | 10.3354/meps146091 | [LAB?][OA] | 5 |
| Healing rate `kc_mult` | Denis et al. 2011 (lesion regen, populations) | (in lab) | [LAB] | 5 |
| Healing rate `kc_mult` | Hall 1997/1998/2001/2015 (Acropora regen series) | (in lab) | [LAB] | 5 |
| Healing rate `kc_mult` | Lirman et al. 2014 (*A. cervicornis* growth) | 10.1371/journal.pone.0107253 | [OA] | 4 |
| Re-symbiosis `iZ` | Bay et al. 2016 (*A. millepora* recovery threshold) | 10.1098/rsos.160322 | [OA] | 5 |
| Re-symbiosis `iZ` | Fitt et al. 2000 (seasonal symbiont density) | 10.4319/lo.2000.45.3.0677 | [OA] | 4 |
| Size at repro `B_mat` | Álvarez-Noriega et al. 2016 (fecundity × morphology) | 10.1002/ecy.1588 | [PW] | 5 |
| Size at repro `B_mat` | Hall & Hughes 1996 (modular reproduction) | 10.2307/2265514 | [PW] | 5 |
| Size at repro `B_mat` | Soong & Lang 1992 (reproductive integration) | (in lab) | [LAB] | 4 |
| Mortality/breakage `m_d` | Madin et al. 2014 (mechanical vulnerability) | 10.1111/ele.12306 | [OA] | 5 |
| Mortality/breakage `m_d` | Highsmith 1982 (fragmentation) | 10.3354/meps007207 | [OA] | 4 |
| Mortality/breakage `m_d` | Lirman 2000a (*Acropora* fragmentation; in lab) | (in lab) | [LAB] | 4 |

## 3. *A. pulchra* model parameters — **library-grounded (2026-06-08)**

Fast–fragile–sensitive branching Acropora. Each value below was extracted from this library
(5 parallel extraction passes) and positioned against the model's Porites/Pocillopora anchors.
The cited, full registry rows (literature_value · units · citation · source PDF · verification)
are in **`~/coral-wound-healing-model/data/phase2_parameters.csv`** (15 `Acropora pulchra` rows).

| Param | Grounded value | vs first-pass | Anchor (Por / Poc) | Key evidence | Confidence |
|---|---|---|---|---|---|
| `v_edge` | **0.045** | = | 0.010 / 0.030 | A. pulchra 13.1–15.8 cm/yr [Yap & Gomez 1984] | **HIGH** ⭐sp. |
| `m_d` | **0.0045** | = | 0.0010 / 0.0030 | A. pulchra 36.4% branch mortality/yr [Yap & Gomez 1984] | **HIGH** ⭐sp. |
| `kc_mult` | **1.2** | ↑ 1.1 | 0.7 / 1.0 | Acropora>Poc>Por [Hall 1997/2001]; Munk ~92%/19d (upper bd) | MED |
| `h` | **0.0002** | = | 0.0004 / 0.00025 | Acropora 94% host-symbiont overlap = top autotroph [Conti-Jerpe 2020] | MED-HIGH rank |
| `E_max` | **0.85** | ↓ 0.9 | 2.0 / 1.0 | A. millepora 6.6 vs Poc 8.1 kJ/g (ratio 0.81) [Schoepf 2013] | HIGH rank |
| `E0` | **0.77** | new | 1.5 / 0.9 | 0.90×E_max (Poc fill-fraction); deep depletion [Leinbach 2021, Mo'orea] | LOW-MED |
| `bT_pulse` | **0.010** | = | 0.0017 / 0.0080 | ED50 ~34–37 °C below congeners [Denis 2024; Cunning 2021]; 58% mort d34 [Berg] | LOW-MED |
| `bleach_z` | **0.15** | = | 0.55 / 0.20 | obligate autotroph; ~29% chl retained at +9 °C [Denis 2024] | LOW-MED |
| `iZ` | **0.004** | = | 0.020 / 0.004 | Mo'orea Acropora rebuild ~4 mo [Thomas 2019, same site] | **MED (best)** |
| `Z0` | **0.0** | new | 0.3 / 0.0 | bleaches white, S:H→0.003 [Bay 2016] | MED |
| `KZ` | **1.0** | new | 1.0 / 1.0 | normalization; A. pulchra 1.345×10⁶ cells/cm² [Anthony 2023] | HIGH |
| `D_E` | **0.014** | ↑ 0.012 | 0.02 / 0.002 | tip porosity 65–80%, connected interior [Roche 2010b] | HIGH porosity ⭐sp. |
| `resid` | **0.07** | ↑ 0.05 | 0.50 / 0.05 | thin branches 0.8–1.8 cm but perforate deep tissue [Roche 2010b] | MED |
| `maxpd` | **1.0** | new | 3.0 / 1.2 | 16 cm branches, 15–25 mm apical tips [Yap & Gomez; Strömgren] | LOW-MED |
| `B_mat` | **4.5** | ↓ 5.0 | 8.0 / 4.0 | A. cervicornis puberty 17 cm branch [Soong 1992]; recruits 5 cm/2 yr [Wallace] | LOW-MED |

⭐sp. = A.-pulchra-specific number; others congener-borrowed (flagged in the registry). **Five
changed from first-pass:** `kc_mult` 1.1→1.2, `E_max` 0.9→0.85, `D_E` 0.012→0.014, `resid` 0.05→0.07,
`B_mat` 5→4.5; plus four newly set (`E0`, `Z0`, `KZ`, `maxpd`). **Weakest-grounded (sweep these):**
`E0`, `maxpd`, `bleach_z`/`bT_pulse` (no A. pulchra loss-rate or ED50 in local PDFs — inferred from
congener thresholds + external Cunning 2024 / Berg 2020 [NB]). **Brainstorming gate:** present this
table for Adrian's biology review before committing the `Apulchra` preset to `ext/species_presets.R`.

## 4. Gaps (genuinely thin in the literature)
- **A. pulchra reproduction / fecundity / size-at-maturity** — only Huang 2009 (histology); rely on congeners (Álvarez-Noriega, Hall & Hughes, Soong & Lang).
- **A. pulchra wounding/regeneration rates** — only Munk 2024 (thesis); rely on Denis 2011/2013, Hall series, Lirman.
- **Acropora symbiont repopulation RATE after bleaching** — none found (only density thresholds). Least-grounded param; borrow from Symbiodiniaceae division-rate lit or fit.

## 4b. Completeness census vs OpenAlex (2026-06-08)

Checked the library against the full published record (OpenAlex). **35 works have "Acropora
pulchra" in the title.** We now hold or track **all relevant, accessible, English ones**; the
remainder are out of scope or non-acquirable:

- **Acquired this round:** Alina et al. 2023 (intertidal microbiome, Mar Biol Res) — the one
  clear English gap the census found; published version pulled via authenticated CDP.
- **Not chased — out of scope:** 3 natural-product **chemistry** papers (Chinese, 2001/2003); 2 **IUCN
  Red List** dataset entries (not papers); duplicate **preprints** of papers already held (Miller 2024,
  Conn 2025).
- **Not chased — non-English / no DOI / low marginal value:** Japanese calcification & photosynthesis
  papers (1995 ×2, 2002) + a 1997 Japanese growth thesis + an Indonesian 2003 growth note — all
  duplicate topics already covered (Yap growth, Comeau calcification, Strømgren light); a 2020 "VOC
  day/night cycling" note and a Tanaka-2008 DOM paper have **no resolvable DOI** in OpenAlex.

**Caveat on "all A. pulchra papers":** title-search captures papers *about* A. pulchra. OpenAlex also
shows **92 abstract-mentions** and **396 full-text mentions** — multi-species studies that merely
include A. pulchra in a species list. The library holds a *curated* subset of those (the model-grounding
+ LTH papers), **not all 396**, which is intentional. So: complete on A.-pulchra-subject papers
(accessible English), deliberately selective on multi-species mentions.

## 5. Acquisition results (2026-06-07)

**37 PDFs in `literature/pdfs/`.** 10 copied from the lab library, 20 via curl (Springer
content/pdf, Nature, PLOS, Frontiers, RSOS, Gigabyte), 7 via EuropePMC render.

**NotebookLM "Acropora pulchra" notebook → 26 sources** (10 prior + 16 new *A. pulchra*-specific
papers ingested). Genus-grounding papers (Section 2, about other corals) were kept local-only to
keep the notebook focused on *A. pulchra*; they belong in the regen notebook if needed there.

### Not obtained (paywalled, no OA route; browser unauthenticated)
| Paper | Why missing | Mitigation |
|---|---|---|
| Cunning 2024, Berg 2020, Comeau 2014, Conetta 2021, Munk 2024 | already **[NB]** (in notebook) | content queryable; local PDF not essential |
| Houlbrèque & Ferrier-Pagès 2009 (Wiley, `h`) | paywall, review article | use Conti-Jerpe 2020 + Radice 2019 (primary, `h`) |
| Álvarez-Noriega 2016 (Wiley) & Hall & Hughes 1996 (JSTOR), `B_mat` | paywall | use Soong 1992 (LAB) for size-at-repro |
| Fitt 2000 (Wiley L&O, `iZ`) | paywall | use Bay 2016 (OA) for symbiont threshold |
| Puisay 2018, Roche 2010a, Tanaka 2006/2009 (Elsevier) | paywall, no PMC | Roche 2010b (OA) covers porosity; lower-rel |
| Huang 2009 (no DOI), Babka 2017 (Bremen thesis) | not locatable / bot-blocked | — |

### Per-paper one-line summaries (acquired PDFs)

**A. pulchra-specific — expression / genomics / microbiome**
- **Conn 2025** — reference genome assembly & annotation for Mo'orea *A. pulchra* (anchors all expression work).
- **⚠ Zhang 2025** — protein + metabolite acclimation to thermal variability in the **macroalga *Halimeda macroloba*** (NOT *A. pulchra*; not a coral). Calcifier thermal-acclimation framing only — *A. pulchra* proteomics remains unstudied.
- **Anthony 2023** — Symbiodiniaceae cellular/phenotypic plasticity under stress.
- **de la Vega 2023** — genome of an *Endozoicomonas* strain isolated from *A. pulchra* (microbiome resource).
- **Lock 2025** — combined Symbiodiniaceae + bacterial microbiome effects on host performance.
- **Miller 2024** — seasonal tissue/mucus microbiome baseline.

**A. pulchra-specific — thermal tolerance / bleaching**
- **Berg 2020** — persistent photosystem damage; 58% mortality by day 34, recovery needs ~3× stress duration ([NB], no local PDF).
- **Shaw 2016** — intraspecific (genotype) variability under warming × acidification.
- **⚠ Denis 2024** — individual thermal-tolerance traits across the GBR in ***A. spathulata*** (NOT *A. pulchra*); CBASS intraspecific-variation framing. Published RSPB version acquired.
- **El-Khaled 2025** — heat *and* cold bleaching vulnerability, incl. *A. pulchra*.
- **Hoadley 2015** — temperature × pCO₂ physiology across genera.
- **Grottoli 2021** — methods/standardization for cross-experiment bleaching comparability.

**A. pulchra-specific — growth / calcification / skeleton**
- **Yap & Gomez 1984** — growth 13.1–15.8 cm/yr; 36.4% branch mortality → `v_edge`, `m_d`.
- **Yap 1985** — growth II, temperature × daylength.
- **Strømgren 1987** — light × growth in intertidal *A. pulchra*.
- **Comeau 2016** — calcification response parameterized vs temp & pCO₂.
- **Roche 2010b** — spatial variation in skeletal porosity → `D_E`, `resid`.
- **Tian 2025** — internal skeletal hydrodynamics.

**A. pulchra-specific — physiology / nutrients / ecology / recovery**
- **Buckingham 2022** — N+P enrichment with skewed N:P stoichiometry.
- **⚠ Almeida 2024** — intertidal vs subtidal recovery, growth, survival in ***Pocillopora* cf. *damicornis*** (Kenya; NOT *A. pulchra*).
- **Ladd 2025** — growth–predation tradeoffs shaping distribution.
- **⚠ Raymundo 2025** — restoration outcomes in a stressful (Guam) environment for ***A. aspera*** outplants (*A. pulchra* is only the background congener at the recipient site).

**Acropora-genus grounding (local PDFs only; for params lacking *A. pulchra* data)**
- **Conti-Jerpe 2020** & **Radice 2019** — heterotrophic vs autotrophic strategy → `h` (low plasticity).
- **Denis 2013** — fast growth trades off against regeneration in *A. muricata* → fast-fragile axis.
- **Schoepf 2013** & **Leinbach 2021** — *Acropora* energy reserves / recovery costs → `E_max`.
- **Denis 2011**, **Hall 1997/1998/2001/2015**, **Meesters 1997**, **Lirman 2014** — lesion/regeneration rates → `kc_mult`.
- **Bay 2016** — symbiont-density recovery threshold → `iZ`.
- **Soong 1992** — reproductive integration / size effects → `B_mat`.
- **Madin 2014** — mechanical vulnerability / colony breakage → `m_d`.
- **Lirman 2000a** — branching-*Acropora* fragmentation dynamics → `m_d`.
