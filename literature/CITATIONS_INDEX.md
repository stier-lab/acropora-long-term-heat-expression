# Citations index — LTH manuscript

> 🗂️ **Citation → PDF index** · Updated 2026-06-12 · Index: [`README.md`](../README.md) · pairs with `manuscript/references.bib` and `manuscript/Manuscript_LTH.md`.

Maps every work cited in `manuscript/Manuscript_LTH.md` to its BibTeX key
(`manuscript/references.bib`), verified DOI, and local PDF (or acquisition status).
**All 20 DOI-bearing references were verified against Crossref on 2026-06-09** and their
BibTeX pulled directly from Crossref (not hand-typed). Use `references.bib` as the
reference list; use this file to find the PDF behind a cite.

Status legend: ✅ local PDF · 📓 NotebookLM corpus only · ⬇ to acquire (needs UCSB
auth — see foot of file) · 📄 method/software (no article PDF).

---

## Mismatch resolution log (2026-06-09)

The earlier inline cites I drafted had several mismatches against the library; all are now resolved:

- **Hughes et al. 2018 → Hughes et al. 2017** — reprosed to the in-library paper (*Global warming
  and recurrent mass bleaching*, Nature, 10.1038/nature21707) which grounds the same claim.
- **Meesters et al. 1994 → Meesters & Bak 1993** — corrected a real citation error: the
  bleaching→regeneration claim belongs to Meesters & Bak 1993 (*Effects of coral bleaching on tissue
  regeneration potential*, 10.3354/meps096189), not the 1994 damage-growth paper.
- **Traylor-Knowles 2016 DOI corrected** in `CITATION_AUDIT.md`: the listed
  10.1186/s12862-016-0791-0 was a wrong-paper typo (a RAD-seq study); correct DOI is
  **10.1007/s00227-016-3011-y**.
- **Acquired the orphan PDFs** that were OA-reachable: Jokiel & Coles 1977 & 1990 (Springer),
  **Cunning et al. 2024** (the thermal-context anchor), Davies 1989, Siebeck 2006 — library 93 → 98.
- **Kept (verified canonical) with an in-library relative** where the exact paper is paywalled:
  Warner 1999 (lib has Warner 1996), Hoegh-Guldberg 1999 (lib has HG 2007), Comeau 2014 (lib has
  Comeau 2016). These prose years are correct; only the local PDF differs — see acquisition list.

---

## Index

| Cite (manuscript) | bib key | DOI | PDF / status |
|---|---|---|---|
| Berg et al. 2020 | `Berg_2020` | 10.1080/17451000.2021.1875245 | 📓 NotebookLM (T&F paywall) |
| Bonesso et al. 2017 | `Bonesso_2017` | 10.7717/peerj.3719 | ✅ Bonesso_2017_SubBleachingImpairsRegen.pdf |
| Comeau et al. 2014 | `Comeau_2014` | 10.1016/j.jembe.2013.12.013 | 📓 NotebookLM (lib has Comeau_2016; Elsevier paywall) |
| Conn et al. 2025 | `Conn_2025` | 10.46471/gigabyte.153 | ✅ Conn_2025_GenomeAssembly.pdf |
| Cunning et al. 2024 | `Cunning_2024` | 10.1007/s00338-024-02577-7 | ✅ Cunning_2024_RapidAcuteHeatToleranceAssays.pdf ★new |
| Davies 1989 | `Spencer_Davies_1989` | 10.1007/BF00428135 | ✅ Davies_1989_BuoyantWeightTechnique.pdf ★new |
| Dixon et al. 2015 | `Dixon_2015` | 10.1126/science.1261224 | ✅ Dixon_2015_GenomicDeterminantsHeatTolerance.pdf |
| Henry & Hart 2005 | `Henry_2005` | 10.1002/iroh.200410759 | ✅ Henry_2005_RegenInjuryResourceAllocation.pdf |
| Highsmith 1982 | `Highsmith_1982` | 10.3354/meps007207 | ⬇ int-res blocks automation |
| Hoegh-Guldberg 1999 | `Hoegh_Guldberg_1999` | 10.1071/MF99078 | ⬇ CSIRO paywall (lib has HoeghGuldberg_2007) |
| Hughes et al. 2017 | `Hughes_2017` | 10.1038/nature21707 | ✅ Hughes_2017_GlobalWarmingMassBleaching.pdf |
| Jokiel & Coles 1977 | `Jokiel_1977` | 10.1007/BF00402312 | ✅ JokielColes_1977_TemperatureMortalityGrowthHawaiianCorals.pdf ★new |
| Jokiel & Coles 1990 | `Jokiel_1990` | 10.1007/BF00265006 | ✅ JokielColes_1990_ResponseIndoPacificElevatedTemperature.pdf ★new |
| Madin et al. 2014 | `Madin_2014` | 10.1111/ele.12306 | ✅ Madin_2014_MechanicalVulnerability.pdf |
| Meesters & Bak 1993 | `Meesters_1993` | 10.3354/meps096189 | ✅ Meesters_1993_BleachingTissueRegen.pdf |
| Shaw et al. 2016 | `Shaw_2016` | 10.1007/s00227-016-2986-8 | ✅ Shaw_2016_IntraspecificVariability.pdf |
| Siebeck et al. 2006 | `Siebeck_2006` | 10.1007/s00338-006-0123-8 | ✅ Siebeck_2006_ColourReferenceCard.pdf ★new |
| Traylor-Knowles et al. 2016 | `Traylor_Knowles_2016` | 10.1007/s00227-016-3011-y | ✅ Traylor-Knowles_2016_WoundHealingTwoTempRegimes.pdf |
| Warner et al. 1999 | `Warner_1999` | 10.1073/pnas.96.14.8007 | ⬇ PNAS needs auth (lib has Warner_1996) |
| Yap & Gomez 1984 | `Yap_1984` | 10.1007/BF00393119 | ✅ YapGomez_1984_Growth.pdf |
| Munk 2024 | `Munk_2024` | — (M.Sc. thesis) | 📓 NotebookLM only; **confirm exact title/institution with author** |
| Jokiel et al. 1978 | `Jokiel_1978` | — (book chapter, UNESCO) | 📄 method (buoyant-weight technique) |
| R Core Team 2025 | `RCoreTeam_2025` | — (software) | 📄 software |

---

## Still to acquire (need UCSB-authenticated browser session; see local setup notes)

These six are DOI-verified and fully cited in `references.bib`; only the local PDF is missing
(each has a close in-library relative or a NotebookLM copy, so the manuscript is not blocked):

| Paper | DOI | Route |
|---|---|---|
| Warner, Fitt & Schmidt 1999 | 10.1073/pnas.96.14.8007 | PNAS (free, but needs a real browser session) |
| Hoegh-Guldberg 1999 | 10.1071/MF99078 | CSIRO Publishing (paywall) |
| Comeau et al. 2014 | 10.1016/j.jembe.2013.12.013 | Elsevier (paywall) |
| Berg et al. 2020 | 10.1080/17451000.2021.1875245 | Taylor & Francis (paywall) |
| Highsmith 1982 | 10.3354/meps007207 | int-res.com (blocks automation; manual save) |
| Munk 2024 (thesis) | — | no online copy located; request from author |

After fetching, drop the PDF into `literature/pdfs/AuthorYear_ShortTitle.pdf`, verify
`head -c5 file.pdf` == `%PDF-`, and update the row above to ✅.
