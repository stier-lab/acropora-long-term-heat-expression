# Library map — paper → *A. pulchra* relevance → model/project connection

> 🗂️ **Paper → relevance → connection map** · Updated 2026-06-12 · Index: [`DOCS_INDEX.md`](../DOCS_INDEX.md) · machine-readable twin: `literature/library_map.csv`.

**What this is.** Every PDF in `pdfs/` mapped to (1) how it relates to *Acropora pulchra* and
(2) the specific wound-healing-model parameter or LTH-project role it supports. The library is
**not** an *A. pulchra*-only bibliography — it is the evidence base for the 3-species wound-healing
model (*Porites*, *Pocillopora*, *A. pulchra*) and the LTH "expression by temperature" project, so
most parameters are grounded on congeners/mechanism papers where *A. pulchra*-specific data are thin.
Machine-readable companion: [`library_map.csv`](library_map.csv). See also `CITATION_AUDIT.md`
(how the grounding papers were found) and `LITERATURE.md` (the *A. pulchra* discovery index).

**Composition of 98 PDFs:** 16 **Core** (*A. pulchra* is the subject) ·
6 **Includes *A. pulchra*** (multi-species) · 76 **Grounding** (genus / other-coral / mechanism).
(The 5 most recently acquired — Jokiel & Coles 1977/1990, Cunning 2024, Davies 1989, Siebeck 2006 —
are methods/grounding citations; see `CITATIONS_INDEX.md`.)

> ⚠ **Species-attribution corrections (2026-06-08).** Four papers were verified against their PDFs
> and reclassified out of the *A. pulchra* groups (rows below flagged `⚠`): **Zhang 2025** =
> macroalga *Halimeda macroloba*; **Almeida 2024** = *Pocillopora* cf. *damicornis*; **Raymundo
> 2025** = *A. aspera*; **Denis 2024** = *A. spathulata*. They stay in the library as grounding/
> framing but are no longer counted as *A. pulchra*-subject. (Counts above already updated.)

---

## Model-parameter cross-reference (which papers ground each parameter)

| Model parameter | Meaning | Grounded by |
|---|---|---|
| `v_edge` | edge growth rate | Lirman et al. 2014, Yap & Gomez 1984, Yap 1985 |
| `m_d` | mortality / breakage | Dias et al. 2018, Lirman 2000a, Madin & Connolly 2006, Madin et al. 2014, Marshall 2000, Meesters et al. 1997, Tunnicliffe 1981, Yap & Gomez 1984 |
| `kc_mult` | wound-healing rate | Bak & Steward-Van Es 1980, Bak 1983, Bonesso et al. 2017, Burmester et al. 2017, Burmester et al. 2018, Counsell et al. 2019, DeFilippo et al. 2016, Denis et al. 2011, Dias et al. 2018, Edmunds & Lenihan 2010, Hall 1997, Hall 1998, Hall 2001, Hall 2015, Lirman 2000b, Meesters & Bak 1993, Meesters et al. 1997, Soong & Chen 2003, Traylor-Knowles et al. 2016 |
| `h` | heterotrophic buffer | Burmester et al. 2018, Conti-Jerpe & Baker 2020, Grottoli et al. 2006, Radice et al. 2019 |
| `E_max` | energy reserves | Denis et al. 2013, Leinbach et al. 2021, Rodrigues & Grottoli 2007, Schoepf et al. 2013 |
| `iZ` | symbiont repopulation | Baker 2003, Bay et al. 2016, Jones et al. 2008, Rouzé et al. 2017, Rouzé et al. 2019, Thomas et al. 2019 |
| `D_E` | perforate translocation | Roche et al. 2010b |
| `resid` | residual reserve / porosity | Roche et al. 2010b |
| `B_mat` | size at maturity | Baird et al. 2009, Carroll et al. 2006, Soong & Lang 1992, Wallace 1985 |
| `bleach_z` | bleaching sensitivity | Shaw et al. 2016 |
| `bT_pulse` | bleaching pulse | Shaw et al. 2016 |
| `Tier-3` | colony energy budget / microfrag | Lock et al. 2022, Rinkevich 1996 |

---

## 1. Core — *A. pulchra* is the study subject (19)
| Paper | Theme | Connection — what it grounds |
|---|---|---|
| ⚠ [Almeida et al. 2024](pdfs/Almeida_2024_IntertidalSubtidalRecovery.pdf) | Recovery | **Pocillopora cf. damicornis, Kenya (NOT A. pulchra)** — intertidal/subtidal recovery-dynamics framing |
| [Anthony et al. 2023](pdfs/Anthony_2023_SymbiontPlasticity.pdf) | Expression/Symbiont | Symbiodiniaceae cellular/phenotypic plasticity under stress |
| [Buckingham et al. 2022](pdfs/Buckingham_2022_NPEnrichment.pdf) | Nutrients | N+P enrichment, skewed N:P |
| [Comeau et al. 2016](pdfs/Comeau_2016_CalcificationTempCO2.pdf) | Calcification | Calcification vs temp & pCO2 → calcification parameterisation |
| [Conn et al. 2025](pdfs/Conn_2025_GenomeAssembly.pdf) | Genomics | Reference genome assembly/annotation (Mo'orea) — underpins all expression work (LTH) |
| [Ladd et al. 2025](pdfs/Ladd_2025_GrowthPredationTradeoffs.pdf) | Ecology | Growth–predation risk tradeoffs constrain A. pulchra distribution |
| [Lock et al. 2025](pdfs/Lock_2025_SymbiontBacterialMicrobiome.pdf) | Microbiome | Symbiont + bacterial microbiome effects on host |
| [Miller & Bentlage 2024](pdfs/Miller_2024_SeasonalMicrobiome.pdf) | Microbiome | Seasonal tissue/mucus microbiome baseline |
| [Nuñez Lendo et al. 2023](pdfs/NunezLendo_2023_IntertidalMicrobialAbundance.pdf) | Microbiome | Low-tide exposure shifts microbial abundance in intertidal A. pulchra |
| ⚠ [Raymundo et al. 2025](pdfs/Raymundo_2025_RestorationGuam.pdf) | Restoration | **A. aspera outplants, Guam (NOT A. pulchra)** — A. pulchra is the background congener at the site |
| [Roche et al. 2010b](pdfs/Roche_2010b_SpatialPorosity.pdf) | Skeleton | Spatial skeletal porosity → D_E, resid (perforate translocation) |
| [Rouzé et al. 2017](pdfs/Rouze_2017_SymbiodiniumMoorea.pdf) | Symbiont/Mo'orea | Symbiodiniaceae of Mo'orea corals incl. A. pulchra → iZ context |
| [Shaw et al. 2016](pdfs/Shaw_2016_IntraspecificVariability.pdf) | Thermal tolerance | Genotype variation under warming×acidification → bleaching sensitivity (bT_pulse, bleach_z) |
| [Strømgren 1987](pdfs/Stromgren_1987_LightGrowth.pdf) | Growth | Light × growth, intertidal → growth response |
| [Tian et al. 2025](pdfs/Tian_2025_SkeletalHydrodynamics.pdf) | Skeleton | Internal skeletal hydrodynamics |
| [Yap & Gomez 1984](pdfs/YapGomez_1984_Growth.pdf) | Growth | Growth 13.1–15.8 cm/yr; 36.4% branch mortality → v_edge, m_d |
| [Yap 1985](pdfs/Yap_1985_GrowthII.pdf) | Growth | Growth II (temp × daylength) → v_edge |
| ⚠ [Zhang et al. 2025](pdfs/Zhang_2025_ProteinMetaboliteAcclim.pdf) | Expression | **Macroalga *Halimeda macroloba* (NOT A. pulchra, not a coral)** — calcifier thermal-acclimation framing only |
| [de la Vega et al. 2023](pdfs/delaVega_2023_EndozoicomonasGenome.pdf) | Microbiome | Endozoicomonas GU-1 genome from A. pulchra — microbiome resource |

## 2. Includes *A. pulchra* among multiple species (7)
| Paper | Theme | Connection — what it grounds |
|---|---|---|
| [Carroll et al. 2006](pdfs/Carroll_2006_AcroporaReproductionMoorea.pdf) | Reproduction/Mo'orea | Sexual reproduction of Acropora at Mo'orea → B_mat (field site) |
| ⚠ [Denis et al. 2024](pdfs/Denis_2024_IndividualThermalTolerance.pdf) | Thermal tolerance | **A. spathulata, GBR (NOT A. pulchra)** — CBASS intraspecific-variation framing |
| [El-Khaled et al. 2025](pdfs/ElKhaled_2025_HeatColdBleaching.pdf) | Thermal tolerance | Heat & cold bleaching vulnerability incl. A. pulchra → bleaching sensitivity |
| [Grottoli et al. 2020](pdfs/Grottoli_2020_BleachingComparability.pdf) | Methods | Standardising bleaching experiments — LTH methods/comparability |
| [Hoadley et al. 2015](pdfs/Hoadley_2015_TempCO2Physiology.pdf) | Thermal physiology | Temp × pCO2 physiology across genera incl. A. pulchra |
| [Rouzé et al. 2019](pdfs/Rouze_2019_SymbiodiniaceaeSignatureMoorea.pdf) | Symbiont/Mo'orea | Quantitative Symbiodiniaceae signature of Mo'orea colonies → iZ |
| [Soong & Chen 2003](pdfs/SoongChen_2003_AcroporaFragmentNursery.pdf) | Wound/regen | Regeneration & growth of Acropora fragments (incl. pulchra) → kc_mult |

## 3. Grounding — genus / other-coral / mechanism papers (67)
*Why they're here: the wound-healing model and LTH thermal-tolerance work need broader coral
literature for parameters and framing that lack A. pulchra-specific data.*

### Wound healing / regeneration → `kc_mult`
| Paper | Theme | Connection — what it grounds |
|---|---|---|
| [Bak & Steward-Van Es 1980](pdfs/Bak_1980_RegenSuperficialDamage.pdf) | Wound healing | Regeneration of superficial damage → kc_mult (foundational) |
| [Bak 1983](pdfs/Bak_1983_NeoplasiaRegenGrowthApalmata.pdf) | Wound healing | Neoplasia/regeneration/growth in A. palmata → kc_mult (foundational) |
| [Bonesso et al. 2017](pdfs/Bonesso_2017_SubBleachingImpairsRegen.pdf) | Wound healing | Sub-bleaching heat impairs recovery/regeneration → kc_mult(temp) |
| [Burmester et al. 2017](pdfs/Burmester_2017_TempSymbiosisLesionRecovery.pdf) | Wound healing | Temperature & symbiosis affect lesion recovery → kc_mult |
| [Burmester et al. 2018](pdfs/Burmester_2018_AutoVsHeteroLesionRecovery.pdf) | Wound healing | Autotrophy vs heterotrophy in lesion recovery → kc_mult, h |
| [Counsell et al. 2019](pdfs/Counsell_2019_ColonySizeDepthWoundRepair.pdf) | Wound healing | Colony size & depth affect wound repair → size-dependent kc_mult |
| [DeFilippo et al. 2016](pdfs/DeFilippo_2016_SurfaceLesionRecoveryAstrangia.pdf) | Wound healing | Surface lesion recovery in Astrangia poculata → kc_mult kinetics |
| [Denis et al. 2011](pdfs/Denis_2011_LesionRegenerationPopulations.pdf) | Wound healing | Lesion regeneration in massive Porites populations → kc_mult |
| [Denis et al. 2013](pdfs/Denis_2013_FastGrowthImpairRegen.pdf) | Wound healing | Fast growth impairs regeneration in A. muricata → fast-fragile axis, E_max |
| [Dias et al. 2018](pdfs/Dias_2018_FragmentationRegenThermalStress.pdf) | Wound healing | Mortality/growth/regen after fragmentation under thermal stress → m_d, kc_mult |
| [Edmunds & Lenihan 2010](pdfs/EdmundsLenihan_2010_SublethalDamageJuvenilePorites.pdf) | Wound healing | Sub-lethal damage to juvenile Porites under temp/flow → kc_mult (Porites) |
| [Hall 1997](pdfs/Hall_1997_InterspecificRegeneration.pdf) | Wound healing | Interspecific differences in injury regeneration → kc_mult |
| [Hall 1998](pdfs/Hall_1998_InjuryRegeneration.pdf) | Wound healing | Injury & regeneration of common reef corals → kc_mult |
| [Hall 2001](pdfs/Hall_2001_AcroporaHyacinthusResponse.pdf) | Wound healing | A. hyacinthus injury response → kc_mult |
| [Hall 2015](pdfs/Hall_2015_LesionRecovery.pdf) | Wound healing | Lesion recovery under low pH → kc_mult |
| [Lirman 2000b](pdfs/Lirman_2000b_LesionRegenAcroporaPalmata.pdf) | Wound healing | Lesion regeneration in branching A. palmata → kc_mult |
| [Meesters & Bak 1993](pdfs/Meesters_1993_BleachingTissueRegen.pdf) | Wound healing | Bleaching reduces tissue-regeneration potential → kc_mult under bleaching |
| [Meesters et al. 1994](pdfs/Meesters_1994_DamageRegenGrowth.pdf) | Wound healing | Damage & regeneration ↔ growth → growth–regeneration tradeoff |
| [Meesters et al. 1997](pdfs/Meesters_1997_PredictingRegeneration.pdf) | Wound healing | Partial mortality vs depth/area; lesion-area model → kc_mult, m_d |
| [Rice et al. 2019](pdfs/Rice_2019_NitrogenCorallivoryRecovery.pdf) | Wound healing | N sources speed recovery from corallivory → browsing-axis recovery |
| [Traylor-Knowles et al. 2016](pdfs/Traylor-Knowles_2016_WoundHealingTwoTempRegimes.pdf) | Wound healing | Wound healing × temperature in Pocillopora damicornis & Acropora hyacinthus (model's 2 existing species) → kc_mult temp-coupling |

### Mechanical / breakage / microfragmentation → `m_d`, Tier-3
| Paper | Theme | Connection — what it grounds |
|---|---|---|
| [Lirman 2000a](pdfs/Lirman_2000a_AcroporaFragmentation.pdf) | Mechanical | Fragmentation in branching A. palmata → m_d |
| [Lock et al. 2022](pdfs/Lock_2022_CaHomeostasisMicrofragmentation.pdf) | Microfragmentation | Ca-homeostasis drives rapid growth after microfragmentation (Porites lobata) → Tier-3 microfrag |
| [Madin & Connolly 2006](pdfs/MadinConnolly_2006_HydrodynamicDisturbances.pdf) | Mechanical | Ecological consequences of hydrodynamic disturbance → m_d |
| [Madin et al. 2014](pdfs/Madin_2014_MechanicalVulnerability.pdf) | Mechanical | Mechanical vulnerability → size-dependent mortality m_d |
| [Marshall 2000](pdfs/Marshall_2000_SkeletalDamageColonyMorphology.pdf) | Mechanical | Skeletal damage vs colony morphology → m_d |
| [Tunnicliffe 1981](pdfs/Tunnicliffe_1981_BreakagePropagationAcervicornis.pdf) | Mechanical | Breakage & propagation of A. cervicornis → m_d, fragmentation |
| [Van de Water et al. 2015](pdfs/VanDeWater_2015_ImmuneResponsePhysicalDamage.pdf) | Mechanical/immune | Elevated temp & coral immune response after physical damage → healing immunity |

### Energetics / reserves / heterotrophy → `h`, `E_max`, Tier-3
| Paper | Theme | Connection — what it grounds |
|---|---|---|
| [Anthony et al. 2007](pdfs/Anthony_2007_BleachingEnergeticsMortality.pdf) | Energetics | Bleaching, energetics & mortality risk → energy-budget mortality |
| [Conti-Jerpe & Baker 2020](pdfs/Conti-Jerpe_2020_TrophicStrategyBleaching.pdf) | Energetics | Trophic strategy & bleaching resistance → h (autotroph-loser) |
| [Grottoli et al. 2006](pdfs/Grottoli_2006_HeterotrophicPlasticity.pdf) | Energetics | Heterotrophic plasticity & resilience in bleached corals → h |
| [Grottoli et al. 2014](pdfs/Grottoli_2014_CumulativeBleachingWinnersLosers.pdf) | Energetics | Cumulative annual bleaching → reserve depletion (winners→losers) |
| [Leinbach et al. 2021](pdfs/Leinbach_2021_RecoveryCosts.pdf) | Energetics | Recovery costs → E_max |
| [Levas et al. 2013](pdfs/Levas_2013_BleachingRecoveryPoritesLobata.pdf) | Energetics | Bleaching/recovery physiology in Porites lobata → energetics |
| [Radice et al. 2019](pdfs/Radice_2019_HeterotrophyFASIA.pdf) | Energetics | Trophic strategy via FA+SIA → h |
| [Rodrigues & Grottoli 2007](pdfs/RodriguesGrottoli_2007_EnergyReservesRecovery.pdf) | Energetics | Energy reserves/metabolism as recovery indicators → E_max |
| [Schoepf et al. 2013](pdfs/Schoepf_2013_AcroporaReserves.pdf) | Energetics | Acropora energy reserves & calcification → E_max |

### Energy allocation → Tier-3 budget
| Paper | Theme | Connection — what it grounds |
|---|---|---|
| [Henry & Hart 2005](pdfs/Henry_2005_RegenInjuryResourceAllocation.pdf) | Energy allocation | Regeneration & resource allocation (review) → energy-allocation framing |
| [Rinkevich 1996](pdfs/Rinkevich_1996_ReproRegenEnergyAllocation.pdf) | Energy allocation | Reproduction vs regeneration compete for energy → Tier-3 energy budget |

### Reproduction / size at maturity → `B_mat`
| Paper | Theme | Connection — what it grounds |
|---|---|---|
| [Baird et al. 2009](pdfs/Baird_2009_ReproductiveBiologyScleractinian.pdf) | Reproduction | Reproductive biology of scleractinians (review) → B_mat framing |
| [Omori et al. 2001](pdfs/Omori_2001_AcroporaFertilizationDrop.pdf) | Reproduction | Drop in Acropora fertilization after bleaching → repro×bleaching |
| [Soong & Lang 1992](pdfs/Soong_1992_ReproductiveIntegration.pdf) | Reproduction | Reproductive integration in reef corals → B_mat |
| [Wallace 1985](pdfs/Wallace_1985_AcroporaReproductionRecruitment.pdf) | Reproduction | Reproduction/recruitment/fragmentation in nine Acropora → B_mat |

### Symbiont identity / shuffling / threshold → `iZ`
| Paper | Theme | Connection — what it grounds |
|---|---|---|
| [Baker 2003](pdfs/Baker_2003_CoralAlgalSymbiosisFlexibility.pdf) | Symbiont | Flexibility/specificity in coral–algal symbiosis (review) → iZ framing |
| [Bay et al. 2016](pdfs/Bay_2016_RecoveryThreshold.pdf) | Symbiont | Symbiont-density recovery threshold → iZ |
| [Berkelmans & van Oppen 2006](pdfs/Berkelmans_2006_ZooxanthellaeThermalTolerance.pdf) | Symbiont | Zooxanthellae role in thermal tolerance → symbiont-mediated tolerance |
| [Jones et al. 2008](pdfs/Jones_2008_EndosymbiontCommunityChange.pdf) | Symbiont | Endosymbiont community change after bleaching → iZ shuffling |
| [LaJeunesse et al. 2018](pdfs/LaJeunesse_2018_SymbiodiniaceaeRevision.pdf) | Symbiont | Systematic revision of Symbiodiniaceae → symbiont taxonomy |

### Thermal tolerance / expression backbone (LTH)
| Paper | Theme | Connection — what it grounds |
|---|---|---|
| [Barshis et al. 2013](pdfs/Barshis_2013_GenomicBasisResilience.pdf) | Expression | Genomic basis for coral resilience → expression of resilience |
| [Bay & Palumbi 2015](pdfs/BayPalumbi_2015_RapidAcclimTranscriptome.pdf) | Expression | Rapid acclimation via transcriptome change → expression plasticity |
| [Cornwell et al. 2021](pdfs/Cornwell_2021_HeatToleranceSymbiontLoadTradeoff.pdf) | Thermal tolerance | Heat tolerance × symbiont load → growth-tolerance tradeoff (A. hyacinthus) |
| [Cunning et al. 2021](pdfs/Cunning_2021_HeatToleranceCensusStaghorn.pdf) | Thermal tolerance | Heat-tolerance census of Florida staghorn Acropora → acute-tolerance variation |
| [Dixon et al. 2015](pdfs/Dixon_2015_GenomicDeterminantsHeatTolerance.pdf) | Expression | Genomic determinants of heat tolerance across latitudes → heritable tolerance |
| [Evensen et al. 2023](pdfs/Evensen_2023_CBASSSystem.pdf) | Thermal/CBASS | CBASS low-cost portable system → LTH assay methods |
| [Fitt et al. 2001](pdfs/Fitt_2001_ThermalToleranceLimitsThresholds.pdf) | Thermal tolerance | Interpretation of thermal tolerance limits & thresholds → bleaching thresholds |
| [Palumbi et al. 2014](pdfs/Palumbi_2014_MechanismsResistance.pdf) | Expression | Mechanisms of reef-coral resistance (A. hyacinthus, Ofu) → expression/acclimation |
| [Parkinson et al. 2018](pdfs/Parkinson_2018_TranscriptionalVariationBiomarker.pdf) | Expression | Transcriptional variation challenges thermal biomarkers → expression variance |
| [Seneca & Palumbi 2010](pdfs/SenecaPalumbi_2010_GeneExpressionNaturalBleaching.pdf) | Expression | Gene expression in a coral undergoing natural bleaching → expression×bleaching |
| [Thomas et al. 2019](pdfs/Thomas_2019_TranscriptomicResilienceShuffling.pdf) | Expression | Transcriptomic resilience + symbiont shuffling → expression + iZ |
| [Voolstra et al. 2020](pdfs/Voolstra_2020_CBASSAcuteHeatAssay.pdf) | Thermal/CBASS | Standardised acute heat-stress assay (CBASS) → LTH assay backbone |
| [Warner et al. 1996](pdfs/Warner_1996_PhotosyntheticEfficiencyZooxanthellae.pdf) | Thermal tolerance | Elevated temp & photosynthetic efficiency of zooxanthellae → photoinhibition |

### Growth (congener) → `v_edge`
| Paper | Theme | Connection — what it grounds |
|---|---|---|
| [Lirman et al. 2014](pdfs/Lirman_2014_AcerGrowth.pdf) | Growth | A. cervicornis growth dynamics → v_edge (congener) |

### Recovery dynamics
| Paper | Theme | Connection — what it grounds |
|---|---|---|
| [Roff et al. 2014](pdfs/Roff_2014_PoritesPhoenixEffectRangiroa.pdf) | Recovery | Porites 'Phoenix effect' recovery, Rangiroa (FP) → recovery dynamics |

### Climate / bleaching context
| Paper | Theme | Connection — what it grounds |
|---|---|---|
| [De'ath et al. 2012](pdfs/Death_2012_GBRCoralDecline.pdf) | Context | 27-year decline of GBR coral cover → context |
| [Hoegh-Guldberg et al. 2007](pdfs/HoeghGuldberg_2007_CoralReefsClimateChange.pdf) | Context | Coral reefs under rapid climate change & OA → intro/discussion framing |
| [Hughes et al. 2017a](pdfs/Hughes_2017_CoralReefsAnthropocene.pdf) | Context | Coral reefs in the Anthropocene → context |
| [Hughes et al. 2017b](pdfs/Hughes_2017_GlobalWarmingMassBleaching.pdf) | Context | Global warming & recurrent mass bleaching → context |

