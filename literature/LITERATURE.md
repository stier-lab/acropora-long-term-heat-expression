# Literature — *Acropora pulchra* (LTH "expression by temperature" + wound-healing model)

Bibliography and synthesis for the LTH heat × wounding project, and for setting *A. pulchra* parameters in the 3-species coral wound-healing model (`~/coral-wound-healing-model`). Merges the former `LITERATURE.md`, `CITATIONS_INDEX.md`, `CITATION_AUDIT.md`, `LIBRARY_MAP.md`, and `KNOWN_UNKNOWN_synthesis.md`. The library holds **101 PDFs** in `pdfs/`; DOIs checked via Crossref/OpenAlex; PDFs named `AuthorYear_ShortTitle.pdf`.

**Species-attribution corrections (checked against source PDFs).** Four papers were misassigned to *A. pulchra* and are flagged ⚠: **Zhang 2025** = the macroalga *Halimeda macroloba* (not a coral); **Almeida 2024** = *Pocillopora* cf. *damicornis* (Kenya); **Denis 2024** = *A. spathulata* (GBR); **Raymundo 2025** = *A. aspera* outplants (Guam). We keep them as background/framing only. Upshot: proteomic and metabolomic thermal acclimation in *A. pulchra* is essentially **unstudied**.

Status legend: `pdfs/…` = local PDF held · **[NB]** = in NotebookLM corpus only · **[⬇]** = still to get (needs UCSB login / paywall / no open-access route) · **[method]** = method or software, no article PDF.

---

# Bibliography

## A. *A. pulchra* — genomics / gene expression / proteomics / microbiome

| Paper | DOI | PDF | Note |
|---|---|---|---|
| Conn et al. 2025 — genome assembly & annotation (Mo'orea) | 10.46471/gigabyte.153 | pdfs/Conn_2025_GenomeAssembly.pdf | reference genome; underpins all expression work (bioRxiv preprint 10.1101/2025.03.27.645822 duplicates) |
| Anthony et al. 2023 — Symbiodiniaceae cellular/phenotypic plasticity | 10.3389/fevo.2023.1288596 | pdfs/Anthony_2023_SymbiontPlasticity.pdf | symbiont expression plasticity; 1.345×10⁶ cells/cm² → `KZ` |
| de la Vega et al. 2023 — *Endozoicomonas* GU-1 genome from *A. pulchra* | 10.1128/mra.01355-22 | pdfs/delaVega_2023_EndozoicomonasGenome.pdf | microbiome resource |
| Lock et al. 2025 — Symbiodiniaceae + bacterial microbiome impacts on host | 10.1111/1462-2920.70175 | pdfs/Lock_2025_SymbiontBacterialMicrobiome.pdf | host performance; cell density drops post-transplant, communities stable |
| Miller & Bentlage 2024 — seasonal tissue/mucus microbiome baseline | 10.7717/peerj.17421 | pdfs/Miller_2024_SeasonalMicrobiome.pdf | microbiome baseline |
| Alina et al. 2023 — low-tide exposure × microbial abundance (intertidal) | 10.1080/17451000.2023.2169464 | pdfs/Alina_2023_IntertidalMicrobialAbundance.pdf | first author Dining Nika Alina (misfiled earlier as Nuñez Lendo) |
| ⚠ Zhang et al. 2025 — protein & metabolite thermal acclimation | 10.3389/fmars.2025.1543591 | pdfs/Zhang_2025_ProteinMetaboliteAcclim.pdf | **macroalga *Halimeda macroloba*, not a coral**; calcifier framing only |

## B. *A. pulchra* — thermal tolerance / bleaching / CBASS

| Paper | DOI | PDF | Note |
|---|---|---|---|
| Cunning et al. 2024 — rapid acute heat-tolerance assays (CBASS) | 10.1007/s00338-024-02577-7 | pdfs/Cunning_2024_RapidAcuteHeatToleranceAssays.pdf | Fv/Fm ED50 34.3–36.6 °C, 20 Mahana/Mo'orea colonies; ED50↔heatwave R=0.74 |
| Berg et al. 2020 — persistent photosystem damage (A. cf. *pulchra*, Guam) | 10.1080/17451000.2021.1875245 | [NB] | 58% mortality d34; photodamage needs ~3× stress to repair |
| Shaw et al. 2016 — intraspecific variability, warming × acidification | 10.1007/s00227-016-2986-8 | pdfs/Shaw_2016_IntraspecificVariability.pdf | genotype variation; growth–tolerance trade-off → `bT_pulse`, `bleach_z` |
| El-Khaled et al. 2025 — heat & cold bleaching vulnerability | 10.1038/s42003-025-09329-5 | pdfs/ElKhaled_2025_HeatColdBleaching.pdf | incl. *A. pulchra* |
| Hoadley et al. 2015 — temp × pCO₂ physiology across genera | 10.1038/srep18371 | pdfs/Hoadley_2015_TempCO2Physiology.pdf | incl. *A. pulchra* |
| Grottoli et al. 2021 — comparability among bleaching experiments | 10.1002/eap.2262 | pdfs/Grottoli_2021_BleachingComparability.pdf | LTH methods/standardization |
| Puisay et al. 2018 — larval thermal resistance/acclimation | 10.1016/j.marenvres.2018.01.005 | [⬇] | larvae; Elsevier paywall |
| Babka 2017 — successive bleaching resilience | 10.60950/bgtl.2017.2.babka.sk | [⬇] | Bremen thesis, bot-blocked |
| ⚠ Denis et al. 2024 — individual thermal-tolerance traits across the GBR | 10.1098/rspb.2024.0587 (preprint 10.1101/2024.01.28.576773) | pdfs/Denis_2024_IndividualThermalTolerance.pdf | **A. spathulata, GBR (not A. pulchra)**; CBASS variation framing |

## C. *A. pulchra* — symbiont community (Mo'orea)

| Paper | DOI | PDF | Note |
|---|---|---|---|
| Rouzé et al. 2017 — Symbiodiniaceae of Mo'orea scleractinians | 10.7717/peerj.2856 | pdfs/Rouze_2017_SymbiodiniumMoorea.pdf | Mo'orea symbiont identity → `iZ` context |
| Rouzé et al. 2019 — quantitative Symbiodiniaceae signature, Mo'orea | 10.1038/s41598-019-44017-5 | pdfs/Rouze_2019_SymbiodiniaceaeSignatureMoorea.pdf | Mo'orea symbiont dynamics → `iZ` |

## D. *A. pulchra* — growth / calcification / skeleton

| Paper | DOI | PDF | Note |
|---|---|---|---|
| Yap & Gomez 1984 — growth of *A. pulchra* | 10.1007/BF00393119 | pdfs/Yap_1984_Growth.pdf | 13.1–15.8 cm/yr; 36.4% branch mortality → `v_edge`, `m_d` |
| Yap 1985 — growth II/III (temp × daylength, transplantation) | 10.1007/BF00539430 | pdfs/Yap_1985_GrowthIII.pdf | warm-season growth suppression → `v_edge` |
| Yap 1981 — growth, regeneration & transplantation (Bolinao) | — (M.Sc. thesis, UP Quezon City) | [⬇] | rare *A. pulchra* regeneration data; no online copy |
| Strømgren 1987 — light × growth, intertidal | 10.1007/BF00302211 | pdfs/Stromgren_1987_LightGrowth.pdf | growth response → `maxpd` |
| Comeau et al. 2014 — irradiance × calcification | 10.1016/j.jembe.2013.12.013 | [NB] | tissue biomass → `E_max`; lib has Comeau 2016 (Elsevier paywall) |
| Comeau et al. 2016 — calcification vs temp & pCO₂ | 10.1007/s00338-016-1425-0 | pdfs/Comeau_2016_CalcificationTempCO2.pdf | thermal optimum ~28 °C; parameterization |
| Roche et al. 2010a — skeletal porosity (X-ray microCT) | 10.1016/j.jembe.2010.10.006 | [⬇] | perforate structure → `D_E`, `resid`; Elsevier paywall |
| Roche et al. 2010b — spatial porosity variation | 10.1007/s00338-010-0679-1 | pdfs/Roche_2010b_SpatialPorosity.pdf | tip porosity 65–80% → `D_E`, `resid` |
| Tian et al. 2025 — internal skeletal hydrodynamics | 10.1016/j.isci.2025.111742 | pdfs/Tian_2025_SkeletalHydrodynamics.pdf | skeletal flow |

## E. *A. pulchra* — physiology / nutrients / reproduction / ecology / recovery

| Paper | DOI | PDF | Note |
|---|---|---|---|
| Conetta 2021 — phenotypic plasticity (URI thesis) | 10.23860/thesis-conetta-dennis-2021 | [NB] | back-reef thickets, restoration workhorse |
| Munk 2024 — host & symbiont physiology during wound regeneration (thesis) | — (M.Sc. thesis) | [NB] | **KEY**: abrasion vs fragmentation × 27.9/29.5 °C, Mo'orea; +1.65 °C robust, biphasic recovery; confirm exact title/institution with author |
| Tanaka et al. 2006 — organic-N translocation/conservation | 10.1016/j.jembe.2006.04.011 | [⬇] | translocation → `D_E`; Elsevier paywall |
| Tanaka et al. 2009 — net DOM release | 10.1016/j.jembe.2009.06.023 | [⬇] | Elsevier paywall |
| Buckingham et al. 2022 — N+P enrichment, skewed N:P | 10.1007/s00338-022-02223-0 | pdfs/Buckingham_2022_NPEnrichment.pdf | nutrient stoichiometry |
| Huang 2009 — oocyte development histology (Sanya) | — (no DOI) | [⬇] | only *A. pulchra* reproduction source found |
| Ladd et al. 2025 — growth–predation tradeoffs, distribution | 10.1038/s41598-025-21028-z | pdfs/Ladd_2025_GrowthPredationTradeoffs.pdf | distribution constraints |
| ⚠ Almeida et al. 2024 — intertidal vs subtidal recovery/growth/survival | 10.1007/s00227-024-04546-8 | pdfs/Almeida_2024_IntertidalSubtidalRecovery.pdf | **Pocillopora cf. damicornis, Kenya (not A. pulchra)**; recovery framing |
| ⚠ Raymundo et al. 2025 — restoration in a stressful environment (Guam) | 10.1016/j.isci.2025.112244 | pdfs/Raymundo_2025_RestorationGuam.pdf | **A. aspera outplants (not A. pulchra)**; pulchra only background congener |

## F. Wound healing / regeneration — grounds `kc_mult`

| Paper | DOI | PDF | Note |
|---|---|---|---|
| Traylor-Knowles et al. 2016 — wound healing, two temperature regimes | 10.1007/s00227-016-3011-y | pdfs/Traylor-Knowles_2016_WoundHealingTwoTempRegimes.pdf | ⭐ *P. damicornis* & *A. hyacinthus* (model's 2 species) → temp-coupling; DOI corrected (was 10.1186/s12862-016-0791-0 typo) |
| Xu et al. 2023 — wound healing & regeneration in *Acropora millepora* | 10.3389/fevo.2022.979278 | pdfs/Xu_2023_WoundHealingRegenAcroporaMillepora.pdf | ⭐⭐ coral wound/regeneration transcriptomics; backbone of the candidate-gene panel (c-Fos, JNK, Wnt, FGF, ADAMTS, Galaxin) — see `docs/rnaseq/` |
| Palmer et al. 2011 — corals use similar immune cells & wound-healing processes | 10.1371/journal.pone.0023992 | pdfs/Palmer_2011_CoralsImmuneCellsWoundHealing.pdf | four-phase wound cellular timeline (immune infiltration → proliferation → maturation) |
| Bonesso et al. 2017 — sub-bleaching SST impairs regeneration | 10.7717/peerj.3719 | pdfs/Bonesso_2017_SubBleachingImpairsRegen.pdf | 32 °C halts apical-tip regrowth (*A. aspera*) → `kc_mult`(temp) |
| Burmester et al. 2017 — temperature & symbiosis × lesion recovery | 10.3354/meps12114 | pdfs/Burmester_2017_TempSymbiosisLesionRecovery.pdf | symbiotic state × temp × healing |
| Burmester et al. 2018 — autotrophy vs heterotrophy in lesion recovery | — | pdfs/Burmester_2018_AutoVsHeteroLesionRecovery.pdf | → `kc_mult`, `h` |
| Dias et al. 2018 — mortality/growth/regen after fragmentation under heat | 10.1016/j.seares.2018.08.008 | pdfs/Dias_2018_FragmentationRegenThermalStress.pdf | regen often rose to 32 °C → `m_d`, `kc_mult` |
| Edmunds & Lenihan 2010 — sub-lethal damage to juvenile *Porites* × temp/flow | 10.1007/s00227-009-1372-1 | pdfs/EdmundsLenihan_2010_SublethalDamageJuvenilePorites.pdf | model species damage × temp |
| Meesters & Bak 1993 — bleaching reduces tissue-regeneration potential | 10.3354/meps096189 | pdfs/Meesters_1993_BleachingTissueRegen.pdf | ⭐ bleaching × regeneration |
| Meesters et al. 1994 — damage & regeneration ↔ growth (*Montastraea*) | 10.3354/meps112119 | pdfs/Meesters_1994_DamageRegenGrowth.pdf | growth↔regeneration trade-off |
| Meesters et al. 1997 — partial mortality; lesion-area model | 10.3354/meps146091 | pdfs/Meesters_1997_PredictingRegeneration.pdf | → `kc_mult`, `m_d` |
| Denis et al. 2011 — lesion regeneration, massive *Porites* populations | — | pdfs/Denis_2011_LesionRegenerationPopulations.pdf | → `kc_mult` |
| Denis et al. 2013 — fast growth impairs regeneration (*A. muricata*) | 10.1371/journal.pone.0072618 | pdfs/Denis_2013_FastGrowthImpairRegen.pdf | fast-fragile axis, `E_max` |
| Hall 1997 — interspecific differences in injury regeneration | — | pdfs/Hall_1997_InterspecificRegeneration.pdf | → `kc_mult` |
| Hall 1998 — injury & regeneration of common reef corals | — | pdfs/Hall_1998_InjuryRegeneration.pdf | → `kc_mult` |
| Hall 2001 — *A. hyacinthus* injury response | — | pdfs/Hall_2001_AcroporaHyacinthusResponse.pdf | → `kc_mult` |
| Hall 2015 — lesion recovery under low pH | — | pdfs/Hall_2015_LesionRecovery.pdf | → `kc_mult` |
| Lirman 2000b — lesion regeneration in branching *A. palmata* | 10.3354/meps197209 | pdfs/Lirman_2000b_LesionRegenAcroporaPalmata.pdf | lesion-area recovery |
| DeFilippo et al. 2016 — surface lesion recovery (*Astrangia poculata*) | 10.1016/j.jembe.2016.03.016 | pdfs/DeFilippo_2016_SurfaceLesionRecoveryAstrangia.pdf | lesion kinetics |
| Counsell et al. 2019 — colony size & depth affect wound repair | 10.1007/s00338-019-01807-7 | pdfs/Counsell_2019_ColonySizeDepthWoundRepair.pdf | size-dependent healing |
| Rice et al. 2019 — N sources speed corallivory recovery + microbiome | 10.7717/peerj.8056 | pdfs/Rice_2019_NitrogenCorallivoryRecovery.pdf | wound healing −66% at 29 °C (browsing axis) |
| Bak 1983 — neoplasia, regeneration & growth (*A. palmata*) | 10.1007/BF00395810 | pdfs/Bak_1983_NeoplasiaRegenGrowthApalmata.pdf | foundational *Acropora* regeneration |
| Bak & Steward-Van Es 1980 — regeneration of superficial damage | — (Bull Mar Sci, no DOI) | pdfs/Bak_1980_RegenSuperficialDamage.pdf | foundational regeneration |
| Henry & Hart 2005 — regeneration & resource allocation (review) | 10.1002/iroh.200410759 | pdfs/Henry_2005_RegenInjuryResourceAllocation.pdf | energy-allocation framing |
| Van de Water et al. 2015 — elevated temp & immune response after damage | 10.1007/s10750-015-2243-z | pdfs/VanDeWater_2015_ImmuneResponsePhysicalDamage.pdf | immune/healing × temp |
| Edmunds & Yarid 2017 — ocean acidification & wound repair in *Porites* | 10.1016/j.jembe.2016.10.001 | [⬇] | ScienceDirect blocks automation |

## G. Mechanical damage / breakage / microfragmentation — grounds `m_d`, Tier-3

| Paper | DOI | PDF | Note |
|---|---|---|---|
| Madin et al. 2014 — mechanical vulnerability | 10.1111/ele.12306 | pdfs/Madin_2014_MechanicalVulnerability.pdf | size-dependent `m_d` |
| Madin & Connolly 2006 — hydrodynamic disturbance consequences | 10.1038/nature05328 | pdfs/Madin_2006_HydrodynamicDisturbances.pdf | → `m_d` |
| Marshall 2000 — skeletal damage vs colony morphology | 10.3354/meps200177 | pdfs/Marshall_2000_SkeletalDamageColonyMorphology.pdf | breakage resistance → `m_d` |
| Tunnicliffe 1981 — breakage & propagation of *A. cervicornis* | 10.1073/pnas.78.4.2427 | pdfs/Tunnicliffe_1981_BreakagePropagationAcervicornis.pdf | fragmentation → `m_d` |
| Highsmith 1982 — reproduction by fragmentation | 10.3354/meps007207 | [⬇] | int-res blocks automation |
| Lirman 2000a — fragmentation in branching *Acropora* | — | pdfs/Lirman_2000a_AcroporaFragmentation.pdf | → `m_d` |
| Lock et al. 2022 — Ca-homeostasis drives growth after microfragmentation | 10.1002/ece3.9345 | pdfs/Lock_2022_CaHomeostasisMicrofragmentation.pdf | ⭐ *Porites lobata* → Tier-3 microfrag |
| Page 2018 — microfragmenting slow-growing massive corals | 10.1016/j.ecoleng.2018.08.017 | [⬇] | ScienceDirect blocks automation |
| Soong & Chen 2003 — regeneration & growth of *Acropora* fragments (incl. pulchra) | 10.1046/j.1526-100x.2003.00100.x | pdfs/Soong_2003_AcroporaFragmentNursery.pdf | fragment regen/growth → `kc_mult` |

## H. Energetics / reserves / heterotrophy — grounds `h`, `E_max`, Tier-3

| Paper | DOI | PDF | Note |
|---|---|---|---|
| Conti-Jerpe & Baker 2020 — trophic strategy & bleaching resistance | 10.1126/sciadv.aaz5443 | pdfs/Conti-Jerpe_2020_TrophicStrategyBleaching.pdf | *Acropora* 94% autotroph → `h` (low plasticity) |
| Houlbrèque & Ferrier-Pagès 2009 — heterotrophy in corals (review) | 10.1111/j.1469-185X.2008.00058.x | [⬇] | Wiley paywall; use Conti-Jerpe + Radice |
| Radice et al. 2019 — trophic strategy via FA + SIA | 10.1371/journal.pone.0222327 | pdfs/Radice_2019_HeterotrophyFASIA.pdf | → `h` |
| Grottoli et al. 2006 — heterotrophic plasticity & resilience | 10.1038/nature04565 | pdfs/Grottoli_2006_HeterotrophicPlasticity.pdf | ⭐ foundational `h` buffer |
| Grottoli et al. 2014 — cumulative annual bleaching, winners→losers | 10.1111/gcb.12658 | pdfs/Grottoli_2014_CumulativeBleachingWinnersLosers.pdf | reserve depletion |
| Schoepf et al. 2013 — *Acropora* energy reserves & calcification | 10.1371/journal.pone.0075049 | pdfs/Schoepf_2013_AcroporaReserves.pdf | 6.6 vs Poc 8.1 kJ/g → `E_max` |
| Leinbach et al. 2021 — recovery costs | 10.1038/s41598-021-02807-w | pdfs/Leinbach_2021_RecoveryCosts.pdf | reserve depletion (Mo'orea) → `E0`, `E_max` |
| Rodrigues & Grottoli 2007 — energy reserves as recovery indicators | 10.4319/lo.2007.52.5.1874 | pdfs/Rodrigues_2007_EnergyReservesRecovery.pdf | → `E_max` |
| Anthony et al. 2007 — bleaching, energetics & mortality risk | 10.4319/lo.2007.52.2.0716 | pdfs/Anthony_2007_BleachingEnergeticsMortality.pdf | energy-budget mortality |
| Anthony et al. 2009 — energetics approach to mortality risk | 10.1111/j.1365-2435.2008.01531.x | [⬇] | companion Anthony 2007 covers ground |
| Levas et al. 2013 — bleaching/recovery physiology (*Porites lobata*) | 10.1371/journal.pone.0063267 | pdfs/Levas_2013_BleachingRecoveryPoritesLobata.pdf | *Porites* recovery |
| Rinkevich 1996 — reproduction vs regeneration compete for energy | 10.3354/meps143297 | pdfs/Rinkevich_1996_ReproRegenEnergyAllocation.pdf | ⭐ Tier-3 energy budget |
| Oren et al. 1997 — oriented ¹⁴C transport (perforate) | 10.3354/meps161117 | [⬇] | translocation → `D_E` |
| Fine & Loya 2002 — bleaching cuts translocation | 10.3354/meps234119 | [⬇] | → `D_E` |
| Swain et al. 2018 — colony integration | 10.3354/meps12445 | [⬇] | → `D_E` |

## I. Reproduction / size at maturity — grounds `B_mat`

| Paper | DOI | PDF | Note |
|---|---|---|---|
| Carroll et al. 2006 — sexual reproduction of *Acropora* at Mo'orea | 10.1007/s00338-005-0057-6 | pdfs/Carroll_2006_AcroporaReproductionMoorea.pdf | ⭐ *Acropora* reproduction AT field site → `B_mat` |
| Wallace 1985 — reproduction/recruitment/fragmentation, nine *Acropora* | 10.1007/BF00392585 | pdfs/Wallace_1985_AcroporaReproductionRecruitment.pdf | fecundity × size → `B_mat` |
| Baird et al. 2009 — reproductive biology of scleractinians (review) | 10.1146/annurev.ecolsys.110308.120220 | pdfs/Baird_2009_ReproductiveBiologyScleractinian.pdf | reproduction framing |
| Soong & Lang 1992 — reproductive integration | — | pdfs/Soong_1992_ReproductiveIntegration.pdf | size effects → `B_mat` |
| Omori et al. 2001 — drop in *Acropora* fertilization after bleaching | 10.4319/lo.2001.46.3.0704 | pdfs/Omori_2001_AcroporaFertilizationDrop.pdf | repro × bleaching |
| Álvarez-Noriega et al. 2016 — fecundity × morphology | 10.1002/ecy.1588 | [⬇] | Wiley paywall; use Soong 1992 |
| Hall & Hughes 1996 — modular reproduction | 10.2307/2265514 | [⬇] | JSTOR paywall; use Soong 1992 |

## J. Symbiont identity / shuffling / recovery threshold — grounds `iZ`

| Paper | DOI | PDF | Note |
|---|---|---|---|
| Bay et al. 2016 — *A. millepora* recovery threshold | 10.1098/rsos.160322 | pdfs/Bay_2016_RecoveryThreshold.pdf | symbiont-density threshold → `iZ`, `Z0` |
| Baker 2003 — flexibility/specificity in coral–algal symbiosis (review) | 10.1146/annurev.ecolsys.34.011802.132417 | pdfs/Baker_2003_CoralAlgalSymbiosisFlexibility.pdf | `iZ` framing |
| Berkelmans & van Oppen 2006 — zooxanthellae role in thermal tolerance | 10.1098/rspb.2006.3567 | pdfs/Berkelmans_2006_ZooxanthellaeThermalTolerance.pdf | ⭐ symbiont-mediated tolerance |
| Jones et al. 2008 — endosymbiont community change after bleaching | 10.1098/rspb.2008.0069 | pdfs/Jones_2008_EndosymbiontCommunityChange.pdf | ⭐ symbiont shuffling → `iZ` |
| Silverstein et al. 2017 — clade D *Symbiodinium* persist at temp extremes | 10.1242/jeb.148239 | [⬇] | thermotolerant symbiont retention |
| LaJeunesse et al. 2018 — systematic revision of Symbiodiniaceae | 10.1016/j.cub.2018.07.008 | pdfs/LaJeunesse_2018_SymbiodiniaceaeRevision.pdf | symbiont taxonomy |
| Fitt et al. 2000 — seasonal symbiont density | 10.4319/lo.2000.45.3.0677 | [⬇] | Wiley L&O paywall; use Bay 2016 |

## K. Thermal-tolerance & gene-expression backbone (LTH "expression by temperature")

| Paper | DOI | PDF | Note |
|---|---|---|---|
| Voolstra et al. 2020 — standardized acute heat-stress assay (CBASS) | 10.1111/gcb.15148 | pdfs/Voolstra_2020_CBASSAcuteHeatAssay.pdf | ⭐⭐ CBASS standard; LTH assay backbone |
| Evensen et al. 2023 — CBASS low-cost portable system | 10.1002/lom3.10555 | pdfs/Evensen_2023_CBASSSystem.pdf | CBASS methods |
| Cunning et al. 2021 — heat-tolerance census, Florida staghorn *Acropora* | 10.1098/rspb.2021.1613 | pdfs/Cunning_2021_HeatToleranceCensusStaghorn.pdf | acute-tolerance variation |
| Barshis et al. 2013 — genomic basis for coral resilience | 10.1073/pnas.1210224110 | pdfs/Barshis_2013_GenomicBasisResilience.pdf | ⭐⭐ expression basis of resilience |
| Cleves et al. 2020 — CRISPR HSF mutants reduce coral thermal tolerance | 10.1073/pnas.1920779117 | pdfs/Cleves_2020_CRISPRHSFThermalTolerance.pdf | causal HSF → thermal tolerance; candidate-gene source for HSF |
| Bay & Palumbi 2015 — rapid acclimation via transcriptome change | 10.1093/gbe/evv085 | pdfs/Bay_2015_RapidAcclimTranscriptome.pdf | expression plasticity |
| Dixon et al. 2015 — genomic determinants of heat tolerance across latitudes | 10.1126/science.1261224 | pdfs/Dixon_2015_GenomicDeterminantsHeatTolerance.pdf | heritable heat tolerance |
| Palumbi et al. 2014 — mechanisms of reef-coral resistance | 10.1126/science.1251336 | pdfs/Palumbi_2014_MechanismsResistance.pdf | ⭐ *A. hyacinthus* expression/acclimation |
| Seneca & Palumbi 2010 — gene expression during natural bleaching | 10.1007/s10126-009-9247-5 | pdfs/Seneca_2010_GeneExpressionNaturalBleaching.pdf | expression × bleaching |
| Parkinson et al. 2018 — transcriptional variation challenges biomarkers | 10.1111/mec.14517 | pdfs/Parkinson_2018_TranscriptionalVariationBiomarker.pdf | expression variance |
| Thomas et al. 2019 — transcriptomic resilience + symbiont shuffling | 10.1111/mec.15143 | pdfs/Thomas_2019_TranscriptomicResilienceShuffling.pdf | expression + `iZ`; Mo'orea rebuild ~4 mo |
| Cornwell et al. 2021 — heat tolerance × symbiont load → growth trade-off | 10.7554/eLife.64790 | pdfs/Cornwell_2021_HeatToleranceSymbiontLoadTradeoff.pdf | ⭐ *A. hyacinthus* growth–tolerance |
| Fitt et al. 2001 — thermal tolerance limits & thresholds | 10.1007/s003380100146 | pdfs/Fitt_2001_ThermalToleranceLimitsThresholds.pdf | foundational bleaching thresholds |
| Warner et al. 1996 — elevated temp & photosynthetic efficiency of zooxanthellae | 10.1046/j.1365-3040.1996.d01-12.x | pdfs/Warner_1996_PhotosyntheticEfficiencyZooxanthellae.pdf | foundational photoinhibition |
| Warner, Fitt & Schmidt 1999 — PSII damage as bleaching determinant | 10.1073/pnas.96.14.8007 | [⬇] | a-priori cite (Fv/Fm decline); PNAS OA, needs browser session; lib has Warner 1996 |

## L. Growth (congener) & recovery dynamics

| Paper | DOI | PDF | Note |
|---|---|---|---|
| Lirman et al. 2014 — *A. cervicornis* growth dynamics | 10.1371/journal.pone.0107253 | pdfs/Lirman_2014_AcerGrowth.pdf | → `v_edge` (congener) |
| Roff et al. 2014 — *Porites* "Phoenix effect" recovery, Rangiroa (FP) | 10.1007/s00227-014-2426-6 | pdfs/Roff_2014_PoritesPhoenixEffectRangiroa.pdf | dramatic FP recovery |

## M. Climate / bleaching context (intro/discussion framing)

| Paper | DOI | PDF | Note |
|---|---|---|---|
| Hoegh-Guldberg et al. 2007 — coral reefs under rapid climate change & OA | 10.1126/science.1152509 | pdfs/HoeghGuldberg_2007_CoralReefsClimateChange.pdf | classic intro cite |
| Hoegh-Guldberg 1999 — climate change, bleaching & the future | 10.1071/MF99078 | [⬇] | a-priori cite (paling, symbiont loss); CSIRO paywall; lib has HG 2007 |
| Hughes et al. 2017a — coral reefs in the Anthropocene | 10.1038/nature22901 | pdfs/Hughes_2017_CoralReefsAnthropocene.pdf | context |
| Hughes et al. 2017b — global warming & recurrent mass bleaching | 10.1038/nature21707 | pdfs/Hughes_2017_GlobalWarmingMassBleaching.pdf | context |
| De'ath et al. 2012 — 27-year decline of GBR coral cover | 10.1073/pnas.1208909109 | pdfs/DeAth_2012_GBRCoralDecline.pdf | context |
| Buddemeier & Fautin 1993 — bleaching as an adaptive mechanism | 10.2307/1312064 | [⬇] | Adaptive Bleaching Hypothesis; JSTOR captcha |
| Jokiel & Coles 1977 — temperature, mortality & growth, Hawaiian corals | 10.1007/BF00402312 | pdfs/Jokiel_1977_TemperatureMortalityGrowthHawaiianCorals.pdf | a-priori cite (heat reduces calcification) |
| Jokiel & Coles 1990 — Indo-Pacific response to elevated temperature | 10.1007/BF00265006 | pdfs/Jokiel_1990_ResponseIndoPacificElevatedTemperature.pdf | a-priori cite (symbiont loss) |

## N. Methods / software

| Paper | DOI | PDF | Note |
|---|---|---|---|
| Davies 1989 — buoyant-weight technique | 10.1007/BF00428135 | pdfs/Davies_1989_BuoyantWeightTechnique.pdf | calcification method |
| Siebeck et al. 2006 — colour reference card for bleaching | 10.1007/s00338-006-0123-8 | pdfs/Siebeck_2006_ColourReferenceCard.pdf | bleaching-score method |
| Jokiel et al. 1978 — buoyant-weight technique (UNESCO) | — (book chapter) | [method] | method |
| R Core Team 2025 — R statistical environment | — (software) | [method] | software |

---

# Model-parameter grounding (*A. pulchra* preset)

A fast-growing, fragile, heat-sensitive branching *Acropora*. Values come from this library, set against the model's *Porites* and *Pocillopora* anchors. Full rows (value · units · citation · source PDF · verification) are in `~/coral-wound-healing-model/data/phase2_parameters.csv` (15 *A. pulchra* rows).

| Param | Value | Anchor (Por / Poc) | Key evidence | Confidence |
|---|---|---|---|---|
| `v_edge` | 0.045 | 0.010 / 0.030 | 13.1–15.8 cm/yr [Yap & Gomez 1984] | HIGH ⭐sp. |
| `m_d` | 0.0045 | 0.0010 / 0.0030 | 36.4% branch mortality/yr [Yap & Gomez 1984] | HIGH ⭐sp. |
| `kc_mult` | 1.2 | 0.7 / 1.0 | Acropora>Poc>Por [Hall 1997/2001]; Munk ~92%/19d | MED |
| `h` | 0.0002 | 0.0004 / 0.00025 | 94% host-symbiont overlap, top autotroph [Conti-Jerpe 2020] | MED-HIGH rank |
| `E_max` | 0.85 | 2.0 / 1.0 | *A. millepora* 6.6 vs Poc 8.1 kJ/g [Schoepf 2013] | HIGH rank |
| `E0` | 0.77 | 1.5 / 0.9 | 0.90×E_max; deep depletion [Leinbach 2021] | LOW-MED |
| `bT_pulse` | 0.010 | 0.0017 / 0.0080 | ED50 ~34–37 °C [Cunning 2024; Denis 2024]; 58% mort d34 [Berg] | LOW-MED |
| `bleach_z` | 0.15 | 0.55 / 0.20 | obligate autotroph [Shaw 2016; Denis 2024] | LOW-MED |
| `iZ` | 0.004 | 0.020 / 0.004 | Mo'orea Acropora rebuild ~4 mo [Thomas 2019] | MED (best) |
| `Z0` | 0.0 | 0.3 / 0.0 | bleaches white, S:H→0.003 [Bay 2016] | MED |
| `KZ` | 1.0 | 1.0 / 1.0 | normalization; 1.345×10⁶ cells/cm² [Anthony 2023] | HIGH |
| `D_E` | 0.014 | 0.02 / 0.002 | tip porosity 65–80%, connected interior [Roche 2010b] | HIGH ⭐sp. |
| `resid` | 0.07 | 0.50 / 0.05 | thin perforate branches [Roche 2010b] | MED |
| `maxpd` | 1.0 | 3.0 / 1.2 | 16 cm branches, 15–25 mm tips [Yap & Gomez; Strømgren] | LOW-MED |
| `B_mat` | 4.5 | 8.0 / 4.0 | *A. cervicornis* puberty 17 cm [Soong 1992]; recruits [Wallace] | LOW-MED |

⭐sp. = value specific to *A. pulchra*; the rest are borrowed from congeners. **Weakest support (sweep these):** `E0`, `maxpd`, `bleach_z`/`bT_pulse` — the local PDFs give no *A. pulchra* loss rate or ED50, so these are inferred from congener thresholds plus Cunning 2024 / Berg 2020. Have a biologist review this table before committing the `Apulchra` preset to `ext/species_presets.R`.

---

# What we know / what we don't

**One-paragraph version.** The literature agrees heat harms *A. pulchra* physiology but disagrees on whether heat blocks wound recovery. It has almost never (i) separated tissue healing from skeletal regeneration in the same coral, or (ii) tied recovery to a genotype's measured thermal tolerance. LTH targets that open core. The disagreement looks like a matter of **dose**: at +1.65 °C *A. pulchra* recovers fine (Munk 2024); at +3 °C / 31 °C tissue heals but skeleton fails to regrow (LTH); at sub-bleaching +5–6 °C skeletal regrowth stops in a congener (Bonesso 2017, *A. aspera*).

**Wounding / healing / regeneration — known.** Branches grow 13.1–15.8 cm/yr; broken branches regrow ~12.6 cm/yr if they survive, but 36.4% die (no warm-season growth) (Yap & Gomez 1984). Recovery is biphasic (Munk 2024): 92% (13/14) of fragments showed ≥2 of 3 signs by day 19, and tissue regrew faster after abrasion than fragmentation (the "phoenix effect"). A restoration workhorse forming back-reef thickets in Mo'orea (Conetta 2021; Soong & Chen 2003). Congener precedent: Bonesso 2017 (*A. aspera*) — 32 °C stopped tip regrowth and shut down wound GFP.

**Wounding — not known.** No fine-scale tissue-closure rates (mm²/d) for *A. pulchra*. Tissue architecture, reserve thickness, and internal translocation distance are unmeasured (the model's weak `D_E`/`resid`). How much heterotrophic feeding compensates during repair is unknown. The molecular basis of healing vs regeneration is pending (LTH RNA-seq).

**Thermal tolerance — known.** CBASS Fv/Fm ED50 = 34.3–36.6 °C for the exact Mahana/Mo'orea population (Cunning 2024). Photodamage takes ~3× the stress duration to repair and is worst during recovery (Berg 2020). Calcification peaks near 28 °C; at 29.8 °C net calcification drops 18–50% (Comeau 2014/2016). The growth–tolerance trade-off is heritable (Shaw 2016). Symbiont basis: a stable, thermotolerant partnership of *Cladocopium* C40 + *Durusdinium* D1 (Anthony 2023; Rouzé 2017/2019; Lock 2025).

**Thermal tolerance — not known.** Almost no data on proteomic/metabolomic thermal acclimation (the earlier "Zhang 2025" source turned out to be a macroalga). No landscape-scale ED50 mapping (it exists for *A. spathulata*, not *A. pulchra*). Whether ED50 predicts anything beyond Fv/Fm — calcification or regeneration under chronic heat — is unknown.

**The heat × wounding interaction — the unresolved core.**
- *Heat impairs recovery:* Meesters & Bak 1993 (bleaching blocks translocation); Bonesso 2017 (32 °C halts skeleton); Rice 2019 (*Pocillopora* Mo'orea, healing −66% at 29 °C); Paradis 2019 (*A. cervicornis*, P:R < 1); Kaufman 2021 (*A. cervicornis*, 35% healed at 31.5 °C vs 99% at 28 °C, donor-history effect).
- *Heat doesn't impair / enhances:* Dias 2018 (regen rose to 32 °C); Burmester 2017 (temperate *Astrangia*, cold is the constraint); Munk 2024 (*A. pulchra* +1.65 °C robust); Traylor-Knowles 2016 (no regime difference, *A. hyacinthus*).
- *How they reconcile (a dose effect):* +1.65 °C is fine → at +3 °C/31 °C tissue still heals but skeleton stops regrowing → at +5–6 °C the skeleton halts. No one has mapped where, for *A. pulchra*, tolerable warming becomes regeneration-blocking heat; the LTH 31 °C treatment sits on that line.
- *Possible mechanisms:* symbiont loss blocks translocation (Fine 2002; Meesters & Bak 1993); negative energy balance (Paradis 2019); the cost of immune, antioxidant, and HSP responses (Madeira 2022; Lock 2022); and a growth-vs-regeneration trade-off (Rice 2019; Denis 2013).

**What LTH adds.** (1) Phase-specificity — heat splits healing from regeneration in the same coral (tissue heals on schedule, but 67% of wounds closed and never regrew skeleton at 31 °C vs 0% at 28 °C; new-corallite Cox HR 0.22). (2) The ED50-to-regeneration link — pending SNP-matching of thickets A/C/D to Cunning's genotyped genets. (3) Mechanism (RNA-seq, pending) and field realism (healing under corallivory and competition as well as heat).

*(Papers cited in the synthesis but not in the local library: Paradis 2019, Kaufman 2021, Madeira 2022 — no DOI/PDF tracked here.)*

---

# Gaps (where the literature is genuinely thin)
- **Reproduction / fecundity / size at maturity** — only Huang 2009 (histology) and Munk-adjacent data; we lean on congeners (Wallace 1985, Carroll 2006, Soong & Lang 1992, Álvarez-Noriega 2016, Hall & Hughes 1996).
- **Wounding / regeneration rates** — only Munk 2024 (thesis) and Yap 1981 (thesis, not obtained); we lean on Denis 2011/2013, the Hall series, Lirman, and Bonesso.
- **Symbiont repopulation rate after bleaching** — none found (only density thresholds; Bay 2016, Jones 2008), making `iZ` the least-supported parameter; borrow from Symbiodiniaceae division-rate studies or fit it.
- **Proteomic/metabolomic thermal acclimation** — no *A. pulchra* data.
