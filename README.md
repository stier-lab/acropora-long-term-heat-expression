# LTH: 15-day heat x wound experiment with RNA-seq samples in *Acropora pulchra*

**October 27 update:** The results summary and analyses include fragment 119's confirmed wounded assignment. Figure 2 now shows side-by-side wound-group boxplots for symbiont density and growth, with standalone panels. The heat-by-wound conclusion is unchanged. See `docs/provenance/summary_refresh_2026-10-27.md` and the rebuilt `docs/team_summary/LTH_results_summary.html`. Molly's Google Doc has not been edited.

**October 13 update:** Rachael's final host SNP clusters and rough symbiont
typing are integrated through `code/06_genetics.R`. Current
sequencing covariates: `output/tables/43_final_snp_rnaseq_covariates.csv`.
All 144 assignments retain their preliminary membership; genetic coordinates
are updated. These are Days 1/3/10 samples, not Day 15. The September `32_*`
outputs are historical. Expression results are still pending.

> Project #17 · Gump Station, Mo'orea, French Polynesia · May–July 2025 · Repo: `stier-lab/regen-2025-acropora-lth-expression` · Updated 2026-10-27

<p align="center">
  <img src="docs/team_summary/img/healing_ambient_vs_heated.gif" width="720"
       alt="Wound-healing time-lapse, ambient 28 °C vs heated 31 °C, Day 0 to Day 15"><br>
  <sub><b>Wound healing over 15 days.</b> <b>Ambient&nbsp;28&nbsp;°C</b> (left, blue) closes the wound <b>and</b> regrows the tip;
  <b>Heated&nbsp;31&nbsp;°C</b> (right, orange) re-covers the wound and rebuilds a tip, but has not budded new corallites by Day 15.</sub>
</p>

A 15-day heat x wound experiment on the branching coral *A. pulchra*: we clipped about 1 cm from the growing tip of half the fragments, held corals at **28 °C** or **31 °C (+3 °C)**, and tracked tissue closure, skeletal regrowth, growth, physiology, and symbionts for 15 days. We also collected tissue for RNA-seq.

**The phenotype result:** all 16 fragments in the separate microscope cohort had tissue covering the wound by Day 1. In the main experiment, new corallites were observed by Day 15 in 12/12 ambient and 4/12 heated wounded fragments. These separate cohorts do not establish a within-fragment sequence, and we do not know whether the remaining heated fragments would have budded later. Final host SNP assignments are integrated for 144 early samples; phenotype associations remain exploratory and expression results are pending.

**Which outputs are current:** `docs/OUTPUTS_STATUS.md` (also explains why files show as changed after a run). **More detail:** plain-language team summary → `docs/team_summary/LTH_results_summary.html` · analysis/diagnostic inventory → `docs/analysis_diagnostic_inventory.md` · next-analysis roadmap → `docs/next_analysis_roadmap.md` · full results → `RESULTS.docx` · figures → `figures/FIGURE_INDEX.md` · data → `data/DATA_DICTIONARY.docx` · docs map → `docs/README.md` · RNA-seq notes → `docs/rnaseq/`.

**Terms used consistently here:**

- **Wounded** = about 1 cm removed from the growing tip on D0 with a band saw.
- **No-wound control** = fragments assigned to the no-wound level. The exact sham/handling procedure still needs to be verified before manuscript wording calls this a "sham" control.
- **Source patch label** = field labels A, C, and D, stored as `thicket` or legacy `genet` columns. Final SNP assignments place all 48 early C samples in one cluster; A and D contain multiple clusters. The Day 15 fragments remain ungenotyped. Source labels are not genetic assignments or confirmed matches to the external Hauru colonies.
- **Healing vs regeneration** = healing is early tissue re-covering of the wound; regeneration is later rebuilding of structures such as the tip and new radial corallites.
- **Species wording** = the repo historically uses *Acropora pulchra* and `apulchra` filenames. Trinity Conn's 2026-08-10 note recommends *Acropora cf. pulchra* for the Mo'orea animal until the taxonomy is settled, because the Mo'orea form may not be the Australian type specimen of *A. pulchra*.

## Quick start

```bash
cd ~/Stier-LTH-expression-by-temperature-2025
Rscript -e 'install.packages(c("tidyverse","here","readxl","janitor","scales","patchwork","broom","broom.mixed","lme4","lmerTest","emmeans","DHARMa","survival","xml2","rmarkdown"))'
Rscript code/run_all.R          # under a minute; rebuilds every number, table and figure the paper uses
```

New to the code? Start with `code/README.md`. Ten files in `code/` make everything the paper uses; scripts 02–05 each end by writing one table of reported numbers to `output/paper/`. Older and exploratory analyses are frozen in `extras/`.

Requires R ≥ 4.3. Everything regenerates from `code/run_all.R`; **never hand-edit numbers**: every number reported in the text is in `output/paper/`. (The older all-statistics table, `20_master_results.csv`, and its manuscript audit are frozen in `extras/`.)

## Design

| Factor | Levels |
|---|---|
| **Temperature** | 28 °C (ambient) vs 31 °C (heated). Ramped 1 °C/day. |
| **Wound** | Wounded (~1 cm clipped from growing tip) vs no-wound control. Applied 7 days after reaching target temp. |
| **Source patch label** | A, C, D — fixed effect (only 3 field source-patch labels; stored as `thicket` or legacy `genet` in some code/tables). |
| **Tank** | 8 total, 4 per temperature (28 °C: 3, 6, 9, 12; 31 °C: 4, 5, 10, 11). Random effect. |
| **Time** | Daily observations; destructive tissue samples D0, D1, D3, D10, D15. |

**n = 208 fragments total.** The main destructive sampling set has 192 fragments; 48 of those were also tracked non-destructively for physiology through Day 15; 144 selected margin libraries form the current RNA-seq design (D1, D3, D10); and a separate 16-fragment photo-only set was photographed for microscope time-lapse imagery. **Site:** Mahana / Tiahura, NW Mo'orea — source patches A (17.49735 °S / 149.91557 °W), C (17.49808 °S / 149.91595 °W), D (17.49726 °S / 149.91581 °W), the same reef as Cunning et al. 2024's CBASS genets. Molly's forwarded Trinity/Ross Hauru table is now stored in `data/external/` as preliminary external context; source patch C's closest coordinate candidate is `Apul-115` (10.5 m away), but exact identity requires SNP comparison.

## Key findings (phenotype half)

Current results and limitations are in `docs/team_summary/LTH_results_summary.html`. `RESULTS.docx` and the local manuscript draft are older documents and have not been reconciled in this refresh. Analysis-specific tables in `output/tables/` provide the numerical record.

1. **Heat reduced coral condition.** Photosynthesis efficiency, color, and symbiont density declined under heat. Mean skeletal mass gain was 34% lower (6.10% versus 4.03%), but the tank-level growth test remains uncertain (p = 0.057).
2. **Wound covering did not imply later regrowth.** Tissue cover was recorded by Day 1 in all 16 microscope fragments. In the separate main cohort, new corallites were observed by Day 15 in 12/12 ambient versus 4/12 heated wounded fragments. This does not mean heat permanently stopped regeneration or that the others never rebuilt skeleton.
3. **Source patch C was less affected by heat.** The clearest evidence is in photosynthesis efficiency, color, and symbiont density. Genetic identity is not established for the endpoint fragments. See `figures/19_genet_dashboard.pdf` and `figures/43_final_snp_structure.pdf`.
4. **Exposure duration matters.** The 31 degrees C treatment is 4.4 degrees below an external short-test mean half-drop temperature of 35.4 degrees C. That external value is not a safe threshold for weeks of exposure or a measured threshold for our fragments.
5. **Trade-offs remain unresolved.** We found no clear growth cost of regeneration. Small samples, source-patch differences, and simultaneous measurement of growth and regrowth limit causal interpretation.

**Models used:** response-specific mixed models with tank effects and, for repeated measurements, fragment effects; daily binary morphology models (exploratory, frozen in `extras/`; not reported); interval-censored and Cox timing models; and multivariate summaries. Source patch is a fixed effect because there are only three patches. Growth is percent skeletal mass change, not areal calcification. Some binary morphology fits remain unstable despite penalization, and timing uncertainty is not fully tank-adjusted. Do not interpret their small p-values as established evidence. See `docs/analysis_diagnostic_inventory.md` for the current checks and gaps.

The separate 16-fragment microscope/photo cohort is cleaned in `code/01_clean_data.R`, tested in `code/04_regeneration.R` and plotted in `code/07_supplement.R`. It is not pooled with the main morphology models because it is wounded-only and includes only source-patch labels A/C. Endpoint (Fisher) and first-day (log-rank) tests for both cohorts are in `output/tables/11f_morphology_dataset_endpoint_tests.csv` and `11f_morphology_dataset_logrank_tests.csv`. The older cross-cohort figures (`11c`–`11e`) are frozen in `extras/figures/`.

## Repository map

| Path | Contents |
|---|---|
| `code/` | Ten files (merged 2026-11-11; guide in `code/README.md`): `00_setup.R`, `01_clean_data.R` through `08_checks.R`, and `run_all.R`. |
| `extras/` | Frozen exploratory analyses, older audits and diagnostics, RNA-seq starter code, Molly's original scripts, the scripts as they were before the merge, and the outputs they made. Not run by `run_all.R`. See `extras/README.md`. |
| `tools/` | `check_output_changes.py`: sorts real output changes from rerun noise. |
| `data/raw/` | Exported from Drive — never hand-edited. Decoded in `data/DATA_DICTIONARY.docx`. |
| `data/processed/` | Cleaned `.rds` the pipeline produces (regenerable). |
| `data/external/` | External heat-tolerance context: Cunning et al. 2024 CBASS ED50 reference plus preliminary Trinity/Ross Hauru colony table from Molly's 2026-09-03 forwarded email. |
| `output/paper/` | One CSV per results section: each number reported in the text, with its test, df, p and source table. Read by the statistics boxes in the team summary. |
| `output/tables/` | Every table the core scripts write, as CSV. |
| `output/` | Also `models/` (fitted `.rds`). Model-diagnostic reports + plots (`output/diagnostics/`, `figures/diagnostics/`) are regenerated by the pipeline and gitignored, not committed. |
| `figures/` | Figures the core scripts make (`.pdf` + `.png`); catalogued in `figures/FIGURE_INDEX.md`. |
| `literature/` | 101 PDFs + `LITERATURE.docx` (bibliography + synthesis). |
| `manuscript/Manuscript_LTH.docx` | Working draft — phenotype Methods + Results. |
| `docs/README.md` | Current documentation map and terminology. |
| `docs/next_analysis_roadmap.md` | Prioritized queue for the next analyses: final genetics, RNA-seq expression, phenotype integration, and supporting checks. |
| `docs/rnaseq/` | RNA-seq design, analysis proposal, preliminary SNP integration, genetic matching, candidate genes (suggestions, not a prescribed pipeline). |
| `docs/team_summary/` | Plain-language results summary (`.Rmd` + HTML) and pulled deck imagery. |
| `RESULTS.docx` | Full results narrative (all responses, source effects, thermal context, §10 limitations). |
| `notes/` | Sequencing plan; `notes/archive/` field notes and superseded plans. |

## Data & Drive

Raw data is exported **from** Google Drive (the project of record: raw Sheets, field notes, manuscript Doc); the repo never writes back. Drive folder `17. LTH_expression_by_temperature_2025` (ID `1sXfnHN-vmSBuwMEfERiYOWDeKRmjFWJP`), reached via the `gws` CLI. When Drive data changes, re-export into `data/raw/` and re-run the pipeline. Photographs live on the Stier Lab NAS (`smb://stier-nas1.eemb.ucsb.edu`); RNA-seq reads are processed at UC Davis. The NCBI BioProject is still TBD.

## RNA-seq (next)

144 selected margin libraries were shipped to UC Davis; expression counts will land in `data/raw/sequencing/` when available. The current RNA-seq analysis design uses D1, D3, and D10 margin samples. Any D0/D15 tissue should be described as collected tissue, not part of the 144-library design, unless the sequencing plan changes. The per-library phenotype covariate table is already built (`output/tables/31_rnaseq_phenotype_covariates.csv`, joined by `Fragment_ID`).

Rachael Bay sent a preliminary SNP cluster/PC file on 2026-09-02 (`data/raw/rnaseq/PRELIM_LTH_genoclusters.csv`). It joins cleanly to all 144 RNA-seq libraries and is integrated by `extras/code/rnaseq/32_prelim_snp_integration.R`, which writes `output/tables/32_prelim_snp_*.csv` and `figures/32_prelim_snp_*.pdf/png`. Treat these outputs as preliminary: the full SNP set should replace this layer for final genotype/kinship models.

Molly forwarded Trinity Conn's Hauru colony table on 2026-09-03. That file is now preserved at `data/external/trinity_conn_hauru_colonies_2026-08-10.xlsx` as preliminary external context; `extras/code/exploratory/26b_trinity_hauru_context.R` writes the cleaned copy, nearest-coordinate candidate table, and map. Do not assign Trinity/Ross ED50 values to source patches A/C/D yet. Use the file to prioritize candidate matches, then compare LTH RNA-seq SNPs to Trinity/Ross whole-genome genotypes.

Molly's follow-ups through 2026-09-10 are handled in `code/04_regeneration.R`. The summary now uses cumulative stage attainment and restricts timing averages to stages reached by everyone at both temperatures. Recovered temperature logs fill Days 0-15 for all eight tanks; see `docs/provenance/molly_followup_2026-09-11.md`.

The ordered analysis queue is in `docs/next_analysis_roadmap.md`. In short: finish the coauthor-summary cleanup, resolve final genetics, run expression genome-wide, then link expression modules or contrasts back to phenotype summaries without assigning destructive RNA-seq samples outcomes they did not directly receive.

## License and funding

**License:** Code MIT · Data CC-BY 4.0 (until publication, then DOI) · Images (separate Zenodo record, in preparation) CC BY-NC 4.0 · Manuscript all rights reserved.
**Funding:** NSF support to the project investigators. Field work at the UC Berkeley Gump Research Station, Mo'orea, French Polynesia.
