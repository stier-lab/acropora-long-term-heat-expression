# extras/: frozen analyses

Frozen 2026-11-11. Nothing here is run by `code/run_all.R`, and the paper doesn't depend on it. The one exception is `output/tables/19_genet_resilience_summary.csv`, which `code/06_genetics.R` reads as a fixed input (three source-level heat-sensitivity scores carried into the RNA-seq library table).

| Folder | What it is |
|---|---|
| `code/exploratory/` | Analyses run to explore, not planned for the paper: physiology PCA, source-patch dashboard, published heat-tolerance context (Cunning et al. 2024; Trinity/Ross Hauru colonies), growth-versus-regeneration trade-offs, alternative stage-timing figures (40, 41). |
| `code/checks/` | Earlier audits and model diagnostics: figure audit, time-series and model-comparison checks, variance partitioning, multiple testing, manuscript audit, analysis inventory, PAM probe-position audit (42), Oct 16 scoring-correction audit (44), fragment 119 label check (46), and diagnostic suites A–H. |
| `code/rnaseq/` | Optional RNA-seq import starter code (21), the September preliminary SNP integration (32, superseded by `code/06_genetics.R`), and the plate-layout figure. |
| `code/archive/` | An earlier Figure 1 draft and Molly's first PAM, growth, color and morphology scripts. |
| `code/analysis_before_merge/` | The analysis and paper scripts exactly as they were before they were merged into `code/01`–`08`. They also made the frozen outputs below that no current script writes: the daily yes/no regeneration models (`12_morph_*`, `12c_*`, `04_morphology_trait_*`), the master results table (`20_master_results*.csv`), the PAM location test (`02b_*`), the milestone-lag tables (`14_milestone_lag_*`), and similar. |
| `data/processed/` | `coral_physio_wide.rds`, a wide table used only by the frozen PCA and manuscript audit. |
| `output/tables/`, `figures/` | Every output those scripts wrote, as last rebuilt before the freeze. The team summary no longer uses any of them (removed 2026-11-11). |

## Running something here

These scripts expect the folder layout from before the merge, so most won't run in place. To rerun one:

```bash
git worktree add ../lth-pre-merge pre-merge-2026-11-11
cd ../lth-pre-merge && Rscript code/run_all.R
```

Then copy the outputs you need back into `extras/`. Don't report p-values from the daily yes/no regeneration models; they change with the computer and package versions (found in a 2026-11-10 rerun check).
