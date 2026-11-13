# extras/: kept exploratory analyses

*Last updated: 2026-11-13*

Nothing here is run by `code/run_all.R`, and the paper doesn't depend on it.
Only two exploratory analyses are kept, because they give outside context the
team may still want:

| Files | What it is |
|---|---|
| `code/exploratory/15_physiology_pca.R`, `figures/15*`, `output/tables/15_*`, `data/processed/coral_physio_wide.rds` | A principal-components summary of the physiology measures (an exploratory overview, not used in the paper). |
| `code/exploratory/26_thermal_context.R`, `26b_trinity_hauru_context.R`, `figures/26*`, `output/tables/26*`, `output/tables/19_genet_resilience_summary.csv` | How our 28 and 31 °C treatments compare with published heat-tolerance data (Cunning et al. 2024) and Trinity/Ross's Hauru colonies. Script 26 reads the retired source-patch resilience table, which is kept only for that reason. Inputs are in `data/external/`. |
| `code/rnaseq/21_rnaseq_import_stub.R`, `plate_layout_figure.R` | Starter code for importing the expression counts, and the sequencing plate-layout figure. |

These scripts expect the folder layout from before the 2026-11-11 merge, so
most won't run in place.

## Removed 2026-11-13

Everything else that used to be here was deleted as deprecated: the scripts
as they were before the merge, earlier audit scripts, Molly's first scripts,
the retired trade-off and resilience analyses, the preliminary SNP
integration, the old stage-timing display options, the daily yes/no
regeneration models (their p-values change between computers), the master
results table, and all their figures and tables. Nothing current used them.

They are all recoverable from these tags (on GitHub):

- `pre-prune-2026-11-13`: the repo just before the deletion.
- `pre-merge-2026-11-11`: the full pre-merge layout, which still runs.
- `pre-refactor-2026-11-11`: the flat layout before that.

```bash
git worktree add ../lth-pre-merge pre-merge-2026-11-11
cd ../lth-pre-merge && Rscript code/run_all.R
```
