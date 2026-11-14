# LTH: heat × wounding in *Acropora* cf. *pulchra*

A 15-day experiment at Gump Station, Mo'orea (May–July 2025). Fragments were held at 28 °C or 31 °C in four tanks per temperature; half had about 1 cm clipped from the growing tip on Day 0. We tracked wound healing and regrowth, photosynthesis efficiency (Fv/Fm), color, symbiont density and skeletal growth, and collected tissue for RNA-seq.

<p align="center">
  <img src="docs/team_summary/img/healing_ambient_vs_heated.gif" width="720"
       alt="Wound healing, 28 °C vs 31 °C, Day 0 to Day 15"><br>
  <sub>Both temperatures re-cover the wound and rebuild a tip; by Day 15 the 31 °C fragment has not budded new corallites.</sub>
</p>

## Main results

Heat at 31 °C did not stop wounds from healing, but it slowed the last step of regrowth and stressed the corals.

<p align="center">
  <img src="figures/39b_stage_timing_cumulative_option.png" width="600"
       alt="Four step plots of the percent of wounded fragments that had reached each regeneration milestone by each day, 28 °C versus 31 °C. The two temperatures track closely for the first three milestones. For radial corallites, 28 °C reaches 100% by Day 15 while 31 °C reaches about 33%."><br>
  <sub><b>Figure 1B.</b> Share of wounded fragments (12 per temperature) that had reached each regeneration milestone. Blue solid = 28 °C, orange dashed = 31 °C. Both temperatures re-formed a tip polyp and corallite; by Day 15 all 12 ambient fragments had budded new side corallites, against 4 of 12 heated ones.</sub>
</p>

<p align="center">
  <img src="figures/45_physiology_figure2.png" width="600"
       alt="Three rows of panels comparing 28 °C and 31 °C. A: photosynthetic efficiency stays near 0.67 at 28 °C but falls to about 0.56 by Day 14 at 31 °C. B: symbiont density holds near 1 million cells per square centimetre at 28 °C but drops to about 0.2 by Day 15 at 31 °C. C: skeletal mass gain is about 16% at 28 °C and 10% at 31 °C, similar for wounded and unwounded fragments."><br>
  <sub><b>Figure 2.</b> Photosynthetic efficiency (A), symbiont density (B) and 15-day skeletal growth (C). Blue = 28 °C, orange = 31 °C; open/white = unwounded, filled/gray = wounded. Heat lowered photosynthetic efficiency (0.67 vs 0.56 on Day 14), symbiont density (0.97 vs 0.21 million cells cm⁻² on Day 15) and growth (16.0% vs 10.2%, p = 0.044). Wounding made little difference to any of the three. Full legend: <code>manuscript/FIGURE_2_LEGEND.md</code>.</sub>
</p>

The full plain-language write-up is the [results summary](https://stier-lab.github.io/acropora-long-term-heat-expression/team_summary/LTH_results_summary.html). Every number comes from `output/paper/*_numbers.csv`.

## Start here

| Need | File |
|---|---|
| Results in plain language | [Results summary (opens in browser)](https://stier-lab.github.io/acropora-long-term-heat-expression/team_summary/LTH_results_summary.html) · source: `docs/team_summary/` |
| What a term means | [`docs/GLOSSARY.md`](docs/GLOSSARY.md) |
| How to rerun or change the analysis | [`docs/HOW_TO_CHANGE.md`](docs/HOW_TO_CHANGE.md) |
| Every number reported in the text | `output/paper/*_numbers.csv` |
| Which figure is which in Molly's draft | `manuscript/MAIN_FIGURE_PLAN.md` |
| How the code is organized | `code/README.md` |
| Which outputs are current | `docs/OUTPUTS_STATUS.md` |
| Model checks | `docs/diagnostics/MODEL_DIAGNOSTICS_2026-11-11.md` |
| Data caveats and collaborator notes | `docs/DATA_NOTES.md` |

## Run it

```bash
Rscript -e 'renv::restore()'   # once: installs the pinned package versions
Rscript code/run_all.R         # about a minute; rebuilds every number, table and figure
```

Package versions are pinned in `renv.lock` (R 4.6.2). Restore them before rerunning: the simulation-based residual checks (DHARMa) give different outlier sets under other versions. Never hand-edit numbers; change the code and rerun. Older and exploratory analyses were removed on 2026-11-13 and 2026-11-14; recover them from the tags `pre-prune-2026-11-13` and `pre-cleanup-2026-11-14`.

## Design

| Factor | Levels |
|---|---|
| Temperature | 28 °C and 31 °C, ramped 1 °C per day. Tanks 3, 6, 9, 12 (28 °C) and 4, 5, 10, 11 (31 °C). |
| Wound | About 1 cm clipped from the tip on Day 0, or not clipped |
| Source patch | A, C and D (Mahana/Tiahura, NW Mo'orea). Stored as `thicket`, or `genet` in older names. These are field labels, not genotypes. |
| Sampling | 192 fragments in the main experiment: 144 sampled on Days 1, 3 and 10 (the RNA-seq libraries) and 48 followed daily, then sampled on Day 15. A separate 16 wounded fragments were photographed under a microscope. |

## Folders

| Path | Contents |
|---|---|
| `code/` | The ten core scripts and `run_all.R` |
| `data/raw/` | Exports from Google Drive; never hand-edited. Codebooks are in `data/metadata/`. |
| `data/external/` | Collaborator and published reference tables (Cunning et al. 2024; Trinity Conn's Hauru colonies) |
| `output/`, `figures/` | Everything the code writes; figures are listed in `figures/FIGURE_INDEX.md` |
| `manuscript/` | Figure map and Figure 2 legend |
| `docs/` | Team summary, diagnostics, data notes, RNA-seq plan, Drive method exports |
| `literature/` | `library_map.csv` (citation, theme and file name for each reference). The PDFs live on Google Drive, not in this public repo (see `literature/README.md`). |
| `tools/` | `check_output_changes.py`, which sorts real output changes from rerun noise |

Google Drive (`17. LTH_expression_by_temperature_2025`, ID `1sXfnHN-vmSBuwMEfERiYOWDeKRmjFWJP`) holds the raw Sheets, field notes and Molly's manuscript draft. Photographs are on the lab NAS (see `docs/DATA_NOTES.md`). RNA-seq processing is at UC Davis (Bay lab).

## Citing, license and funding

To cite this code and data, use GitHub's "Cite this repository" button (from `CITATION.cff`).

Code MIT · data CC-BY 4.0 · images CC BY-NC 4.0 · manuscript all rights reserved. NSF support to the project investigators; fieldwork at the UC Berkeley Gump Research Station.
