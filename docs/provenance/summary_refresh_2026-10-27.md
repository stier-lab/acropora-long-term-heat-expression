# Results-summary refresh, 2026-10-27

## Confirmed correction

Molly's October 19 email in `Re: updated LTH summary HTML` (Gmail message `1a15556cc398d8aa`, thread `1a0639048ef4f103`) confirms fragment 119 is wounded and asks for side-by-side boxplots for buoyant-weight growth and symbiont density, plus individual panels.

Read-only verification of the live metadata sheet on October 27 confirmed `id=119`, `wound=yes`, `treatment=31`, `tank=10`, `biopsy_day=15`. Sheet: `1_Tdh7oQlCiH_g3TfWPqHoj76gEJ1Lr09hTBMuyVMuUI`, tab `metadata`.

Archived raw exports are preserved. `correct_wound_labels()` in `code/00_setup.R` applies the single confirmed correction in metadata, PAM, color, growth, symbiont, and paired tip/base imports. The main morphology sheet already identifies this fragment as wounded. Code 46 validates its assignment in all processed cohorts and verifies 12 fragments in each temperature-by-wound endpoint group. No response values are changed.

## Live score checks

Both scoring sheets were checked read-only on October 27. All 769 rows and all cell values in the main morphology sheet match the October 16 import when empty trailing cells are treated consistently; byte serialization differs only in empty trailing cells. The microscope-sheet export is byte-identical, SHA256 `3df738f75bdc5ccde36692b6acd7abcac390ae4fd8f82de383aa5aa2ed4c5f2c`. No new score corrections were required.

## Analysis and presentation

The full pipeline was rerun with the corrected label. Baseline tests are preserved in `data/metadata/119_before_correction_anova_2026-10-16.csv`; code 46 writes before/after heat-by-wound tests to `output/tables/46_wound_label_model_comparison.csv`. Raw measurements and temperature assignments are unchanged. Morphology attainment counts are unaffected.

Figure 2 retains repeated-measurement PAM lines but now shows symbiont density and growth as side-by-side wound-group boxplots. Each panel is also exported separately. These are distributions, not paired individual samples or tank-adjusted confidence intervals.

The correction does not reverse the heat-by-wound conclusion. In the canonical models, the heat-by-wound-by-time p-value changed from 0.235 to 0.207 for PAM and from 0.143 to 0.112 for log symbiont density; color remains 0.182. The growth heat-by-wound p-value changed from 0.908 to 0.997. Lower-order interactions in time-dependent models describe the reference time and must not be read as overall wound-effect tests. These checks do not establish equivalence between wound groups.

Summary changes distinguish the last observed regrowth stage from complete regeneration, early repair features from directly observed wound covering, and lack of clear wound effects from equivalence. Source-score values are read from analysis outputs during rendering. The summary includes adjusted trade-off uncertainty and remaining diagnostic limitations rather than treating unstable binary-model p-values as evidence.

Molly's Google Doc is unchanged. No email, Drive upload, or Git push is part of this refresh.

## Verification

All analysis stages in the original runner completed in 6.24 minutes. Editing the runner while Rscript was still reading it caused a final session-info parsing error after those stages finished. The saved runner was subsequently parsed successfully; the new code 46 checks, diagnostic inventory, session-info write, and summary render were rerun in a clean process and exited successfully. Logs: `output/logs/summary_refresh_2026-10-27.log` and `output/logs/summary_refresh_2026-10-27_render.log`.

All six processed datasets identify fragment 119 as wounded. All four endpoint temperature-by-wound groups contain 12 fragments. Before/after model tests are recorded in code 46's comparison table. The rebuilt HTML contains 24 embedded images and two videos; the embedded Figure 2 bytes match the newly generated PNG exactly. The combined figure was visually inspected for layout and legends. HTML checks were static, not browser-layout validation. The local Word manuscript audit still flags four stale phenotype numbers; that document and Molly's Google Doc were not rewritten.
