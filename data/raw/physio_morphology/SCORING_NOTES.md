# Morphology scoring notes — *A. pulchra* wound healing (LTH expression-by-temperature 2025)

Provenance for the binary wound-healing trait columns in `data.csv` /
`physio_characterization_log_-_complete.xlsx`, plus a documented data-quality
issue. Compiled from the lab's own records:

- **`READ.ME` tab** of `physio_characterization_log_-_complete.xlsx` (column definitions below)
- **Google Doc "A. Pulchra Healing General Observations by Day"** (the day-by-day
  healing timeline below) — Stier Lab Drive
- The live master sheet **"physio characterization log - complete"** (Google Sheets)

> Note: the Drive files `MORPHOLOGY_RUBRIC.md` / `MORPHOLOGY_RELIABILITY.md`
> describe a **different** study (the Tiahura massive-*Porites* colony-shape blitz:
> `bumpy` / `nobby` / `columnar` …) and do **not** apply to these traits.

## Trait definitions (from the spreadsheet `READ.ME` tab)

All traits are scored `yes`/`no` by visual assessment from the daily color-card
photo; blank/`NA` = not scored that visit.

| Column | Definition |
|---|---|
| `polyps_out` | The majority of polyps were out and relaxed. |
| `hole_in_center` | A hole was present in the center of the wound. |
| `polyp_in_hole` | A polyp was present in the hole. |
| `wound_smoothed` | The wounded tissue showed signs of skeletal remodeling / smoothing. |
| `pigment_over_wound` | The wound area appeared pigmented. |
| `tip_exist` | The axial corallite was present. |
| `tip_extension` | The axial corallite extended vertically. |
| `new_corallites_on_tip` | Newly budding corallites were present along the axial corallite. |
| `algae_on_wound` | Algae was present on the wound area. |

Design columns: `day` (post-wound day, D0 = wounding), `treatment` (28/31 °C),
`tank`, `thicket` (genet, pre-genotype alphabetical label), `id` (coral),
`wounded` (yes = wounded, no = control).

## Healing timeline (unheated corals; from the field-notes doc)

The healing milestones are **sequential**, not simultaneous:

| Day | Observation |
|---|---|
| 0–1 | Blunt wound; coarse skeleton, no tissue. |
| 2 | Skeletal smoothing begins; **axial polyp hole begins to be defined**. |
| 3 | Smoothing continues; axial polyp hole becomes more defined. |
| 4 | Obvious axial polyp hole present **with axial polyp present**; skeleton rounding out. |
| 6 | Definitive axial hole with polyp and smoothed skeleton; ~50% chance of tip/axial corallite beginning. |
| 8 | Greater chance of axial corallite existing, possible extension. |
| 10 | Corallite extension continues. |
| 12 | New corallites begin forming on the axial corallite "like branches on a tree". |
| 13 | New corallites continue; polyps may be present in these corallites. |

So `hole_in_center` is expected to appear ~D2–3 and `polyp_in_hole` ~D4 — the
polyp should **lag** the hole by roughly 1–2 days. The "axial polyp hole" is
where the axial corallite (the branch's growth axis) regenerates, so
`polyp_in_hole` is an **early regeneration** milestone, not just wound closure.

## ⚠️ KNOWN DATA-QUALITY ISSUE: `polyp_in_hole` duplicates `hole_in_center`

In both `data.csv` and the live master Google Sheet, the `polyp_in_hole` column
is **byte-for-byte identical** to `hole_in_center` — same `yes`/`no`/blank/`NA`
value in every row (290 `yes`, 94 `no`, 384 `NA`; perfect diagonal cross-tab).
Evidence this is a data-entry duplication, not real data:

- Not a CSV-export artifact (identical in the source `.xlsx`) and not a formula
  (the `.xlsx` cells store literal values, 0 formula cells in those columns).
- Specific to this pair — `polyp_in_hole` differs from its other neighbor
  `wound_smoothed` (28 cells), so it is not a whole-sheet column shift.
- Biologically impossible given the timeline above: for **all 24** wounded
  corals the two traits switch to `yes` on the **same visit** (zero lag), and
  there is **no** `hole=yes, polyp=no` row — i.e. no "hole present but empty"
  interval, which the D2→D4 sequence requires.

**Consequence in the pipeline:** the genuine `polyp_in_hole` observations were
never recorded digitally (not even in the master sheet), so they cannot be
recovered from the data. Until the column is corrected, the analysis code
(`04_physio_morphology.R` and downstream `12_models.R`, `14_morphology_kaplan.R`,
`sensitivity/28_multiple_testing.R`) **detects the duplication at runtime and
excludes `polyp_in_hole`** so it is not modelled, plotted, or counted twice in
the multiple-testing family. The guard removes itself automatically once the two
columns are no longer identical.

**To restore the trait:** re-score `polyp_in_hole` (yes/no, "was a polyp present
in the central hole?") for each wounded coral × day from the daily color-card
photos, and write the corrected values into `data.csv` /
`physio_characterization_log_-_complete.xlsx`. Worth doing because it is a
regeneration-phase milestone (relevant to the heat-impairs-regeneration result).

A ready-to-fill template is provided: **`polyp_in_hole_RESCORING_TEMPLATE.csv`**
(384 rows = 24 wounded corals × D0–D15). Each row carries the photo `date`,
`id`, `treatment`/`tank`/`thicket`, and `hole_in_center_REFERENCE`; open that
day's color-card photo and enter `yes`/`no` in `polyp_in_hole_RESCORED`. Note:
`polyp_in_hole` should be a subset of `hole_in_center` (no polyp without a hole),
so `polyp = yes` is only plausible where the reference column is `yes`.
