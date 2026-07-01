# Morphology scoring notes — *A. pulchra* wound healing (LTH expression-by-temperature 2025)

Where the binary wound-healing trait columns in `data.csv` /
`physio_characterization_log_-_complete.xlsx` come from, plus a known data-quality
issue. Compiled from the lab's records:

- **`READ.ME` tab** of `physio_characterization_log_-_complete.xlsx` (column definitions below)
- **Google Doc "A. Pulchra Healing General Observations by Day"** (the day-by-day
  healing timeline below) — Stier Lab Drive
- The live master sheet **"physio characterization log - complete"** (Google Sheets)

> The Drive files `MORPHOLOGY_RUBRIC.md` / `MORPHOLOGY_RELIABILITY.md`
> describe a **different** study (the Tiahura massive-*Porites* colony-shape survey:
> `bumpy` / `nobby` / `columnar` …) and do **not** apply to these traits.

## Trait definitions (from the spreadsheet `READ.ME` tab)

All traits are scored `yes`/`no` by eye from the daily color-card photo.
Blank/`NA` = not scored that visit.

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

The healing milestones happen **in sequence**:

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

## `hole_in_center` and `polyp_in_hole` are one trait → `axial_polyp_formation`

In both `data.csv` and the live master Google Sheet, `polyp_in_hole` is
**byte-for-byte identical** to `hole_in_center` — same `yes`/`no`/blank/`NA`
in every row (290 `yes`, 94 `no`, 384 `NA`; perfect diagonal cross-tab), switching
to `yes` on the same visit for all 24 wounded corals.

**This is not a data-entry error — it is how the trait was scored.** M.
Brzezinski (pers. comm., 2026) confirmed the central "hole" *is* the axial polyp
hole, which forms around the regenerating axial polyp, so the two go together (a
hole without a polyp was seen only once or twice): a single observable. On Molly's
recommendation, the two are **combined and renamed `axial_polyp_formation`** (the
axial corallite/calyx + polyp structure).

**In the pipeline:** `04_physio_morphology.R` keeps both original columns (for
provenance), checks they are still identical, and creates
`axial_polyp_formation = hole_in_center`. All downstream analysis
(`12_models.R`, `14_morphology_kaplan.R`, `sensitivity/28_multiple_testing.R`,
diagnostics) uses `axial_polyp_formation`, so the trait is modelled, plotted, and
entered into the BH multiple-testing family exactly once. No re-scoring is needed.

> Phase note: Table 1 of the manuscript places this trait in the
> **tissue-healing** phase. This fits the biphasic framework, which scores
> *regeneration* at the branch tip via `tip_exist` / `tip_extension` /
> `new_corallites_on_tip` — even though the axial polyp hole is where the axial
> corallite eventually regenerates.
