# Data notes

Facts about the data that came from Molly and collaborators and that the code alone doesn't show. Updated 2026-11-12.

## Source sheets (Google Drive, read-only; exported to `data/raw/`)

| Data | Sheet ID |
|---|---|
| Main-experiment morphology scores | `1_6CCmYHGpocI7q9XDzaJHya5NaR-SqArgd24aOSpuoM` |
| Microscope-experiment scores | `1rB5VFJqQov0ZXsPjRKnjwVDOtiPxgbfMMDm9l_9Y-Q4` |
| Fragment metadata | `1_Tdh7oQlCiH_g3TfWPqHoj76gEJ1Lr09hTBMuyVMuUI` |
| YSI daily log | `11MrOICcoIp-lOJ9HiHcB9xWzf-_MgeMLQbt79Bmbb84` |

## Scoring and labels

- **October 16 corrections (Molly).** Field-sheet corrections to tip presence, tip extension and pigment. The earlier exports are kept as `data/raw/*/data_before_2026-10-16.csv`.
- **Fragment 119 is wounded** (Molly, October 19). `correct_wound_labels()` in `code/00_setup.R` applies this to every data set.
- **Budding can be scored before tip extension** (Molly's note on fragment 115). Stages are not forced into a strict order.
- **Pigment over the wound** was scored only on Days 3, 4, 9, 12 and 15 for wounded fragments. Other days are missing, not absent.
- **Fragments 115 and 121 are unresolved.** 121 has budding on Day 15 with no tip recorded that day; 115 budded on Day 14, before extension on Day 15. The color-card photos don't settle either. Still to ask Molly which observation supports 121's Day 15 budding. Without both fragments, budding is 11/11 ambient vs 3/11 heated, so the result doesn't depend on them.
- **Day 4:** every fragment was dipped to treat Acropora-eating flatworms, and all polyps were retracted that day.

## Water

- YSI temperature was recorded in °F (converted to °C), oxygen in mg/L and % saturation, salinity unitless. The pH scale (NBS or total) is not recorded; Molly has been asked.
- No YSI readings on June 4 or 7. Readings were taken between late morning and evening.
- Molly recovered the late-study Apex temperature logs on September 10, 2026; all eight tanks are covered through Day 15.

## Genetics and symbionts (Rachael Bay, October 12, 2026)

- Host clusters in `FINAL_LTH_genoclusters_syms.csv` are final. All 48 patch C samples are one cluster; A and D each span several.
- Symbiont typing is rough: reads mapped to four reference genomes. Group A is *Symbiodinium*; the other letters are Rachael's reference labels and haven't been confirmed as genera. Read shares are not cell fractions.
- Every sequencing label contains "wound", even for unwounded controls. Take wound status from the metadata, never from file names.
- The `total.ab` column is undocumented, so it isn't used.
- No Day 15 genotypes and no expression counts yet; Rachael and Shreya are rerunning the expression processing.

## External context (Trinity Conn, August 10, 2026, forwarded by Molly)

- Trinity recommends calling the Mo'orea animal *Acropora* cf. *pulchra* until the taxonomy is settled.
- Her preliminary Hauru colony table is in `data/external/`. Patch C sits 10.5 m from colony `Apul-115`, but Hauru clones span thickets, so identity needs a SNP comparison against her whole-genome data. Don't assign her ED50 values to our patches before that.

## Photographs (lab NAS, `smb://stier-nas1.eemb.ucsb.edu/Stier_Lab`)

- Color-card photos: `Moorea 2025/Moorea 2025/LTH 2025/A. pulchra/Color Card Photos/` (date folders 05302025 to 06182025).
- Microscope photos: `Moorea 2025/Moorea 2025/LTH 2025/A. pulchra/Ambient (28C)/` and `Heated (31C)/`, day folders D0–D15 (.tif and .vsi). These are the separate 16-fragment experiment, so they can't verify main-experiment scores.
