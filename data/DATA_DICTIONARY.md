# Data Dictionary — LTH Expression by Temperature (Project #17, 2025)

**Experiment:** Heat (28 °C vs 31 °C) × wounding (unwounded vs wounded) factorial on the
coral *Acropora pulchra*, Mo'orea (Gump Station), 2025. Three genets/thickets (a, c, d),
8 experimental tanks, ~208 individually tracked fragments. Each fragment carries a unique
integer `id` that links every data stream (PAM, color, morphology, growth, symbiont/chl,
wax, worm, RNA-seq).

**This file is the canonical guide to the data folders.** It exists so the collaborator
leading the gene-expression analysis (Shreya) can join RNA-seq libraries to phenotype
covariates by `id` without reverse-engineering the spreadsheets.

## Folder map

| Folder | Contents |
|---|---|
| `data/raw/<stream>/` | As-exported Google Drive sheets (`.xlsx` originals + flattened `data.csv` / `Sheet1.csv`). **Treat as read-only; never edit.** Some sheets store live Excel **formulas** (e.g. `=L2/100`); the ingest scripts recompute these from the underlying numeric columns. |
| `data/processed/*.rds` | Cleaned, typed, analysis-ready R objects produced by `code/0N_*.R`. One `.rds` per stream. These are what analyses load. |
| `data/external/` | Outside reference data (Cunning et al. 2024 *A. pulchra* CBASS ED50). |
| `data/metadata/*_codebook.csv` | Per-sheet column codebooks written by the original field team (source for much of this dictionary). |
| `data/raw/plate_layout/` | RNA-seq plate maps and the 144-library selected-sample design (see §RNA-seq). |

Each raw stream has a `code/0N_<name>.R` ingest script whose header `# Purpose/Input/Output`
block is authoritative. Reading the script header + this section together tells you exactly
what each cleaned `.rds` contains.

---

## Shared coding conventions (read this first)

These conventions are identical across **every** stream. `code/00_setup.R` defines paths and
palettes; each ingest script harmonizes types the same way.

| Variable | Raw form(s) | Cleaned form | Meaning |
|---|---|---|---|
| `id` | integer (1–301) | `integer` | **The master key.** Unique per fragment; links all streams (incl. RNA-seq via `Fragment_ID`). 208 unique ids. |
| `treatment` | `28` / `31` (numeric, °C) | `factor("28C","31C")` | Chronic temperature treatment. `28C` = ambient, `31C` = +3 °C heated. |
| `wound` | `no` / `yes` (also `U`/`W` in plate files; col named `wounded` in some raw sheets) | `factor("no","yes")` | Wounded fragments received a standardized lesion on D0. |
| `thicket` (= genet) | `a` / `c` / `d` (col is `"thicket "` with a trailing space) | `chr` lowercase | Source colony / genotype. Labels are **arbitrary, pre-genotyping** — not matched to Cunning's numbered genets. Counts: a=72, c=72, d=64. |
| `tank` | 3,4,5,6,9,10,11,12 | `integer` | Physical tank. **See tank→treatment map below.** |
| `day` / `biopsy_day` | integer | `integer` | Experiment day. **D0 = wounding day.** Negative days (e.g. −4, −1) are pre-wound baseline measurements. Destructive biopsies at D1, D3, D10, D15. |
| `sample` | `gene` / `gene/physio` (`physio/gene`) / `photo` | `chr` | Fragment's intended use: `gene`=transcriptomics only (144), `gene/physio`=transcriptomics + destructive physiology (48), `photo`=microscope imaging (16). |
| `species` | `a. pulchra` | `chr` | Constant; *Acropora pulchra* throughout. |

**Tank → treatment map (hard-coded in scripts; memorize this):**

| Treatment | Tanks |
|---|---|
| **28C** (ambient) | 3, 6, 9, 12 |
| **31C** (heated) | 4, 5, 10, 11 |

**Experiment timeline (from metadata):** collection/acclimation start ≈ 2025-05-18;
heat ramp 2025-05-28 → 31 °C hold reached 2025-05-31 (+1 °C/day); D0 wounding 2025-06-04;
biopsies span 2025-06-05 (D1) → 2025-06-19 (D15). (A second cohort with later dates appears
in the metadata — collection 2025-06-17, ramp ending 2025-06-29 — so do not assume a single
calendar offset; use `biopsy_day` / `day`, not raw dates, as the time axis.)

**Repeated vs terminal measures.** PAM, color, and morphology are **repeated** on the same
`id` across days (longitudinal). Symbiont/chl, wax SA, and buoyant weight are **terminal**
(one row per fragment at its biopsy day). Join longitudinal streams on `(id, day)`;
terminal streams on `id` alone.

---

## Datasets

### `coral_metadata` — master fragment lookup
- **Raw:** `data/raw/metadata/metadata.csv` (one row per fragment; ~800 trailing empty Google-export rows are dropped on load).
- **Script:** `code/01_load_clean_metadata.R` → **`data/processed/coral_metadata.rds`** (208 × 23; also written as `coral_metadata.csv`).
- **Role:** the canonical per-fragment table every other script joins against by `id`.
- **Key columns:** `id`, `species`, `thicket`, `sample`, `wound`, `tank`, `treatment` (factor 28C/31C; `treatment_c` keeps numeric), `biopsy_day`, `biopsy_date`, `collection_date`, `coord_lat`/`coord_long` (collection site per thicket), `acclimation_start_date`, `heat_ramp_start_date`, `heat_ramp_end_date`, `biopsy_notes`, `collected_in_liquid_nitrogen`, `in_dna_rna_shield`, `sub_samples_taken` (which sub-sections — tip/wound/middle/far — went into RNA/DNA shield), `percent_growth_bw`, `calculated_sa` (cm², wax-curve SA), `chlorophyll_ug_cm2`, `zooxanthellae_cells_cm2`.
- **DATA-QUALITY NOTE:** `chlorophyll_ug_cm2` **and** `zooxanthellae_cells_cm2` are **entirely empty (0/208)** in this table — chl-a/zoox were never merged back into the master sheet. Use `symbiont_chl_clean` for symbiont densities (see below). Chlorophyll-a appears to be **not yet measured / pending** (consistent with the project's "awaiting chl-a + RNA-seq" status).
- `biopsy_day` is `NA` for the 16 `photo` fragments.

### `pam_clean` — photosynthetic efficiency (Fv/Fm)
- **Raw:** `data/raw/pam/PAM_data.csv`. **Script:** `code/02_pam_analysis.R` → **`pam_clean.rds`** (336 × 8).
- Repeated measures of Diving-PAM dark-adapted yield. `fv_fm` was a spreadsheet formula in some rows (`=L2/100`) and is recomputed from `F` (F0), `M` (Fm), `Y` (yield), `E` (rel. ETR).
- Raw `location` = `top` (near wound) / `bottom` (~1 cm away); aggregated in cleaning.
- **Columns:** `date`, `day`, `treatment`, `tank`, `thicket`, `wound`, `id`, `fv_fm`.
- Join to metadata on `id`; longitudinal on `(id, day)`.

### `color_clean` — color-card pigmentation (bleaching proxy)
- **Raw:** `data/raw/color_card/data.csv`. **Script:** `code/03_color_card_analysis.R` → **`color_clean.rds`** (336 × 15).
- Siebeck D-scale scores (`D1`–`D6`). `color_num` = numeric D-scale; split scores (`D3/D4`) averaged to the midpoint (3.5). Lower = paler/more bleached.
- **Columns:** `species`, `date`, `day`, `treatment`, `tank`, `thicket`, `sample`, `id`, `wound`, `health_status`, `color` (text D-scale), `paling` (yes/no), `hole_at_center`, `notes`, `color_num`.

### `physio_clean` — morphological wound-healing traits
- **Raw:** `data/raw/physio_morphology/data.csv`. **Script:** `code/04_physio_morphology.R` → **`physio_clean.rds`** (768 × 21).
- 9 binary (yes/no → 1/0) morphology traits scored over time, mainly for wounded corals: `polyps_out`, `hole_in_center`, `polyp_in_hole`, `wound_smoothed`, `pigment_over_wound`, `tip_exist`, `tip_extension`, `new_corallites_on_tip`, `algae_on_wound`.
- Also `health_status`, `disease_status`, plus the shared keys (`id`, `day`, `treatment`, `tank`, `thicket`, `wound`, `sample`).

### `buoyant_weight_clean` — calcification / growth
- **Raw:** `data/raw/buoyant_weight/data.csv` (stores Excel formulas; dry-mass conversion recomputed via Davies 1989 / Jokiel 1978). **Script:** `code/05_buoyant_weight.R` → **`buoyant_weight_clean.rds`** (48 × 39).
- **48 growth fragments only**, all terminal at `biopsy_day == 15` (one clean 15-day window).
- **Primary metric:** `areal_calc` = areal calcification rate (mg CaCO₃ cm⁻² d⁻¹), normalized to day-15 wax SA (`sa_cm2`). Robustness metrics also present: `pct_growth` (% mass change), `sgr` (% d⁻¹), `g_per_day`, `delta_g`.
- Surface area (`sa_cm2`) comes from `wax_clean`.

### `wax_clean` — surface area calibration
- **Raw:** `data/raw/wax_dipping/data.csv` + `Standard_curve.csv` (14-point cylinder curve). **Script:** `code/07_wax_dipping.R` → **`wax_clean.rds`** (192 × 9).
- Per-fragment surface area from wax mass (`sa_curve_cm2`, primary) cross-checked against caliper geometry (`sa_caliper_cm2`).
- **Columns:** `id`, `treatment`, `biopsy_day`, `thicket`, `wound`, `dry_g`, `wax_g`, `sa_caliper_cm2`, `sa_curve_cm2`. The denominator for `buoyant_weight_clean$areal_calc` and `symbiont_chl_clean`.

### `symbiont_chl_clean` — symbiont density (+ chl-a slot)
- **Raw:** `data/raw/symbiont_counts/Raw_counts.csv` (4 hemocytometer quadrant counts Q1–Q4 per fragment) + `metadata_ordered_merge.csv` (SA, slurry volume). **Script:** `code/06_symbiont_chl.R` → **`symbiont_chl_clean.rds`** (192 × 9).
- **Columns:** `id`, `treatment`, `wound`, `biopsy_day`, `thicket`, `tank`, `sa_cm2`, `cells_per_cm2` (symbiont density, **populated 192/192**), `chlorophyll_ug_cm2`.
- **DATA-QUALITY NOTE:** `chlorophyll_ug_cm2` here is **also empty (0/192)** — chl-a is pending. The chl-a column is plumbed through the pipeline but contains no values yet.

### `worm_clean` — flatworm contamination check
- **Raw:** `data/raw/worm_presence/Sheet1.csv` (3 dates: 06/07, 06/08, 06/12). **Script:** `code/10_worms.R` → **`worm_clean.rds`** (576 × 9).
- *Acropora*-eating flatworm presence per fragment per date (`n_worms`, `present` 0/1). A QC/confound check on tank-level contamination, not a response variable.

### `ysi_clean` — daily water-chemistry spot checks
- **Raw:** `data/raw/ysi/Sheet1.csv`. **Script:** `code/09_ysi_water_chem.R` → **`ysi_clean.rds`** (72 × 14).
- Daily per-tank YSI: `temp_c` (converted from recorded °F), `do_pct`, `do_mgl`, `sal`, `ph`. `treatment` derived from the tank map.

### `apex_temperature` / `apex_temperature_daily` — continuous tank temperature
- **Raw:** `data/raw/apex/datalog*.xml` (Neptune Apex controller logs). **Script:** `code/08_apex_temperature.R` → **`apex_temperature.rds`** (hourly, 62k rows) and **`apex_temperature_daily.rds`** (daily).
- **Columns:** `datetime`/`date`, `probe`, `value_mean`, `value_sd`, `n`. Tank water temps are probes **`Temp1`–`Temp12`** (probe `TempN` = tank N); the file also logs many non-temperature probes (heaters, pumps, pH, ORP) — filter to `Temp1`–`Temp12` for tank temperature. This is the record of the heat ramp and 31 °C steady state.

### `coral_physio_wide` — endpoint physiology summary (derived)
- **Script output (no single raw file):** **`coral_physio_wide.rds`** (48 × 8) — one row per growth fragment with terminal/endpoint values: `pam_end`, `color_end`, `growth_areal`, `zoox_end`, keyed by `id` + design factors. Convenience table for multivariate/endpoint analyses.

### External — `cunning2024_apulchra_ed50.csv`
- Per-genet acute thermal-tolerance ED50 (°C) for *A. pulchra* from Mahana, Mo'orea (CBASS assay; Cunning et al. 2024 *Coral Reefs*). Used by `code/26_thermal_context.R` to place the chronic 28/31 °C treatments on a calibrated tolerance axis. **Caveat:** acute (18 h) ≠ chronic (weeks); Cunning's numbered genets are **not** matched to LTH thickets a/c/d. See `data/external/README.md`.

---

## RNA-seq sample metadata → phenotype join (for the gene-expression analysis)

The transcriptomics design and library→fragment mapping live in `data/raw/plate_layout/`.
The counts matrix and a final `sample_metadata.csv` are **not yet in the repo** (the ingest
scaffold is `code/21_rnaseq_stub.R`, which is intentionally a stub — the gene-expression
pipeline is yours to design).

**Files:**

| File | What it is |
|---|---|
| `Selected_Samples.csv` | **The 144-library sequencing design.** One row per library. |
| `Plate_1.csv`, `Plate_2.csv` | 96-well plate maps (`Well` A1–H12) with `Role` ∈ {sample, control, unused}. Plate_1 = 72 samples + 8 controls + 16 unused; together the two plates hold the 144 selected libraries. |
| `LTH_PlateLayout_with_IDs.xlsx` | Original Excel workbook the CSVs were exported from. |

**The 144-library design** is a fully crossed, **margin-only** subset:

> 3 genets (A,C,D) × 2 wound (U,W) × 2 temp (28 °C,31 °C) × 3 days (D1,D3,D10) × 4 tanks
> = 36 cells × 4 replicates (one per tank) = **144 libraries.**

Verified from `Selected_Samples.csv`: every (Temp × Wound × Day × Genotype) cell has exactly
4 libraries, one from each of the 4 tanks in that temperature. Note RNA-seq uses days
**D1, D3, D10** (not D15) and only `sample == "gene"` fragments.

**How to join libraries to phenotype covariates:**

1. Each library row carries a `Sample_ID` of the form
   **`{temp}_T{tank}_{genet}_D{day}_{wound}_id{ID}`** — e.g. `28C_T3_A_D1_W_id209`.
   `Selected_Samples.csv` additionally exposes `Fragment_ID` = the integer **`id`** directly.
2. **`Fragment_ID` (= `id`) is the join key to everything else.** Join your library/sample
   table to `coral_metadata.rds` on `id`, then optionally to `pam_clean`, `color_clean`,
   `physio_clean`, `symbiont_chl_clean`, etc. (terminal streams on `id`; longitudinal on
   `(id, day)`).
3. The design factors are also embedded in `Sample_ID` and as explicit columns
   (`Temp_C`, `Day`, `Wound` U/W, `Tank`, `Genotype` A/C/D, `WoundOrder`), so you can rebuild
   the model matrix directly from the plate file if `id`-level metadata is incomplete.
4. **Factor-coding heads-up:** plate files use `U`/`W` for wound and uppercase `A`/`C`/`D`
   for genet, whereas the cleaned phenotype `.rds` use `no`/`yes` and lowercase `a`/`c`/`d`.
   Harmonize on join (the stub `code/21_rnaseq_stub.R` sets `wound = factor(c("U","W"))`,
   `treatment = factor(c("28C","31C"))`, `day = factor(c(1,3,10))`).

**Suggested DE model factors:** `~ treatment * wound * day` with genet (`thicket`) and `tank`
as blocking terms — matching the phenotype models in `code/02`–`code/04`.
