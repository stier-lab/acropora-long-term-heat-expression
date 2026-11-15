# Raw Data Provenance

Last audited: 2026-09-03 · Updated 2026-11-15

Canonical Drive project folder:
`17. LTH_expression_by_temperature_2025`
`https://drive.google.com/drive/folders/1sXfnHN-vmSBuwMEfERiYOWDeKRmjFWJP`

This folder contains local raw exports from the Google Drive project data
folder plus provenance-only methods notes needed to interpret them. The detailed
source-by-source map is in `data/metadata/source_inventory.csv`. Each raw data
folder has a matching codebook in `data/metadata/*_codebook.csv`.

## Analysis Inputs

- `metadata`: master experimental metadata. This is the design spine used by
  the analysis.
- `plate_layout`: planned RNA-seq plate layout and selected sample manifest.
  Used for RNA-seq covariate/handoff tables.
- `rnaseq`: preliminary Bay lab SNP clustering for the 144 RNA-seq libraries.
  Used only by the exploratory preliminary-SNP integration script.
- `pam`: repeated PAM Fv/Fm measurements for the physiology subset.
- `color_card`: repeated color-card/health observations for the physiology
  subset.
- `physio_morphology`: gross daily morphology/healing characterization for the
  physiology subset.
- `microscope_physio`: separate photo-only microscope cohort, analyzed in its
  own branch of the pipeline.
- `buoyant_weight`: growth/calcification measurements and calculation sheets.
- `wax_dipping`: surface-area measurements and standard curve used to normalize
  destructive physiology.
- `symbiont_counts`: hemocytometer counts and metadata used to calculate
  symbiont density per surface area.
- `ysi`: daily handheld water-quality spot checks.
- `apex`: continuous Apex XML temperature logs.
- `worm_presence`: AEFW/worm surveillance checks.

## Provenance-Only Imports

- `shipping`: sample handling, box, and preservation metadata. Imported and
  codebooked for RNA-seq handoff provenance; the current pipeline does not read
  these sheets directly.

## Archived 2026-11-15

These were in the repo but nothing used them. They are still on Drive and in git
tag `pre-cleanup-2026-11-15`:

- `daily_health_log/`: a Drive sheet whose observation columns (date, health,
  polyps, algae, tip, notes) were never filled in.
- `wax_dipping_standard_curve/`: a second copy of the standard curve. The code
  reads `wax_dipping/Standard_curve.csv`, which has the same numbers plus Molly's
  note on which fit to use.
- `printable_data_sheets/`: blank printable YSI, PAM, color-card and daily-health
  field sheets.
- `chlorophyll/`: an example chlorophyll calculation workbook. The chl-a assay was
  planned but never run for this project.
- `buoyant_weight_archive/`: an older working copy of the buoyant-weight
  calculations. The growth calculation now uses `buoyant_weight/Buoyant_Weight_-_complete.xlsx`.
- `buoyant_weight/calcuations.csv` and `water_temp_table.csv`: CSV copies of two
  tabs of that workbook. The code reads the workbook itself.
- `symbiont_counts/Copy_of_calculations.csv` and `Ordered_averages.csv`:
  spreadsheet intermediates. The code recomputes density from `Raw_counts.csv`.

## Important Conventions

- Raw temperature treatments appear as `28` and `31` in many sheets; analysis
  code standardizes them to `28C` and `31C`.
- Tank-to-temperature plumbing is fixed: tanks 3, 6, 9, and 12 are ambient;
  tanks 4, 5, 10, and 11 are heated.
- Source patch labels are pre-genotype-analysis field labels: `a`, `c`, and
  `d` for the main experiment; the microscope photo cohort includes `a` and
  `c`. Some code and exported tables still call these `thicket` or legacy
  `genet` values. Rachael Bay's preliminary 2026-09-02 SNP clusters suggest
  source patch C is coherent, while A and D include multiple clusters and share
  some clusters. Treat those clusters as preliminary until the full SNP set is
  delivered.
- Some Google Sheet exports preserve printable page headers or blank trailing
  rows. Analysis scripts filter those rows explicitly when a sheet is used.
- `physio_morphology/data.csv` and its companion `.xlsx` workbook were refreshed
  from the current Drive sheet on 2026-09-01. The prior local CSV export had
  extra blank grid rows; after removing blank IDs and ignoring date-format
  class, non-date values matched the current Drive export.
- `buoyant_weight/Buoyant_Weight_-_complete.xlsx` was refreshed from the
  authenticated Drive export on 2026-09-01; the analysis still reads the
  tab-level CSV exports in the same folder.
- Microscope photo scoring is not pooled with gross morphology. It has a
  different cohort, tanks, source patches, and photo-only purpose.
- The microscope `pigment_over_wound` column is scored only through day 7 in
  the current sheet; blanks after day 7 mean not scored, not no pigmentation.
- Raw color-card and microscope image files live on the lab NAS, not in Drive.
  The NAS location is in `docs/DATA_NOTES.md` (Photographs); the older
  path-and-count log was removed on 2026-11-12 and is in git history
  (`git show 281f4f9:docs/provenance/nas_image_sources.md`).
