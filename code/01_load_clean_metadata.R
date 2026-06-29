# =============================================================================
# Purpose: Load + tidy the master metadata (one row per coral fragment) and
#          export a canonical lookup that all other scripts can join against.
#
# What & why: this is the experiment's master table. Every fragment that went
#   into the LTH heat × wound × genet experiment has one row here recording its
#   identity (id), which field colony / genet it came from (thicket = a/c/d), its
#   temperature treatment (28 °C ambient vs 31 °C heated), whether it was
#   wounded (apical tip clipped) or sham, which tank it was in, and its
#   biopsy schedule. Standardising these identifiers and factor levels ONCE here
#   means every downstream response script (PAM, colour, growth, symbionts, …)
#   can join back to a single canonical lookup instead of re-parsing the raw
#   field spreadsheet. Consistent factor codings here make the models in later
#   scripts comparable.
# Input:   data/raw/metadata/metadata.csv
# Output:  data/processed/coral_metadata.rds
#          output/tables/01_metadata_summary.csv  (sample-size tally)
# =============================================================================

source(here::here("code", "00_setup.R"))

# ---- Load ------------------------------------------------------------------
# Read the raw field spreadsheet and trim it to real records. guess_max = 5000
# makes readr inspect many rows before deciding each column's type (so a column
# that is blank for the first few hundred rows isn't mis-typed). clean_names()
# snake_cases the headers; the final filter drops the ~800 empty trailing rows
# Google Sheets exports by requiring both species and id present.
meta_raw <- suppressWarnings(
  read_csv(file.path(DATA_RAW, "metadata", "metadata.csv"),
           show_col_types = FALSE,
           guess_max = 5000)
) |>
  janitor::clean_names() |>
  filter(!is.na(species) & !is.na(id))

# ---- Tidy ------------------------------------------------------------------
# Coerce every column to a predictable type and encode the experimental design
# as factors. The CSV reads many fields as text, and consistent factor levels
# here guarantee the contrasts (set in 00_setup.R) behave the same way in every
# downstream model.
# Drop trailing whitespace in `thicket ` column; harmonize types
meta <- meta_raw |>
  # The header is "thicket " with a trailing space; matches() finds it
  # regardless so we get a `thicket` column. thicket = the genet / field colony
  # of origin (a, c, d) — modelled as a FIXED effect later (only 3 levels).
  rename(thicket = matches("^thicket")) |>
  mutate(
    species        = str_squish(species),            # collapse stray whitespace
    thicket        = str_to_lower(str_squish(thicket)),  # force a/c/d lowercase to match palettes
    id             = as.integer(id),                 # fragment ID — the join key
    sample         = str_squish(sample),
    wound          = factor(wound, levels = c("no", "yes")),  # no = sham, yes = apical tip clipped
    tank           = as.integer(tank),               # physical tank (1 of 8); maps to treatment
    treatment_c    = as.integer(treatment),          # raw numeric temp code (28 or 31)
    # Relabel the numeric temp as an ordered 2-level factor: 28C ambient is the
    # first level (the control), 31C heated second. Downstream figures/models
    # rely on this "28C, 31C" ordering.
    treatment      = factor(treatment_c, levels = c(28, 31),
                            labels = c("28C", "31C")),
    biopsy_day     = as.integer(biopsy_day),         # destructive sampling day (0/1/3/10/15)
    collection_date= as_date(collection_date),
    biopsy_date    = as_date(biopsy_date),
    coord_lat      = as.numeric(coord_lat),
    coord_long     = as.numeric(coord_long),
    # Rename the awkward auto-cleaned "_cm_2" columns to readable "_cm2" names.
    chlorophyll_ug_cm2 = as.numeric(chlorophyll_ug_cm_2),
    zooxanthellae_cells_cm2 = as.numeric(zooxanthellae_cells_cm_2),
    percent_growth_bw = as.numeric(percent_growth_bw),
    calculated_sa  = as.numeric(calculated_sa)
  ) |>
  # Drop the intermediate numeric temp code and the original "_cm_2" duplicates
  # now that we've kept tidy-named copies.
  select(-treatment_c, -ends_with("_cm_2"))

# Sanity: expected 192 destructive + 16 microscope + a small set of replicates
# (print the row count so the run log flags any load problem early).
n_total <- nrow(meta)
message(glue::glue("Loaded metadata: {n_total} rows"))

# ---- Summary tally ---------------------------------------------------------
# Cross-tabulate the full factorial design (temperature × wound × genet) to
# show the sample sizes per cell — a check that the design is balanced.
summary_tbl <- meta |>
  count(treatment, wound, thicket, name = "n") |>
  arrange(treatment, wound, thicket)

# ---- Save ------------------------------------------------------------------
# Write the canonical lookup as .rds (type-preserving — what other scripts read)
# and .csv (human-readable), plus the sample-size tally table.
dir.create(DATA_PROC, recursive = TRUE, showWarnings = FALSE)
dir.create(TBL_DIR,   recursive = TRUE, showWarnings = FALSE)
saveRDS(meta, file.path(DATA_PROC, "coral_metadata.rds"))
write_csv(summary_tbl, file.path(TBL_DIR, "01_metadata_summary.csv"))
write_csv(meta,        file.path(DATA_PROC, "coral_metadata.csv"))

message("Wrote coral_metadata.rds, coral_metadata.csv, 01_metadata_summary.csv")
print(summary_tbl)
