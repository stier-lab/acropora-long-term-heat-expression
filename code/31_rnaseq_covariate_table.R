# =============================================================================
# Purpose: HANDOFF CONVENIENCE for the RNA-seq analysis (Shreya / Bay lab) —
#          build one tidy table with a row per RNA-seq library, joined to the
#          phenotype covariates, so the gene-expression data can be modeled
#          against organismal traits without re-deriving the joins.
#
#          Coding is harmonized to the phenotype convention (data dictionary):
#          wound U/W -> no/yes; genotype A/C/D -> lowercase a/c/d; Temp_C -> 28C/31C.
#          Join key: Fragment_ID == id.
#
#          Covariate types:
#            - design        : treatment, wound, genet, tank, day, plate (per library)
#            - per-fragment   : symbiont density at that biopsy (where measured)
#            - per-genet      : resilience score / PCA displacement / rank
#                               (genet-level; identical for all libraries of a genet)
#          (Per-coral regeneration outcome is NOT attached: RNA-seq fragments were
#          destructively biopsied at D1/D3/D10, before the regeneration trajectory.)
#
# Input:   data/raw/plate_layout/Selected_Samples.csv
#          data/processed/symbiont_chl_clean.rds
#          output/tables/19_genet_resilience_summary.csv
# Output:  output/tables/31_rnaseq_phenotype_covariates.csv      (harmonized)
#          output/tables/31_rnaseq_library_lookup_raw.csv        (raw, un-recoded)
# =============================================================================

source(here::here("code", "00_setup.R"))

# ---- RNA-seq library rows, as they appear in the source plate layout ---------
sel_raw <- read_csv(file.path(DATA_RAW, "plate_layout", "Selected_Samples.csv"),
                    show_col_types = FALSE) |>
  filter(SampleType == "gene")

# ---- RAW, un-recoded lookup --------------------------------------------------
# library_id <-> fragment_id <-> the original design fields, verbatim from
# Selected_Samples.csv. We impose NO recoding here (no 28C/31C, no a/c/d, no
# no/yes) so Shreya can set her own factor levels and reference categories. The
# harmonized table below is the convenience join; this is the source of truth
# for "what was the raw design value for this library."
raw_lookup <- sel_raw |>
  transmute(library_id  = Sample_ID,
            fragment_id = Fragment_ID,
            Temp_C, Day, Wound, Tank, Genotype, Plate, WoundOrder)
write_csv(raw_lookup, file.path(TBL_DIR, "31_rnaseq_library_lookup_raw.csv"))

# ---- RNA-seq library list (144 'gene' libraries), harmonized to phenotype coding
libs <- sel_raw |>
  transmute(
    library_id = Sample_ID,
    id         = as.integer(Fragment_ID),
    treatment  = paste0(as.integer(Temp_C), "C"),     # 28C / 31C
    wound      = if_else(Wound == "W", "yes", "no"),  # U/W -> no/yes
    genet      = str_to_lower(Genotype),              # A/C/D -> a/c/d
    tank       = as.integer(Tank),
    day        = as.integer(Day),
    plate      = as.integer(Plate)
  )

# ---- Per-fragment symbiont density (destructive biopsy; one row per coral) ---
sym <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
  distinct(id, .keep_all = TRUE) |>
  transmute(id, symbiont_cells_per_cm2 = cells_per_cm2)

# ---- Per-genet resilience covariates (genet-level; from the dashboard) -------
res <- read_csv(file.path(TBL_DIR, "19_genet_resilience_summary.csv"),
                show_col_types = FALSE) |>
  transmute(genet                       = thicket,
            genet_mean_heat_sensitivity = round(mean_sensitivity, 3),
            genet_pca_displacement      = round(pca_displacement, 3),
            genet_resilience_rank       = rank_overall)

covariates <- libs |>
  left_join(sym, by = "id") |>
  left_join(res, by = "genet") |>
  arrange(plate, library_id)

write_csv(covariates, file.path(TBL_DIR, "31_rnaseq_phenotype_covariates.csv"))

cat("\n=== RNA-seq phenotype covariate table ===\n")
cat(sprintf("Libraries: %d  (design cells: %d temp x %d wound x %d genet x %d day)\n",
            nrow(covariates),
            n_distinct(covariates$treatment), n_distinct(covariates$wound),
            n_distinct(covariates$genet), n_distinct(covariates$day)))
cat(sprintf("With a symbiont-density value: %d / %d\n",
            sum(!is.na(covariates$symbiont_cells_per_cm2)), nrow(covariates)))
cat("Wrote output/tables/31_rnaseq_phenotype_covariates.csv (harmonized)\n")
cat("Wrote output/tables/31_rnaseq_library_lookup_raw.csv (raw, un-recoded)\n")
