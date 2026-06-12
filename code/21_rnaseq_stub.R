# =============================================================================
# Purpose: OPTIONAL starting scaffold for importing the RNA-seq counts — NOT a
#          prescribed pipeline. The gene-expression analysis is Shreya's (Bay
#          lab) to design and lead; this just shows one way to read counts +
#          sample metadata into the repo conventions, if useful. Replace freely.
#
#          Not in code/_run_all.R. Optional suggestions live in
#          docs/for_shreya/ (analysis_proposal, gene_expression_integration_map).
# =============================================================================

source(here::here("code", "00_setup.R"))

# ---- Expected inputs (NOT YET PRESENT) -------------------------------------
# data/raw/sequencing/counts.csv         — gene × library count matrix
# data/raw/sequencing/sample_metadata.csv — library → (treatment, wound, day,
#                                                       tank, genet, plate, well)
# data/raw/plate_layout/Plate_1.csv,
# data/raw/plate_layout/Plate_2.csv      — already present

if (!file.exists(file.path(DATA_RAW, "sequencing", "counts.csv"))) {
  message("RNA-seq counts not yet present. This script is a stub — see NEXT_STEPS.md.")
  quit(save = "no", status = 0)
}

# ---- Import example only: read counts + sample metadata into the repo ------
# NOTE: this is deliberately just an *import* example. We do NOT fit a DE model
# here — the differential-expression design (factors, normalization, fixed/
# random structure, tool) is Shreya's (Bay lab) to specify. There is no single
# pre-decided model in this repo; the open design *questions* (not a formula)
# are written out in docs/for_shreya/analysis_proposal.md §3. Picking the model
# here would just hard-code one of those choices, so we don't.

counts   <- read_csv(file.path(DATA_RAW, "sequencing", "counts.csv"),
                     show_col_types = FALSE)
sample_md <- read_csv(file.path(DATA_RAW, "sequencing", "sample_metadata.csv"),
                      show_col_types = FALSE) |>
  janitor::clean_names()

# A per-library phenotype covariate table, already harmonized to the phenotype
# coding, is at output/tables/31_rnaseq_phenotype_covariates.csv (built by
# code/31). A RAW, un-recoded design lookup is alongside it
# (31_rnaseq_library_lookup_raw.csv) so factor levels / reference are yours to set.

# Basic sanity check that libraries line up between the two files.
stopifnot(all(sample_md$library_id %in% colnames(counts)))

cat(sprintf("Imported %d genes x %d libraries; %d sample-metadata rows.\n",
            nrow(counts), ncol(counts) - 1L, nrow(sample_md)))
cat("Next (yours to design): differential expression, module analysis, GO,\n")
cat("expression x phenotype, SNP calling — see docs/for_shreya/.\n")
