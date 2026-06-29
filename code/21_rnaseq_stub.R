# =============================================================================
# Purpose: OPTIONAL starting scaffold for importing the RNA-seq counts — NOT a
#          prescribed pipeline. The gene-expression analysis is Shreya's (Bay
#          lab) to design and lead; this shows one way to read counts +
#          sample metadata into the repo conventions, if useful. Replace freely.
#
#          Not in code/_run_all.R. Optional suggestions live in
#          docs/for_shreya/ (analysis_proposal, gene_expression_integration_map).
#
# What & why: this is a PLACEHOLDER. The RNA-seq count matrix does not exist in
#   the repo yet (it is Shreya's to generate), so there is nothing to analyze.
#   The script deliberately stops early — and does nothing — until the counts
#   arrive. We keep it in the repo so that (a) the README's "pending work" has a
#   concrete home, and (b) whoever picks up the expression analysis has a tested
#   example of how to read counts + sample metadata using this project's path
#   helpers (DATA_RAW, etc.) and column conventions. It is intentionally NOT a
#   differential-expression model: choosing the DE design here would pre-empt a
#   decision that belongs to the lead author (see the NOTE below).
#
# Input:   data/raw/sequencing/counts.csv          (NOT YET PRESENT — see below)
#          data/raw/sequencing/sample_metadata.csv (NOT YET PRESENT — see below)
# Output:  none (import sanity check only; prints a summary to the console)
# =============================================================================

# 00_setup.R loads packages and defines shared paths (DATA_RAW, ...). Even a stub
# sources it so the path helpers below resolve correctly.
source(here::here("code", "00_setup.R"))

# ---- Expected inputs (NOT YET PRESENT) -------------------------------------
# data/raw/sequencing/counts.csv         — gene × library count matrix
# data/raw/sequencing/sample_metadata.csv — library → (treatment, wound, day,
#                                                       tank, genet, plate, well)
# data/raw/plate_layout/Plate_1.csv,
# data/raw/plate_layout/Plate_2.csv      — already present

# ---- Early exit if the counts are not here yet -----------------------------
# Stub behavior: if the count matrix is absent (the normal state of
# the repo today), report it and quit cleanly with status 0 (a SUCCESS exit, so it
# never breaks _run_all.R or a CI run). Everything below this block only runs
# once Shreya has dropped real count data into data/raw/sequencing/.
if (!file.exists(file.path(DATA_RAW, "sequencing", "counts.csv"))) {
  message("RNA-seq counts not yet present. This script is a stub — see README.md (Pending work).")
  quit(save = "no", status = 0)
}

# ---- Import example only: read counts + sample metadata into the repo ------
# NOTE: this is deliberately just an *import* example. We do NOT fit a DE model
# here — the differential-expression design (factors, normalization, fixed/
# random structure, tool) is Shreya's (Bay lab) to specify. There is no single
# pre-decided model in this repo; the open design *questions* (not a formula)
# are written out in docs/for_shreya/analysis_proposal.md §3. Picking the model
# here would just hard-code one of those choices, so we don't.

# Read the gene x library count matrix (genes in rows, one column per library).
counts   <- read_csv(file.path(DATA_RAW, "sequencing", "counts.csv"),
                     show_col_types = FALSE)
# Read the per-library metadata; clean_names() standardises headers to snake_case
# (e.g. "Library ID" -> library_id) so the join/check below can rely on the name.
sample_md <- read_csv(file.path(DATA_RAW, "sequencing", "sample_metadata.csv"),
                      show_col_types = FALSE) |>
  janitor::clean_names()

# A per-library phenotype covariate table, already harmonized to the phenotype
# coding, is at output/tables/31_rnaseq_phenotype_covariates.csv (built by
# code/31). A RAW, un-recoded design lookup is alongside it
# (31_rnaseq_library_lookup_raw.csv) so factor levels / reference categories can be set by the analyst.

# Basic sanity check that libraries line up between the two files: every library
# named in the metadata must appear as a column in the count matrix, or the two
# files are out of sync and any downstream join would silently mis-pair samples.
# stopifnot() halts with an error if the condition is FALSE.
stopifnot(all(sample_md$library_id %in% colnames(counts)))

# Report what we read. ncol(counts) - 1L subtracts the leading gene-ID column so
# the count is the number of actual sample libraries.
cat(sprintf("Imported %d genes x %d libraries; %d sample-metadata rows.\n",
            nrow(counts), ncol(counts) - 1L, nrow(sample_md)))
cat("Next (yours to design): differential expression, module analysis, GO,\n")
cat("expression x phenotype, SNP calling — see docs/for_shreya/.\n")
