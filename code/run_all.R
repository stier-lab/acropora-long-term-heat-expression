# =============================================================================
# run_all.R: rebuild every number, table and figure the paper uses
#
# Runs the eight numbered scripts in order. Each runs in its own empty
# workspace, so no leftover objects carry from one script into the next;
# each starts by loading 00_setup.R (packages, paths, colours, random seed).
# Run R with UTF-8 text encoding (the Mac default) so symbols such as ° and
# × are written correctly. Before the first run on a new machine, run
# `renv::restore()` once: it installs the exact package versions recorded in
# renv.lock.
#
# Afterwards, `python3 tools/check_output_changes.py` separates real changes
# from files that change on every run (PDF timestamps). Scripts in extras/ are
# kept for the record, no longer updated, and are not run here.
#
# Sections
#   1. Scripts to run
#   2. Run them in order
#   3. Session log
#
# Reads:  everything the numbered scripts read (data/raw, data/processed, ...)
# Writes: everything the numbered scripts write, plus output/session_info.txt
# Run:    Rscript code/run_all.R
# =============================================================================


# ---- 1. Scripts to run -----------------------------------------------------
# Work from the project folder so the file paths below, written from there,
# point to the right place.
project_root <- here::here()
setwd(project_root)

# Order matters: 01 writes the cleaned data that every later script reads.
scripts <- c(
  "01_clean_data.R",      # raw spreadsheets -> data/processed/*.rds
  "02_tanks_water.R",     # tank temperature and water quality
  "03_physiology.R",      # Fv/Fm, color, symbionts, growth; Figure 2
  "04_regeneration.R",    # healing and regeneration stages; Figure 1B
  "05_source_patches.R",  # does the heat effect differ among source patches?
  "06_genetics.R",        # final genetic clusters; RNA-seq library table
  "07_supplement.R",      # supplement-only results
  "08_checks.R"           # data validation and model checks
)


# ---- 2. Run them in order --------------------------------------------------
# new.env(parent = globalenv()) gives each script an empty workspace; the
# timer prints how long each one took.
t0 <- Sys.time()
for (s in scripts) {
  cat("\n=== Running ", s, " ===\n", sep = "")
  tic <- Sys.time()
  source(file.path("code", s), local = new.env(parent = globalenv()))
  cat("  done in", round(as.numeric(Sys.time() - tic, units = "secs"), 1), "s\n")
}
cat("\nPipeline complete in", round(as.numeric(Sys.time() - t0, units = "mins"), 2), "min.\n")
cat("Numbers for the results text are in output/paper/.\n")


# ---- 3. Session log --------------------------------------------------------
# Record the R and package versions that made these outputs; trimws() strips
# trailing spaces so git does not show changes where nothing real changed.
writeLines(trimws(capture.output(sessionInfo()), which = "right"), "output/session_info.txt")
