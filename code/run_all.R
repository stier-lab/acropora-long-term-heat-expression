# =============================================================================
# run_all.R: rebuild every number, table and figure the paper uses
#
#   Rscript code/run_all.R
#
# Runs the eight numbered scripts in order. Each runs in its own workspace, so
# nothing leaks from one script into the next; each starts by loading
# 00_setup.R (packages, paths, colors, random seed). Run R in a UTF-8 locale
# (the Mac default) so symbols such as ° and × are written correctly.
#
# Afterwards, `python3 tools/check_output_changes.py` separates real changes
# from files that change on every run (PDF timestamps).
#
# Scripts in extras/ are frozen and are not run here.
# =============================================================================

project_root <- here::here()
setwd(project_root)

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

t0 <- Sys.time()
for (s in scripts) {
  cat("\n=== Running ", s, " ===\n", sep = "")
  tic <- Sys.time()
  source(file.path("code", s), local = new.env(parent = globalenv()))
  cat("  done in", round(as.numeric(Sys.time() - tic, units = "secs"), 1), "s\n")
}
cat("\nPipeline complete in", round(as.numeric(Sys.time() - t0, units = "mins"), 2), "min.\n")
cat("Numbers for the results text are in output/paper/.\n")

writeLines(trimws(capture.output(sessionInfo()), which = "right"), "output/session_info.txt")
