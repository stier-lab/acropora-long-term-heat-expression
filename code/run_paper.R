# =============================================================================
# Purpose: Rebuild the numbers and figures used in the paper's results.
#          Run from the project root:  Rscript code/run_paper.R
#
# What you get (output/paper/):
#   tanks_water_numbers.csv     tank temperature, salinity, oxygen, pH
#   physiology_numbers.csv      Fv/Fm, color, symbionts, growth (Figure 2)
#   regeneration_numbers.csv    regeneration stages and timing (Figure 1)
#   source_patches_numbers.csv  differences among source patches A, C and D
# Each row is one number used in the text, with its test and the table it came
# from. code/README.md explains every script; run_all.R rebuilds everything.
# =============================================================================

project_root <- here::here()
setwd(project_root)

scripts <- c(
  "analysis/01_metadata.R",
  "analysis/02_photosynthesis.R",
  "analysis/03_color.R",
  "analysis/04_regeneration_scores.R",
  "analysis/07_surface_area.R",
  "analysis/05_growth.R",
  "analysis/06_symbionts.R",
  "analysis/08_temperature.R",
  "analysis/09_water_quality.R",
  "analysis/11_microscope_experiment.R",
  "analysis/11c_stage_counts.R",
  "analysis/12_models.R",
  "analysis/13_source_patch_tests.R",
  "analysis/14_stage_timing_models.R",
  "analysis/36_stage_timing_tables.R",
  "analysis/39_stage_timing_figures.R",
  "analysis/45_figure2_physiology.R",
  "paper/tanks_water.R",
  "paper/physiology.R",
  "paper/regeneration.R",
  "paper/source_patches.R"
)

t0 <- Sys.time()
for (s in scripts) {
  cat("\n=== Running ", s, " ===\n", sep = "")
  source(file.path("code", s))
}
cat("\nPaper outputs rebuilt in", round(as.numeric(Sys.time() - t0, units = "mins"), 2),
    "min. See output/paper/.\n")
