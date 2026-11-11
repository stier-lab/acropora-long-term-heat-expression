# =============================================================================
# Purpose: Rebuild everything: data cleaning, the paper results, exploratory
#          analyses, checks, and the master results table.
#          Run from the project root:  Rscript code/run_all.R   (about 5 min)
#
# For writing the paper you only need code/run_paper.R, which is faster.
# Run this full version before sharing results or rebuilding the team summary.
#
# The order matters: later scripts read files written by earlier ones. Examples:
#   - analysis/07_surface_area.R runs before 05_growth.R (growth needs surface areas)
#   - rnaseq/31 runs after exploratory/19 (it joins the source-patch summary)
#   - rnaseq/32 (preliminary SNP clusters, superseded by analysis/43) still runs,
#     because exploratory/37 and analysis/36 read its joinability table
#   - checks/30_manuscript_audit.R runs last; it only warns, never stops the run
# =============================================================================

project_root <- here::here()
setwd(project_root)

scripts <- c(
  # Data cleaning and the analyses behind the paper
  "analysis/01_metadata.R",
  "analysis/02_photosynthesis.R",
  "analysis/03_color.R",
  "analysis/04_regeneration_scores.R",
  "analysis/07_surface_area.R",
  "analysis/05_growth.R",
  "analysis/06_symbionts.R",
  "analysis/08_temperature.R",
  "analysis/09_water_quality.R",
  "checks/10_flatworms.R",
  "analysis/11_microscope_experiment.R",
  "analysis/11c_stage_counts.R",
  "analysis/12_models.R",
  "analysis/13_source_patch_tests.R",
  "analysis/14_stage_timing_models.R",
  # Exploratory analyses and checks the master table draws on
  "exploratory/15_physiology_pca.R",
  "checks/17_figure_audit.R",
  "checks/18_data_validation.R",
  "exploratory/19_source_patch_dashboard.R",
  "checks/22_sensitivity_flagged.R",
  "checks/23_timeseries_diagnostics.R",
  "checks/24_headline_model_comparison.R",
  "exploratory/26_thermal_context.R",
  "exploratory/26b_trinity_hauru_context.R",
  "checks/27_variance_partitioning.R",
  "checks/28_multiple_testing.R",
  "exploratory/29_morphology_prob_contrasts.R",
  "checks/diagnostics/A_continuous_lmm.R",
  "checks/diagnostics/B_morphology_glmm.R",
  "checks/diagnostics/C_cox_diagnostics.R",
  "checks/diagnostics/D_pca_lrt.R",
  "checks/diagnostics/E_design_alignment.R",
  "checks/diagnostics/F_model_reproducibility.R",
  "checks/diagnostics/G_diagnostic_plots.R",
  "analysis/20_master_results_table.R",
  "checks/diagnostics/H_spreadsheet_coverage.R",
  # RNA-seq handoff and trade-off screens
  "rnaseq/31_rnaseq_covariate_table.R",
  "rnaseq/32_prelim_snp_integration.R",
  "exploratory/33_growth_tradeoff_screen.R",
  "exploratory/34_regeneration_tradeoff_screen.R",
  "exploratory/35_tradeoff_summary.R",
  "exploratory/37_tradeoff_formal_tests.R",
  # Stage timing, genetics, Figure 2 and remaining checks
  "analysis/36_stage_timing_tables.R",
  "checks/25_model_diagnostic_coverage.R",
  "analysis/39_stage_timing_figures.R",
  "exploratory/40_stage_mean_among_reachers.R",
  "exploratory/41_stage_progression_figure.R",
  "checks/42_pam_location_and_score_audit.R",
  "checks/44_corrected_scoring_audit.R",
  "analysis/43_genetic_clusters.R",
  "analysis/45_figure2_physiology.R",
  "checks/46_wound_label_validation.R",
  # Numbers and figures for each results section (output/paper/)
  "paper/tanks_water.R",
  "paper/physiology.R",
  "paper/regeneration.R",
  "paper/source_patches.R",
  "checks/38_analysis_inventory.R",
  "checks/30_manuscript_audit.R"
)

t0 <- Sys.time()
for (s in scripts) {
  cat("\n=== Running ", s, " ===\n", sep = "")
  tic <- Sys.time()
  source(file.path("code", s))
  cat("  done in", round(as.numeric(Sys.time() - tic, units = "secs"), 1), "s\n")
}
cat("\nFull pipeline complete in", round(as.numeric(Sys.time() - t0, units = "mins"), 2), "min.\n")

# Record the R and package versions used for this run.
dir.create("output", showWarnings = FALSE)
writeLines(trimws(capture.output(sessionInfo()), which = "right"), "output/session_info.txt")
