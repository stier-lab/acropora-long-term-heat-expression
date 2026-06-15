# =============================================================================
# Purpose: Run the full LTH analysis pipeline in dependency order.
#          Source this once from the project root: `Rscript code/_run_all.R`
#          Total runtime ~3 minutes (Apex XML parse dominates first pass).
# =============================================================================

# Source from project root (resolve via here::here, which uses .Rproj/.here)
project_root <- here::here()
setwd(project_root)

scripts <- c(
  "01_load_clean_metadata.R",
  "02_pam_analysis.R",
  "03_color_card_analysis.R",
  "04_physio_morphology.R",
  "07_wax_dipping.R",          # must run before 05 — 05 reads wax_clean.rds (surface areas)
  "05_buoyant_weight.R",
  "06_symbiont_chl.R",
  "08_apex_temperature.R",
  "09_ysi_water_chem.R",
  "10_worms.R",
  "11_combined_figure.R",
  "12_extended_stats.R",
  "12b_color_clmm_robustness.R",
  "12c_morph_blme.R",
  "13_genet_interaction.R",
  "14_morphology_kaplan.R",
  "15_multivariate.R",
  "16_main_figure.R",
  "17_figure_audit.R",
  "18_data_validation.R",
  "19_genet_dashboard.R",
  "22_sensitivity_flagged.R",
  "23_timeseries_diagnostics.R",
  "24_headline_model_comparison.R",
  "26_thermal_context.R",
  "27_variance_partitioning.R",
  "28_multiple_testing.R",
  "29_morphology_prob_contrasts.R",
  "diagnostics/A_continuous_lmm.R",
  "diagnostics/B_morphology_glmm.R",
  "diagnostics/C_cox_diagnostics.R",
  "diagnostics/D_pca_lrt.R",
  "diagnostics/E_design_alignment.R",
  "diagnostics/F_model_reproducibility.R",
  "diagnostics/G_diagnostic_plots.R",
  "20_master_results_table.R",
  "diagnostics/H_spreadsheet_coverage.R",
  "31_rnaseq_covariate_table.R",    # handoff: per-library phenotype covariates for the RNA-seq
  "25_model_diagnostic_coverage.R",
  "30_manuscript_audit.R"          # advisory phenotype reproducibility check — warns (never fails) if phenotype numbers drift
)

t0 <- Sys.time()
for (s in scripts) {
  cat("\n============================================================\n",
      "Running:", s, "\n",
      "============================================================\n", sep = "")
  tic <- Sys.time()
  source(file.path("code", s))
  toc <- Sys.time()
  cat("  -> done in", round(as.numeric(toc - tic, units = "secs"), 1), "s\n")
}
t1 <- Sys.time()
cat("\n\nFull pipeline complete in",
    round(as.numeric(t1 - t0, units = "mins"), 2), "min.\n")

# Save session info for reproducibility provenance
dir.create("output", showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), "output/session_info.txt")
