# =============================================================================
# Purpose: Numbers for the source-patch results (patches A, C and D).
# Inputs:  output/tables/13_genet_anova.csv             (analysis/13_source_patch_tests.R)
#          output/tables/12_genet_treatment_effects.csv (analysis/12_models.R)
# Output:  output/paper/source_patches_numbers.csv
# Methods: each physiology model was compared with a version without the terms
#   for source patch x temperature (including their interactions with wound
#   status and day) using likelihood-ratio tests. Ambient-minus-heated
#   differences for each patch are model estimates at Day 14, averaged over
#   wound status. Source patches are field labels, not confirmed genotypes.
# =============================================================================
source(here::here("code", "00_setup.R"))

lrt <- readr::read_csv(file.path(TBL_DIR, "13_genet_anova.csv"), show_col_types = FALSE)
eff <- readr::read_csv(file.path(TBL_DIR, "12_genet_treatment_effects.csv"), show_col_types = FALSE)

labels <- c(pam_fvfm = "Fv/Fm", color_dscale = "Color", log_zoox = "Symbionts",
            growth_pct = "Growth")
rows <- list()
for (r in names(labels)) {
  d <- lrt[lrt$response == r, ]
  if (nrow(d) != 1) next
  rows[[length(rows) + 1]] <- paper_row(
    paste0("patch_x_temp_", r), paste0(labels[[r]], ": does the heat effect differ among source patches?"),
    test = "likelihood-ratio test", statistic = d$lrt_chisq, df1 = d$lrt_df, p_value = d$lrt_p,
    n = d$n_obs, source = "output/tables/13_genet_anova.csv")
}
for (r in c("pam_fvfm", "color_dscale")) {
  for (p in c("a", "c", "d")) {
    d <- eff[eff$response == r & eff$thicket == p & eff$contrast == "28C - 31C", ]
    if (nrow(d) == 0) next
    rows[[length(rows) + 1]] <- paper_row(
      paste0("heat_drop_", r, "_patch_", toupper(p)),
      paste0(labels[[r]], ": ambient minus heated on Day 14, patch ", toupper(p), " (mean over wound status)"),
      estimate = mean(d$estimate), test = "model estimate (emmeans)",
      source = "output/tables/12_genet_treatment_effects.csv")
  }
}
write_paper(rows, "source_patches")
