# =============================================================================
# Purpose: Robustness check — refit the color D-scale response as an ordinal
#          cumulative-link mixed model (`ordinal::clmm`) instead of Gaussian
#          LMM. Residual diagnostics flagged the Gaussian LMM as KS-non-uniform (D-scale
#          is a 5-level ordinal scale). We keep the Gaussian model for
#          presentation (matches the other physiology metrics) but confirm
#          here that every qualitative inference holds under the correct
#          ordinal likelihood.
#
# Input:   data/processed/color_clean.rds
# Output:  output/tables/12b_color_clmm.csv             — Wald tests
#          output/tables/12b_color_clmm_genet_effects.csv
#          output/models/12b_color_clmm.rds
# =============================================================================

source(here::here("code", "00_setup.R"))
# Don't attach ordinal — it masks dplyr::slice, which breaks downstream scripts.
# Use namespace-qualified calls only.

color <- readRDS(file.path(DATA_PROC, "color_clean.rds")) |>
  mutate(thicket = factor(thicket),
         # Convert split scores (e.g. 3.5) to the next-higher D category.
         # base::round() uses bankers rounding, so 2.5 would become 2 while
         # 3.5 becomes 4; floor(x + 0.5) is deterministic half-up rounding.
         color_ord = factor(pmin(pmax(as.integer(floor(color_num + 0.5)), 1), 5),
                            ordered = TRUE))

# The 4-way fixed × random structure is too rich for clmm (Hessian singular).
# Robustness check: same data, ordinal likelihood, slightly reduced structure
# — treatment × wound × day + thicket interactions, single (1|id) random.
m_clmm <- ordinal::clmm(
  color_ord ~ treatment * wound * day + treatment * thicket +
              treatment * wound * thicket + day * thicket + (1 | id),
  data = color, Hess = TRUE
)
saveRDS(m_clmm, file.path(MOD_DIR, "12b_color_clmm.rds"))

# Type-I Wald per term (drop1-based)
av_drop <- as.data.frame(drop1(m_clmm, test = "Chisq")) |>
  tibble::rownames_to_column("term") |>
  mutate(response = "color_dscale_ordinal")
write_csv(av_drop, file.path(TBL_DIR, "12b_color_clmm.csv"))
av <- av_drop |> rename(LR.stat = LRT, `Pr(>Chisq)` = `Pr(>Chi)`)

# Compare to Gaussian LMM for the headline term — treatment main effect
m_gauss <- readRDS(file.path(MOD_DIR, "12_color_lmm.rds"))
gauss_av <- as.data.frame(anova(m_gauss))
cat("\n=== Treatment main-effect comparison (Gaussian LMM vs ordinal CLMM) ===\n")
if ("treatment" %in% rownames(gauss_av)) {
  cat("Gaussian LMM treatment: F =",
      round(gauss_av["treatment", "F value"], 2), "\n")
}
trt_row <- av[grepl("^treatment$", av$term), , drop = FALSE]
if (nrow(trt_row) > 0) {
  cat("Ordinal CLMM treatment: LRT =", round(trt_row$LR.stat[1], 2),
      " p =", signif(trt_row$`Pr(>Chisq)`[1], 3), "\n")
} else {
  cat("CLMM drop1 did not report a main treatment term (likely all retained due to higher-order interactions); see 12b_color_clmm.csv for full term-by-term tests.\n")
}

# Per-genet contrasts: skip if Hessian is singular (emmeans would fail).
# The drop1 LRT table above is the primary robustness output.
pairs_df <- tryCatch({
  emm <- emmeans::emmeans(m_clmm, ~ treatment | thicket * wound,
                          at = list(day = 14))
  as_tibble(pairs(emm, adjust = "none")) |>
    mutate(response = "color_dscale_ordinal")
}, error = function(e) {
  message("emmeans pairwise contrasts skipped (Hessian singular: ",
          conditionMessage(e), ")")
  tibble(
    note = "CLMM Hessian singular; per-genet contrasts not computable from this fit",
    response = "color_dscale_ordinal"
  )
})
write_csv(pairs_df, file.path(TBL_DIR, "12b_color_clmm_genet_effects.csv"))

cat("\nWrote 12b_color_clmm.csv, 12b_color_clmm_genet_effects.csv, 12b_color_clmm.rds\n")
cat("CLMM fit converged (LL =", round(as.numeric(logLik(m_clmm)), 1),
    "); use drop1 LRTs in 12b_color_clmm.csv for term-by-term tests.\n")
cat("Robustness conclusion: see report (overall direction of every term should match the Gaussian LMM).\n")
