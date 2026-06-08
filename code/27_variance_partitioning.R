# =============================================================================
# Purpose: Variance partitioning / intraclass correlation (ICC) for every
#          mixed model — how much variance sits in tank vs colony (id) vs
#          residual. Complements the marginal/conditional R² already reported
#          (script 12) by attributing the random-effect variance explicitly.
#          (Parallels the latent-scale ICC reporting in the wound-type repos.)
#
#          Gaussian LMMs: ICC_grp = var_grp / (sum of all variances).
#          Binomial GLMMs: LATENT-scale ICC with residual variance = pi^2/3
#          (Nakagawa & Schielzeth 2010).
# Input:   output/models/12_{pam,color,zoox}_lmm.rds,
#          output/models/12c_morph_*_blme.rds
# Output:  output/tables/27_variance_partitioning.csv
# =============================================================================

source(here::here("code", "00_setup.R"))

PI2_3 <- pi^2 / 3

icc_from_model <- function(path, label, scale = c("gaussian", "latent")) {
  scale <- match.arg(scale)
  m <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  vc <- as.data.frame(lme4::VarCorr(m))
  vc <- vc[is.na(vc$var2), c("grp", "vcov")]   # variance terms only
  resid <- if (scale == "gaussian") {
    vc$vcov[vc$grp == "Residual"]
  } else PI2_3
  ranef_rows <- vc[vc$grp != "Residual", , drop = FALSE]
  total <- sum(ranef_rows$vcov) + resid
  bind_rows(
    tibble(model = label, component = ranef_rows$grp,
           variance = ranef_rows$vcov, icc = ranef_rows$vcov / total),
    tibble(model = label, component = "Residual",
           variance = resid, icc = resid / total)
  )
}

rows <- bind_rows(
  icc_from_model(file.path(MOD_DIR, "12_pam_lmm.rds"),   "PAM Fv/Fm (LMM)",        "gaussian"),
  icc_from_model(file.path(MOD_DIR, "12_color_lmm.rds"), "Color D-scale (LMM)",    "gaussian"),
  icc_from_model(file.path(MOD_DIR, "12_zoox_lmm.rds"),  "log symbionts (LMM)",    "gaussian")
)

# Morphology binomial GLMMs (blme): latent-scale ICC of the tank random effect
morph_files <- list.files(MOD_DIR, pattern = "^12c_morph_.*_blme\\.rds$",
                          full.names = TRUE)
morph_icc <- map_dfr(morph_files, function(f) {
  trait <- sub("^12c_morph_(.*)_blme\\.rds$", "\\1", basename(f))
  icc_from_model(f, sprintf("morph: %s (GLMM)", trait), "latent")
})

out <- bind_rows(rows, morph_icc) |>
  mutate(variance = round(variance, 6), icc = round(icc, 4))
write_csv(out, file.path(TBL_DIR, "27_variance_partitioning.csv"))

cat("\n=== Variance partitioning / ICC ===\n")
print(as.data.frame(out))
cat("\nInterpretation: ICC = fraction of (latent) variance attributable to each",
    "grouping factor. High tank ICC would flag tank effects; high id ICC =",
    "consistent among-colony differences.\n")
cat("Wrote output/tables/27_variance_partitioning.csv\n")
