# =============================================================================
# Purpose: Variance partitioning / intraclass correlation (ICC) for every
#          mixed model — how much variance sits in tank vs colony (id) vs
#          residual. Complements the marginal/conditional R² already reported
#          (script 12) by attributing the random-effect variance explicitly.
#
#          Gaussian LMMs: ICC_grp = var_grp / (sum of all variances).
#          Binomial GLMMs: LATENT-scale ICC with residual variance = pi^2/3
#          (Nakagawa & Schielzeth 2010).
#
# What & why: every model in this study has random effects — tank (the unit
#   heated) and colony/id (repeated measures on the same coral). The ICC
#   (intraclass correlation) gives the fraction of the random-effect variance
#   attributable to each grouping factor. Two reasons to report it: (1) a high
#   tank ICC would be a concern — it would mean tank identity, not treatment,
#   drives the response (pseudoreplication risk); a low tank ICC indicates the
#   treatment effect isn't a tank artifact. (2) a high colony/id ICC means corals
#   are consistently different from one another (a genotype/individual signal
#   worth noting). For the Gaussian (continuous) LMMs the residual variance is
#   estimated directly. For binomial GLMMs there is no single residual variance on
#   the 0/1 scale, so we use the standard logistic latent-scale approach: fix the
#   residual at pi^2/3 (the variance of the standard logistic distribution) and
#   compute the ICC on that underlying continuous scale (Nakagawa & Schielzeth 2010).
# Input:   output/models/12_{pam,color,zoox}_lmm.rds,
#          output/models/12c_morph_*_blme.rds
# Output:  output/tables/27_variance_partitioning.csv
# =============================================================================

# 00_setup.R loads packages and shared paths (MOD_DIR, TBL_DIR, ...).
source(here::here("code", "00_setup.R"))

# Latent-scale residual variance for a logit-link binomial model: the variance
# of the standard logistic distribution. Used as the fixed "residual" below.
PI2_3 <- pi^2 / 3

# Helper: read one fitted model and return a tidy ICC table for it. `scale`
# picks how the residual variance is obtained (estimated vs fixed logistic).
icc_from_model <- function(path, label, scale = c("gaussian", "latent")) {
  scale <- match.arg(scale)                         # validate the scale argument
  m <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(m)) return(NULL)                      # skip silently if model absent
  vc <- as.data.frame(lme4::VarCorr(m))             # variance/covariance components
  vc <- vc[is.na(vc$var2), c("grp", "vcov")]   # variance terms only (drop covariances)
  # Gaussian: residual is the model's own estimate; latent (binomial): fixed pi^2/3.
  resid <- if (scale == "gaussian") {
    vc$vcov[vc$grp == "Residual"]
  } else PI2_3
  ranef_rows <- vc[vc$grp != "Residual", , drop = FALSE]   # the grouping factors
  total <- sum(ranef_rows$vcov) + resid             # denominator = all variance
  # Stack one row per grouping factor + one residual row; icc = share of total.
  bind_rows(
    tibble(model = label, component = ranef_rows$grp,
           variance = ranef_rows$vcov, icc = ranef_rows$vcov / total),
    tibble(model = label, component = "Residual",
           variance = resid, icc = resid / total)
  )
}

# ---- Continuous (Gaussian LMM) responses ----------------------------------
# The three main physiology models, each partitioned with the estimated residual.
rows <- bind_rows(
  icc_from_model(file.path(MOD_DIR, "12_pam_lmm.rds"),   "PAM Fv/Fm (LMM)",        "gaussian"),
  icc_from_model(file.path(MOD_DIR, "12_color_lmm.rds"), "Color D-scale (LMM)",    "gaussian"),
  icc_from_model(file.path(MOD_DIR, "12_zoox_lmm.rds"),  "log symbionts (LMM)",    "gaussian")
)

# ---- Binary (morphology GLMM) responses -----------------------------------
# Morphology binomial GLMMs (blme): latent-scale ICC of the tank random effect.
# Glob every saved trait model and partition each on the logistic latent scale.
morph_files <- list.files(MOD_DIR, pattern = "^12c_morph_.*_blme\\.rds$",
                          full.names = TRUE)
morph_icc <- map_dfr(morph_files, function(f) {
  trait <- sub("^12c_morph_(.*)_blme\\.rds$", "\\1", basename(f))   # pull trait name from filename
  icc_from_model(f, sprintf("morph: %s (GLMM)", trait), "latent")
})

# ---- Combine + write ------------------------------------------------------
out <- bind_rows(rows, morph_icc) |>
  mutate(variance = round(variance, 6), icc = round(icc, 4))
write_csv(out, file.path(TBL_DIR, "27_variance_partitioning.csv"))

cat("\n=== Variance partitioning / ICC ===\n")
print(as.data.frame(out))
cat("\nInterpretation: ICC = fraction of (latent) variance attributable to each",
    "grouping factor. High tank ICC would flag tank effects; high id ICC =",
    "consistent among-colony differences.\n")
cat("Wrote output/tables/27_variance_partitioning.csv\n")
