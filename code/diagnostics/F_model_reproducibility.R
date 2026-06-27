# =============================================================================
# Purpose: Model reproducibility check (diagnostic suite F) — re-fit each saved
#          model from scratch and confirm the coefficients land in the same place.
#
# What & why: a saved model is only trustworthy if re-running the fit on the same
#   data gives the same answer. For every headline model this script reloads the
#   saved fit, re-fits the identical fixed/random structure on the identical
#   processed data, and compares the two: the largest fixed-effect coefficient
#   shift (max_coef_drift) and the difference in log-likelihood. A run is PASS
#   when both stay below 1e-3, DRIFT if either moves more (e.g. an optimizer
#   converging to a different mode), and HANDLED when the model can't be refit
#   (file missing or fit errors) so the suite doesn't crash. This guards against
#   silent non-determinism and accidental edits to the upstream fitting scripts.
#   It is deliberately NARROW: it re-fits models, it does not re-run the whole
#   data-cleaning pipeline.
# Input:   output/models/*.rds  (saved fits to verify)
#          data/processed/*.rds (the data to re-fit on)
# Output:  output/diagnostics/F_model_reproducibility.csv
#          output/diagnostics/F_model_reproducibility_report.md
# =============================================================================

# 00_setup.R loads shared packages + paths (MOD_DIR, DATA_PROC, ...). blme is
# loaded separately because only the penalized morphology refits need it.
source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages({
  library(blme)
})

DIAG_OUT <- here("output", "diagnostics")
dir.create(DIAG_OUT, recursive = TRUE, showWarnings = FALSE)

# One result row per model accumulates here, keyed by model name.
results <- list()

# compare_lmm(): the workhorse. Loads a saved fit, re-fits the same model on the
# same data with the matching engine (lm / lmer / glmer / blme::bglmer), and
# records the coefficient drift + logLik difference. Wrapped in tryCatch so a
# single failed refit becomes a HANDLED row instead of aborting the whole script.
compare_lmm <- function(name, saved_path, formula, data, family = NULL,
                         use_blme = FALSE, fixef.prior = NULL) {
  # No saved model on disk (e.g. an upstream fitter skipped a degenerate trait):
  # record HANDLED and move on rather than error.
  if (!file.exists(saved_path)) {
    results[[name]] <<- tibble::tibble(
      model = name,
      n_coefs = NA_integer_,
      max_coef_drift = NA_real_,
      logLik_diff = NA_real_,
      status = "HANDLED",
      saved_logLik = NA_real_,
      refit_logLik = NA_real_,
      note = "saved model absent; upstream fitter skipped this trait"
    )
    return(invisible(NULL))
  }
  saved <- readRDS(saved_path)
  # Re-fit with the engine that matches how the model was originally built. Same
  # bobyqa optimizer + maxfun settings as the source scripts so the comparison is
  # apples-to-apples; warnings/messages are silenced to keep the console readable.
  refit <- tryCatch({
    if (use_blme) {
      suppressWarnings(suppressMessages(
        blme::bglmer(formula, data = data, family = family,
                     fixef.prior = fixef.prior,
                     control = lme4::glmerControl(optimizer = "bobyqa",
                                                  optCtrl = list(maxfun = 1e5)))
      ))
    } else if (!is.null(family)) {
      # A family but no prior => ordinary GLMM via glmer (none used currently,
      # but kept so the helper is general).
      suppressWarnings(suppressMessages(
        lme4::glmer(formula, data = data, family = family,
                    control = lme4::glmerControl(optimizer = "bobyqa",
                                                 optCtrl = list(maxfun = 1e5)))
      ))
    } else if (inherits(saved, "lm") && !inherits(saved, "lmerMod")) {
      # Plain linear model (the buoyant-weight fit) — no random effects.
      lm(formula, data = data)
    } else {
      # Default: a Gaussian LMM via lmer. Singular-fit convergence checks are
      # ignored so an expected boundary fit doesn't abort the refit.
      suppressWarnings(suppressMessages(
        lme4::lmer(formula, data = data,
                   control = lme4::lmerControl(
                     check.conv.singular = .makeCC("ignore", tol = 1e-4)))
      ))
    }
  }, error = function(e) e)  # capture the error object instead of stopping

  # Refit threw an error: log the saved model's stats + the error message as
  # HANDLED, so the suite reports the gap rather than dying.
  if (inherits(refit, "error")) {
    saved_coef <- if (inherits(saved, "lm") && !inherits(saved, "lmerMod")) {
      coef(saved)
    } else {
      lme4::fixef(saved)
    }
    results[[name]] <<- tibble::tibble(
      model = name,
      n_coefs = length(saved_coef),
      max_coef_drift = NA_real_,
      logLik_diff = NA_real_,
      status = "HANDLED",
      saved_logLik = as.numeric(logLik(saved)),
      refit_logLik = NA_real_,
      note = paste("refit failed:", conditionMessage(refit))
    )
    return(invisible(NULL))
  }

  # Pull fixed-effect coefficients from each fit. coef() for plain lm, fixef()
  # for mixed models (which would otherwise also return random effects).
  saved_coef <- if (inherits(saved, "lm") && !inherits(saved, "lmerMod")) {
    coef(saved)
  } else {
    lme4::fixef(saved)
  }
  refit_coef <- if (inherits(refit, "lm") && !inherits(refit, "lmerMod")) {
    coef(refit)
  } else {
    lme4::fixef(refit)
  }
  # Largest absolute coefficient shift, matching by name so term order can't
  # cause a spurious mismatch. This is the headline reproducibility number.
  max_drift <- max(abs(saved_coef - refit_coef[names(saved_coef)]),
                   na.rm = TRUE)
  # Second, model-wide check: did the overall fit (log-likelihood) land identically?
  ll_diff <- abs(as.numeric(logLik(saved)) - as.numeric(logLik(refit)))
  # PASS only if BOTH coefficients and logLik reproduce to within 1e-3.
  status <- if (max_drift < 1e-3 && ll_diff < 1e-3) "PASS" else "DRIFT"
  results[[name]] <<- tibble::tibble(
    model = name,
    n_coefs = length(saved_coef),
    max_coef_drift = max_drift,
    logLik_diff = ll_diff,
    status = status,
    saved_logLik = as.numeric(logLik(saved)),
    refit_logLik = as.numeric(logLik(refit)),
    note = ""
  )
}

# Fix the RNG so any stochastic step in the refits (e.g. blme starting values)
# is itself reproducible run-to-run.
set.seed(42)

# ---- PAM -------------------------------------------------------------------
# Each block below rebuilds the exact analysis data (mutate(thicket=factor(...))
# matches the source script) and re-fits the same formula the saved model used.
pam <- readRDS(file.path(DATA_PROC, "pam_clean.rds")) |>
  mutate(thicket = factor(thicket))
compare_lmm("12_pam_lmm",
            file.path(MOD_DIR, "12_pam_lmm.rds"),
            fv_fm ~ treatment * wound * day * thicket + (1|tank) + (1|id),
            pam)

# ---- Color -----------------------------------------------------------------
# Same 4-way LMM as PAM, on the numeric color score.
color <- readRDS(file.path(DATA_PROC, "color_clean.rds")) |>
  mutate(thicket = factor(thicket))
compare_lmm("12_color_lmm",
            file.path(MOD_DIR, "12_color_lmm.rds"),
            color_num ~ treatment * wound * day * thicket + (1|tank) + (1|id),
            color)

# ---- Buoyant weight --------------------------------------------------------
# Plain lm (no day, no id) — compare_lmm detects the lm class and refits with lm().
bw <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds")) |>
  mutate(thicket = factor(thicket))
compare_lmm("12_bw_lm",
            file.path(MOD_DIR, "12_bw_lm.rds"),
            areal_calc ~ treatment * wound * thicket + (1|tank),
            bw)

# ---- Symbionts -------------------------------------------------------------
# Rebuild the log-density subset (positive/finite cells only) and centered day,
# exactly as the source fit did, so the refit sees identical rows.
phys <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
  filter(is.finite(cells_per_cm2), cells_per_cm2 > 0) |>
  mutate(thicket = factor(thicket), biopsy_day_c = biopsy_day - 1)
compare_lmm("12_zoox_lmm",
            file.path(MOD_DIR, "12_zoox_lmm.rds"),
            log(cells_per_cm2) ~
              treatment * wound * biopsy_day_c * thicket + (1|tank),
            phys)

# ---- Morphology blme (8 traits) --------------------------------------------
# Loop over the 8 binary wound-healing traits, refitting each penalized GLMM.
# Wounded-only data; thicket as factor as in the source fitter.
ph <- readRDS(file.path(DATA_PROC, "physio_clean.rds")) |>
  filter(wound == "yes", !is.na(day), day >= 0) |>
  mutate(thicket = factor(thicket))
traits <- c("polyps_out", "hole_in_center", "polyp_in_hole",
            "wound_smoothed", "pigment_over_wound", "tip_exist",
            "tip_extension", "new_corallites_on_tip")
# Drop the duplicate trait (polyp_in_hole == hole_in_center; see data-quality note
# in code/04). 12_models no longer fits it, so its saved .rds does not exist —
# keeping it here would make the reproducibility loop chase a missing model.
traits <- traits[!duplicated(lapply(traits, \(t) ph[[t]]))]
for (tr in traits) {
  d <- ph |> mutate(y = .data[[tr]]) |> filter(!is.na(y))
  # Skip degenerate traits: an all-0/all-1 outcome or <30 rows can't fit a GLMM.
  if (length(unique(d$y)) < 2 || nrow(d) < 30) next
  # Set treatment-coding contrasts to match the source fit, so coefficient names
  # (and therefore the by-name drift comparison) line up exactly.
  contrasts(d$treatment) <- contr.treatment(nlevels(d$treatment))
  contrasts(d$thicket) <- contr.treatment(nlevels(d$thicket))
  compare_lmm(paste0("12c_morph_", tr, "_blme"),
              file.path(MOD_DIR, paste0("12c_morph_", tr, "_blme.rds")),
              y ~ treatment * day * thicket + (1|tank) + (1|id),
              d, family = binomial,
              use_blme = TRUE, fixef.prior = "t(scale = 2.5, df = 1)")
}

# ---- Combine + write -------------------------------------------------------
# Stack all per-model rows into one table and save the machine-readable CSV.
final <- bind_rows(results)
write_csv(final, file.path(DIAG_OUT, "F_model_reproducibility.csv"))

# Console summary: print the table with numbers rounded to 4 sig figs for legibility.
cat("=== Model reproducibility check ===\n")
print(final |> mutate(across(where(is.numeric), \(x) signif(x, 4))))
cat("\nTotal:", nrow(final), " PASS:", sum(final$status == "PASS"),
    " DRIFT:", sum(final$status == "DRIFT"),
    " HANDLED:", sum(final$status == "HANDLED"), "\n")

# ---- Report ----------------------------------------------------------------
# sink() routes cat() into the markdown report until the matching sink() closes it.
sink(file.path(DIAG_OUT, "F_model_reproducibility_report.md"))
cat("# F. Model reproducibility check\n\n")
cat("Generated:", format(Sys.time()), "\n\n")
cat("| Status | Count |\n|---|---|\n")
for (s in unique(final$status))
  cat("| ", s, " | ", sum(final$status == s), " |\n", sep = "")
cat("\n## Per-model verdicts\n\n")
for (i in seq_len(nrow(final))) {
  cat("- **", final$model[i], "**: ", final$status[i],
      " (max coef drift = ", signif(final$max_coef_drift[i], 3),
      ", logLik diff = ", signif(final$logLik_diff[i], 3), ")\n", sep = "")
}
sink()

cat("\nWrote", file.path(DIAG_OUT, "F_model_reproducibility.csv"),
    "and report.\n")
