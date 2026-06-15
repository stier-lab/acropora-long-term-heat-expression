# =============================================================================
# Model reproducibility.
#
# For each saved RDS in output/models/, re-source the model from script 12 /
# 12b / 12c / 14 and verify coefficients reproduce within 1e-3.
#
# This is a *narrow* reproducibility check — not a full re-pipeline. It
# refits the same fixed/random structure on the same processed data and
# compares fixed-effect estimates + logLik.
# =============================================================================

source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages({
  library(blme)
})

DIAG_OUT <- here("output", "diagnostics")
dir.create(DIAG_OUT, recursive = TRUE, showWarnings = FALSE)

results <- list()

compare_lmm <- function(name, saved_path, formula, data, family = NULL,
                         use_blme = FALSE, fixef.prior = NULL) {
  saved <- readRDS(saved_path)
  if (use_blme) {
    refit <- suppressWarnings(suppressMessages(
      blme::bglmer(formula, data = data, family = family,
                   fixef.prior = fixef.prior,
                   control = lme4::glmerControl(optimizer = "bobyqa",
                                                optCtrl = list(maxfun = 1e5)))
    ))
  } else if (!is.null(family)) {
    refit <- suppressWarnings(suppressMessages(
      lme4::glmer(formula, data = data, family = family,
                  control = lme4::glmerControl(optimizer = "bobyqa",
                                               optCtrl = list(maxfun = 1e5)))
    ))
  } else if (inherits(saved, "lm") && !inherits(saved, "lmerMod")) {
    refit <- lm(formula, data = data)
  } else {
    refit <- suppressWarnings(suppressMessages(
      lme4::lmer(formula, data = data,
                 control = lme4::lmerControl(
                   check.conv.singular = .makeCC("ignore", tol = 1e-4)))
    ))
  }
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
  max_drift <- max(abs(saved_coef - refit_coef[names(saved_coef)]),
                   na.rm = TRUE)
  ll_diff <- abs(as.numeric(logLik(saved)) - as.numeric(logLik(refit)))
  status <- if (max_drift < 1e-3 && ll_diff < 1e-3) "PASS" else "DRIFT"
  results[[name]] <<- tibble::tibble(
    model = name,
    n_coefs = length(saved_coef),
    max_coef_drift = max_drift,
    logLik_diff = ll_diff,
    status = status,
    saved_logLik = as.numeric(logLik(saved)),
    refit_logLik = as.numeric(logLik(refit))
  )
}

set.seed(42)

# PAM
pam <- readRDS(file.path(DATA_PROC, "pam_clean.rds")) |>
  mutate(thicket = factor(thicket))
compare_lmm("12_pam_lmm",
            file.path(MOD_DIR, "12_pam_lmm.rds"),
            fv_fm ~ treatment * wound * day * thicket + (1|tank) + (1|id),
            pam)

# Color
color <- readRDS(file.path(DATA_PROC, "color_clean.rds")) |>
  mutate(thicket = factor(thicket))
compare_lmm("12_color_lmm",
            file.path(MOD_DIR, "12_color_lmm.rds"),
            color_num ~ treatment * wound * day * thicket + (1|tank) + (1|id),
            color)

# Buoyant weight
bw <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds")) |>
  mutate(thicket = factor(thicket))
compare_lmm("12_bw_lm",
            file.path(MOD_DIR, "12_bw_lm.rds"),
            areal_calc ~ treatment * wound * thicket,
            bw)

# Symbionts
phys <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
  filter(is.finite(cells_per_cm2), cells_per_cm2 > 0) |>
  mutate(thicket = factor(thicket), biopsy_day_c = biopsy_day - 1)
compare_lmm("12_zoox_lmm",
            file.path(MOD_DIR, "12_zoox_lmm.rds"),
            log(cells_per_cm2) ~
              treatment * wound * biopsy_day_c * thicket + (1|tank),
            phys)

# Morphology blme (8 traits)
ph <- readRDS(file.path(DATA_PROC, "physio_clean.rds")) |>
  filter(wound == "yes", !is.na(day), day >= 0) |>
  mutate(thicket = factor(thicket))
traits <- c("polyps_out", "hole_in_center", "polyp_in_hole",
            "wound_smoothed", "pigment_over_wound", "tip_exist",
            "tip_extension", "new_corallites_on_tip")
for (tr in traits) {
  d <- ph |> mutate(y = .data[[tr]]) |> filter(!is.na(y))
  if (length(unique(d$y)) < 2 || nrow(d) < 30) next
  contrasts(d$treatment) <- contr.treatment(nlevels(d$treatment))
  contrasts(d$thicket) <- contr.treatment(nlevels(d$thicket))
  compare_lmm(paste0("12c_morph_", tr, "_blme"),
              file.path(MOD_DIR, paste0("12c_morph_", tr, "_blme.rds")),
              y ~ treatment * day * thicket + (1|tank),
              d, family = binomial,
              use_blme = TRUE, fixef.prior = "t(scale = 2.5, df = 1)")
}

# Combine + write
final <- bind_rows(results)
write_csv(final, file.path(DIAG_OUT, "F_model_reproducibility.csv"))

cat("=== Model reproducibility check ===\n")
print(final |> mutate(across(where(is.numeric), \(x) signif(x, 4))))
cat("\nTotal:", nrow(final), " PASS:", sum(final$status == "PASS"),
    " DRIFT:", sum(final$status == "DRIFT"), "\n")

# Report
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
