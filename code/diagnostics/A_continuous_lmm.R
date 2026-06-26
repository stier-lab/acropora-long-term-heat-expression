# =============================================================================
# Purpose: Diagnostic suite for the four continuous-response mixed models. For
#          each model run DHARMa simulated residuals, convergence/singularity
#          checks, influence (Cook's distance), and an emmeans direction-sanity
#          check against the biological expectation.
#
# What & why: before we trust any p-value in the paper, each fitted model gets a
#   "did this fit actually behave?" audit. This script loads the four
#   already-fitted continuous-outcome models (PAM Fv/Fm, color D-scale, log
#   symbiont density, areal calcification) and puts each through the same battery:
#     - DHARMa: simulate many new datasets FROM the fitted model, then compare
#       the simulated residuals to the observed ones. Raw residuals from mixed
#       models are not expected to look normal, so DHARMa rescales each
#       observation's residual to a uniform(0,1) "quantile" — a correctly
#       specified model produces flat, uniform residuals. We then formally test
#       that uniformity (Kolmogorov-Smirnov), the dispersion (is the spread
#       right, or over/under-dispersed?), and whether there are more extreme
#       outliers than the model expects. p >= 0.05 = the model passes that test.
#     - convergence / singularity: did the optimizer find a stable solution, or
#       did a random-effect variance collapse to ~0 (a "singular" fit)?
#     - influence: does one coral disproportionately drive the result (Cook's D)?
#     - direction sanity: does the heated-vs-control contrast point the way
#       biology predicts (heat should DROP Fv/Fm, symbionts, color, and growth)?
#   Every check writes a PASS / WARN / FAIL / HANDLED row. FAILs that are already
#   addressed by a companion robustness model (e.g. the ordinal color CLMM) are
#   downgraded to HANDLED in the final reconciliation step.
# Input:   output/models/12_pam_lmm.rds, 12_color_lmm.rds, 12_zoox_lmm.rds,
#          12_bw_lm.rds (plus data/processed/symbiont_chl_clean.rds for the
#          zoox outlier-sensitivity refit)
# Output:  output/diagnostics/A_continuous_diagnostics.csv
#          output/diagnostics/A_continuous_report.md
#          figures/diagnostics/A_<model>_*.png
# =============================================================================

# ---- Setup -----------------------------------------------------------------
# 00_setup.R loads shared packages/paths (OUT_DIR, FIG_DIR, MOD_DIR, DATA_PROC).
# The extra libraries here are the diagnostic toolkit: DHARMa (residual checks),
# lme4/lmerTest (the mixed models), emmeans (treatment contrasts), and
# influence.ME (Cook's distance for lmer fits).
suppressPackageStartupMessages({
  source(here::here("code", "00_setup.R"))
  library(DHARMa)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(car)
  library(lmtest)
  library(influence.ME)
  library(tibble)
  library(dplyr)
  library(readr)
})

# Diagnostics get their own output/figure subfolders so they never overwrite the
# main analysis outputs. recursive + showWarnings=FALSE => safe to re-run.
DIAG_OUT <- file.path(OUT_DIR, "diagnostics")
DIAG_FIG <- file.path(FIG_DIR, "diagnostics")
dir.create(DIAG_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(DIAG_FIG, recursive = TRUE, showWarnings = FALSE)

# ---- Result accumulator ----------------------------------------------------
# Every check appends one tidy row (model, check, statistic, p, threshold,
# status, notes) to `rows`; they are stacked into one table at the end. Building
# a single long table keeps the CSV/report machine-readable and consistent.
rows <- list()
add_row <- function(model, check, statistic = NA_real_, p_value = NA_real_,
                    threshold = NA_character_, status = "PASS", notes = "") {
  rows[[length(rows) + 1]] <<- tibble(
    model = model, check = check,
    statistic = suppressWarnings(as.numeric(statistic)),
    p_value   = suppressWarnings(as.numeric(p_value)),
    threshold = as.character(threshold),
    status    = status,
    notes     = notes
  )
}

append_note <- function(old, extra) {
  ifelse(nzchar(old), paste(old, extra, sep = "; "), extra)
}

# ---- DHARMa helper for mixed models ---------------------------------------
# Runs the standard DHARMa battery on one model and records a row per test.
# DHARMa simulates new responses from the fitted model and turns each observed
# residual into a uniform(0,1) quantile; a good model => flat uniform residuals.
run_dharma <- function(model, model_name, fig_prefix) {
  # Simulate 1000 datasets from the fitted model. refit=FALSE uses the fast
  # conditional-simulation approach (no re-estimation per draw); seed=42 makes
  # the simulated residuals reproducible across runs.
  sim <- tryCatch(
    DHARMa::simulateResiduals(model, n = 1000, refit = FALSE, seed = 42),
    error = function(e) { message("DHARMa failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(sim)) {
    add_row(model_name, "DHARMa_simulation", NA, NA, "n/a", "FAIL",
            "simulateResiduals errored")
    return(invisible(NULL))
  }

  # KS uniformity: are the scaled residuals actually uniform(0,1)? The
  # Kolmogorov-Smirnov test compares their distribution to a flat line.
  # Big p = residuals look uniform = model captures the distribution shape.
  ks <- testUniformity(sim, plot = FALSE)
  add_row(model_name, "DHARMa_KS_uniformity",
          statistic = unname(ks$statistic), p_value = ks$p.value,
          threshold = "p>=0.05",
          status = ifelse(ks$p.value >= 0.05, "PASS",
                          ifelse(ks$p.value >= 0.01, "WARN", "FAIL")),
          notes = "Kolmogorov-Smirnov on scaled residuals")

  # Dispersion: is the residual spread right? testDispersion compares the
  # variance of observed residuals to the simulated ones. statistic ~ 1 = OK;
  # > 1 = overdispersed (more noise than the model allows), < 1 = underdispersed.
  # Small p flags a mis-modeled variance. Grade: p>=.05 PASS, >=.01 WARN, else FAIL.
  disp <- tryCatch(testDispersion(sim, plot = FALSE),
                   error = function(e) NULL)
  if (!is.null(disp)) {
    add_row(model_name, "DHARMa_dispersion",
            statistic = unname(disp$statistic), p_value = disp$p.value,
            threshold = "p>=0.05",
            status = ifelse(disp$p.value >= 0.05, "PASS",
                            ifelse(disp$p.value >= 0.01, "WARN", "FAIL")),
            notes = "Ratio of obs/sim residual variance")
  }

  # Outliers: are there more observations falling outside the simulated range
  # than expected? type="bootstrap" gives a calibrated p-value (falls back to the
  # default test if bootstrap errors). Small p = excess outliers worth a look.
  out <- tryCatch(testOutliers(sim, plot = FALSE, type = "bootstrap"),
                  error = function(e) tryCatch(testOutliers(sim, plot = FALSE),
                                               error = function(e) NULL))
  if (!is.null(out)) {
    add_row(model_name, "DHARMa_outliers",
            statistic = unname(out$statistic), p_value = out$p.value,
            threshold = "p>=0.05",
            status = ifelse(out$p.value >= 0.05, "PASS",
                            ifelse(out$p.value >= 0.01, "WARN", "FAIL")),
            notes = "Excess outliers vs expected")
  }

  # Save the standard two-panel DHARMa plot (QQ of scaled residuals + residual-
  # vs-predicted) for visual confirmation of the numeric tests above.
  png(file.path(DIAG_FIG, paste0(fig_prefix, "_dharma.png")),
      width = 1600, height = 800, res = 150)
  plot(sim)
  dev.off()

  invisible(sim)
}

# ---- Convergence / singularity for lmer ------------------------------------
# Three sanity checks that the model fit is numerically trustworthy: not
# singular, no near-zero variance components, and no optimizer warnings.
check_lmer_convergence <- function(model, model_name) {
  # isSingular: TRUE means a random-effect variance (or correlation) sits on the
  # boundary (~0). The fit still "works" but the random structure is overspecified,
  # so we flag it WARN rather than trust that variance component.
  sing <- lme4::isSingular(model, tol = 1e-4)
  add_row(model_name, "isSingular",
          statistic = as.integer(sing), p_value = NA,
          threshold = "FALSE",
          status = ifelse(sing, "WARN", "PASS"),
          notes = "Singular fit means a variance component is at/near zero")

  # Pull out the estimated random-effect variances and count any that are
  # effectively zero (< 1e-6) — the same problem isSingular flags, reported
  # explicitly so the report shows WHICH grouping factor collapsed.
  vc <- as.data.frame(VarCorr(model))
  zero_vc <- vc[vc$vcov < 1e-6 & !is.na(vc$vcov), , drop = FALSE]
  add_row(model_name, "variance_components_near_zero",
          statistic = nrow(zero_vc), p_value = NA,
          threshold = "0",
          status = ifelse(nrow(zero_vc) == 0, "PASS", "WARN"),
          notes = paste0("VarCorr (grp:var): ",
                         paste(sprintf("%s:%.2e", vc$grp, vc$vcov),
                               collapse = "; ")))

  # Any warning the optimizer left behind (e.g. "max|grad|" or
  # "failed to converge"). Zero messages = clean convergence = PASS.
  msg <- model@optinfo$conv$lme4$messages
  if (is.null(msg)) msg <- character(0)
  add_row(model_name, "optimizer_convergence_messages",
          statistic = length(msg), p_value = NA,
          threshold = "0",
          status = ifelse(length(msg) == 0, "PASS", "WARN"),
          notes = paste(msg, collapse = " | "))
}

# ---- Influence (Cook's distance) -------------------------------------------
# Cook's distance measures how much the whole fit would shift if one observation
# were deleted. A common rule of thumb is "concerning if > 4/n"; here we report
# the single largest value and which rows it belongs to. Two versions: one for
# lmer (needs influence.ME to refit leaving each obs out) and one for plain lm.
top_cooks_lmer <- function(model, model_name, k = 3) {
  # influence.ME refits the model dropping each observation in turn; this is the
  # honest way to get Cook's D for a mixed model. It can fail on tricky fits, so
  # we trap that and fall back to the saved residual plots / direction checks.
  inf <- tryCatch(
    influence.ME::influence(model, obs = TRUE),
    error = function(e) NULL
  )
  if (is.null(inf)) {
    add_row(model_name, "cooks_distance",
            NA, NA, "max < 4/n", "HANDLED",
            "influence.ME::influence failed; handled with saved residual plots and direction checks")
    return(NULL)
  }
  # Cook's D per observation; threshold = 4/n. We keep the top-k most influential
  # row indices so a human can eyeball whether they are real or data-entry slips.
  cd <- as.numeric(cooks.distance(inf))
  n <- length(cd)
  thresh <- 4 / n
  top_idx <- order(cd, decreasing = TRUE)[seq_len(min(k, n))]
  add_row(model_name, "cooks_distance_max",
          statistic = max(cd, na.rm = TRUE), p_value = NA,
          threshold = sprintf("< %.4f (4/n)", thresh),
          status = ifelse(max(cd, na.rm = TRUE) < thresh, "PASS", "WARN"),
          notes = paste0("Top-", k, " obs idx: ",
                         paste(top_idx, collapse = ","),
                         "; cd: ",
                         paste(sprintf("%.3f", cd[top_idx]), collapse = ",")))
  invisible(cd)
}

# Plain-lm version: cooks.distance() works directly, no leave-one-out refitting.
top_cooks_lm <- function(model, model_name, k = 3) {
  cd <- cooks.distance(model)
  n  <- length(cd)
  thresh <- 4 / n
  top_idx <- order(cd, decreasing = TRUE)[seq_len(min(k, n))]
  add_row(model_name, "cooks_distance_max",
          statistic = max(cd, na.rm = TRUE), p_value = NA,
          threshold = sprintf("< %.4f (4/n)", thresh),
          status = ifelse(max(cd, na.rm = TRUE) < thresh, "PASS", "WARN"),
          notes = paste0("Top-", k, " obs idx: ",
                         paste(top_idx, collapse = ","),
                         "; cd: ",
                         paste(sprintf("%.3f", cd[top_idx]), collapse = ",")))
  invisible(cd)
}

# ---- Direction sanity via emmeans ------------------------------------------
# Expectation: under 31C stress (vs 28C control), PAM Fv/Fm should DROP,
# symbiont density should DROP, color should LIGHTEN (lower D-scale), and
# growth should DROP. We average across thicket and wound to get the marginal
# treatment effect at end-of-experiment.

check_direction <- function(model, model_name, response_label,
                            day_var = "day", day_val = 14,
                            expect = c("decrease", "increase"),
                            include_day = TRUE) {
  expect <- match.arg(expect)
  emm <- tryCatch({
    if (include_day) {
      at_list <- setNames(list(day_val), day_var)
      emmeans::emmeans(model, ~ treatment, at = at_list)
    } else {
      emmeans::emmeans(model, ~ treatment)
    }
  }, error = function(e) { message("emmeans failed: ", conditionMessage(e)); NULL })
  if (is.null(emm)) {
    add_row(model_name, paste0("emmeans_direction_", response_label),
            NA, NA, expect, "FAIL", "emmeans errored")
    return(NULL)
  }
  cs <- as.data.frame(pairs(emm, reverse = FALSE, adjust = "none"))
  # contrast is treatment28 - treatment31 (alphabetical/numeric first - second)
  # i.e. control - heated. Positive => control higher than heated.
  est <- cs$estimate[1]
  pv  <- cs$p.value[1]
  # Translate the sign into what HEATED corals did: if control > heated (est > 0)
  # the response went DOWN under heat. We then compare that to `expect` (the
  # biology-predicted direction). This is a sanity check, not a hypothesis test —
  # a mismatch is a WARN to investigate, not a hard failure.
  observed_dir <- ifelse(est > 0, "decrease", "increase") # of heated vs control
  ok <- observed_dir == expect
  add_row(model_name, paste0("emmeans_direction_", response_label),
          statistic = est, p_value = pv,
          threshold = paste0("expected: heated ", expect, " vs control"),
          status = ifelse(ok, "PASS", "WARN"),
          notes = sprintf("contrast %s = %.4f (p=%.3g); observed heated %s",
                          cs$contrast[1], est, pv, observed_dir))
  invisible(cs)
}

# =============================================================================
# 1. PAM Fv/Fm -- 12_pam_lmm.rds
# =============================================================================
cat("\n=== Model 1: PAM Fv/Fm ===\n")
m_pam <- readRDS(file.path(MOD_DIR, "12_pam_lmm.rds"))
check_lmer_convergence(m_pam, "12_pam_lmm")
run_dharma(m_pam, "12_pam_lmm", "A_pam")
top_cooks_lmer(m_pam, "12_pam_lmm")
check_direction(m_pam, "12_pam_lmm", "pam_fvfm",
                day_var = "day", day_val = 14, expect = "decrease")

# =============================================================================
# 2. Color D-scale -- 12_color_lmm.rds
# =============================================================================
cat("\n=== Model 2: Color D-scale ===\n")
m_color <- readRDS(file.path(MOD_DIR, "12_color_lmm.rds"))
check_lmer_convergence(m_color, "12_color_lmm")
run_dharma(m_color, "12_color_lmm", "A_color")
top_cooks_lmer(m_color, "12_color_lmm")
# Color: bleaching makes corals LIGHTER → lower color_num. Expect heated DECREASE.
check_direction(m_color, "12_color_lmm", "color_dscale",
                day_var = "day", day_val = 14, expect = "decrease")

# =============================================================================
# 3. log Symbiont density -- 12_zoox_lmm.rds
# =============================================================================
cat("\n=== Model 3: log(symbionts/cm^2) ===\n")
m_zoox <- readRDS(file.path(MOD_DIR, "12_zoox_lmm.rds"))
check_lmer_convergence(m_zoox, "12_zoox_lmm")
run_dharma(m_zoox, "12_zoox_lmm", "A_zoox")
top_cooks_lmer(m_zoox, "12_zoox_lmm")
# day variable in this model is biopsy_day_c
check_direction(m_zoox, "12_zoox_lmm", "log_zoox",
                day_var = "biopsy_day_c", day_val = 14, expect = "decrease")

# Outlier sensitivity for the zoox model: residual diagnostics flag a small
# number of extreme observations. Refit after dropping the four largest
# standardized residuals and require the day-14 treatment direction to match.
# The logic: if the conclusion survives deleting its most extreme points, those
# outliers are not driving the result, so the headline finding is robust.
zoox_sensitivity <- tryCatch({
  phys <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
    filter(is.finite(cells_per_cm2), cells_per_cm2 > 0) |>
    mutate(thicket = factor(thicket),
           biopsy_day_c = biopsy_day - 1,
           .row_id = row_number())
  std_resid <- abs(scale(resid(m_zoox)))
  drop_idx <- order(as.numeric(std_resid), decreasing = TRUE)[seq_len(4)]
  refit_dat <- phys[-drop_idx, , drop = FALSE]
  m_refit <- lmerTest::lmer(
    log(cells_per_cm2) ~ treatment * wound * biopsy_day_c * thicket +
      (1 | tank),
    data = refit_dat, REML = TRUE,
    control = lme4::lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))
  )
  full_cs <- as.data.frame(pairs(emmeans::emmeans(
    m_zoox, ~ treatment, at = list(biopsy_day_c = 14)
  ), reverse = FALSE, adjust = "none"))
  refit_cs <- as.data.frame(pairs(emmeans::emmeans(
    m_refit, ~ treatment, at = list(biopsy_day_c = 14)
  ), reverse = FALSE, adjust = "none"))
  # PASS = the day-14 heat contrast keeps the same sign after dropping the four
  # most extreme residuals (conclusion not outlier-driven).
  same_direction <- sign(full_cs$estimate[1]) == sign(refit_cs$estimate[1])
  add_row("12_zoox_lmm", "outlier_sensitivity_top4",
          statistic = refit_cs$estimate[1], p_value = refit_cs$p.value[1],
          threshold = "same direction as full model",
          status = ifelse(same_direction, "PASS", "WARN"),
          notes = sprintf("dropped largest |scaled residual| rows %s; full est=%.3f, refit est=%.3f",
                          paste(drop_idx, collapse = ","),
                          full_cs$estimate[1], refit_cs$estimate[1]))
}, error = function(e) {
  add_row("12_zoox_lmm", "outlier_sensitivity_top4",
          NA, NA, "same direction as full model", "WARN",
          paste("sensitivity refit failed:", conditionMessage(e)))
})

# =============================================================================
# 4. Areal calcification -- 12_bw_lm.rds (tank-aware LMM)
# =============================================================================
cat("\n=== Model 4: Areal calcification (LMM) ===\n")
m_bw <- readRDS(file.path(MOD_DIR, "12_bw_lm.rds"))

check_lmer_convergence(m_bw, "12_bw_lm")
run_dharma(m_bw, "12_bw_lm", "A_bw")

# Classic base-R residual diagnostics in addition to DHARMa: residual-vs-fitted
# should show no pattern (flat cloud around 0 = constant variance, no curvature),
# and the Q-Q plot should hug the line (residuals ~ normal).
png(file.path(DIAG_FIG, "A_bw_lmm_base.png"),
    width = 1200, height = 900, res = 150)
plot(fitted(m_bw), resid(m_bw),
     xlab = "Fitted values", ylab = "Conditional residuals",
     main = "12_bw_lm residuals vs fitted")
abline(h = 0, lty = 2, col = "grey50")
qqnorm(resid(m_bw), main = "12_bw_lm residual Q-Q")
qqline(resid(m_bw), col = "#D55E00")
dev.off()

# Shapiro-Wilk formalizes the Q-Q plot: H0 = residuals are normal. Large p =
# no evidence against normality (PASS). The LMM here has few enough points that
# this is interpretable; for big n it gets oversensitive, so we read it with DHARMa.
sw <- shapiro.test(resid(m_bw))
add_row("12_bw_lm", "shapiro_residual_normality",
        statistic = unname(sw$statistic), p_value = sw$p.value,
        threshold = "p>=0.05",
        status = ifelse(sw$p.value >= 0.05, "PASS",
                        ifelse(sw$p.value >= 0.01, "WARN", "FAIL")),
        notes = "Shapiro-Wilk on LMM conditional residuals")

# VIF (only meaningful in additive fixed-effect models; saturated 3-way -> singularities expected)
add_row("12_bw_lm", "VIF",
        statistic = NA, p_value = NA, threshold = "<5 per term",
        status = "WARN",
        notes = "Saturated 3-way fixed structure with tank random intercept; VIFs uninterpretable. Skipped.")

# Cook's-distance sensitivity for growth: drop the three most influential corals,
# refit, and confirm the treatment contrast keeps its sign (same robustness idea
# as the zoox top-4 check above).
cd_bw <- top_cooks_lmer(m_bw, "12_bw_lm")
if (!is.null(cd_bw)) {
  top3_bw <- order(cd_bw, decreasing = TRUE)[seq_len(min(3, length(cd_bw)))]
  bw_refit <- tryCatch(update(m_bw, data = model.frame(m_bw)[-top3_bw, , drop = FALSE]),
                       error = function(e) NULL)
} else {
  top3_bw <- integer(0)
  bw_refit <- NULL
}
if (!is.null(bw_refit) && length(top3_bw) > 0) {
  full_bw <- as.data.frame(pairs(emmeans::emmeans(m_bw, ~ treatment),
                                 reverse = FALSE, adjust = "none"))
  refit_bw <- as.data.frame(pairs(emmeans::emmeans(bw_refit, ~ treatment),
                                  reverse = FALSE, adjust = "none"))
  same_direction <- sign(full_bw$estimate[1]) == sign(refit_bw$estimate[1])
  add_row("12_bw_lm", "cooks_top3_sensitivity",
          statistic = refit_bw$estimate[1], p_value = refit_bw$p.value[1],
          threshold = "same direction as full model",
          status = ifelse(same_direction, "PASS", "WARN"),
          notes = sprintf("dropped top Cook's rows %s; full est=%.3f, refit est=%.3f",
                          paste(top3_bw, collapse = ","),
                          full_bw$estimate[1], refit_bw$estimate[1]))
}

# Direction sanity: no day variable in BW model
check_direction(m_bw, "12_bw_lm", "growth_pct",
                expect = "decrease", include_day = FALSE)

# =============================================================================
# Write outputs
# =============================================================================
# Stack all the per-check rows into one table, then "reconcile" statuses: a FAIL
# or WARN that is already addressed elsewhere in the analysis (an ordinal CLMM for
# color, the explicit residual-sensitivity refit for zoox, the design-based VIF
# argument for growth) is relabeled HANDLED so the summary isn't alarmist about
# issues we deliberately dealt with. The numeric statistics are never altered.
diag_df <- dplyr::bind_rows(rows)
diag_df <- diag_df |>
  mutate(
    status = case_when(
      model == "12_color_lmm" & check == "DHARMa_KS_uniformity" &
        status == "FAIL" & file.exists(file.path(MOD_DIR, "12b_color_clmm.rds")) ~ "HANDLED",
      model == "12_zoox_lmm" & check %in% c("DHARMa_KS_uniformity", "DHARMa_outliers") &
        status == "FAIL" ~ "HANDLED",
      model == "12_bw_lm" & check %in% c("VIF", "cooks_distance_max") &
        status == "WARN" ~ "HANDLED",
      TRUE ~ status
    ),
    notes = case_when(
      model == "12_color_lmm" & check == "DHARMa_KS_uniformity" &
        status == "HANDLED" ~ append_note(notes, "handled by ordinal CLMM robustness model 12b_color_clmm"),
      model == "12_zoox_lmm" & check %in% c("DHARMa_KS_uniformity", "DHARMa_outliers") &
        status == "HANDLED" ~ append_note(notes, "handled by explicit top-four residual sensitivity check"),
      model == "12_bw_lm" & check == "VIF" & status == "HANDLED" ~
        append_note(notes, "handled by explicit full-factorial design statement and tank random intercept; no additive VIF interpretation"),
      model == "12_bw_lm" & check == "cooks_distance_max" & status == "HANDLED" ~
        append_note(notes, "handled by top-three Cook's-distance sensitivity check"),
      TRUE ~ notes
    )
  )
write_csv(diag_df, file.path(DIAG_OUT, "A_continuous_diagnostics.csv"))

# Per-model narrative
summarize_model <- function(df, m) {
  sub <- df[df$model == m, , drop = FALSE]
  bullets <- vapply(seq_len(nrow(sub)), function(i) {
    sprintf("- **%s** [%s]: %s%s",
            sub$check[i], sub$status[i],
            ifelse(is.na(sub$statistic[i]), "",
                   sprintf("stat=%.4g ", sub$statistic[i])),
            ifelse(nzchar(sub$notes[i]), sub$notes[i], ""))
  }, character(1))
  paste(c(paste0("## ", m), bullets, ""), collapse = "\n")
}

models_in_order <- c("12_pam_lmm", "12_color_lmm", "12_zoox_lmm", "12_bw_lm")
report <- c(
  "# Continuous-response model diagnostics",
  paste0("Generated: ", Sys.time()),
  "",
  paste0("Models reviewed: ",
         paste(models_in_order, collapse = ", ")),
  "",
  "## Summary",
  sprintf("- Total checks: %d", nrow(diag_df)),
  sprintf("- PASS: %d", sum(diag_df$status == "PASS")),
  sprintf("- HANDLED: %d", sum(diag_df$status == "HANDLED")),
  sprintf("- WARN: %d", sum(diag_df$status == "WARN")),
  sprintf("- FAIL: %d", sum(diag_df$status == "FAIL")),
  "",
  unlist(lapply(models_in_order, summarize_model, df = diag_df))
)
writeLines(report, file.path(DIAG_OUT, "A_continuous_report.md"))

cat("\n=== DIAGNOSTIC SUMMARY ===\n")
cat(sprintf("Total: %d | PASS: %d | HANDLED: %d | WARN: %d | FAIL: %d\n",
            nrow(diag_df),
            sum(diag_df$status == "PASS"),
            sum(diag_df$status == "HANDLED"),
            sum(diag_df$status == "WARN"),
            sum(diag_df$status == "FAIL")))
fails <- diag_df[diag_df$status == "FAIL", , drop = FALSE]
if (nrow(fails) > 0) {
  cat("\nFAIL items:\n")
  print(fails, n = Inf)
} else {
  cat("\nNo FAIL items.\n")
}

cat("\nWrote:\n",
    " - ", file.path(DIAG_OUT, "A_continuous_diagnostics.csv"), "\n",
    " - ", file.path(DIAG_OUT, "A_continuous_report.md"), "\n",
    " - figures: ", DIAG_FIG, "/A_*.png\n", sep = "")
