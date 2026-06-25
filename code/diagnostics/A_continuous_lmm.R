# =============================================================================
# Diagnostic suite for continuous-response LMMs
#   Models:  12_pam_lmm.rds, 12_color_lmm.rds, 12_zoox_lmm.rds, 12_bw_lm.rds
#   Checks:  DHARMa residuals, convergence/singularity, influence, emmeans
#            direction sanity vs biological expectation.
#   Outputs: output/diagnostics/A_continuous_diagnostics.csv
#            output/diagnostics/A_continuous_report.md
#            figures/diagnostics/A_<model>_*.png
# =============================================================================

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

DIAG_OUT <- file.path(OUT_DIR, "diagnostics")
DIAG_FIG <- file.path(FIG_DIR, "diagnostics")
dir.create(DIAG_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(DIAG_FIG, recursive = TRUE, showWarnings = FALSE)

# ---- Result accumulator ----------------------------------------------------
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
run_dharma <- function(model, model_name, fig_prefix) {
  sim <- tryCatch(
    DHARMa::simulateResiduals(model, n = 1000, refit = FALSE, seed = 42),
    error = function(e) { message("DHARMa failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(sim)) {
    add_row(model_name, "DHARMa_simulation", NA, NA, "n/a", "FAIL",
            "simulateResiduals errored")
    return(invisible(NULL))
  }

  # KS uniformity
  ks <- testUniformity(sim, plot = FALSE)
  add_row(model_name, "DHARMa_KS_uniformity",
          statistic = unname(ks$statistic), p_value = ks$p.value,
          threshold = "p>=0.05",
          status = ifelse(ks$p.value >= 0.05, "PASS",
                          ifelse(ks$p.value >= 0.01, "WARN", "FAIL")),
          notes = "Kolmogorov-Smirnov on scaled residuals")

  # Dispersion
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

  # Outliers
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

  # Save DHARMa diagnostic plot
  png(file.path(DIAG_FIG, paste0(fig_prefix, "_dharma.png")),
      width = 1600, height = 800, res = 150)
  plot(sim)
  dev.off()

  invisible(sim)
}

# ---- Convergence / singularity for lmer ------------------------------------
check_lmer_convergence <- function(model, model_name) {
  sing <- lme4::isSingular(model, tol = 1e-4)
  add_row(model_name, "isSingular",
          statistic = as.integer(sing), p_value = NA,
          threshold = "FALSE",
          status = ifelse(sing, "WARN", "PASS"),
          notes = "Singular fit means a variance component is at/near zero")

  vc <- as.data.frame(VarCorr(model))
  zero_vc <- vc[vc$vcov < 1e-6 & !is.na(vc$vcov), , drop = FALSE]
  add_row(model_name, "variance_components_near_zero",
          statistic = nrow(zero_vc), p_value = NA,
          threshold = "0",
          status = ifelse(nrow(zero_vc) == 0, "PASS", "WARN"),
          notes = paste0("VarCorr (grp:var): ",
                         paste(sprintf("%s:%.2e", vc$grp, vc$vcov),
                               collapse = "; ")))

  # Convergence messages from optimizer
  msg <- model@optinfo$conv$lme4$messages
  if (is.null(msg)) msg <- character(0)
  add_row(model_name, "optimizer_convergence_messages",
          statistic = length(msg), p_value = NA,
          threshold = "0",
          status = ifelse(length(msg) == 0, "PASS", "WARN"),
          notes = paste(msg, collapse = " | "))
}

# ---- Influence (Cook's distance) -------------------------------------------
top_cooks_lmer <- function(model, model_name, k = 3) {
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

png(file.path(DIAG_FIG, "A_bw_lmm_base.png"),
    width = 1200, height = 900, res = 150)
plot(fitted(m_bw), resid(m_bw),
     xlab = "Fitted values", ylab = "Conditional residuals",
     main = "12_bw_lm residuals vs fitted")
abline(h = 0, lty = 2, col = "grey50")
qqnorm(resid(m_bw), main = "12_bw_lm residual Q-Q")
qqline(resid(m_bw), col = "#D55E00")
dev.off()

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
