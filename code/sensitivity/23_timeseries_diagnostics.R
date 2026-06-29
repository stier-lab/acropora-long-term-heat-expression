# =============================================================================
# Purpose: Time-series-specific diagnostics for the repeated-measures responses.
#          The primary LMMs (script 12) model `day` as a LINEAR fixed effect
#          with random intercepts only. For data where the SAME coral is
#          measured repeatedly over time (PAM Fv/Fm and color: ~7 obs/coral),
#          two assumptions are untested by the generic DHARMa diagnostics and
#          are the ones that matter most for time series:
#
#            1. TEMPORAL AUTOCORRELATION of within-coral residuals. If present,
#               treatment×day p-values from the intercept-only model are
#               anti-conservative. Tested by refitting with an AR(1) correlation
#               structure (nlme::lme, corAR1) and an LRT vs the no-correlation
#               model; we also check whether the treatment×day conclusion holds
#               under AR(1). (Molly flagged autocorrelation in her notes.)
#
#            2. LINEARITY OF TIME. A linear `day` term assumes a constant rate
#               of change. We compare a linear vs quadratic time trend (and its
#               interaction with treatment) by LRT/AIC.
#
#          Also reported: random-slope vs random-intercept comparison
#          ((day | id) vs (1 | id)) for the repeated-measures responses.
#
#          Symbiont density is CROSS-SECTIONAL (1 destructive obs/coral, a
#          different cohort at each biopsy day), so within-coral autocorrelation
#          does not apply — only the linearity-of-time check is run for it.
#
# What & why: the manuscript's main models keep time simple — a straight line in
#   `day` with one random intercept per coral. That is easy to interpret, but for
#   repeated measurements of the same coral it can be wrong in two ways that
#   matter: (1) measurements close in time are correlated (autocorrelation),
#   which makes p-values look more significant than they should; and (2) the
#   trajectory may curve rather than rise/fall at a constant rate. This script
#   tests both assumptions and checks whether the treatment×time conclusion still
#   holds under the more careful model. Nothing here replaces the primary models;
#   it documents that they are defensible.
# Input:   data/processed/{pam_clean,color_clean,symbiont_chl_clean}.rds
# Output:  output/tables/23_timeseries_diagnostics.csv
#          output/diagnostics/I_timeseries_report.md
#          figures/diagnostics/I_<response>_acf.png
# =============================================================================

# 00_setup.R loads packages and shared paths (DATA_PROC, TBL_DIR, FIG_DIR, ...).
source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages({ library(nlme) })   # nlme::lme handles corAR1

# ---- Results collector -----------------------------------------------------
# Every diagnostic appends one tidy row via add(); bind_rows(rows) at the end
# assembles the full table. Storing a "conclusion" string with each row lets the
# markdown report present plain-language conclusions.
rows <- list()
add <- function(response, check, statistic, df, p_value, conclusion, detail = "") {
  rows[[length(rows) + 1]] <<- tibble(
    response, check,
    statistic = round(as.numeric(statistic), 4),
    df = as.character(df),
    p_value = if (is.na(p_value)) NA_real_ else signif(as.numeric(p_value), 4),
    conclusion, detail
  )
}

# Where the ACF diagnostic PNGs go.
DIAG_DIR <- file.path(FIG_DIR, "diagnostics")
dir.create(DIAG_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Repeated-measures diagnostics (PAM, color)
# ---------------------------------------------------------------------------
# For one repeated-measures response: run all three diagnostics + the ACF plot.
# Generic column names (.y response, .t time) let the same code serve PAM and color.
ts_repeated <- function(data, response, time = "day", label) {
  data <- data |>
    mutate(thicket = factor(thicket), tank = factor(tank), id = factor(id)) |>
    rename(.y = all_of(response), .t = all_of(time)) |>
    filter(!is.na(.y), !is.na(.t)) |>
    arrange(tank, id, .t)               # AR(1) needs rows in time order within coral

  # Full fixed-effects structure mirrors the primary model (all interactions).
  fixed <- .y ~ treatment * wound * .t * thicket

  # --- 1. AR(1) temporal autocorrelation -----------------------------------
  # Fit twice with nlme::lme: a baseline assuming independent within-coral
  # residuals, then the same model with an AR(1) correlation structure (residuals
  # one timepoint apart correlate by phi, two apart by phi^2, ...). Nested in
  # tank/id. method = "ML" (not REML) so the two are comparable by likelihood.
  base <- tryCatch(
    nlme::lme(fixed, random = ~ 1 | tank/id, data = data, method = "ML",
              control = nlme::lmeControl(opt = "optim", maxIter = 200,
                                         msMaxIter = 200)),
    error = function(e) { message(label, " base lme failed: ", conditionMessage(e)); NULL })
  ar1 <- if (!is.null(base)) tryCatch(
    update(base, correlation = nlme::corAR1(form = ~ .t | tank/id)),
    error = function(e) { message(label, " AR1 lme failed: ", conditionMessage(e)); NULL }) else NULL

  if (!is.null(base) && !is.null(ar1)) {
    # Likelihood-ratio test: is the AR(1) model a significantly better fit? A
    # small p means the autocorrelation is real and the simpler model understated
    # uncertainty. phi is the estimated lag-1 correlation (0 = none, ->1 = strong).
    lr  <- anova(base, ar1)
    phi <- tryCatch(coef(ar1$modelStruct$corStruct, unconstrained = FALSE)[[1]],
                    error = function(e) NA_real_)
    p_ac <- lr$`p-value`[2]
    add(label, "AR(1) temporal autocorrelation (LRT vs no-corr)",
        lr$L.Ratio[2], lr$df[2] - lr$df[1], p_ac,
        if (!is.na(p_ac) && p_ac < 0.05)
          sprintf("autocorrelation PRESENT (phi=%.2f); refit AR(1)", phi)
        else "no significant residual autocorrelation",
        sprintf("phi=%.3f; AIC base=%.1f AR1=%.1f", phi, AIC(base), AIC(ar1)))

    # Does the treatment×day conclusion survive AR(1)? Pull the treatment×time
    # coefficient's t/p from each fit; the regex finds the interaction row whose
    # name contains both "treatment" and the time term (.t).
    get_td <- function(m) {
      tt <- summary(m)$tTable
      r  <- grep("treatment.*:.t$|:.t$", rownames(tt))
      r  <- r[grepl("treatment", rownames(tt)[r])][1]
      if (is.na(r)) return(c(NA, NA))
      tt[r, c("t-value", "p-value")]
    }
    td_b <- get_td(base); td_a <- get_td(ar1)
    # Robustness check: flag only if significance flips (crosses 0.05) between the
    # baseline and AR(1) fits — that would mean the conclusion depended on ignoring
    # autocorrelation.
    add(label, "treatment×time robustness to AR(1)",
        td_a[1], NA, td_a[2],
        sprintf("treatment×time p: base=%.3g, AR(1)=%.3g — %s",
                td_b[2], td_a[2],
                if (is.na(td_a[2])) "n/a" else
                if ((td_b[2] < 0.05) == (td_a[2] < 0.05)) "conclusion unchanged"
                else "CONCLUSION CHANGES under AR(1)"),
        "")
  }

  # --- 2. Random slope vs intercept ----------------------------------------
  # Should each coral get its OWN time slope, or just its own intercept? Compare
  # (1 | id) vs (1 + .t | id) by LRT. Built via as.formula() strings so the
  # random part can be swapped while the fixed part stays identical. (The first
  # m_int fit is immediately overwritten by the explicit string-built version
  # below — only the explicit ones are compared.) The singular-fit check is set
  # to "ignore" because these rich random structures often hit boundaries.
  m_int <- tryCatch(lme4::lmer(
    update(fixed, . ~ . - .t + .t),   # keep formula; randoms differ below
    data = data, REML = FALSE,
    control = lme4::lmerControl(check.conv.singular = .makeCC("ignore", 1e-4))),
    error = function(e) NULL)
  # build explicitly to control random terms
  f_chr <- "(.y) ~ treatment * wound * .t * thicket"
  m_int <- tryCatch(lme4::lmer(
    as.formula(paste(f_chr, "+ (1 | tank) + (1 | id)")),
    data = data, REML = FALSE,
    control = lme4::lmerControl(check.conv.singular = .makeCC("ignore", 1e-4))),
    error = function(e) NULL)
  m_slp <- tryCatch(lme4::lmer(
    as.formula(paste(f_chr, "+ (1 | tank) + (1 + .t | id)")),
    data = data, REML = FALSE,
    control = lme4::lmerControl(check.conv.singular = .makeCC("ignore", 1e-4))),
    error = function(e) NULL)
  if (!is.null(m_int) && !is.null(m_slp)) {
    # Small p => corals really do follow different trajectories and the random
    # slope is warranted; otherwise the intercept-only model is adequate.
    lr <- anova(m_int, m_slp)
    p_rs <- lr$`Pr(>Chisq)`[2]
    add(label, "random slope of time by coral (LRT vs intercept-only)",
        lr$Chisq[2], lr$Df[2], p_rs,
        if (!is.na(p_rs) && p_rs < 0.05)
          "random slope improves fit (individual trajectories differ)"
        else "random intercept sufficient",
        sprintf("AIC int=%.1f slope=%.1f", AIC(m_int), AIC(m_slp)))
  }

  # --- 3. Linearity of time (linear vs quadratic) --------------------------
  # Delegate to the shared helper (also used for the cross-sectional symbionts).
  ts_linearity(data, ".y", ".t", label, has_id = TRUE)

  # --- ACF figure ----------------------------------------------------------
  # Visual companion to test 1: autocorrelation-function plots of the normalized
  # residuals, baseline vs AR(1). Bars beyond the dashed bounds = leftover
  # autocorrelation; a good AR(1) fit should flatten them.
  if (!is.null(base)) {
    png(file.path(DIAG_DIR, paste0("I_", gsub("[^a-z]", "", tolower(label)), "_acf.png")),
        width = 1100, height = 500, res = 150)
    op <- par(mfrow = c(1, 2))
    r <- residuals(base, type = "normalized")
    acf(r, main = paste0(label, ": residual ACF (no-corr model)"), na.action = na.pass)
    if (!is.null(ar1)) {
      r2 <- residuals(ar1, type = "normalized")
      acf(r2, main = paste0(label, ": residual ACF (AR(1) model)"), na.action = na.pass)
    }
    par(op); dev.off()
  }
}

# ---------------------------------------------------------------------------
# Linearity-of-time check (works for repeated or cross-sectional)
# ---------------------------------------------------------------------------
# Fit a linear-in-time model and a quadratic-in-time model, then LRT them.
# has_id toggles the random structure: repeated measures get (1|tank)+(1|id);
# cross-sectional data (symbionts: one obs/coral) get (1|tank) only — there is
# no within-coral replication to support an id random effect.
ts_linearity <- function(data, ycol, tcol, label, has_id) {
  d <- data |> rename(.y = all_of(ycol), .t = all_of(tcol)) |>
    filter(!is.na(.y), !is.na(.t)) |>
    mutate(thicket = factor(thicket), tank = factor(tank))
  if (has_id) d$id <- factor(d$id)
  rand <- if (has_id) "+ (1 | tank) + (1 | id)" else "+ (1 | tank)"
  m_lin  <- tryCatch(lme4::lmer(
    as.formula(paste("(.y) ~ treatment * .t + wound + thicket", rand)),
    data = d, REML = FALSE,
    control = lme4::lmerControl(check.conv.singular = .makeCC("ignore", 1e-4))),
    error = function(e) NULL)
  # poly(.t, 2) adds a curvature (squared-time) term and its treatment interaction.
  m_quad <- tryCatch(lme4::lmer(
    as.formula(paste("(.y) ~ treatment * poly(.t, 2) + wound + thicket", rand)),
    data = d, REML = FALSE,
    control = lme4::lmerControl(check.conv.singular = .makeCC("ignore", 1e-4))),
    error = function(e) NULL)
  if (!is.null(m_lin) && !is.null(m_quad)) {
    # Small p => the quadratic fits significantly better, i.e. the trajectory
    # curves and a single linear slope is only an approximation.
    lr <- anova(m_lin, m_quad)
    p_nl <- lr$`Pr(>Chisq)`[2]
    add(label, "nonlinearity of time (quadratic vs linear, LRT)",
        lr$Chisq[2], lr$Df[2], p_nl,
        if (!is.na(p_nl) && p_nl < 0.05)
          "trajectory is NONLINEAR — linear slope is an approximation"
        else "linear time adequate",
        sprintf("AIC linear=%.1f quad=%.1f", AIC(m_lin), AIC(m_quad)))
  }
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
cat("=== Time-series diagnostics ===\n")

# Both repeated-measures responses get all checks (AR(1) + slopes + curvature).
# day >= 0 so the residual diagnostics describe the same post-wounding model that
# 12_models fits (the pre-treatment baseline is excluded there).
pam <- readRDS(file.path(DATA_PROC, "pam_clean.rds")) |> filter(day >= 0)
ts_repeated(pam, "fv_fm", "day", "PAM Fv/Fm")

color <- readRDS(file.path(DATA_PROC, "color_clean.rds")) |> filter(day >= 0)
ts_repeated(color, "color_num", "day", "Color (D-scale)")

# Symbiont density: cross-sectional (1 obs/coral) — linearity only. Log-transform
# the cell counts (right-skewed, multiplicative) before the curvature check;
# time here is biopsy_day (which cohort was sacrificed), not repeated sampling.
zoox <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
  filter(is.finite(cells_per_cm2), cells_per_cm2 > 0) |>
  mutate(.logz = log(cells_per_cm2))
ts_linearity(zoox |> rename(yy = .logz), "yy", "biopsy_day",
             "log symbionts (cross-sectional)", has_id = FALSE)

# ---- Collect + write -------------------------------------------------------
out <- bind_rows(rows)
write_csv(out, file.path(TBL_DIR, "23_timeseries_diagnostics.csv"))

# ---- Report ----------------------------------------------------------------
report <- c(
  "# Time-series diagnostics (repeated-measures responses)",
  paste0("Generated by code/sensitivity/23_timeseries_diagnostics.R. PAM Fv/Fm and color ",
         "are true repeated measures (~7 obs/coral); symbiont density is ",
         "cross-sectional (1 obs/coral)."),
  "",
  "| Response | Check | Statistic | df | p | Conclusion |",
  "|---|---|---|---|---|---|"
)
report <- c(report, out |>
  mutate(line = sprintf("| %s | %s | %s | %s | %s | %s |",
                        response, check,
                        ifelse(is.na(statistic), "—", format(statistic)),
                        ifelse(is.na(df), "—", df),
                        ifelse(is.na(p_value), "—", format(p_value)),
                        conclusion)) |>
  pull(line))
writeLines(report, file.path(here::here("output", "diagnostics"),
                             "I_timeseries_report.md"))

cat("\n"); print(as.data.frame(out))
cat("\nWrote output/tables/23_timeseries_diagnostics.csv,",
    "output/diagnostics/I_timeseries_report.md, figures/diagnostics/I_*_acf.png\n")
