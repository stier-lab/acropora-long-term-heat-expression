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
# Input:   data/processed/{pam_clean,color_clean,symbiont_chl_clean}.rds
# Output:  output/tables/23_timeseries_diagnostics.csv
#          output/diagnostics/I_timeseries_report.md
#          figures/diagnostics/I_<response>_acf.png
# =============================================================================

source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages({ library(nlme) })

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

DIAG_DIR <- file.path(FIG_DIR, "diagnostics")
dir.create(DIAG_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Repeated-measures diagnostics (PAM, color)
# ---------------------------------------------------------------------------
ts_repeated <- function(data, response, time = "day", label) {
  data <- data |>
    mutate(thicket = factor(thicket), tank = factor(tank), id = factor(id)) |>
    rename(.y = all_of(response), .t = all_of(time)) |>
    filter(!is.na(.y), !is.na(.t)) |>
    arrange(tank, id, .t)

  fixed <- .y ~ treatment * wound * .t * thicket

  # --- 1. AR(1) temporal autocorrelation -----------------------------------
  base <- tryCatch(
    nlme::lme(fixed, random = ~ 1 | tank/id, data = data, method = "ML",
              control = nlme::lmeControl(opt = "optim", maxIter = 200,
                                         msMaxIter = 200)),
    error = function(e) { message(label, " base lme failed: ", conditionMessage(e)); NULL })
  ar1 <- if (!is.null(base)) tryCatch(
    update(base, correlation = nlme::corAR1(form = ~ .t | tank/id)),
    error = function(e) { message(label, " AR1 lme failed: ", conditionMessage(e)); NULL }) else NULL

  if (!is.null(base) && !is.null(ar1)) {
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

    # Does the treatment×day conclusion survive AR(1)?
    get_td <- function(m) {
      tt <- summary(m)$tTable
      r  <- grep("treatment.*:.t$|:.t$", rownames(tt))
      r  <- r[grepl("treatment", rownames(tt)[r])][1]
      if (is.na(r)) return(c(NA, NA))
      tt[r, c("t-value", "p-value")]
    }
    td_b <- get_td(base); td_a <- get_td(ar1)
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
  ts_linearity(data, ".y", ".t", label, has_id = TRUE)

  # --- ACF figure ----------------------------------------------------------
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
  m_quad <- tryCatch(lme4::lmer(
    as.formula(paste("(.y) ~ treatment * poly(.t, 2) + wound + thicket", rand)),
    data = d, REML = FALSE,
    control = lme4::lmerControl(check.conv.singular = .makeCC("ignore", 1e-4))),
    error = function(e) NULL)
  if (!is.null(m_lin) && !is.null(m_quad)) {
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

pam <- readRDS(file.path(DATA_PROC, "pam_clean.rds"))
ts_repeated(pam, "fv_fm", "day", "PAM Fv/Fm")

color <- readRDS(file.path(DATA_PROC, "color_clean.rds"))
ts_repeated(color, "color_num", "day", "Color (D-scale)")

# Symbiont density: cross-sectional (1 obs/coral) — linearity only
zoox <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
  filter(is.finite(cells_per_cm2), cells_per_cm2 > 0) |>
  mutate(.logz = log(cells_per_cm2))
ts_linearity(zoox |> rename(yy = .logz), "yy", "biopsy_day",
             "log symbionts (cross-sectional)", has_id = FALSE)

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
