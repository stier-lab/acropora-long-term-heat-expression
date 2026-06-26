# =============================================================================
# Purpose: Decide whether to upgrade the primary PAM/color models. The
#          time-series diagnostics (script 23) showed the linear / random-
#          intercept models are a simplification (nonlinear trajectories,
#          random slopes improve fit, AR(1) autocorrelation in PAM). Here we
#          fit the upgraded "maximum-rigor" specification and compare it
#          head-to-head with the current primary model on the quantities the
#          manuscript actually reports:
#            - model fit (AIC from ML fits, because fixed effects differ)
#            - significance of the treatment × time interaction
#            - the day-14 heat effect (28C - 31C marginal contrast, with CI)
#          If the upgrade does not change the effect size or the conclusion,
#          the simpler model is retained and this is documented; if it does,
#          the upgrade is adopted.
#
#          Models compared (PAM and color):
#            current : resp ~ treatment*wound*day*thicket + (1|tank)+(1|id)
#            upgraded: resp ~ treatment*wound*poly(day,2)*thicket
#                             + (1+day|id) + (1|tank)        [random slopes + quad time]
#            ar1     : (PAM only) nlme lme, linear day, random slope, corAR1
#
# What & why: script 23 showed the simple models leave a little on the table
#   (curved trajectories, coral-specific slopes, autocorrelation). The natural
#   next question — "so should we switch to the fancier model?" — only matters if
#   the answer the manuscript reports actually changes. So here we fit both the
#   current and an upgraded "maximum-rigor" model and compare them on what the
#   paper claims: the size and significance of the day-14 heat effect and the
#   treatment×time interaction. If the upgrade barely moves the number and keeps
#   the same conclusion, we keep the simpler, more interpretable model and cite
#   this as robustness; if it changes the story, we adopt it. A decision aid, not
#   a hypothesis test.
# Output:  output/tables/24_headline_model_comparison.csv
#          output/diagnostics/J_headline_model_comparison.md
# =============================================================================

# 00_setup.R loads packages and shared paths (DATA_PROC, TBL_DIR, ...).
source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages({ library(nlme); library(emmeans) })

# Solver settings shared by every lmer fit: tolerate singular fits (the rich
# random structures push boundaries) and give bobyqa a high iteration budget.
ctrl <- lme4::lmerControl(check.conv.singular = .makeCC("ignore", 1e-4),
                          optimizer = "bobyqa",
                          optCtrl = list(maxfun = 2e5))

rows <- list()
add <- function(...) rows[[length(rows) + 1]] <<- tibble(...)   # one model -> one row

# ---- The quantity the manuscript reports -----------------------------------
# Estimated-marginal-means contrast: ambient minus heated (28C - 31C) Fv/Fm (or
# colour) at day 14, averaged over wound and genet. emmeans does the heavy
# lifting of evaluating the fitted model at day = 14 and forming the difference;
# we return the estimate, its 95% CI, and the p-value. The *.limit args raise
# emmeans' safety caps so it won't refuse the (large) df calculation.
day14_effect <- function(m, dayval = 14) {
  e <- tryCatch(
    emmeans::emmeans(m, ~ treatment, at = list(day = dayval), lmerTest.limit = 1e5,
                     pbkrtest.limit = 1e5),
    error = function(err) NULL)
  if (is.null(e)) return(c(est = NA, lo = NA, hi = NA, p = NA))
  p <- as.data.frame(pairs(e))            # the contrast estimate + p
  ci <- as.data.frame(confint(pairs(e)))  # its confidence interval
  c(est = p$estimate[1], lo = ci$lower.CL[1], hi = ci$upper.CL[1],
    p = p$p.value[1])
}

# Fit current + upgraded (+ optional AR(1)) models for one response and record
# each one's AIC, treatment×time p, and day-14 effect into the results table.
compare_response <- function(data, response, label, do_ar1 = FALSE) {
  data <- data |>
    mutate(thicket = factor(thicket), tank = factor(tank), id = factor(id)) |>
    rename(.y = all_of(response)) |>
    filter(!is.na(.y), !is.na(day))

  # --- Current primary model ---
  # What the manuscript actually uses: linear day, random intercepts only.
  m_cur <- lme4::lmer(.y ~ treatment * wound * day * thicket +
                        (1 | tank) + (1 | id),
                      data = data, REML = FALSE, control = ctrl)

  # --- Upgraded: random slopes + quadratic time ---
  # The "maximum-rigor" alternative: curvature in time (poly(day,2)) and a
  # per-coral time slope ((1 + day | id)) — exactly the features script 23 flagged.
  m_up <- lme4::lmer(.y ~ treatment * wound * poly(day, 2, raw = TRUE) * thicket +
                       (1 + day | id) + (1 | tank),
                     data = data, REML = FALSE, control = ctrl)

  # Treatment×time significance (Type-II). For the current model that's the
  # single treatment:day term; for the upgraded model there are two interaction
  # terms (linear + quadratic), so we take the smallest p (strongest evidence).
  a_cur <- car::Anova(m_cur, type = 2)
  a_up  <- car::Anova(m_up,  type = 2)
  p_cur <- a_cur[grep("^treatment:day$", rownames(a_cur)), "Pr(>Chisq)"]
  # upgraded: any treatment:poly(day) term
  up_terms <- grep("treatment:poly", rownames(a_up), value = TRUE)
  p_up <- if (length(up_terms)) min(a_up[up_terms, "Pr(>Chisq)"], na.rm = TRUE) else NA

  # The headline number from each model, for direct comparison.
  e_cur <- day14_effect(m_cur)
  e_up  <- day14_effect(m_up)

  add(response = label, model = "current (linear, rand-intercept)",
      AIC = round(AIC(m_cur), 1),
      trt_time_p = signif(p_cur, 3),
      day14_effect = round(e_cur["est"], 4),
      day14_CI = sprintf("[%.3f, %.3f]", e_cur["lo"], e_cur["hi"]),
      day14_p = signif(e_cur["p"], 3))
  add(response = label, model = "upgraded (quad time + rand-slope)",
      AIC = round(AIC(m_up), 1),
      trt_time_p = signif(p_up, 3),
      day14_effect = round(e_up["est"], 4),
      day14_CI = sprintf("[%.3f, %.3f]", e_up["lo"], e_up["hi"]),
      day14_p = signif(e_up["p"], 3))

  # --- AR(1) (PAM only): nlme, linear day, random slope ---
  # A third contender for PAM, where script 23 found real autocorrelation: an
  # nlme model that adds the AR(1) correlation structure on top of a random
  # time slope. corAR1(form = ~ day | id) correlates successive within-coral
  # residuals; returnObject = TRUE keeps a usable fit even if it stops short of
  # full convergence.
  if (do_ar1) {
    m_ar1 <- tryCatch(
      nlme::lme(.y ~ treatment * wound * day * thicket,
                random = ~ 1 + day | id,
                correlation = nlme::corAR1(form = ~ day | id),
                data = data, method = "ML",
                control = nlme::lmeControl(opt = "optim", maxIter = 300,
                                           msMaxIter = 300, returnObject = TRUE)),
      error = function(e) { message(label, " AR1 lme failed: ", conditionMessage(e)); NULL })
    if (!is.null(m_ar1)) {
      # nlme reports per-coefficient t/p (no Anova), so grab the treatment:day row.
      tt <- summary(m_ar1)$tTable
      r <- grep("treatment.*:day$", rownames(tt)); r <- r[1]
      e_ar <- day14_effect(m_ar1)
      add(response = label, model = "AR(1) + rand-slope (nlme)",
          AIC = round(AIC(m_ar1), 1),
          trt_time_p = if (!is.na(r)) signif(tt[r, "p-value"], 3) else NA,
          day14_effect = round(e_ar["est"], 4),
          day14_CI = sprintf("[%.3f, %.3f]", e_ar["lo"], e_ar["hi"]),
          day14_p = signif(e_ar["p"], 3))
    }
  }
}

# ---- Run -------------------------------------------------------------------
# PAM gets the AR(1) contender (autocorrelation was significant); colour does not.
cat("=== Headline model comparison ===\n")
pam <- readRDS(file.path(DATA_PROC, "pam_clean.rds"))
compare_response(pam, "fv_fm", "PAM Fv/Fm", do_ar1 = TRUE)

color <- readRDS(file.path(DATA_PROC, "color_clean.rds"))
compare_response(color, "color_num", "Color (D-scale)", do_ar1 = FALSE)

out <- bind_rows(rows)
write_csv(out, file.path(TBL_DIR, "24_headline_model_comparison.csv"))

# ---- Verdict --------------------------------------------------------------
# Turn the comparison into a yes/no recommendation per response. "worth_it" is
# TRUE only if the upgrade moves the day-14 effect by more than 15% OR flips its
# significance — i.e. the upgrade changes the manuscript's story. Otherwise the
# simpler model stands and this is logged as a robustness check. pct_diff guards
# against divide-by-zero when the current effect is essentially zero.
verdict <- out |>
  group_by(response) |>
  summarise(
    cur_eff = day14_effect[model == "current (linear, rand-intercept)"][1],
    up_eff  = day14_effect[grepl("upgraded", model)][1],
    cur_sig = day14_p[model == "current (linear, rand-intercept)"][1] < 0.05,
    up_sig  = day14_p[grepl("upgraded", model)][1] < 0.05,
    .groups = "drop") |>
  mutate(
    pct_diff = if_else(abs(cur_eff) < 1e-8, NA_real_,
                       round(100 * (up_eff - cur_eff) / cur_eff, 1)),
    same_conclusion = cur_sig == up_sig,
    worth_it = coalesce(abs(pct_diff) > 15, FALSE) | !same_conclusion
  )

report <- c(
  "# Headline model upgrade — is it worth it?",
  "",
  "Compares the current primary PAM/color models (linear `day`, random",
  "intercepts) with an upgraded specification (quadratic time + random slopes;",
  "AR(1) for PAM) on the day-14 heat effect (28C - 31C, averaged over wound",
  "and genet) and the treatment x time interaction. Generated by",
  "code/sensitivity/24_headline_model_comparison.R.",
  "",
  "| Response | Model | AIC | trt×time p | day-14 effect | 95% CI | day-14 p |",
  "|---|---|---|---|---|---|---|"
)
report <- c(report, out |>
  mutate(line = sprintf("| %s | %s | %s | %s | %s | %s | %s |",
                        response, model, AIC, trt_time_p, day14_effect,
                        day14_CI, day14_p)) |> pull(line))
report <- c(report, "", "## Verdict", "")
report <- c(report, verdict |>
  mutate(line = sprintf("- **%s**: day-14 effect %.3f (current) vs %.3f (upgraded), %+.1f%%; conclusion %s. Upgrade %s.",
                        response, cur_eff, up_eff, pct_diff,
                        ifelse(same_conclusion, "unchanged", "CHANGES"),
                        ifelse(worth_it, "WORTH ADOPTING", "not necessary (document as robustness)"))) |>
  pull(line))
writeLines(report, file.path(here::here("output", "diagnostics"),
                             "J_headline_model_comparison.md"))

cat("\n"); print(as.data.frame(out))
cat("\n=== Verdict ===\n"); print(as.data.frame(verdict))
cat("\nWrote output/tables/24_headline_model_comparison.csv,",
    "output/diagnostics/J_headline_model_comparison.md\n")
