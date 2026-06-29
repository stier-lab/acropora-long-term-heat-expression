# =============================================================================
# Purpose: Time-to-onset analysis for each wound-healing morphological trait.
#          For each wounded coral, compute the day when each binary trait
#          first switches from 0 -> 1 ("first observed event" day). Because
#          morphology was scored at discrete visits, also construct an
#          interval-censored endpoint bounded by the previous scored day and
#          the first scored day with trait expression.
#
#          Three Cox tests per trait (nested for LRT comparison):
#            (a) baseline:    Surv ~ treatment
#            (b) +genet:      Surv ~ treatment + thicket
#            (c) interaction: Surv ~ treatment * thicket
#          LRT (a→b) tests main effect of genet on hazard; LRT (b→c) tests
#          whether the temperature effect depends on genet.
#
# What & why: the question is WHEN each wound-healing milestone first appears, not
#   only whether it ever did. Two simpler approaches fail. A t-test on "day
#   reached" drops every coral that never reaches the milestone (at 31 °C ~67%
#   close their wound but never regenerate a tip), which biases the result toward
#   the fast healers. A "% reached by day 15" discards the timing and the fact
#   that some corals were observed only briefly. Survival analysis handles
#   CENSORING: corals known to reach the milestone only AFTER their last
#   observation (right-censored) or BETWEEN two scoring visits (interval-
#   censored). Each coral contributes the information it has, so survival analysis
#   is appropriate when many corals never reach the endpoint.
#
#   Three complementary survival tools appear below, in order of rigor:
#     1. Interval-censored Weibull AFT (survreg, the PRIMARY inferential test):
#        morphology was scored only at discrete visits (D0/D1/D3/D10/D15), so the
#        exact day a trait switched on is unknown — only that it happened in the
#        window between the last "0" visit and the first "1" visit. The AFT
#        (Accelerated Failure Time) model uses that [lower, upper] interval
#        directly. Its effect is a TIME RATIO: time-ratio > 1 means 31 °C corals
#        take LONGER to reach the milestone (the program is slowed/blocked).
#     2. Kaplan-Meier (KM) curves: the descriptive, assumption-light picture —
#        the cumulative fraction of corals that have reached the milestone by
#        each day, one curve per group. Used for figures; no p-value by itself.
#     3. Cox proportional-hazards (coxph): models the HAZARD (instantaneous rate
#        of reaching the milestone). Its effect is a HAZARD RATIO: HR < 1 means
#        31 °C corals reach the milestone at a slower rate. KM/Cox here use a
#        "first-observed-day" approximation (they pin the event to the visit day,
#        ignoring the interval); they support the figures and provide a check,
#        while the interval-censored AFT is the model used for inference. Cox
#        assumes the HR is constant over time (the "proportional hazards"
#        assumption); cox.zph checks this for every model by testing whether the
#        scaled Schoenfeld residuals trend with time. A small p flags a PH
#        violation, refit with a time-varying coefficient (see end of script).
#   Time ratio (AFT) and hazard ratio (Cox) point the same way but are not the
#   same number: a time ratio describes how much longer it takes; a hazard ratio
#   describes the rate (how much less likely per unit time).
#
# Input:   data/processed/physio_clean.rds
# Output:  output/tables/14_km_event_summary.csv         — per-trait/genet/treatment
#          output/tables/14_cox_hazard_ratios.csv        — HR for 31 vs 28 (first-observed approximation)
#          output/tables/14_interval_survreg.csv         — interval-censored Weibull AFT tests
#          output/tables/14_cox_genet_LRT.csv            — does genet × treatment improve fit?
#          figures/14_morphology_KM.{pdf,png}            — KM by treatment
#          figures/14b_morphology_KM_by_genet.{pdf,png}  — KM by treatment × genet
# =============================================================================

# 00_setup.R loads tidyverse/here, shared paths (DATA_PROC, TBL_DIR, FIG_DIR),
# theme_pub(), and save_fig(). survival provides Surv/survreg/coxph/cox.zph;
# survminer provides helpers for survival plots.
source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages({
  library(survival)
  library(survminer)
})

# ---- Load + set factor levels ---------------------------------------------
# Set the reference level: 28C is the baseline, so every model's
# "treatment31C" term reads as "31 °C relative to 28 °C". thicket = coral genet
# (A/C/D). contr.treatment makes the contrasts explicit dummy coding rather than
# relying on R's session default (keeps results reproducible across machines).
ph <- readRDS(file.path(DATA_PROC, "physio_clean.rds")) |>
  mutate(
    treatment = factor(treatment, levels = c("28C", "31C")),
    wound = factor(wound, levels = c("no", "yes")),
    thicket = factor(thicket)
  )
contrasts(ph$treatment) <- contr.treatment(nlevels(ph$treatment))

# ---- Define the wound-healing milestones ----------------------------------
# Each trait is a binary 0/1 scored at each visit. They span the healing
# sequence from wound closure (early) to skeletal regeneration (late). The
# contrast is wound_smoothed (closure) vs new_corallites_on_tip
# (regeneration) — heat blocks the latter but not the former.
# Trait -> (facet label, plain-language event gloss). Single source of truth: the
# display labels and the CSV glosses both derive from this, so they can never
# drift out of order. The gloss is carried into the output tables so a reader of
# the CSV knows what "reaching the milestone" means. Several traits can switch on
# and off between visits (non-monotonic), which is why we anchor on the
# FIRST observed "1" rather than the final state.
# hole_in_center + polyp_in_hole are one observable, combined upstream (code/04)
# into axial_polyp_formation (M. Brzezinski pers. comm.: the central hole is the
# axial polyp hole, scored together); it enters here as a single milestone.
trait_meta <- tibble::tribble(
  ~trait,                  ~label,                  ~event_interpretation,
  "axial_polyp_formation", "Axial polyp formation", "first observed axial polyp / central hole",
  "wound_smoothed",        "Wound smoothed",        "first observed smoothed wound surface",
  "pigment_over_wound",    "Pigment over wound",    "first observed pigment over wound; expression can be non-monotonic",
  "tip_exist",             "Tip exists",            "first observed visible tip; expression can be non-monotonic",
  "tip_extension",         "Tip extension",         "first observed tip extension; expression can be non-monotonic",
  "new_corallites_on_tip", "New corallites on tip", "first observed new corallites on tip"
)
traits <- trait_meta$trait
trait_interpretation <- trait_meta |> select(trait, event_interpretation)

# ---- Turn repeated 0/1 scores into one survival record per coral -----------
# This implements the censoring logic. For each (coral, trait) we collapse
# the visit-by-visit 0/1 history into:
#   event_day   = first day the trait was seen as 1 (used by KM/Cox); if never
#                 seen, the last day we observed the coral (a right-censored time)
#   event_lower = last day still scored 0 before the first 1 (0 if the very first
#                 visit was already 1) — the EARLIEST the event could have happened
#   event_upper = first day scored 1 — the LATEST it could have happened; Inf if
#                 it never happened (right-censored: event is somewhere after the
#                 last visit, possibly never)
#   event       = 1 if the trait was ever observed, else 0 (the censoring flag)
# The [event_lower, event_upper] pair makes the survreg fit INTERVAL-censored: the
# event is known only to fall inside that window, not the exact day. Only
# wound == "yes" corals can heal a wound, so unwounded controls are excluded here.
compute_events <- function(d, trait) {
  d |>
    filter(wound == "yes", !is.na(day), day >= 0) |>
    mutate(y = .data[[trait]]) |>
    group_by(id, treatment, tank, thicket) |>
    arrange(day, .by_group = TRUE) |>   # order visits so "first 1" is meaningful
    summarise(
      # which(y == 1)[1] = index of the first visit at which the trait was 1.
      # NA means it never switched on -> right-censor at the last observed day.
      event_day = {
        first1 <- which(y == 1)[1]
        if (is.na(first1)) max(day, na.rm = TRUE) else day[first1]
      },
      # lower bound of the interval: the visit just before the first 1 (so we
      # know the event happened after this day). 0 if it was on at visit 1.
      event_lower = {
        first1 <- which(y == 1)[1]
        if (is.na(first1)) max(day, na.rm = TRUE) else if (first1 == 1) 0 else day[first1 - 1]
      },
      # upper bound of the interval: the first day we saw the 1. Inf = never seen
      # = right-censored (survreg reads [lower, Inf] as "event after lower").
      event_upper = {
        first1 <- which(y == 1)[1]
        if (is.na(first1)) Inf else day[first1]
      },
      event = as.integer(any(y == 1, na.rm = TRUE)),  # 1 = observed, 0 = censored
      .groups = "drop"
    ) |>
    mutate(trait = trait)
}
# Run the collapse for every trait and stack the results (one row per coral×trait)
events <- map_dfr(traits, ~ compute_events(ph, .x))
# Re-assert factor levels/contrasts on the summarised table (group_by/summarise
# can drop them) so every downstream model again treats 28C as the reference.
events <- events |>
  mutate(
    treatment = factor(treatment, levels = c("28C", "31C")),
    thicket = factor(thicket)
  ) |>
  left_join(trait_interpretation, by = "trait")
contrasts(events$treatment) <- contr.treatment(nlevels(events$treatment))

# ---- KM summary tables ----------------------------------------------------
# Descriptive median time-to-event per trait × genet × treatment. Medians/IQRs
# are computed only over corals that reached the milestone (event == 1);
# a median of a censored sample is undefined, hence the NA guards.
# Sample sizes for the figure subtitles, computed from the data (one record per
# coral × trait, so collapse to distinct corals first).
n_km_per_trt  <- events |> distinct(id, treatment) |> count(treatment) |>
  pull(n) |> (\(x) if (length(unique(x)) == 1) unique(x) else paste(range(x), collapse = "–"))()
n_km_per_cell <- events |> distinct(id, treatment, thicket) |>
  count(treatment, thicket) |>
  pull(n) |> (\(x) if (length(unique(x)) == 1) unique(x) else paste(range(x), collapse = "–"))()

km_summary <- events |>
  group_by(trait, treatment, thicket) |>
  summarise(
    n_corals    = n(),
    n_events    = sum(event),
    median_day  = if (sum(event) > 0) median(event_day[event == 1]) else NA_real_,
    iqr_low     = if (sum(event) > 0) quantile(event_day[event == 1], 0.25, na.rm = TRUE) else NA_real_,
    iqr_high    = if (sum(event) > 0) quantile(event_day[event == 1], 0.75, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )
write_csv(km_summary, file.path(TBL_DIR, "14_km_event_summary.csv"))

# ---- Interval-censored survival models (PRIMARY inference) -----------------
# Primary inferential model. survreg fits a Weibull Accelerated Failure Time (AFT)
# model: it models the (log) time-to-event directly, so its treatment effect is
# a TIME RATIO. Surv(lower, upper, type = "interval2") feeds in the
# [event_lower, event_upper] windows from above, with Inf upper = right-censored.
# We adjust for thicket (genet) as a fixed effect. Guard: need >=5 events and
# both treatments present, or the model is not estimable.
fit_interval_survreg <- function(tr) {
  d <- events |> filter(trait == tr)
  if (sum(d$event) < 5 || length(unique(d$treatment)) < 2) return(NULL)
  fit <- tryCatch(
    survreg(Surv(event_lower, event_upper, type = "interval2") ~ treatment + thicket,
            data = d, dist = "weibull"),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(tibble(
      trait = tr, term = "treatment", n = nrow(d), n_event = sum(d$event),
      time_ratio_31_vs28 = NA_real_, ratio_lo = NA_real_, ratio_hi = NA_real_,
      z = NA_real_, p = NA_real_, model = "interval-censored Weibull AFT",
      event_interpretation = trait_interpretation$event_interpretation[
        match(tr, trait_interpretation$trait)
      ],
      note = "survreg failed or treatment coefficient unavailable"
    ))
  }
  # Pull the treatment coefficient row by name (its exact label can vary with
  # contrast coding, hence the grep fallback).
  tab <- summary(fit)$table
  trt_row <- if ("treatment31C" %in% rownames(tab)) {
    "treatment31C"
  } else {
    grep("^treatment", rownames(tab), value = TRUE)[1]
  }
  if (is.na(trt_row)) {
    return(tibble(
      trait = tr, term = "treatment", n = nrow(d), n_event = sum(d$event),
      time_ratio_31_vs28 = NA_real_, ratio_lo = NA_real_, ratio_hi = NA_real_,
      z = NA_real_, p = NA_real_, model = "interval-censored Weibull AFT",
      event_interpretation = trait_interpretation$event_interpretation[
        match(tr, trait_interpretation$trait)
      ],
      note = "treatment coefficient unavailable"
    ))
  }
  # AFT coefficients are on the log-time scale, so exp() turns them into a TIME
  # RATIO and exp(est ± 1.96·SE) gives its 95% CI. time_ratio > 1 => 31 °C corals
  # take proportionally longer to reach the milestone (delayed/blocked program);
  # a CI excluding 1 is the "significant" case.
  est <- tab[trt_row, "Value"]
  se <- tab[trt_row, "Std. Error"]
  tibble(
    trait = tr,
    term = trt_row,
    n = nrow(d),
    n_event = sum(d$event),
    time_ratio_31_vs28 = exp(est),
    ratio_lo = exp(est - 1.96 * se),
    ratio_hi = exp(est + 1.96 * se),
    z = tab[trt_row, "z"],
    p = tab[trt_row, "p"],
    model = "interval-censored Weibull AFT",
    event_interpretation = trait_interpretation$event_interpretation[
      match(tr, trait_interpretation$trait)
    ],
    note = "time ratio >1 means later first expression at 31C"
  )
}
# Fit the primary AFT for every trait and collect into one table.
interval_survreg <- map_dfr(traits, fit_interval_survreg)
write_csv(interval_survreg, file.path(TBL_DIR, "14_interval_survreg.csv"))

# ---- Inter-milestone lag: wound closure -> regeneration -------------------
# The result is "heat impairs regeneration, not closure." Quantify it
# per coral: the LAG (days) from achieving wound closure
# (`wound_smoothed`) to forming new skeleton at the tip
# (`new_corallites_on_tip`). Corals that close but never regenerate within the
# experiment are right-censored (no lag) and counted separately — heat is
# expected to inflate that censored fraction.
# Reshape to one row per coral with a column per trait (event_<trait> and
# event_day_<trait>) so we can subtract one milestone's day from another's.
ev_wide <- events |>
  select(id, treatment, thicket, trait, event, event_day) |>
  pivot_wider(names_from = trait, values_from = c(event, event_day))

lag_pairs <- list(
  c(from = "wound_smoothed", to = "new_corallites_on_tip"),  # closure -> regeneration (primary)
  c(from = "wound_smoothed", to = "tip_extension"),          # closure -> tip extension
  c(from = "axial_polyp_formation", to = "new_corallites_on_tip")  # axial polyp -> regeneration
)

lag_long <- map_dfr(lag_pairs, function(p) {
  fe <- paste0("event_", p["from"]);  fd <- paste0("event_day_", p["from"])
  te <- paste0("event_", p["to"]);    td <- paste0("event_day_", p["to"])
  ev_wide |>
    transmute(
      id, treatment, thicket,
      pair       = sprintf("%s -> %s", p["from"], p["to"]),
      reached_from = .data[[fe]] == 1,
      reached_both = .data[[fe]] == 1 & .data[[te]] == 1,
      # lag only defined when both milestones reached and order is sensible
      lag_days   = ifelse(.data[[fe]] == 1 & .data[[te]] == 1,
                          .data[[td]] - .data[[fd]], NA_real_)
    )
})

# Summary by pair × treatment: median/IQR lag among corals reaching both, plus
# the fraction that closed but never regenerated (censored regeneration).
lag_summary <- lag_long |>
  group_by(pair, treatment) |>
  summarise(
    n_closed          = sum(reached_from, na.rm = TRUE),
    n_reached_both    = sum(reached_both, na.rm = TRUE),
    pct_closed_no_regen = round(100 * (sum(reached_from, na.rm = TRUE) -
                                       sum(reached_both, na.rm = TRUE)) /
                                pmax(sum(reached_from, na.rm = TRUE), 1), 1),
    median_lag        = if (sum(reached_both, na.rm = TRUE) > 0)
                          median(lag_days, na.rm = TRUE) else NA_real_,
    iqr_low           = if (sum(reached_both, na.rm = TRUE) > 0)
                          quantile(lag_days, 0.25, na.rm = TRUE) else NA_real_,
    iqr_high          = if (sum(reached_both, na.rm = TRUE) > 0)
                          quantile(lag_days, 0.75, na.rm = TRUE) else NA_real_,
    mean_lag          = mean(lag_days, na.rm = TRUE),
    sd_lag            = sd(lag_days, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(lag_summary, file.path(TBL_DIR, "14_milestone_lag_summary.csv"))

# Wilcoxon rank-sum (non-parametric, robust to the small skewed lag distribution)
# on the primary closure->regeneration lag, 28 vs 31, among corals that reached
# both milestones. This is a secondary, conditional comparison: it sees only the
# corals that regenerated, so it understates the heat effect (the corals most
# affected by heat never regenerate and are absent here). Interpret it alongside
# the pct_closed_no_regen column, which captures that missing fraction.
prim <- lag_long |> filter(pair == "wound_smoothed -> new_corallites_on_tip",
                           reached_both)
lag_test <- if (length(unique(prim$treatment)) == 2 &&
                all(table(prim$treatment) >= 2)) {
  w <- suppressWarnings(wilcox.test(lag_days ~ treatment, data = prim))
  tibble(pair = "wound_smoothed -> new_corallites_on_tip",
         test = "Wilcoxon rank-sum (28C vs 31C lag)",
         W = unname(w$statistic), p_value = w$p.value,
         n_28 = sum(prim$treatment == "28C"),
         n_31 = sum(prim$treatment == "31C"))
} else {
  tibble(pair = "wound_smoothed -> new_corallites_on_tip",
         test = "Wilcoxon rank-sum (28C vs 31C lag)",
         W = NA_real_, p_value = NA_real_,
         n_28 = sum(prim$treatment == "28C"),
         n_31 = sum(prim$treatment == "31C"),
         note = "too few corals reached both milestones in one treatment")
}
write_csv(lag_test, file.path(TBL_DIR, "14_milestone_lag_test.csv"))

cat("\n=== Inter-milestone lag (closure -> regeneration) ===\n")
print(as.data.frame(lag_summary |>
  mutate(across(where(is.numeric), \(x) round(x, 2)))))
cat("\nPrimary lag test (closure -> new corallites, 28C vs 31C):\n")
print(as.data.frame(lag_test))

# ---- Cox PH per trait (overall + per genet) -------------------------------
# Supporting (not primary) analysis: Cox models the HAZARD (rate of reaching the
# milestone) and reports a HAZARD RATIO. It uses event_day (the first-observed
# day), so it ignores the interval — an approximation that pairs with
# the KM figures. strata(thicket) lets each genet keep its own baseline hazard
# shape while sharing one treatment effect (i.e. control for genet without
# assuming the genets share a baseline hazard shape). HR < 1 = 31 °C reaches the
# milestone at a slower rate. Guard: >=5 events or the fit is unstable.
fit_cox_overall <- function(tr) {
  d <- events |> filter(trait == tr)
  if (sum(d$event) < 5) return(NULL)
  fit <- tryCatch(
    coxph(Surv(event_day, event) ~ treatment + strata(thicket), data = d),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  tibble(
    trait      = tr,
    scope      = "overall first-observed approximation (strata=thicket)",
    n          = s$n,
    n_event    = s$nevent,
    HR_31_vs28 = s$conf.int[1, "exp(coef)"],
    HR_lo      = s$conf.int[1, "lower .95"],
    HR_hi      = s$conf.int[1, "upper .95"],
    z          = s$coefficients[1, "z"],
    p          = s$coefficients[1, "Pr(>|z|)"]
  )
}

# Same Cox model, but fit separately within each genet (no strata, since each
# fit is already one genet) to see whether the heat effect differs by genotype.
# Sparse cells (<3 events) return NA rather than an unreliable estimate.
fit_cox_per_genet <- function(tr) {
  events |>
    filter(trait == tr) |>
    group_by(thicket) |>
    group_modify(\(d, k) {
      if (sum(d$event) < 3) {
        return(tibble(HR_31_vs28 = NA_real_, HR_lo = NA_real_,
                      HR_hi = NA_real_, z = NA_real_, p = NA_real_,
                      n = nrow(d), n_event = sum(d$event)))
      }
      fit <- tryCatch(
        coxph(Surv(event_day, event) ~ treatment, data = d),
        error = function(e) NULL
      )
      if (is.null(fit)) {
        return(tibble(HR_31_vs28 = NA_real_, HR_lo = NA_real_,
                      HR_hi = NA_real_, z = NA_real_, p = NA_real_,
                      n = nrow(d), n_event = sum(d$event)))
      }
      s <- summary(fit)
      tibble(HR_31_vs28 = s$conf.int[1, "exp(coef)"],
             HR_lo      = s$conf.int[1, "lower .95"],
             HR_hi      = s$conf.int[1, "upper .95"],
             z          = s$coefficients[1, "z"],
             p          = s$coefficients[1, "Pr(>|z|)"],
             n          = s$n, n_event = s$nevent)
    }) |>
    ungroup() |>
    mutate(trait = tr, scope = paste0("genet=", thicket)) |>
    select(trait, scope, n, n_event, HR_31_vs28, HR_lo, HR_hi, z, p)
}

# Run both Cox variants across all traits and stack into one HR table.
cox_overall   <- map_dfr(traits, fit_cox_overall)
cox_per_genet <- map_dfr(traits, fit_cox_per_genet)
cox_results   <- bind_rows(cox_overall, cox_per_genet) |>
  mutate(time_scale = "first-observed-day approximation")
write_csv(cox_results, file.path(TBL_DIR, "14_cox_hazard_ratios.csv"))

# ---- Proportional-hazards diagnostics for EVERY overall Cox model ----------
# Cox's assumption: the hazard ratio is CONSTANT over time (proportional
# hazards). cox.zph tests this by checking whether the scaled Schoenfeld
# residuals (per-event-time deviations of the covariate from its risk-set mean)
# TREND with time. No trend -> flat residuals -> PH holds; a slope -> the effect
# changes over time -> PH violated. A large chisq / small p flags a violation.
# We run and save the test + plot for every trait (ph_ok = p >= 0.05), not just
# the one that fails, so the assumption is documented for each model.
DIAG_DIR <- file.path(FIG_DIR, "diagnostics")
dir.create(DIAG_DIR, recursive = TRUE, showWarnings = FALSE)
ph_rows <- list()
for (tr in traits) {
  d <- events |> filter(trait == tr)
  if (sum(d$event) < 5) next
  fit <- tryCatch(coxph(Surv(event_day, event) ~ treatment + strata(thicket),
                        data = d), error = function(e) NULL)
  if (is.null(fit)) next
  zph <- tryCatch(cox.zph(fit), error = function(e) NULL)  # the Schoenfeld PH test
  if (is.null(zph)) next
  ph_rows[[tr]] <- tibble(
    trait      = tr,
    n_event    = fit$nevent,
    zph_chisq  = zph$table["treatment", "chisq"],
    zph_p      = zph$table["treatment", "p"],
    ph_ok      = zph$table["treatment", "p"] >= 0.05
  )
  ttl <- sprintf("PH check: %s (zph p=%.3f)", tr, zph$table["treatment", "p"])
  tt  <- as.numeric(zph$time); yy <- as.numeric(zph$y)
  png(file.path(DIAG_DIR, paste0("14_cox_ph_", tr, ".png")),
      width = 800, height = 500, res = 130)
  # plot.cox.zph's smoothing spline needs >4 distinct event times; with fewer
  # it silently skips the panel. Use the spline plot when there are enough
  # distinct times, otherwise a manual scaled-Schoenfeld scatter (always valid).
  if (length(unique(tt)) > 4) {
    ok <- tryCatch({ plot(zph, resid = TRUE, main = ttl); TRUE },
                   error = function(e) FALSE)
  } else ok <- FALSE
  if (!ok) {
    plot(tt, yy, xlab = "Time (days)", ylab = "Scaled Schoenfeld residual",
         main = ttl, pch = 19, col = "grey40",
         ylim = range(c(0, yy), na.rm = TRUE))
    abline(h = 0, lty = 2, col = "grey60")
    if (length(unique(tt)) > 2)
      try(lines(lowess(tt, yy), col = "#D55E00", lwd = 2), silent = TRUE)
  }
  dev.off()
}
cox_ph <- bind_rows(ph_rows)
write_csv(cox_ph, file.path(TBL_DIR, "14_cox_ph_tests.csv"))
cat("\n=== Proportional-hazards (cox.zph) tests, all overall Cox models ===\n")
print(as.data.frame(cox_ph |> mutate(across(where(is.numeric), \(x) round(x, 4)))))

# ---- LRT: does genet × treatment improve over treatment + genet? ----------
# Likelihood-ratio tests on NESTED Cox models. A bigger model can only fit at
# least as well, so the LRT asks whether the extra terms improve fit MORE than
# expected by chance (chi-square on the degrees of freedom added):
#   base  = treatment only
#   add   = treatment + genet      -> LRT(base, add)  tests a genet main effect
#   inter = treatment * genet      -> LRT(add, inter) tests whether the heat
#                                     effect DEPENDS on genet (the interaction)
# Need >=10 events for the interaction model to be estimable.
fit_lrt <- function(tr) {
  d <- events |> filter(trait == tr)
  if (sum(d$event) < 10) return(NULL)
  m_base  <- tryCatch(coxph(Surv(event_day, event) ~ treatment, data = d),
                      error = function(e) NULL)
  m_add   <- tryCatch(coxph(Surv(event_day, event) ~ treatment + thicket, data = d),
                      error = function(e) NULL)
  m_inter <- tryCatch(coxph(Surv(event_day, event) ~ treatment * thicket, data = d),
                      error = function(e) NULL)
  if (any(sapply(list(m_base, m_add, m_inter), is.null))) return(NULL)
  lrt_g  <- anova(m_base, m_add)    # genet main effect
  lrt_gi <- anova(m_add, m_inter)   # genet × treatment interaction
  tibble(
    trait               = tr,
    chisq_genet_main    = lrt_g$Chisq[2],
    df_genet_main       = lrt_g$Df[2],
    p_genet_main        = lrt_g$`Pr(>|Chi|)`[2],
    chisq_genet_trt_int = lrt_gi$Chisq[2],
    df_genet_trt_int    = lrt_gi$Df[2],
    p_genet_trt_int     = lrt_gi$`Pr(>|Chi|)`[2]
  )
}
cox_lrt <- map_dfr(traits, fit_lrt)
write_csv(cox_lrt, file.path(TBL_DIR, "14_cox_genet_LRT.csv"))

# ---- KM facet figure ------------------------------------------------------
# Build the Kaplan-Meier curves manually (rather than survfit) so they integrate
# with ggplot. The KM estimator is a product over each event day of
# (1 - events/at-risk): at each day, the fraction still "surviving" (i.e. not yet
# expressing the trait) is multiplied down. We plot 1 - survival = the
# cumulative % of corals that HAVE reached the milestone by each day. The leading
# c(0, ...) / c(1, ...) seeds the curve at day 0 with 100% not-yet-expressed.
km_curves <- events |>
  group_by(trait, treatment) |>
  group_modify(\(d, k) {
    d <- dplyr::arrange(d, event_day)
    times <- sort(unique(d$event_day))         # distinct days something happened
    # Build a step at every distinct day, NA-safe
    n_at_risk <- sapply(times, \(t) sum(d$event_day >= t))         # still eligible
    n_events  <- sapply(times, \(t) sum(d$event_day == t & d$event == 1))  # events that day
    haz       <- n_events / pmax(n_at_risk, 1)  # day-specific hazard (pmax avoids /0)
    surv      <- cumprod(1 - haz)               # KM survival = running product
    out_day   <- c(0, times)
    out_surv  <- c(1, surv)
    tibble(day = out_day,
           surv = out_surv,
           cum_event = 1 - out_surv)            # what we actually plot
  }) |>
  ungroup() |>
  mutate(trait = factor(trait,
                        levels = trait_meta$trait,
                        labels = trait_meta$label))

# geom_step draws the KM staircase (the curve only changes on event days);
# one panel per trait, blue = 28 °C, orange = 31 °C.
p_km <- ggplot(km_curves, aes(day, cum_event,
                              colour = treatment, group = treatment)) +
  geom_step(linewidth = 0.7) +
  facet_wrap(~ trait, ncol = 4) +
  scale_colour_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                      name = "Temperature") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(x = "Days after wounding",
       y = "Cumulative % corals expressing trait",
       title = "Heat delays the regenerative-tip program but not wound closure",
       subtitle = sprintf("Kaplan-Meier curves for each wound-healing milestone (n = %s corals per treatment)",
                          n_km_per_trt)) +
  theme_pub(9)

save_fig(p_km, "14_morphology_KM", width = 200, height = 130)

# ---- Genet-resolved KM figure --------------------------------------------
# Identical KM construction, but split a third way by genet so each
# treatment × genet cell gets its own curve (line colour = temperature,
# linetype = genet). Reveals whether all three genotypes show the same pattern.
km_curves_genet <- events |>
  group_by(trait, treatment, thicket) |>
  group_modify(\(d, k) {
    d <- dplyr::arrange(d, event_day)
    times <- sort(unique(d$event_day))
    n_at_risk <- sapply(times, \(t) sum(d$event_day >= t))
    n_events  <- sapply(times, \(t) sum(d$event_day == t & d$event == 1))
    haz       <- n_events / pmax(n_at_risk, 1)
    surv      <- cumprod(1 - haz)
    tibble(day = c(0, times), cum_event = 1 - c(1, surv))
  }) |>
  ungroup() |>
  mutate(trait = factor(trait,
                        levels = trait_meta$trait,
                        labels = trait_meta$label))

p_km_genet <- ggplot(km_curves_genet,
                      aes(day, cum_event, colour = treatment,
                          linetype = thicket,
                          group = interaction(treatment, thicket))) +
  geom_step(linewidth = 0.55, alpha = 0.85) +
  facet_wrap(~ trait, ncol = 4) +
  scale_colour_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                      name = "Temperature") +
  scale_linetype_manual(values = c(a = "solid", c = "22", d = "44"),
                        name = "Genet") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(x = "Days after wounding",
       y = "Cumulative % corals expressing trait",
       title = "Genet-specific healing trajectories",
       subtitle = sprintf("Kaplan-Meier curves by treatment × genet (n = %s per cell)",
                          n_km_per_cell)) +
  theme_pub(9)

save_fig(p_km_genet, "14b_morphology_KM_by_genet", width = 210, height = 135)

# ---- Time-varying-coefficient refit for the one PH violation -------------
# Follow-up for the single model that failed the PH check above. When PH is
# violated, a single constant HR is misleading, so we let the treatment effect
# CHANGE with time: tt() supplies a time transform, here multiplying the 31 °C
# indicator by log(t + 1), so the effective HR is exp(coef · log(t + 1)) — it can
# grow or shrink as the experiment proceeds rather than being pinned at one value.
# This is a robustness check on that one cell, not a primary result.
pig_c <- events |>
  filter(trait == "pigment_over_wound", thicket == "c", event_day > 0)

cox_tt_pig_c <- if (sum(pig_c$event) >= 3 && length(unique(pig_c$treatment)) > 1) {
  fit_tt <- tryCatch(
    coxph(Surv(event_day, event) ~ tt(treatment),
          data = pig_c,
          tt   = function(x, t, ...) (as.integer(x == "31C")) * log(t + 1)),
    error = function(e) { message("tt() fit failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(fit_tt)) {
    tibble(trait="pigment_over_wound", scope="genet=c (tt log(t+1))",
           n=nrow(pig_c), n_event=sum(pig_c$event),
           coef=NA_real_, se=NA_real_, z=NA_real_, p=NA_real_,
           note="tt() fit failed")
  } else {
    s <- summary(fit_tt)
    tibble(trait="pigment_over_wound",
           scope="genet=c (tt log(t+1))",
           n=s$n, n_event=s$nevent,
           coef=s$coefficients[1,"coef"],
           se=s$coefficients[1,"se(coef)"],
           z=s$coefficients[1,"z"],
           p=s$coefficients[1,"Pr(>|z|)"],
           note="time-varying coefficient on treatment; HR is exp(coef * log(t+1))")
  }
} else {
  tibble(trait="pigment_over_wound", scope="genet=c (tt log(t+1))",
         n=nrow(pig_c), n_event=sum(pig_c$event),
         coef=NA_real_, se=NA_real_, z=NA_real_, p=NA_real_,
         note="insufficient events")
}
write_csv(cox_tt_pig_c, file.path(TBL_DIR, "14c_cox_tt_pigment_genetC.csv"))

cat("\n=== Cox PH hazard ratios (31 °C vs 28 °C) — overall + per genet ===\n")
print(cox_results |> mutate(across(where(is.numeric), \(x) round(x, 3))))

cat("\n=== LRT: does adding genet × treatment improve Cox fit? ===\n")
print(cox_lrt |> mutate(across(where(is.numeric), \(x) round(x, 3))))

cat("\n=== Time-varying-coefficient refit for pigment_over_wound, genet c ===\n")
print(cox_tt_pig_c |> mutate(across(where(is.numeric), \(x) round(x, 3))))
