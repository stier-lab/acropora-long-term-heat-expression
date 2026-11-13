# =============================================================================
# 04_regeneration.R: healing and regeneration stages (Figure 1B)
#
# Asks whether heat (31 °C vs 28 °C ambient) changes whether, and how fast,
# wounded fragments heal and regrow. Only wounded fragments are used. Main
# experiment: 12 per temperature, scored daily Days 0-15. Microscope
# experiment: 8 per temperature, analysed separately (only it recorded living
# tissue, the coenosarc, covering the cut). Stages in order: healing
# (coenosarc covers the wound), then 1 axial polyp forms, 2 wound surface
# remodels, 3 axial corallite/tip forms, 4 tip extends, 5 new radial
# corallites bud. "Source patch" = the thicket (A/C/D) a fragment came from.
#
# Sections
#   1. Day 15 counts and log-rank tests (both experiments)
#   2. Timing models (main experiment)
#   3. Stage table and stage figures
#   4. Figure 1B: cumulative percent reached
#   5. Numbers for the results text
#
# Reads:  data/processed/physio_clean.rds,
#         data/processed/microscope_physio_clean.rds
# Writes: output/tables/11f_morphology_dataset_{endpoint,logrank}_tests.csv,
#         14_event_intervals.csv, 14_interval_survreg.csv,
#         14_interval_survreg_distributions.csv, 14_cox_hazard_ratios.csv,
#         14_cox_genet_LRT.csv, 36_stage_timing.csv,
#         36_stage_timing_mean_iqr.csv, 39_stage_cumulative_percent.csv
#         figures/39b_stage_timing_cumulative_option (Figure 1B),
#         39d_stage_cumulative_all_stages, 36_stage_timing,
#         36b_stage_violin, 14_morphology_KM
#         output/paper/regeneration_numbers.csv
# Run:    sourced by code/run_all.R (or: Rscript code/run_all.R)
# =============================================================================

source(here::here("code", "00_setup.R"))

suppressPackageStartupMessages(library(survival))


# ---- 1. Day 15 counts and log-rank tests (both experiments) ----------------
# Two simple tests per trait, each experiment tested on its own: (i) Fisher's
# exact test on how many fragments showed the trait on the last scored day,
# and (ii) a log-rank test on the first day each trait was seen. Fisher's
# exact test: an exact test for a 2 x 2 table, safe with small counts.
# Log-rank test: compares two curves of "time until the event" between groups.

## 1a. Load data and describe the traits -------------------------------------
ph <- readRDS(file.path(DATA_PROC, "physio_clean.rds"))
micro <- readRDS(file.path(DATA_PROC, "microscope_physio_clean.rds"))

# One row per scored trait: display labels, the healing phase it belongs to,
# and which experiment(s) scored it (main / micro columns).
trait_meta <- tibble::tribble(
  ~trait,                  ~short_label,       ~label,                    ~phase,              ~main, ~micro,
  "tissue_over_wound",     "Tissue cover",     "Tissue over wound",       "Wound covering",    FALSE, TRUE,
  "polyps_out",            "Polyps out",       "Polyps out",             "General condition", TRUE,  FALSE,
  "axial_polyp_formation", "Axial polyp",      "Axial polyp formation",  "Wound closure",     TRUE,  TRUE,
  "wound_smoothed",        "Wound smoothed",   "Wound smoothed",         "Wound closure",     TRUE,  TRUE,
  "pigment_over_wound",    "Pigment",          "Pigment over wound",     "Wound covering",    TRUE,  TRUE,
  "tip_exist",             "Tip exists",       "Tip exists",             "Regeneration",      TRUE,  TRUE,
  "tip_extension",         "Tip extension",    "Tip extension",          "Regeneration",      TRUE,  TRUE,
  "new_corallites_on_tip", "New corallites",   "New corallites on tip",  "Regeneration",      TRUE,  TRUE,
  "algae_on_wound",        "Algae",            "Algae on wound",         "Fouling",           TRUE,  FALSE
) |>
  mutate(
    trait_order = row_number(),
    # factor(levels = x) keeps the table's row order in plots and tables
    trait = factor(trait, levels = trait),
    short_label = factor(short_label, levels = short_label),
    phase = factor(phase, levels = c("General condition", "Wound covering",
                                     "Wound closure", "Regeneration", "Fouling")),
    phase_short = dplyr::recode(
      as.character(phase),
      `General condition` = "Condition",
      `Wound covering` = "Covering",
      `Wound closure` = "Closure",
      Regeneration = "Regeneration",
      Fouling = "Fouling"
    )
  )

## 1b. Reshape both experiments to one row per fragment x day x trait --------
# Keeps wounded fragments only. dataset_id makes fragment IDs unique across
# the two experiments (the same id number can occur in both).
prep_long <- function(dat, dataset, dataset_label, trait_cols) {
  dat |>
    filter(as.character(wound) == "yes") |>
    select(id, date, day, treatment, tank, thicket, all_of(trait_cols)) |>
    pivot_longer(all_of(trait_cols), names_to = "trait",
                 values_to = "expressed") |>
    mutate(
      dataset = dataset,
      dataset_label = dataset_label,
      dataset_short = if_else(dataset == "physio_morphology", "Main", "Microscope"),
      dataset_id = paste(dataset, id, sep = "_"),
      treatment = factor(as.character(treatment), levels = c("28C", "31C")),
      expressed = as.integer(expressed)
    ) |>
    left_join(trait_meta, by = "trait") |>
    filter(!is.na(expressed))
}

main_traits <- trait_meta |> filter(main) |> pull(trait) |> as.character()
micro_traits <- trait_meta |> filter(micro) |> pull(trait) |> as.character()

long <- bind_rows(
  prep_long(ph, "physio_morphology", "Main morphology", main_traits),
  prep_long(micro, "microscope_physio", "Microscope photos", micro_traits)
) |>
  mutate(
    dataset_label = factor(dataset_label,
                           levels = c("Main morphology", "Microscope photos")),
    facet_event = paste(dataset_short, phase_short, sep = "\n")
  )

## 1c. Summarise each fragment's record for each trait -----------------------
# One row per fragment x trait: was the trait ever seen, the first day it
# was seen, and whether it was present on the last scored day.
events_by_coral <- long |>
  arrange(dataset, trait, treatment, dataset_id, day) |>
  group_by(dataset, dataset_label, dataset_id, id, treatment, thicket, tank,
           trait, short_label, label, trait_order, phase, phase_short,
           facet_event) |>
  summarise(
    n_scored_days = n(),
    first_scored_day = min(day),
    last_scored_day = max(day),
    event = any(expressed == 1),
    first_yes_day = if (any(expressed == 1)) min(day[expressed == 1]) else NA_real_,
    endpoint_day = max(day),
    endpoint_yes = expressed[which.max(day)],
    .groups = "drop"
  )

## 1d. Fisher's exact test on the last scored day ----------------------------
# For each trait: 2 x 2 table of temperature x present/absent on the last day.
# group_modify() runs the function on each trait's rows (.x) and stacks the
# one-row results. The odds ratio adds 0.5 to every cell so that a zero count
# does not give an infinite ratio.
endpoint_tests <- events_by_coral |>
  group_by(dataset, dataset_label, trait, short_label, label, trait_order,
           phase, phase_short) |>
  group_modify(\(.x, .y) {
    d <- .x |> filter(!is.na(endpoint_yes),
                      treatment %in% c("28C", "31C"))
    n_28 <- sum(d$treatment == "28C")
    n_31 <- sum(d$treatment == "31C")
    yes_28 <- sum(d$endpoint_yes[d$treatment == "28C"] == 1)
    yes_31 <- sum(d$endpoint_yes[d$treatment == "31C"] == 1)
    no_28 <- n_28 - yes_28
    no_31 <- n_31 - yes_31
    endpoint_day <- if (nrow(d) > 0) max(d$endpoint_day, na.rm = TRUE) else NA_real_

    # No fragments in one temperature: nothing to compare, return NAs
    if (n_28 == 0 || n_31 == 0) {
      return(tibble(
        endpoint_day = endpoint_day,
        n_28 = n_28, yes_28 = yes_28, prop_28 = NA_real_,
        n_31 = n_31, yes_31 = yes_31, prop_31 = NA_real_,
        risk_diff_31_minus_28 = NA_real_,
        odds_ratio_31_vs_28 = NA_real_,
        fisher_p = NA_real_,
        test_status = "not_estimable"
      ))
    }

    mat <- matrix(c(yes_28, no_28, yes_31, no_31),
                  nrow = 2, byrow = TRUE,
                  dimnames = list(treatment = c("28C", "31C"),
                                  endpoint = c("yes", "no")))
    # tryCatch(): if the test errors, return NULL instead of stopping the script
    ft <- tryCatch(fisher.test(mat), error = function(e) NULL)
    odds_31 <- (yes_31 + 0.5) / (no_31 + 0.5)
    odds_28 <- (yes_28 + 0.5) / (no_28 + 0.5)

    tibble(
      endpoint_day = endpoint_day,
      n_28 = n_28,
      yes_28 = yes_28,
      prop_28 = yes_28 / n_28,
      n_31 = n_31,
      yes_31 = yes_31,
      prop_31 = yes_31 / n_31,
      risk_diff_31_minus_28 = yes_31 / n_31 - yes_28 / n_28,
      odds_ratio_31_vs_28 = odds_31 / odds_28,
      fisher_p = if (is.null(ft)) NA_real_ else ft$p.value,
      test_status = if (is.null(ft)) "not_estimable" else "computed"
    )
  }) |>
  ungroup() |>
  mutate(
    test = "Fisher exact test at last scored day",
    test_scope = "cohorts tested separately"
  ) |>
  arrange(dataset, trait_order)

write_csv(endpoint_tests,
          file.path(TBL_DIR, "11f_morphology_dataset_endpoint_tests.csv"))

## 1e. Log-rank test on the first day each trait was seen --------------------
# event_time = first day seen, or the last scored day for fragments that never
# showed it (these are "censored": we only know it had not happened yet).
# This treats the first-seen day as the exact event day; section 2 uses the
# better interval-censored models for the main experiment.
logrank_tests <- events_by_coral |>
  mutate(event_time = if_else(event, first_yes_day, last_scored_day)) |>
  group_by(dataset, dataset_label, trait, short_label, label, trait_order,
           phase, phase_short) |>
  group_modify(\(.x, .y) {
    d <- .x |>
      filter(!is.na(event_time), treatment %in% c("28C", "31C")) |>
      mutate(treatment = factor(as.character(treatment),
                                levels = c("28C", "31C")),
             event_int = as.integer(event))
    n_28 <- sum(d$treatment == "28C")
    n_31 <- sum(d$treatment == "31C")
    n_event_28 <- sum(d$event_int[d$treatment == "28C"])
    n_event_31 <- sum(d$event_int[d$treatment == "31C"])

    # Missing a temperature, or no events at all: test not possible
    if (n_28 == 0 || n_31 == 0 || sum(d$event_int) == 0) {
      return(tibble(
        n_28 = n_28, n_31 = n_31,
        n_event_28 = n_event_28, n_event_31 = n_event_31,
        logrank_chisq = NA_real_, logrank_df = NA_real_,
        logrank_p = NA_real_,
        test_status = "not_estimable"
      ))
    }

    # survdiff() is the log-rank test; Surv(time, status) marks censoring
    fit <- tryCatch(
      survdiff(Surv(event_time, event_int) ~ treatment, data = d),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(tibble(
        n_28 = n_28, n_31 = n_31,
        n_event_28 = n_event_28, n_event_31 = n_event_31,
        logrank_chisq = NA_real_, logrank_df = NA_real_,
        logrank_p = NA_real_,
        test_status = "not_estimable"
      ))
    }

    # df = number of groups - 1; p from the chi-square distribution
    df <- length(fit$n) - 1
    tibble(
      n_28 = n_28,
      n_31 = n_31,
      n_event_28 = n_event_28,
      n_event_31 = n_event_31,
      logrank_chisq = unname(fit$chisq),
      logrank_df = df,
      logrank_p = pchisq(fit$chisq, df = df, lower.tail = FALSE),
      test_status = "computed"
    )
  }) |>
  ungroup() |>
  mutate(
    test = "Log-rank test on first-observed day",
    test_scope = paste(
      "cohorts tested separately;",
      "main-experiment timing results use the interval-censored models (04_regeneration.R, section 2)"
    )
  ) |>
  arrange(dataset, trait_order)

write_csv(logrank_tests,
          file.path(TBL_DIR, "11f_morphology_dataset_logrank_tests.csv"))


# ---- 2. Timing models (main experiment) ------------------------------------
# How much does heat speed up or slow down each stage? A stage is "reached" on
# the first day it was scored yes. Fragments that never reached it by Day 15
# stay in as censored (not reached yet), not dropped. Three model types:
#   - interval-censored AFT models (the main timing result),
#   - Cox models stratified by source patch (hazard ratios),
#   - likelihood-ratio tests for source patch x temperature.
# Kaplan-Meier curves are drawn as a supporting figure.

## 2a. Reload data and list the timing traits --------------------------------
ph <- readRDS(file.path(DATA_PROC, "physio_clean.rds")) |>
  mutate(
    treatment = factor(treatment, levels = c("28C", "31C")),
    wound = factor(wound, levels = c("no", "yes")),
    thicket = factor(thicket)
  )

# Treatment contrasts: 28C is the reference, so model coefficients are 31C vs 28C
contrasts(ph$treatment) <- contr.treatment(nlevels(ph$treatment))

# Note: this replaces the section 1 trait_meta with the six timing traits
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

## 2b. Event intervals: when could each stage have started? ------------------
# Interval-censored: we only know the stage began between the last "no" day
# (event_lower) and the first "yes" day (event_upper). Not reached by the end:
# lower = last scored day, upper = Inf (right-censored). Seen on the first
# visit: lower = 0.
compute_events <- function(d, trait) {
  d |>
    filter(wound == "yes", !is.na(day), day >= 0) |>
    mutate(y = .data[[trait]]) |>
    filter(!is.na(y)) |>
    group_by(id, treatment, tank, thicket) |>
    arrange(day, .by_group = TRUE) |>   # order visits so "first 1" is meaningful
    summarise(
      event_day = {
        first1 <- which(y == 1)[1]
        if (is.na(first1)) max(day, na.rm = TRUE) else day[first1]
      },
      event_lower = {
        first1 <- which(y == 1)[1]
        if (is.na(first1)) max(day, na.rm = TRUE) else if (first1 == 1) 0 else
          day[first1 - 1]
      },
      event_upper = {
        first1 <- which(y == 1)[1]
        if (is.na(first1)) Inf else day[first1]
      },
      event = as.integer(any(y == 1, na.rm = TRUE)),  # 1 = observed, 0 = censored
      .groups = "drop"
    ) |>
    mutate(trait = trait)
}

events <- map_dfr(traits, ~ compute_events(ph, .x))
events <- events |>
  mutate(
    treatment = factor(treatment, levels = c("28C", "31C")),
    thicket = factor(thicket)
  ) |>
  left_join(trait_interpretation, by = "trait")

# Guarantee no stage was already present on Day 0 (the day of clipping). That
# would give an event time of exactly 0, which the log-scale timing models
# cannot fit; the fit would fail and silently return a row of NAs.
stopifnot(all(events$event_upper > 0))

write_csv(events, file.path(TBL_DIR, "14_event_intervals.csv"))

contrasts(events$treatment) <- contr.treatment(nlevels(events$treatment))


## 2c. Interval-censored log-logistic AFT models -----------------------------
# AFT (accelerated failure time) model: models log(time to the stage), so
# exp(coefficient) is a time ratio. Time ratio > 1 = heated fragments took
# longer. Source patch (thicket) is included as a fixed effect.
#
# Distribution: log-logistic (chosen 2026-11-11). Compared by AIC across
# Weibull, lognormal, log-logistic and exponential
# (14_interval_survreg_distributions.csv), it is the best fit or within 0.6
# AIC of it for every reported stage. It does not converge for pigment over
# the wound (not reported), which falls back to lognormal; the `model` column
# says which was used. Weibull, used until 2026-11-11, fitted worst for tip
# formation and tip extension (dAIC 22, 12) and made tip extension look slower
# under heat.
AFT_DIST <- "loglogistic"
AFT_FALLBACK <- "lognormal"

fit_interval_survreg <- function(tr) {
  d <- events |> filter(trait == tr)
  # Skip stages with fewer than 5 events or only one temperature
  if (sum(d$event) < 5 || length(unique(d$treatment)) < 2) return(NULL)

  # A warning (e.g. non-convergence) counts as failure, like an error
  fit_dist <- function(dist) tryCatch(
    survreg(Surv(event_lower, event_upper, type = "interval2") ~ treatment + thicket,
            data = d, dist = dist),
    error = function(e) NULL, warning = function(w) NULL
  )

  # Try log-logistic first; if it fails, use the lognormal fallback
  used <- AFT_DIST
  fit <- fit_dist(AFT_DIST)
  if (is.null(fit)) { used <- AFT_FALLBACK; fit <- fit_dist(AFT_FALLBACK) }

  # Both failed: return a row of NAs with a note
  if (is.null(fit)) {
    return(tibble(
      trait = tr, term = "treatment", n = nrow(d), n_event = sum(d$event),
      time_ratio_31_vs28 = NA_real_, ratio_lo = NA_real_, ratio_hi = NA_real_,
      z = NA_real_, p = NA_real_, model = paste("interval-censored", used, "AFT"),
      event_interpretation = trait_interpretation$event_interpretation[
        match(tr, trait_interpretation$trait)
      ],
      note = "survreg failed or treatment coefficient unavailable"
    ))
  }

  # Find the 31C-vs-28C coefficient row in the model summary
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
      z = NA_real_, p = NA_real_, model = paste("interval-censored", used, "AFT"),
      event_interpretation = trait_interpretation$event_interpretation[
        match(tr, trait_interpretation$trait)
      ],
      note = "treatment coefficient unavailable"
    ))
  }

  # Coefficient is on the log-time scale: exp() gives the time ratio, and
  # exp(estimate +/- 1.96 SE) gives its 95% confidence interval
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
    model = paste("interval-censored", used, "AFT"),
    event_interpretation = trait_interpretation$event_interpretation[
      match(tr, trait_interpretation$trait)
    ],
    note = "time ratio >1 means later first expression at 31C"
  )
}

interval_survreg <- map_dfr(traits, fit_interval_survreg)
write_csv(interval_survreg, file.path(TBL_DIR, "14_interval_survreg.csv"))

## 2d. Distribution check: AIC for each candidate distribution ---------------
# Fits the same model under four distributions and records AIC (lower = better
# fit) and the time ratio, so the choice above can be checked. delta_AIC =
# AIC minus the best AIC for that stage; `used` marks the one reported.
aft_dists <- map_dfr(traits, function(tr) {
  d <- events |> filter(trait == tr)
  map_dfr(c("weibull", "lognormal", "loglogistic", "exponential"), function(dist) {
    f <- tryCatch(survreg(Surv(event_lower, event_upper, type = "interval2") ~ treatment + thicket,
                          data = d, dist = dist), error = function(e) NULL)
    if (is.null(f)) return(tibble(trait = tr, dist = dist, converged = FALSE))
    tab <- summary(f)$table; r <- grep("^treatment", rownames(tab), value = TRUE)[1]
    tibble(trait = tr, dist = dist, converged = TRUE, AIC = AIC(f),
           time_ratio_31_vs28 = exp(tab[r, "Value"]), p = tab[r, "p"])
  })
}) |>
  group_by(trait) |>
  mutate(delta_AIC = AIC - min(AIC, na.rm = TRUE)) |>
  ungroup() |>
  left_join(
    interval_survreg |>
      transmute(trait, used_dist = sub("interval-censored (.*) AFT", "\\1", model)),
    by = "trait"
  ) |>
  mutate(used = dist == used_dist) |>
  select(-used_dist) |>
  mutate(across(where(is.double), \(x) signif(x, 4)))

write_csv(aft_dists, file.path(TBL_DIR, "14_interval_survreg_distributions.csv"))

## 2e. Cox models: hazard ratio for heat -------------------------------------
# Cox proportional-hazards model: hazard ratio (HR) > 1 = heated fragments
# reached the stage sooner. strata(thicket) lets each source patch have its own
# baseline timing, so the heat effect is compared within patches. Uses the
# first-seen day as the event time (an approximation; the AFT models above
# handle the interval properly).
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

# Same Cox model fitted separately within each source patch ("genet")
fit_cox_per_genet <- function(tr) {
  events |>
    filter(trait == tr) |>
    group_by(thicket) |>
    group_modify(\(d, k) {
      # Too few events in this patch: return NAs
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

cox_overall   <- map_dfr(traits, fit_cox_overall)
cox_per_genet <- map_dfr(traits, fit_cox_per_genet)
cox_results   <- bind_rows(cox_overall, cox_per_genet) |>
  mutate(time_scale = "first-observed-day approximation")

write_csv(cox_results, file.path(TBL_DIR, "14_cox_hazard_ratios.csv"))

## 2f. Likelihood-ratio tests for source patch and patch x temperature -------
# Likelihood-ratio test (LRT): compares nested models; a small p means the
# larger model fits meaningfully better. Base (temperature only) vs +patch
# tests the patch main effect; +patch vs patch x temperature tests whether the
# heat effect differs among patches. Needs at least 10 events.
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

## 2g. Kaplan-Meier curves (supporting figure) -------------------------------
# Kaplan-Meier: the share of fragments that have NOT yet reached a stage,
# updated each day an event happens. Computed by hand here; we plot
# 1 - survival = cumulative share that reached the stage.
km_curves <- events |>
  group_by(trait, treatment) |>
  group_modify(\(d, k) {
    d <- dplyr::arrange(d, event_day)
    times <- sort(unique(d$event_day))         # distinct days something happened
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

p_km <- ggplot(km_curves, aes(day, cum_event,
                              colour = treatment, group = treatment)) +
  geom_step(linewidth = 0.7) +
  facet_wrap(~ trait, ncol = 4) +
  scale_colour_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                      name = "Temperature") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(x = "Days after wounding",
       y = "Cumulative % of fragments that reached the stage") +
  theme_pub(9)

save_fig(p_km, "14_morphology_KM", width = 200, height = 130)


# ---- 3. Stage table and stage figures --------------------------------------
# Builds the stage table: the first day each fragment reached each of six
# stages, then counts and first-day summaries per temperature (averages only
# where every fragment reached the stage). Stage 1 (tissue cover) comes from
# the microscope experiment; stages 2-6 from the main experiment. (In the
# manuscript these are "healing" and stages 1-5; see section 4.)

## 3a. Define the six stages and which data each comes from ------------------
physio <- readRDS(file.path(DATA_PROC, "physio_clean.rds"))
micro <- readRDS(file.path(DATA_PROC, "microscope_physio_clean.rds"))

stage_map <- tibble::tribble(
  ~stage_order, ~phase,          ~stage,                                ~dataset,              ~dataset_label,        ~trait,                  ~measurement_note,
  1L,           "Healing",       "Coenosarc barrier re-established",    "microscope_physio",   "Microscope photos",   "tissue_over_wound",     "microscope-only coenosarc/tissue-cover score",
  2L,           "Regeneration",  "Axial polyp forms",                  "physio_morphology",   "Main morphology",     "axial_polyp_formation", "gross-photo central axial-polyp/hole score",
  3L,           "Regeneration",  "Wound surface remodels",             "physio_morphology",   "Main morphology",     "wound_smoothed",        "gross-photo smoothed wound-surface score",
  4L,           "Regeneration",  "Axial corallite or tip forms",        "physio_morphology",   "Main morphology",     "tip_exist",             "gross-photo visible regenerated-tip score",
  5L,           "Regeneration",  "Axial corallite or tip extends",      "physio_morphology",   "Main morphology",     "tip_extension",         "gross-photo tip-extension score",
  6L,           "Regeneration",  "New radial corallites bud",           "physio_morphology",   "Main morphology",     "new_corallites_on_tip", "gross-photo new-corallites-on-tip score"
) |>
  mutate(
    phase = factor(phase, levels = c("Healing", "Regeneration")),
    stage = factor(stage, levels = stage)
  )

## 3b. First day each fragment reached each stage ----------------------------
# One row per fragment x stage, Days 0-15 only, wounded fragments only.
compute_stage_events <- function(dat, stage_row) {
  tr <- stage_row$trait
  dat |>
    filter(as.character(wound) == "yes", !is.na(day), day >= 0, day <= 15) |>
    select(id, treatment, tank, thicket, day, y = all_of(tr)) |>
    filter(!is.na(y)) |>
    mutate(treatment = factor(as.character(treatment), levels = c("28C", "31C"))) |>
    group_by(id, treatment, tank, thicket) |>
    arrange(day, .by_group = TRUE) |>
    summarise(
      first_scored_day = min(day),
      last_scored_day = max(day),
      event = any(y == 1, na.rm = TRUE),
      first_yes_day = if (any(y == 1, na.rm = TRUE)) {
        min(day[y == 1], na.rm = TRUE)
      } else NA_real_,
      .groups = "drop"
    ) |>
    mutate(
      stage_order = stage_row$stage_order,
      phase = stage_row$phase,
      stage = stage_row$stage,
      dataset = stage_row$dataset,
      dataset_label = stage_row$dataset_label,
      trait = stage_row$trait,
      measurement_note = stage_row$measurement_note
    )
}

# pmap_dfr() runs the function once per row of stage_map (columns become
# arguments) and stacks the results. Fragments not reached are placed at
# plot_day 16.35, the "not reached" column of the dot plot.
stage_events <- purrr::pmap_dfr(stage_map, function(stage_order, phase, stage,
                                                    dataset, dataset_label,
                                                    trait, measurement_note) {
  stage_row <- tibble(stage_order = stage_order, phase = phase, stage = stage,
                      dataset = dataset, dataset_label = dataset_label,
                      trait = trait, measurement_note = measurement_note)
  dat <- if (dataset == "microscope_physio") micro else physio
  compute_stage_events(dat, stage_row)
}) |>
  mutate(
    stage = factor(stage, levels = stage_map$stage),
    phase = factor(phase, levels = c("Healing", "Regeneration")),
    treatment = factor(treatment, levels = c("28C", "31C")),
    event_status = if_else(event, "Reached", "Not reached by final day"),
    plot_day = if_else(event, first_yes_day, 16.35),
    stage_label = paste0(as.character(phase), ": ", as.character(stage)),
    stage_label = factor(stage_label,
                         levels = rev(paste0(stage_map$phase, ": ",
                                             stage_map$stage)))
  ) |>
  arrange(stage_order, treatment, id)

write_csv(stage_events,
          file.path(TBL_DIR, "36_stage_timing.csv"))

## 3c. Counts and first-day summaries per stage x temperature ----------------
# Mean/median/IQR of the first day are kept only when every fragment at both
# temperatures reached the stage; otherwise they would describe only the
# faster fragments and are set to NA.
stage_summary <- stage_events |>
  group_by(stage_order, phase, stage, dataset_label, trait, measurement_note,
           treatment) |>
  summarise(
    n_corals = n(),
    n_reached = sum(event),
    n_not_reached = n_corals - n_reached,
    pct_reached = round(100 * n_reached / n_corals, 1),
    mean_first_day = if (n_reached > 0) {
      round(mean(first_yes_day[event], na.rm = TRUE), 2)
    } else NA_real_,
    median_first_day = if (n_reached > 0) {
      round(median(first_yes_day[event]), 2)
    } else NA_real_,
    q25_first_day = if (n_reached > 0) {
      round(unname(quantile(first_yes_day[event], 0.25)), 2)
    } else NA_real_,
    q75_first_day = if (n_reached > 0) {
      round(unname(quantile(first_yes_day[event], 0.75)), 2)
    } else NA_real_,
    iqr_first_day = if (n_reached > 0) {
      round(IQR(first_yes_day[event]), 2)
    } else NA_real_,
    .groups = "drop"
  ) |>
  group_by(stage_order) |>
  mutate(complete_both_temperatures = all(n_reached == n_corals),
         across(c(mean_first_day, median_first_day, q25_first_day,
                  q75_first_day, iqr_first_day),
                ~ if_else(complete_both_temperatures, .x, NA_real_))) |>
  ungroup() |>
  arrange(stage_order, treatment)

# Reader-friendly version of the same table
stage_summary_display <- stage_summary |>
  transmute(
    stage_order,
    phase = as.character(phase),
    stage = as.character(stage),
    dataset = dataset_label,
    treatment = dplyr::recode(as.character(treatment),
                              `28C` = "28 °C", `31C` = "31 °C"),
    n_reached,
    n_scored = n_corals,
    n_not_reached,
    pct_reached,
    complete_both_temperatures,
    mean_day = mean_first_day,
    median_day = median_first_day,
    iqr_days = if_else(
      complete_both_temperatures,
      sprintf("%.2f-%.2f", q25_first_day, q75_first_day),
      NA_character_
    )
  )

write_csv(stage_summary_display,
          file.path(TBL_DIR, "36_stage_timing_mean_iqr.csv"))

## 3d. Dot plot of first day per fragment (figure 36) ------------------------
# One row per stage; 28 °C dots sit just below the row line and 31 °C just
# above (-/+ 0.18). "n reached / n scored" labels go to the right.
label_data <- stage_summary |>
  mutate(
    stage_label = paste0(as.character(phase), ": ", as.character(stage)),
    stage_label = factor(stage_label,
                         levels = rev(paste0(stage_map$phase, ": ",
                                             stage_map$stage))),
    reached_label = paste0(n_reached, "/", n_corals),
    y_plot = as.numeric(stage_label) +
      if_else(treatment == "28C", -0.18, 0.18)
  )

stage_levels <- levels(stage_events$stage_label)

# Fragments with the same stage, temperature and day would overlap, so stack
# them in a small grid (up to 4 per row) around the true position.
# %% = remainder, %/% = integer division.
stage_events_plot <- stage_events |>
  group_by(stage_label, treatment, event_status, plot_day) |>
  arrange(id, .by_group = TRUE) |>
  mutate(
    n_in_cell = n(),
    dot_index = row_number(),
    n_cols = pmin(n_in_cell, 4L),
    stack_col = ((dot_index - 1L) %% n_cols) + 1L,
    stack_row = ((dot_index - 1L) %/% n_cols) + 1L,
    n_rows = max(stack_row),
    x_offset = if_else(n_in_cell == 1L, 0,
                       (stack_col - (n_cols + 1) / 2) * 0.11),
    y_offset = if_else(n_in_cell == 1L, 0,
                       (stack_row - (n_rows + 1) / 2) * 0.08),
    y_plot = as.numeric(stage_label) +
      if_else(treatment == "28C", -0.18, 0.18) + y_offset,
    x_plot = plot_day + x_offset
  ) |>
  ungroup()

# Invisible points (NA x/y) that exist only to build the shape key
status_key <- tibble(
  x = NA_real_,
  y = NA_real_,
  event_status = factor(c("Reached", "Not reached by final day"),
                        levels = c("Reached", "Not reached by final day"))
)

# Grey band on the right = the "not reached" column; dotted line = Day 15
p_stage <- ggplot() +
  annotate("rect", xmin = 15.75, xmax = 16.85, ymin = -Inf, ymax = Inf,
           fill = "grey94", colour = NA) +
  geom_vline(xintercept = 15, linetype = "dotted",
             linewidth = 0.35, colour = "grey55") +
  geom_point(
    data = filter(stage_events_plot, event_status == "Reached"),
    aes(x_plot, y_plot, fill = treatment),
    shape = 21, colour = "white", size = 1.75, stroke = 0.25, alpha = 0.95
  ) +
  geom_point(
    data = filter(stage_events_plot, event_status == "Not reached by final day"),
    aes(x_plot, y_plot, colour = treatment),
    shape = 4, size = 1.95, stroke = 0.8, alpha = 0.95
  ) +
  geom_point(
    data = status_key,
    aes(x, y, shape = event_status),
    colour = "#444444", fill = "#444444", size = 2.8, stroke = 0.8,
    inherit.aes = FALSE, na.rm = TRUE
  ) +
  geom_text(
    data = label_data,
    aes(x = 17.15, y = y_plot, label = reached_label,
        colour = treatment),
    inherit.aes = FALSE,
    size = 2.45,
    hjust = 0,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = PAL_TEMP, name = NULL) +
  scale_colour_manual(values = PAL_TEMP, guide = "none") +
  scale_shape_manual(values = c(`Reached` = 21,
                                `Not reached by final day` = 4),
                     name = NULL) +
  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(shape = 21, colour = "white", size = 3)
    ),
    shape = guide_legend(
      order = 2,
      override.aes = list(colour = "#444444", fill = "#444444", size = 3)
    )
  ) +
  scale_x_continuous(
    breaks = c(0, 1, 3, 5, 10, 15, 16.35),
    labels = c("0", "1", "3", "5", "10", "15", "not reached"),
    limits = c(-0.2, 18.0),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  scale_y_continuous(
    breaks = seq_along(stage_levels),
    labels = stage_levels,
    limits = c(0.5, length(stage_levels) + 0.5),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = "First day the stage was seen",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +   # clip = "off" lets the n labels sit outside
  theme_pub(8) +
  theme(
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 40, 8, 8, "pt")
  )

save_fig(p_stage, "36_stage_timing", width = 190, height = 110)

## 3e. Violin plot of first day per stage (figure 36b) -----------------------
# Violins show the spread of first days among fragments that reached each
# stage; boxplots and mean diamonds are added only for stages every fragment
# reached. Fragments not reached are crosses in the grey band at the top, with
# the % not reached printed above.
stage_axis_levels <- paste0("Stage ", seq_len(nrow(stage_map)))
stage_axis_labels <- c(
  "1\nCoenosarc\ncover",
  "2\nMain polyp\nforms",
  "3\nWound surface\nremodels",
  "4\nTip\nforms",
  "5\nTip\nextends",
  "6\nNew corallites\nbud"
)
dodge_width <- 0.72   # side-by-side spacing of the two temperatures

stage_events_violin <- stage_events |>
  mutate(
    stage_axis = factor(paste0("Stage ", stage_order),
                        levels = stage_axis_levels),
    event_status = factor(if_else(event, "Reached",
                                  "Not reached by Day 15"),
                          levels = c("Reached",
                                     "Not reached by Day 15"))
  )

stage_summary_violin <- stage_summary |>
  mutate(
    stage_axis = factor(paste0("Stage ", stage_order),
                        levels = stage_axis_levels),
    pct_not_reached_label = if_else(
      n_not_reached > 0,
      paste0(round(100 * n_not_reached / n_corals), "%"),
      NA_character_
    )
  )

# Jittered points use fixed seeds (36, 37) so the figure is reproducible
p_stage_violin <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 16.05, ymax = Inf,
           fill = "grey94", colour = NA) +
  geom_hline(yintercept = 16.05, linetype = "dotted",
             linewidth = 0.4, colour = "grey50") +
  geom_violin(
    data = filter(stage_events_violin, event),
    aes(x = stage_axis, y = first_yes_day, fill = treatment,
        group = interaction(stage_axis, treatment)),
    position = position_dodge(width = dodge_width),
    scale = "count",
    trim = TRUE,
    adjust = 0.85,
    alpha = 0.72,
    colour = "#202020",
    linewidth = 0.35,
    width = 0.68
  ) +
  geom_boxplot(
    data = filter(stage_events_violin, event, stage_order %in%
                    stage_summary$stage_order[stage_summary$complete_both_temperatures]),
    aes(x = stage_axis, y = first_yes_day,
        group = interaction(stage_axis, treatment)),
    position = position_dodge(width = dodge_width),
    width = 0.12,
    fill = NA,
    colour = "#202020",
    linewidth = 0.38,
    outlier.shape = NA
  ) +
  geom_point(
    data = filter(stage_summary_violin, complete_both_temperatures),
    aes(x = stage_axis, y = mean_first_day,
        group = treatment),
    position = position_dodge(width = dodge_width),
    inherit.aes = FALSE,
    shape = 23,
    size = 2,
    stroke = 0.35,
    colour = "#202020",
    fill = "white"
  ) +
  geom_point(
    data = filter(stage_events_violin, event),
    aes(x = stage_axis, y = first_yes_day, fill = treatment),
    position = position_jitterdodge(
      jitter.width = 0.07, jitter.height = 0,
      dodge.width = dodge_width, seed = 36
    ),
    shape = 21,
    colour = "white",
    size = 1.55,
    stroke = 0.25,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  geom_point(
    data = filter(stage_events_violin, !event),
    aes(x = stage_axis, y = 17.05, colour = treatment),
    position = position_jitterdodge(
      jitter.width = 0.07, jitter.height = 0.06,
      dodge.width = dodge_width, seed = 37
    ),
    shape = 4,
    size = 2.05,
    stroke = 0.85,
    show.legend = FALSE
  ) +
  geom_text(
    data = filter(stage_summary_violin, !is.na(pct_not_reached_label)),
    aes(x = stage_axis, y = 17.72, label = pct_not_reached_label,
        colour = treatment, group = treatment),
    position = position_dodge(width = dodge_width),
    inherit.aes = FALSE,
    size = 2.65,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = PAL_TEMP,
    labels = c(`28C` = "28 °C", `31C` = "31 °C"),
    name = "Temperature"
  ) +
  scale_colour_manual(
    values = PAL_TEMP,
    labels = c(`28C` = "28 °C", `31C` = "31 °C"),
    guide = "none"
  ) +
  scale_x_discrete(labels = stage_axis_labels, drop = FALSE) +
  scale_y_continuous(
    breaks = c(1, 3, 5, 7, 9, 11, 13, 15, 17.05),
    labels = c("1", "3", "5", "7", "9", "11", "13", "15",
               "Not\nreached"),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = "Healing or regeneration stage",
    y = "Day first seen"
  ) +
  coord_cartesian(ylim = c(0.5, 18.15), clip = "on") +
  theme_pub(8) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(lineheight = 0.9, margin = margin(t = 4)),
    plot.margin = margin(8, 12, 8, 8, "pt")
  )

save_fig(p_stage_violin, "36b_stage_violin",
         width = 190, height = 130)


# ---- 4. Figure 1B: cumulative percent reached ------------------------------
# For each day, the percent of fragments that had reached each stage by then.
# Every fragment stays in the denominator, including those that had not
# reached a stage by Day 15. Stage names and numbers follow the manuscript
# (healing, then regeneration stages 1-5). Healing (coenosarc cover) comes
# from the microscope experiment (8 per temperature) and appears only in the
# all-stage version (39d); stages 1-5 come from the main experiment (12
# wounded fragments per temperature). Those sample sizes go in the caption.

## 4a. Read the stage table and attach manuscript stage names ----------------
STAGE_LABELS <- c("Healing: coenosarc\ncovers wound", "1. Axial polyp\nforms",
                  "2. Wound surface\nremodels", "3. Axial corallite\nforms",
                  "4. Axial corallite\nextends", "5. Radial corallites\nform")

stage_events <- readr::read_csv(
  file.path(TBL_DIR, "36_stage_timing.csv"),
  show_col_types = FALSE
) |>
  mutate(
    treatment = factor(treatment, levels = c("28C", "31C")),
    stage_order = as.integer(stage_order),
    stage_short = factor(
      dplyr::recode(
        as.character(stage),
        `Coenosarc barrier re-established` = STAGE_LABELS[1],
        `Axial polyp forms` = STAGE_LABELS[2],
        `Wound surface remodels` = STAGE_LABELS[3],
        `Axial corallite or tip forms` = STAGE_LABELS[4],
        `Axial corallite or tip extends` = STAGE_LABELS[5],
        `New radial corallites bud` = STAGE_LABELS[6]
      ),
      levels = STAGE_LABELS
    ),
    event_status = factor(if_else(event, "Reached", "Not reached by Day 15"),
                          levels = c("Reached", "Not reached by Day 15"))
  )

## 4b. Cumulative percent reached on each scoring day ------------------------
score_days <- 0:15

# Guarantee: one row per stage x fragment, and every fragment that never
# reached a stage was still being scored on Day 15 (so it is truly "not
# reached", not lost early).
stopifnot(!anyDuplicated(select(stage_events, stage_order, id)),
          all(stage_events$event | stage_events$last_scored_day == 15))

# Build every stage x temperature x day combination, attach all fragments to
# each day (relationship = "many-to-many" says the row duplication is
# intended), and count those whose first day is on or before that day.
cumulative <- tidyr::expand_grid(
  stage_order = sort(unique(stage_events$stage_order)),
  treatment = levels(stage_events$treatment),
  day = score_days
) |>
  left_join(distinct(stage_events, stage_order, stage_short),
            by = "stage_order") |>
  left_join(
    stage_events |>
      select(stage_order, treatment, id, event, first_yes_day),
    by = c("stage_order", "treatment"),
    relationship = "many-to-many"
  ) |>
  group_by(stage_order, stage_short, treatment, day) |>
  summarise(
    n_scored = n(),
    n_reached = sum(event & first_yes_day <= day, na.rm = TRUE),
    pct_reached = 100 * n_reached / n_scored,
    .groups = "drop"
  ) |>
  mutate(treatment = factor(treatment, levels = c("28C", "31C")))

write_csv(cumulative, file.path(TBL_DIR, "39_stage_cumulative_percent.csv"))

# Guarantee: never more fragments reached than were scored
stopifnot(all(cumulative$n_reached <= cumulative$n_scored))

## 4c. Plot the cumulative curves --------------------------------------------
# Step lines (the percent jumps on scoring days). Colour, line type and point
# shape all mark temperature, so the figure still reads in greyscale.
plot_cumulative <- function(dat, ncol) {
  ggplot(dat, aes(day, pct_reached, colour = treatment,
                  linetype = treatment, shape = treatment)) +
    geom_step(linewidth = 0.75, direction = "hv") +
    geom_point(size = 1.5) +
    facet_wrap(~ stage_short, ncol = ncol) +
    scale_colour_manual(values = PAL_TEMP, name = "Temperature",
                        labels = c(`28C` = "28 °C", `31C` = "31 °C")) +
    scale_linetype_manual(values = c(`28C` = "solid", `31C` = "dashed"),
                          name = "Temperature",
                          labels = c(`28C` = "28 °C", `31C` = "31 °C")) +
    scale_shape_manual(values = c(`28C` = 16, `31C` = 17),
                       name = "Temperature",
                       labels = c(`28C` = "28 °C", `31C` = "31 °C")) +
    scale_x_continuous(breaks = c(0, 1, 3, 5, 10, 15),
                       limits = c(0, 15)) +
    scale_y_continuous(labels = \(x) paste0(x, "%"),
                       limits = c(0, 100),
                       breaks = c(0, 25, 50, 75, 100)) +
    labs(x = "Days after wounding", y = "Fragments that had reached the stage (%)") +
    theme_pub(10)
}

# Figure 1B in the working draft: the four regeneration milestones (axial
# polyp, axial corallite forms, axial corallite extends, radial corallites),
# all from the main experiment. Stage numbers match Supplementary figure 3.
FIG1B_STAGES <- STAGE_LABELS[c(2, 4, 5, 6)]
p_cumulative <- plot_cumulative(filter(cumulative, stage_short %in% FIG1B_STAGES),
                                ncol = 2)
save_fig(p_cumulative, "39b_stage_timing_cumulative_option",
         width = 150, height = 135)

# All six stages, including healing and wound-surface remodelling (supplement).
p_cumulative_all <- plot_cumulative(cumulative, ncol = 3)
save_fig(p_cumulative_all, "39d_stage_cumulative_all_stages",
         width = 200, height = 145)


# ---- 5. Numbers for the results text ---------------------------------------
# Collects every number quoted in the results text into one table
# (output/paper/regeneration_numbers.csv). paper_row() (from 00_setup.R) makes
# one row; rows are appended to the list `rows` and written at the end. Most
# numbers are re-read from the tables written above so the text always matches
# them. In the labels, "estimate" and "statistic" are reused as two slots
# (e.g. ambient count, heated count); each label says what they hold.

## 5a. Re-read the result tables ---------------------------------------------
endpoint <- readr::read_csv(
  file.path(TBL_DIR, "11f_morphology_dataset_endpoint_tests.csv"),
  show_col_types = FALSE
)
logrank  <- readr::read_csv(
  file.path(TBL_DIR, "11f_morphology_dataset_logrank_tests.csv"),
  show_col_types = FALSE
)
aft      <- readr::read_csv(file.path(TBL_DIR, "14_interval_survreg.csv"),
                            show_col_types = FALSE)
cox      <- readr::read_csv(file.path(TBL_DIR, "14_cox_hazard_ratios.csv"),
                            show_col_types = FALSE)
cox_gx   <- readr::read_csv(file.path(TBL_DIR, "14_cox_genet_LRT.csv"),
                            show_col_types = FALSE)
stages   <- readr::read_csv(file.path(TBL_DIR, "36_stage_timing_mean_iqr.csv"),
                            show_col_types = FALSE)

# Plain-language names for the five main-experiment stages (replaces the
# section 2 `traits` vector)
traits <- c(axial_polyp_formation = "Axial polyp forms",
            wound_smoothed = "Cut surface smooths over",
            tip_exist = "New tip forms",
            tip_extension = "Tip extends",
            new_corallites_on_tip = "New radial corallites bud")
rows <- list()

## 5b. Day 15 counts, stages reached and log-rank ----------------------------
# Stage present on the last scored day (Fisher's exact test), main experiment.
# stopifnot(): stop if there is not exactly one matching row.
for (tr in names(traits)) {
  d <- endpoint[endpoint$dataset == "physio_morphology" & endpoint$trait == tr, ]
  stopifnot(nrow(d) == 1)
  rows[[length(rows) + 1]] <- paper_row(
    paste0(tr, "_present_day15"),
    paste0(traits[[tr]], ": fragments with the stage on Day 15 (estimate = ambient, statistic = heated; n = per temperature)"),
    estimate = d$yes_28, statistic = d$yes_31, n = d$n_28,
    test = "Fisher's exact test", p_value = d$fisher_p,
    source = "output/tables/11f_morphology_dataset_endpoint_tests.csv")
}

# Stage reached at any point by Day 15, main experiment (stage_order 2-6)
stage_trait <- c(`2` = "axial_polyp_formation", `3` = "wound_smoothed",
                 `4` = "tip_exist",
                 `5` = "tip_extension", `6` = "new_corallites_on_tip")
for (so in names(stage_trait)) {
  d28 <- stages[stages$stage_order == as.integer(so) & grepl("^28", stages$treatment), ]
  d31 <- stages[stages$stage_order == as.integer(so) & grepl("^31", stages$treatment), ]
  if (nrow(d28) != 1 || nrow(d31) != 1) next
  tr <- stage_trait[[so]]
  rows[[length(rows) + 1]] <- paper_row(
    paste0(tr, "_reached_by_day15"),
    paste0(traits[[tr]], ": fragments that reached the stage by Day 15 (estimate = ambient, statistic = heated)"),
    estimate = d28$n_reached, statistic = d31$n_reached, n = d28$n_scored,
    source = "output/tables/36_stage_timing_mean_iqr.csv")
}

# Log-rank test for budding, main experiment
d <- logrank[logrank$dataset == "physio_morphology" &
               logrank$trait == "new_corallites_on_tip", ]
rows[[length(rows) + 1]] <- paper_row(
  "budding_logrank", "Budding: log-rank test", test = "log-rank",
  statistic = d$logrank_chisq, df1 = d$logrank_df, p_value = d$logrank_p,
  source = "output/tables/11f_morphology_dataset_logrank_tests.csv")

## 5c. Timing models: time ratios, hazard ratio, patch x temperature ---------
# Time to each stage (interval-censored log-logistic; time ratio > 1 = slower
# under heat)
for (tr in names(traits)) {
  d <- aft[aft$trait == tr, ]
  if (nrow(d) != 1) next
  rows[[length(rows) + 1]] <- paper_row(
    paste0(tr, "_time_ratio"), paste0(traits[[tr]], ": time ratio, heated / ambient"),
    estimate = d$time_ratio_31_vs28, ci_low = d$ratio_lo, ci_high = d$ratio_hi,
    test = paste0(sub(" AFT$", "", d$model), " model"),
    statistic = d$z, p_value = d$p, n = d$n,
    source = "output/tables/14_interval_survreg.csv")
}

# Budding hazard ratio from the overall (patch-stratified) Cox model
d <- cox[cox$trait == "new_corallites_on_tip" & grepl("^overall", cox$scope), ]
stopifnot(nrow(d) == 1)
rows[[length(rows) + 1]] <- paper_row(
  "budding_hazard_ratio", "Budding: hazard ratio, heated / ambient (Cox)",
  estimate = d$HR_31_vs28, ci_low = d$HR_lo, ci_high = d$HR_hi,
  test = "Cox model, stratified by source patch",
  statistic = d$z, p_value = d$p, n = d$n,
  source = "output/tables/14_cox_hazard_ratios.csv")

# Budding source patch x temperature LRT (unstable; see the label)
d <- cox_gx[cox_gx$trait == "new_corallites_on_tip", ]
rows[[length(rows) + 1]] <- paper_row(
  "budding_patch_x_temp", "Budding: source patch x temperature (unstable: no heated patch A fragment budded, so the model cannot estimate it; report the counts below)",
  test = "likelihood-ratio test",
  statistic = d$chisq_genet_trt_int, df1 = d$df_genet_trt_int,
  p_value = d$p_genet_trt_int, source = "output/tables/14_cox_genet_LRT.csv")

## 5d. Budding tank by tank: exact permutation test --------------------------
# Temperature was applied to tanks (4 per temperature), so tanks, not
# fragments, are the true replicates. Compare the proportion of fragments that
# budded per tank, 4 vs 4 tanks. Exact permutation test: try all 70 ways
# (choose(8, 4)) of labelling 4 of the 8 tanks "ambient", and ask how often
# the difference is at least as large as the real one. The 1e-12 guards
# against rounding error in the >= comparison.
ph <- readRDS(file.path(DATA_PROC, "physio_clean.rds"))
by_tank <- ph |>
  dplyr::filter(wound == "yes") |>
  dplyr::group_by(treatment, tank, id) |>
  dplyr::summarise(budded = any(new_corallites_on_tip %in% c("yes", 1, TRUE)),
                   .groups = "drop") |>
  dplyr::group_by(treatment, tank) |>
  dplyr::summarise(n = dplyr::n(), k = sum(budded), prop = k / n, .groups = "drop")

props <- setNames(by_tank$prop, by_tank$tank)
amb <- as.character(by_tank$tank[by_tank$treatment == "28C"])
diff_for <- function(set) mean(props[set]) - mean(props[setdiff(names(props), set)])
all_sets <- utils::combn(names(props), length(amb), simplify = FALSE)
obs <- diff_for(amb)
perm_p <- mean(vapply(all_sets, function(s) abs(diff_for(s)) >= abs(obs) - 1e-12,
                      logical(1)))

# Range of budded fragments per tank, for each temperature
for (tr in c("28C", "31C")) {
  lab <- if (tr == "28C") "ambient" else "heated"
  d <- by_tank[by_tank$treatment == tr, ]
  rows[[length(rows) + 1]] <- paper_row(
    paste0("budding_per_tank_", lab),
    paste0("Budding by Day 15, ", lab, " tanks: wounded fragments that budded per tank (estimate = fewest, statistic = most)"),
    estimate = min(d$k), statistic = max(d$k), n = d$n[1], df1 = nrow(d),
    source = "data/processed/physio_clean.rds")
}

rows[[length(rows) + 1]] <- paper_row(
  "budding_tank_permutation", "Budding: exact permutation test of tank proportions, 4 vs 4 tanks",
  estimate = obs, test = "exact permutation (70 assignments); smallest possible p = 2/70",
  p_value = perm_p, n = length(props), source = "data/processed/physio_clean.rds")

## 5e. Budding by source patch -----------------------------------------------
# Counts behind the patch x temperature pattern
by_patch <- ph |>
  dplyr::filter(wound == "yes") |>
  dplyr::group_by(treatment, thicket, id) |>
  dplyr::summarise(budded = any(new_corallites_on_tip %in% c("yes", 1, TRUE)),
                   .groups = "drop") |>
  dplyr::group_by(treatment, thicket) |>
  dplyr::summarise(n = dplyr::n(), k = sum(budded), .groups = "drop")

for (i in seq_len(nrow(by_patch))) {
  lab <- if (by_patch$treatment[i] == "28C") "ambient" else "heated"
  rows[[length(rows) + 1]] <- paper_row(
    paste0("budding_", lab, "_patch_", toupper(by_patch$thicket[i])),
    paste0("Budding by Day 15, ", lab, ", source patch ", toupper(by_patch$thicket[i]),
           ": wounded fragments that budded"),
    estimate = by_patch$k[i], n = by_patch$n[i],
    source = "data/processed/physio_clean.rds")
}

## 5f. Numbers for the Figure 1 text -----------------------------------------
# Share reached by a given day, and the range of first days among fragments
# that reached the stage. Each spec is (stage_order, day, key); the label
# prints stage_order - 1 to match the manuscript stage numbers.
fig1 <- read_csv(file.path(TBL_DIR, "39_stage_cumulative_percent.csv"),
                 show_col_types = FALSE)
pct_at <- function(stage, d, trt) fig1 |>
  dplyr::filter(stage_order == stage, day == d, treatment == trt)
for (spec in list(list(5, 8, "tip_extension_pct_day8"),
                  list(6, 12, "budding_pct_day12"))) {
  for (trt in c("28C", "31C")) {
    r <- pct_at(spec[[1]], spec[[2]], trt)
    lab <- if (trt == "28C") "ambient" else "heated"
    rows[[length(rows) + 1]] <- paper_row(
      paste0(spec[[3]], "_", lab),
      sprintf("Wounded fragments that had reached stage %d by Day %d, %s (%%)",
              spec[[1]] - 1, spec[[2]], lab),
      estimate = r$pct_reached, statistic = r$n_reached, n = r$n_scored,
      source = "output/tables/39_stage_cumulative_percent.csv")
  }
}

# Earliest and latest first day for tip extension (5) and budding (6)
st <- read_csv(file.path(TBL_DIR, "36_stage_timing.csv"), show_col_types = FALSE)
for (s in c(5, 6)) for (trt in c("28C", "31C")) {
  r <- st |> dplyr::filter(stage_order == s, treatment == trt, event)
  key <- c(`5` = "tip_extension", `6` = "budding")[as.character(s)]
  lab <- if (trt == "28C") "ambient" else "heated"
  rows[[length(rows) + 1]] <- paper_row(
    paste0(key, "_first_day_range_", lab),
    paste0("First day seen, among fragments that reached it, ", lab, " (estimate = earliest, statistic = latest)"),
    estimate = min(r$first_yes_day), statistic = max(r$first_yes_day), n = nrow(r),
    source = "output/tables/36_stage_timing.csv")
}

## 5g. Was budding delayed or stalled under heat? (open question) ------------
# At 28 C every fragment budded within a few days of its tip extending. Count
# heated fragments whose tip had been extending for longer than the slowest
# ambient fragment needed, and that still had not budded by Day 15.
# pivot_wider() gives one row per fragment with ext/bud columns side by side.
ext_bud <- st |>
  dplyr::filter(stage_order %in% c(5, 6)) |>
  dplyr::mutate(s = ifelse(stage_order == 5, "ext", "bud")) |>
  dplyr::select(id, treatment, s, event, first_yes_day) |>
  tidyr::pivot_wider(names_from = s, values_from = c(event, first_yes_day))

# Days from tip extension to budding, ambient fragments that did both
lag_amb <- ext_bud |>
  dplyr::filter(treatment == "28C", event_ext, event_bud) |>
  dplyr::mutate(lag = first_yes_day_bud - first_yes_day_ext) |>
  dplyr::pull(lag)

stalled <- ext_bud |>
  dplyr::filter(treatment == "31C", event_ext, !event_bud,
                15 - first_yes_day_ext > max(lag_amb))

rows[[length(rows) + 1]] <- paper_row(
  "bud_lag_ambient", "28 C: days from tip extension to budding (estimate = shortest, statistic = longest)",
  estimate = min(lag_amb), statistic = max(lag_amb), n = length(lag_amb),
  source = "output/tables/36_stage_timing.csv")
rows[[length(rows) + 1]] <- paper_row(
  "bud_stalled_heated", "31 C: fragments with an extending tip for longer than any 28 C fragment needed, but no budding by Day 15",
  estimate = nrow(stalled), n = sum(ext_bud$treatment == "31C"),
  source = "output/tables/36_stage_timing.csv")

## 5h. Microscope experiment: tissue cover and budding -----------------------
# Tissue cover: first day living tissue covered the cut, per fragment
mic <- readRDS(file.path(DATA_PROC, "microscope_physio_clean.rds"))
first_cover <- mic |>
  dplyr::filter(tissue_over_wound %in% c("yes", 1, TRUE)) |>
  dplyr::group_by(id) |>
  dplyr::summarise(first_day = min(day), .groups = "drop")

rows[[length(rows) + 1]] <- paper_row(
  "tissue_cover_by_day1", "Microscope experiment: fragments with tissue covering the cut by Day 1",
  estimate = sum(first_cover$first_day <= 1), n = dplyr::n_distinct(mic$id),
  source = "data/processed/microscope_physio_clean.rds")

# Budding in the microscope experiment, a second cohort (open question)
for (trt in c("28C", "31C")) {
  b <- mic |>
    dplyr::filter(treatment == trt, wound == "yes") |>
    dplyr::group_by(id) |>
    dplyr::summarise(bud = any(new_corallites_on_tip %in% c("yes", 1, TRUE)),
                     .groups = "drop")
  rows[[length(rows) + 1]] <- paper_row(
    paste0("micro_budding_day15_", if (trt == "28C") "ambient" else "heated"),
    paste("Microscope experiment: wounded fragments that budded by Day 15,", trt),
    estimate = sum(b$bud), n = nrow(b),
    source = "data/processed/microscope_physio_clean.rds")
}

# Write all rows to output/paper/regeneration_numbers.csv
write_paper(rows, "regeneration")
