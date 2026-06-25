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
# Input:   data/processed/physio_clean.rds
# Output:  output/tables/14_km_event_summary.csv         — per-trait/genet/treatment
#          output/tables/14_cox_hazard_ratios.csv        — HR for 31 vs 28 (first-observed approximation)
#          output/tables/14_interval_survreg.csv         — interval-censored Weibull AFT tests
#          output/tables/14_cox_genet_LRT.csv            — does genet × treatment improve fit?
#          figures/14_morphology_KM.{pdf,png}            — KM by treatment
#          figures/14b_morphology_KM_by_genet.{pdf,png}  — KM by treatment × genet
# =============================================================================

source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages({
  library(survival)
  library(survminer)
})

ph <- readRDS(file.path(DATA_PROC, "physio_clean.rds")) |>
  mutate(
    treatment = factor(treatment, levels = c("28C", "31C")),
    wound = factor(wound, levels = c("no", "yes")),
    thicket = factor(thicket)
  )
contrasts(ph$treatment) <- contr.treatment(nlevels(ph$treatment))

traits <- c("hole_in_center", "polyp_in_hole", "wound_smoothed",
            "pigment_over_wound", "tip_exist", "tip_extension",
            "new_corallites_on_tip")

trait_interpretation <- tibble(
  trait = traits,
  event_interpretation = c(
    "first observed hole closure",
    "first observed polyp within wound",
    "first observed smoothed wound surface",
    "first observed pigment over wound; expression can be non-monotonic",
    "first observed visible tip; expression can be non-monotonic",
    "first observed tip extension; expression can be non-monotonic",
    "first observed new corallites on tip"
  )
)

# For each (coral, trait) compute first observed event day plus interval bounds:
#   event_lower = previous observed day before first 1
#   event_upper = first observed day with value 1
# Right-censor at last observation if never observed.
compute_events <- function(d, trait) {
  d |>
    filter(wound == "yes", !is.na(day), day >= 0) |>
    mutate(y = .data[[trait]]) |>
    group_by(id, treatment, tank, thicket) |>
    arrange(day, .by_group = TRUE) |>
    summarise(
      event_day = {
        first1 <- which(y == 1)[1]
        if (is.na(first1)) max(day, na.rm = TRUE) else day[first1]
      },
      event_lower = {
        first1 <- which(y == 1)[1]
        if (is.na(first1)) max(day, na.rm = TRUE) else if (first1 == 1) 0 else day[first1 - 1]
      },
      event_upper = {
        first1 <- which(y == 1)[1]
        if (is.na(first1)) Inf else day[first1]
      },
      event = as.integer(any(y == 1, na.rm = TRUE)),
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
contrasts(events$treatment) <- contr.treatment(nlevels(events$treatment))

# ---- KM summary tables ----------------------------------------------------
# Per trait × genet × treatment median time-to-event
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

# ---- Interval-censored survival models ------------------------------------
# These are the inferential survival tests for discretely scored morphology.
# Cox/KM outputs below are retained as first-observed-day summaries/figures.
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
interval_survreg <- map_dfr(traits, fit_interval_survreg)
write_csv(interval_survreg, file.path(TBL_DIR, "14_interval_survreg.csv"))

# ---- Inter-milestone lag: wound closure -> regeneration -------------------
# The headline result is "heat impairs regeneration, not closure." Quantify it
# directly per coral: the LAG (days) from achieving wound closure
# (`wound_smoothed`) to forming new skeleton at the tip
# (`new_corallites_on_tip`). Corals that close but never regenerate within the
# experiment are right-censored (no lag) and counted separately — heat is
# expected to inflate that censored fraction.
ev_wide <- events |>
  select(id, treatment, thicket, trait, event, event_day) |>
  pivot_wider(names_from = trait, values_from = c(event, event_day))

lag_pairs <- list(
  c(from = "wound_smoothed", to = "new_corallites_on_tip"),  # closure -> regeneration (primary)
  c(from = "wound_smoothed", to = "tip_extension"),          # closure -> tip extension
  c(from = "hole_in_center", to = "new_corallites_on_tip")   # initial closure -> regeneration
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

# Wilcoxon test of the primary closure->regeneration lag (28 vs 31) among
# corals that reached both, where estimable.
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

cox_overall   <- map_dfr(traits, fit_cox_overall)
cox_per_genet <- map_dfr(traits, fit_cox_per_genet)
cox_results   <- bind_rows(cox_overall, cox_per_genet) |>
  mutate(time_scale = "first-observed-day approximation")
write_csv(cox_results, file.path(TBL_DIR, "14_cox_hazard_ratios.csv"))

# ---- Proportional-hazards diagnostics for EVERY overall Cox model ----------
# cox.zph (Schoenfeld residual test + plot) for all traits, not just the one
# that violated PH. Documents the assumption visually for each model.
DIAG_DIR <- file.path(FIG_DIR, "diagnostics")
dir.create(DIAG_DIR, recursive = TRUE, showWarnings = FALSE)
ph_rows <- list()
for (tr in traits) {
  d <- events |> filter(trait == tr)
  if (sum(d$event) < 5) next
  fit <- tryCatch(coxph(Surv(event_day, event) ~ treatment + strata(thicket),
                        data = d), error = function(e) NULL)
  if (is.null(fit)) next
  zph <- tryCatch(cox.zph(fit), error = function(e) NULL)
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
  lrt_g  <- anova(m_base, m_add)
  lrt_gi <- anova(m_add, m_inter)
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
# Build per-trait, per-treatment survival curves manually for ggplot
km_curves <- events |>
  group_by(trait, treatment) |>
  group_modify(\(d, k) {
    d <- dplyr::arrange(d, event_day)
    times <- sort(unique(d$event_day))
    # Build a step at every distinct day, NA-safe
    n_at_risk <- sapply(times, \(t) sum(d$event_day >= t))
    n_events  <- sapply(times, \(t) sum(d$event_day == t & d$event == 1))
    haz       <- n_events / pmax(n_at_risk, 1)
    surv      <- cumprod(1 - haz)
    out_day   <- c(0, times)
    out_surv  <- c(1, surv)
    tibble(day = out_day,
           surv = out_surv,
           cum_event = 1 - out_surv)
  }) |>
  ungroup() |>
  mutate(trait = factor(trait,
                        levels = traits,
                        labels = c("Hole in center", "Polyp in hole",
                                   "Wound smoothed", "Pigment over wound",
                                   "Tip exists", "Tip extension",
                                   "New corallites on tip")))

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
       subtitle = "Kaplan-Meier curves for each wound-healing milestone (n = 24 corals × treatment)") +
  theme_pub(9)

save_fig(p_km, "14_morphology_KM", width = 200, height = 130)

# ---- Genet-resolved KM figure --------------------------------------------
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
                        levels = traits,
                        labels = c("Hole in center", "Polyp in hole",
                                   "Wound smoothed", "Pigment over wound",
                                   "Tip exists", "Tip extension",
                                   "New corallites on tip")))

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
       subtitle = "Kaplan-Meier curves by treatment × genet (n ≈ 8 per cell)") +
  theme_pub(9)

save_fig(p_km_genet, "14b_morphology_KM_by_genet", width = 210, height = 135)

# ---- Time-varying-coefficient refit for the one PH violation -------------
# Cox diagnostics flagged pigment_over_wound in genet c (cox.zph chisq = 6.59,
# p = 0.0103, n_event = 5) as violating PH. Refit with a log(t+1)
# time-varying coefficient on treatment so the HR is allowed to evolve.
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
