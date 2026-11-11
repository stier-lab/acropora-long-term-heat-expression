# =============================================================================
# 04_regeneration.R: healing and regeneration stages (Figure 1B)
#
# Only wounded fragments are used. Main experiment: 12 per temperature, scored
# daily Days 0-15. Microscope experiment: 8 per temperature, analysed
# separately (only it recorded living tissue covering the cut).
#
#   1 Day 15 counts: Fisher's exact test per stage; log-rank test on the first
#     day each stage was seen. Each experiment tested on its own.
#   2 Timing models (main experiment): interval-censored Weibull models (the
#     stage began between the last "no" and the first "yes"), Cox models
#     stratified by source patch, and a likelihood-ratio test for source patch
#     x temperature. Kaplan-Meier curves as a supporting figure.
#   3 Stage table: first day each fragment reached each of six stages; counts
#     and first-day summaries (averages only where every fragment reached it).
#   4 Figure 1B: cumulative percent of fragments that reached each stage by day.
#   5 Numbers for the results text, including the tank-level budding check.
#
# Writes: output/tables/11f_morphology_dataset_{endpoint,logrank}_tests.csv,
#         14_event_intervals.csv, 14_interval_survreg.csv, 14_cox_hazard_ratios.csv,
#         14_cox_genet_LRT.csv, 36_stage_timing.csv, 36_stage_timing_mean_iqr.csv,
#         39_stage_cumulative_percent.csv
#         figures/39b_stage_timing_cumulative_option (Figure 1B), 36_stage_timing,
#         36b_stage_violin, 14_morphology_KM
#         output/paper/regeneration_numbers.csv
# =============================================================================

source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages(library(survival))

# ---- 1. Day 15 counts and log-rank tests (both experiments) ----------------
ph <- readRDS(file.path(DATA_PROC, "physio_clean.rds"))
micro <- readRDS(file.path(DATA_PROC, "microscope_physio_clean.rds"))
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
    facet_trait = paste(dataset_short, short_label, sep = ": "),
    facet_event = paste(dataset_short, phase_short, sep = "\n")
  )
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
    test_scope = "paired with figures 11c-11e; cohorts tested separately"
  ) |>
  arrange(dataset, trait_order)
write_csv(endpoint_tests,
          file.path(TBL_DIR, "11f_morphology_dataset_endpoint_tests.csv"))
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
    if (n_28 == 0 || n_31 == 0 || sum(d$event_int) == 0) {
      return(tibble(
        n_28 = n_28, n_31 = n_31,
        n_event_28 = n_event_28, n_event_31 = n_event_31,
        logrank_chisq = NA_real_, logrank_df = NA_real_,
        logrank_p = NA_real_,
        test_status = "not_estimable"
      ))
    }
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
      "paired with figures 11c-11e; cohorts tested separately;",
      "main-cohort primary inference remains script 14 interval-censored AFT"
    )
  ) |>
  arrange(dataset, trait_order)
write_csv(logrank_tests,
          file.path(TBL_DIR, "11f_morphology_dataset_logrank_tests.csv"))

# ---- 2. Timing models (main experiment) ------------------------------------
# A stage is "reached" on the first day it was scored yes. Fragments that never
# reached it by Day 15 stay in as censored (not reached yet), not dropped.
ph <- readRDS(file.path(DATA_PROC, "physio_clean.rds")) |>
  mutate(
    treatment = factor(treatment, levels = c("28C", "31C")),
    wound = factor(wound, levels = c("no", "yes")),
    thicket = factor(thicket)
  )
contrasts(ph$treatment) <- contr.treatment(nlevels(ph$treatment))
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
        if (is.na(first1)) max(day, na.rm = TRUE) else if (first1 == 1) 0 else day[first1 - 1]
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
write_csv(events, file.path(TBL_DIR, "14_event_intervals.csv"))
contrasts(events$treatment) <- contr.treatment(nlevels(events$treatment))
n_km_per_trt  <- events |> distinct(id, treatment) |> count(treatment) |>
  pull(n) |> (\(x) if (length(unique(x)) == 1) unique(x) else paste(range(x), collapse = "–"))()
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
       y = "Cumulative % corals expressing trait",
       title = "Heat delays the regenerative-tip program but not wound closure",
       subtitle = sprintf("Kaplan-Meier curves for each wound-healing milestone (n = %s corals per treatment)",
                          n_km_per_trt)) +
  theme_pub(9)
save_fig(p_km, "14_morphology_KM", width = 200, height = 130)

# ---- 3. Stage table and stage figures --------------------------------------
# Stage 1 (tissue cover) comes from the microscope experiment; stages 2-6 from
# the main experiment.
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
                ~if_else(complete_both_temperatures, .x, NA_real_))) |>
  ungroup() |>
  arrange(stage_order, treatment)
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
status_key <- tibble(
  x = NA_real_,
  y = NA_real_,
  event_status = factor(c("Reached", "Not reached by final day"),
                        levels = c("Reached", "Not reached by final day"))
)
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
    y = NULL,
    title = "Heat mostly affects the late regeneration stage",
    subtitle = "Small stacks separate corals that reached the same stage on the same day. X marks corals that never reached the stage.",
    caption = "Each point is one wounded coral. Right labels show reached/total corals.\nThe healing row uses microscope photos; regeneration rows use the main morphology data set."
  ) +
  coord_cartesian(clip = "off") +
  theme_pub(8) +
  theme(
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(8, 40, 8, 8, "pt")
  )
save_fig(p_stage, "36_stage_timing", width = 190, height = 110)
stage_axis_levels <- paste0("Stage ", seq_len(nrow(stage_map)))
stage_axis_labels <- c(
  "1\nCoenosarc\ncover",
  "2\nMain polyp\nforms",
  "3\nWound surface\nremodels",
  "4\nTip\nforms",
  "5\nTip\nextends",
  "6\nNew corallites\nbud"
)
dodge_width <- 0.72
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
    y = "Day first seen",
    title = "Heat reduced how many corals reached the final regrowth stage",
    subtitle = "Violin width scales with the number of corals that reached each stage; X marks corals still short of a stage by Day 15.",
    caption = "Boxes and means appear only where all fragments at both temperatures reached the stage.\nPoints and violin shapes for incomplete stages describe only fragments that reached them; X marks retain the others."
  ) +
  coord_cartesian(ylim = c(0.5, 18.15), clip = "on") +
  theme_pub(8) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(lineheight = 0.9, margin = margin(t = 4)),
    plot.margin = margin(8, 12, 8, 8, "pt")
  )
save_fig(p_stage_violin, "36b_stage_violin",
         width = 190, height = 130)

# ---- 4. Figure 1B: cumulative percent reached ------------------------------
# Every fragment stays in the denominator, including those that had not
# reached a stage by Day 15.
stage_events <- readr::read_csv(
  file.path(TBL_DIR, "36_stage_timing.csv"),
  show_col_types = FALSE
) |>
  mutate(
    treatment = factor(treatment, levels = c("28C", "31C")),
    treatment_label = factor(
      dplyr::recode(as.character(treatment),
                    `28C` = "28 °C", `31C` = "31 °C"),
      levels = c("28 °C", "31 °C")
    ),
    stage_order = as.integer(stage_order),
    stage_short = factor(
      dplyr::recode(
        as.character(stage),
        `Coenosarc barrier re-established` = "Living tissue\ncovers wound",
        `Axial polyp forms` = "Main polyp\nforms",
        `Wound surface remodels` = "Wound surface\nremodels",
        `Axial corallite or tip forms` = "Tip\nforms",
        `Axial corallite or tip extends` = "Tip\nextends",
        `New radial corallites bud` = "New skeletal\ncups bud"
      ),
      levels = c("Living tissue\ncovers wound", "Main polyp\nforms",
                 "Wound surface\nremodels", "Tip\nforms",
                 "Tip\nextends", "New skeletal\ncups bud")
    ),
    event_status = factor(if_else(event, "Reached", "Not reached by Day 15"),
                          levels = c("Reached", "Not reached by Day 15"))
  )
score_days <- 0:15
stopifnot(!anyDuplicated(select(stage_events, stage_order, id)),
          all(stage_events$event | stage_events$last_scored_day == 15))
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
stopifnot(all(cumulative$n_reached <= cumulative$n_scored))
p_cumulative <- ggplot(cumulative, aes(day, pct_reached, colour = treatment,
                                      linetype = treatment, shape = treatment)) +
  geom_step(linewidth = 0.75, direction = "hv") +
  geom_point(size = 1.5) +
  facet_wrap(~ stage_short, ncol = 3) +
  scale_colour_manual(values = PAL_TEMP, name = "Temperature",
                      labels = c(`28C` = "28 °C", `31C` = "31 °C")) +
  scale_linetype_manual(values = c(`28C` = "solid", `31C` = "dashed"),
                        name = "Temperature", labels = c(`28C` = "28 °C", `31C` = "31 °C")) +
  scale_shape_manual(values = c(`28C` = 16, `31C` = 17),
                     name = "Temperature", labels = c(`28C` = "28 °C", `31C` = "31 °C")) +
  scale_x_continuous(breaks = c(0, 1, 3, 5, 10, 15),
                     limits = c(0, 15)) +
  scale_y_continuous(labels = \(x) paste0(x, "%"),
                     limits = c(0, 100),
                     breaks = c(0, 25, 50, 75, 100)) +
  labs(
    title = "Fewer heated fragments reached the final regrowth step by Day 15",
    subtitle = "All fragments stay in the denominator, including those that had not reached a stage",
    x = "Day",
    y = "Fragments first seen at this stage by each day",
    caption = "Tissue cover: 8 fragments per temperature in a separate photo experiment, all covered on Day 1.\nOther panels: 12 fragments per temperature in the main experiment."
  ) +
  theme_pub(10) +
  theme(legend.position = "bottom")
save_fig(p_cumulative, "39b_stage_timing_cumulative_option",
         width = 200, height = 145)

# ---- 5. Numbers for the results text ---------------------------------------
endpoint <- readr::read_csv(file.path(TBL_DIR, "11f_morphology_dataset_endpoint_tests.csv"), show_col_types = FALSE)
logrank  <- readr::read_csv(file.path(TBL_DIR, "11f_morphology_dataset_logrank_tests.csv"), show_col_types = FALSE)
aft      <- readr::read_csv(file.path(TBL_DIR, "14_interval_survreg.csv"), show_col_types = FALSE)
cox      <- readr::read_csv(file.path(TBL_DIR, "14_cox_hazard_ratios.csv"), show_col_types = FALSE)
cox_gx   <- readr::read_csv(file.path(TBL_DIR, "14_cox_genet_LRT.csv"), show_col_types = FALSE)
stages   <- readr::read_csv(file.path(TBL_DIR, "36_stage_timing_mean_iqr.csv"), show_col_types = FALSE)

traits <- c(axial_polyp_formation = "Axial polyp forms",
            wound_smoothed = "Cut surface smooths over",
            tip_exist = "New tip forms",
            tip_extension = "Tip extends",
            new_corallites_on_tip = "New radial corallites bud")
rows <- list()

# Stage present on the last scored day (Fisher's exact test), main experiment
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

# Stage reached at any point by Day 15, main experiment
stage_trait <- c(`2` = "axial_polyp_formation", `3` = "wound_smoothed", `4` = "tip_exist",
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

d <- logrank[logrank$dataset == "physio_morphology" & logrank$trait == "new_corallites_on_tip", ]
rows[[length(rows) + 1]] <- paper_row(
  "budding_logrank", "Budding: log-rank test", test = "log-rank",
  statistic = d$logrank_chisq, df1 = d$logrank_df, p_value = d$logrank_p,
  source = "output/tables/11f_morphology_dataset_logrank_tests.csv")

# Time to each stage (interval-censored Weibull; time ratio > 1 = slower under heat)
for (tr in names(traits)) {
  d <- aft[aft$trait == tr, ]
  if (nrow(d) != 1) next
  rows[[length(rows) + 1]] <- paper_row(
    paste0(tr, "_time_ratio"), paste0(traits[[tr]], ": time ratio, heated / ambient"),
    estimate = d$time_ratio_31_vs28, ci_low = d$ratio_lo, ci_high = d$ratio_hi,
    test = "interval-censored Weibull model", statistic = d$z, p_value = d$p, n = d$n,
    source = "output/tables/14_interval_survreg.csv")
}
d <- cox[cox$trait == "new_corallites_on_tip" & grepl("^overall", cox$scope), ]
stopifnot(nrow(d) == 1)
rows[[length(rows) + 1]] <- paper_row(
  "budding_hazard_ratio", "Budding: hazard ratio, heated / ambient (Cox)",
  estimate = d$HR_31_vs28, ci_low = d$HR_lo, ci_high = d$HR_hi,
  test = "Cox model, stratified by source patch", statistic = d$z, p_value = d$p, n = d$n,
  source = "output/tables/14_cox_hazard_ratios.csv")
d <- cox_gx[cox_gx$trait == "new_corallites_on_tip", ]
rows[[length(rows) + 1]] <- paper_row(
  "budding_patch_x_temp", "Budding: source patch x temperature (exploratory; 4 fragments per patch and temperature)",
  test = "likelihood-ratio test", statistic = d$chisq_genet_trt_int, df1 = d$df_genet_trt_int,
  p_value = d$p_genet_trt_int, source = "output/tables/14_cox_genet_LRT.csv")

# Budding tank by tank (temperature was applied to tanks)
ph <- readRDS(file.path(DATA_PROC, "physio_clean.rds"))
by_tank <- ph |>
  dplyr::filter(wound == "yes") |>
  dplyr::group_by(treatment, tank, id) |>
  dplyr::summarise(budded = any(new_corallites_on_tip %in% c("yes", 1, TRUE)), .groups = "drop") |>
  dplyr::group_by(treatment, tank) |>
  dplyr::summarise(n = dplyr::n(), k = sum(budded), prop = k / n, .groups = "drop")
props <- setNames(by_tank$prop, by_tank$tank)
amb <- as.character(by_tank$tank[by_tank$treatment == "28C"])
diff_for <- function(set) mean(props[set]) - mean(props[setdiff(names(props), set)])
all_sets <- utils::combn(names(props), length(amb), simplify = FALSE)
obs <- diff_for(amb)
perm_p <- mean(vapply(all_sets, function(s) abs(diff_for(s)) >= abs(obs) - 1e-12, logical(1)))
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

# Tissue cover, separate microscope experiment
mic <- readRDS(file.path(DATA_PROC, "microscope_physio_clean.rds"))
first_cover <- mic |>
  dplyr::filter(tissue_over_wound %in% c("yes", 1, TRUE)) |>
  dplyr::group_by(id) |>
  dplyr::summarise(first_day = min(day), .groups = "drop")
rows[[length(rows) + 1]] <- paper_row(
  "tissue_cover_by_day1", "Microscope experiment: fragments with tissue covering the cut by Day 1",
  estimate = sum(first_cover$first_day <= 1), n = dplyr::n_distinct(mic$id),
  source = "data/processed/microscope_physio_clean.rds")

write_paper(rows, "regeneration")
