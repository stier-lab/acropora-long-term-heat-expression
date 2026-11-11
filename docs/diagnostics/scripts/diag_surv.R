suppressMessages(source(here::here("code", "00_setup.R"))); suppressMessages(library(survival))
ev <- read_csv(file.path(TBL_DIR, "14_event_intervals.csv"), show_col_types = FALSE) |>
  mutate(treatment = factor(treatment, c("28C", "31C")), thicket = factor(thicket))
contrasts(ev$treatment) <- contr.treatment(2)
cat("\n--- Events / censored by trait and temperature ---\n")
print(as.data.frame(ev |> group_by(trait, treatment) |> summarise(n = n(), events = sum(event), censored = n - events,
  interval_censored = sum(event & event_lower > 0 & event_upper > event_lower, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = treatment, values_from = c(events, censored, n, interval_censored))))
cat("\n--- Zero-event or all-event cells by trait x temperature x patch (separation risk) ---\n")
print(as.data.frame(ev |> group_by(trait, treatment, thicket) |> summarise(k = sum(event), n = n(), .groups = "drop") |>
  filter(k == 0 | k == n) |> group_by(trait) |> summarise(cells = paste0(treatment, "-", thicket, ":", k, "/", n, collapse = " "), .groups = "drop")))
cat("\n--- PH tests (08_checks.R) ---\n"); print(as.data.frame(read_csv(file.path(TBL_DIR, "14_cox_ph_tests.csv"), show_col_types = FALSE)))
cat("\n--- AFT: distribution choice and tank clustering (time ratio heated/ambient) ---\n")
rows <- list()
for (tr in unique(ev$trait)) {
  d <- ev |> filter(trait == tr)
  for (dist in c("weibull", "lognormal", "loglogistic", "exponential")) {
    f <- tryCatch(survreg(Surv(event_lower, event_upper, type = "interval2") ~ treatment + thicket, data = d, dist = dist), error = function(e) NULL)
    if (is.null(f)) next
    fc <- tryCatch(survreg(Surv(event_lower, event_upper, type = "interval2") ~ treatment + thicket + cluster(tank), data = d, dist = dist, robust = TRUE), error = function(e) NULL)
    b <- coef(f)[grep("^treatment", names(coef(f)))]; se <- sqrt(vcov(f)[grep("^treatment", rownames(vcov(f))), grep("^treatment", rownames(vcov(f)))])
    sec <- if (!is.null(fc)) sqrt(vcov(fc)[grep("^treatment", rownames(vcov(fc))), grep("^treatment", rownames(vcov(fc)))]) else NA
    rows[[length(rows) + 1]] <- tibble(trait = tr, dist, AIC = AIC(f), time_ratio = exp(b), p = 2 * pnorm(-abs(b / se)),
                                       p_tank_robust = 2 * pnorm(-abs(b / sec)), scale = f$scale)
  }
}
aft <- bind_rows(rows) |> group_by(trait) |> mutate(dAIC = AIC - min(AIC)) |> ungroup()
print(as.data.frame(aft |> select(-AIC) |> mutate(across(where(is.double), \(x) signif(x, 3)))))
cat("\n--- How the reported AFT handles intervals: compare with 14_interval_survreg.csv ---\n")
print(as.data.frame(read_csv(file.path(TBL_DIR, "14_interval_survreg.csv"), show_col_types = FALSE) |> select(any_of(c("trait","time_ratio_31_vs28","p","n","scale")))))
cat("\n--- Cox (stratified by patch): naive vs tank-robust SE ---\n")
for (tr in unique(ev$trait)) {
  d <- ev |> filter(trait == tr)
  f <- suppressWarnings(coxph(Surv(event_day, event) ~ treatment + strata(thicket), data = d))
  fc <- suppressWarnings(coxph(Surv(event_day, event) ~ treatment + strata(thicket) + cluster(tank), data = d))
  s <- summary(fc)$coefficients
  cat(sprintf("%-24s HR=%.3f naive p=%.3g  tank-robust p=%.3g  events=%d\n", tr, exp(coef(f)), summary(f)$coefficients[, "Pr(>|z|)"], s[, "Pr(>|z|)"], sum(d$event)))
}
cat("\n--- Ties: share of events sharing a scoring day ---\n")
print(as.data.frame(ev |> filter(event == 1) |> group_by(trait) |> summarise(events = n(), distinct_days = n_distinct(event_day), .groups = "drop")))
