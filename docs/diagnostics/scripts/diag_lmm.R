# Scratch diagnostics for the four physiology LMMs (stats-diagnostics skill audit)
suppressMessages(source(here::here("code", "00_setup.R")))
suppressMessages(library(lmerTest))
set.seed(42)
ctrl <- lme4::lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))

pam   <- readRDS(file.path(DATA_PROC, "pam_clean.rds")) |> mutate(thicket = factor(thicket)) |> filter(day >= 0)
color <- readRDS(file.path(DATA_PROC, "color_clean.rds")) |> mutate(thicket = factor(thicket)) |> filter(day >= 0)
bw    <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds")) |> mutate(thicket = factor(thicket)) |> filter(is.finite(pct_growth))
sym   <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
  filter(is.finite(cells_per_cm2), cells_per_cm2 > 0) |>
  mutate(thicket = factor(thicket), biopsy_day_c = biopsy_day - 1, y = log(cells_per_cm2))

specs <- list(
  pam_fvfm   = list(d = pam,   f = fv_fm ~ treatment * wound * day * thicket + (1 | tank) + (1 | id),
                    keys = c("treatment:day", "treatment:wound:day", "wound"), time = "day"),
  color      = list(d = color, f = color_num ~ treatment * wound * day * thicket + (1 | tank) + (1 | id),
                    keys = c("treatment:day", "treatment:wound:day"), time = "day"),
  symbionts  = list(d = sym,   f = y ~ treatment * wound * biopsy_day_c * thicket + (1 | tank),
                    keys = c("treatment:biopsy_day_c", "treatment:wound:biopsy_day_c"), time = NA),
  growth     = list(d = bw,    f = pct_growth ~ treatment * wound * thicket + (1 | tank),
                    keys = c("treatment", "wound", "treatment:wound"), time = NA)
)
fit <- function(f, d) lmerTest::lmer(f, data = d, REML = TRUE, control = ctrl)
keyp <- function(m, keys) { a <- as.data.frame(anova(m, type = 3)); a[keys, c("F value", "Pr(>F)")] }

out <- list(); loto <- list(); alt <- list()
for (nm in names(specs)) {
  s <- specs[[nm]]; d <- s$d; m <- fit(s$f, d)
  vc <- as.data.frame(lme4::VarCorr(m))
  r <- residuals(m); fv <- fitted(m)
  lev <- car::leveneTest(r ~ d$treatment)
  sd_ratio <- sd(r[d$treatment == "31C"]) / sd(r[d$treatment == "28C"])
  spear <- suppressWarnings(cor.test(abs(r), fv, method = "spearman"))
  ac <- NA_real_
  if (!is.na(s$time)) {
    dd <- d |> mutate(r = r) |> arrange(id, .data[[s$time]]) |> group_by(id) |> mutate(r_lag = dplyr::lag(r)) |> ungroup()
    ac <- cor(dd$r, dd$r_lag, use = "complete.obs")
  }
  kp <- keyp(m, s$keys)
  out[[nm]] <- tibble(model = nm, n_obs = nobs(m), n_frag = n_distinct(d$id), n_tank = n_distinct(d$tank),
    rank_deficient_cols = ncol(model.matrix(m)) - Matrix::rankMatrix(model.matrix(m))[1],
    singular = lme4::isSingular(m),
    var_tank = vc$vcov[vc$grp == "tank"], var_id = if (any(vc$grp == "id")) vc$vcov[vc$grp == "id"] else NA, var_resid = vc$vcov[vc$grp == "Residual"],
    resid_sd_ratio_heated_ambient = sd_ratio, levene_p = lev$`Pr(>F)`[1],
    absresid_vs_fitted_rho = unname(spear$estimate), absresid_vs_fitted_p = spear$p.value,
    lag1_resid_autocorr_within_fragment = ac,
    shapiro_p = shapiro.test(r)$p.value)
  # Leave one tank out
  for (tk in sort(unique(d$tank))) {
    m2 <- suppressMessages(fit(s$f, d[d$tank != tk, ]))
    k2 <- keyp(m2, s$keys)
    loto[[length(loto) + 1]] <- tibble(model = nm, dropped_tank = tk, term = s$keys, F = k2[, 1], p = k2[, 2])
  }
  # Alternative structures for repeated measures
  if (!is.na(s$time)) {
    f_slope_id   <- update(s$f, . ~ . - (1 | id) + (1 + day | id))
    f_slope_tank <- update(s$f, . ~ . - (1 | tank) + (1 + day | tank))
    f_dayfac     <- update(s$f, . ~ . - day - treatment:day - wound:day - day:thicket - treatment:wound:day - treatment:day:thicket - wound:day:thicket - treatment:wound:day:thicket)
    for (lab in c("random day slope per fragment", "random day slope per tank")) {
      ff <- if (lab == "random day slope per fragment") f_slope_id else f_slope_tank
      mm <- tryCatch(suppressMessages(fit(ff, d)), error = function(e) NULL)
      if (is.null(mm)) next
      k3 <- keyp(mm, s$keys)
      alt[[length(alt) + 1]] <- tibble(model = nm, structure = lab, term = s$keys, F = k3[, 1], p = k3[, 2],
                                       singular = lme4::isSingular(mm), AIC_vs_base = AIC(refitML(mm)) - AIC(refitML(m)))
    }
    # Curvature: day as a factor, test treatment x day
    dfac <- d |> mutate(dayf = factor(day))
    mf <- suppressMessages(fit(update(s$f, . ~ . + treatment:dayf - treatment:day), dfac))
    mq <- suppressMessages(fit(update(s$f, . ~ . + treatment:I(day^2) + I(day^2)), d))
    aq <- as.data.frame(anova(mq, type = 3))
    alt[[length(alt) + 1]] <- tibble(model = nm, structure = "add day^2 and temperature x day^2", term = "treatment:I(day^2)",
                                     F = aq["treatment:I(day^2)", "F value"], p = aq["treatment:I(day^2)", "Pr(>F)"], singular = lme4::isSingular(mq),
                                     AIC_vs_base = AIC(refitML(mq)) - AIC(refitML(m)))
  }
}
res <- bind_rows(out); lo <- bind_rows(loto); al <- bind_rows(alt)
options(width = 220, pillar.sigfig = 3)
print(as.data.frame(res |> mutate(across(where(is.double), \(x) signif(x, 3)))))
cat("\n--- Leave one tank out: range of p for key terms ---\n")
print(as.data.frame(lo |> group_by(model, term) |> summarise(F_min = min(F), F_max = max(F), p_min = min(p), p_max = max(p),
  worst_tank = dropped_tank[which.max(p)], .groups = "drop") |> mutate(across(where(is.double), \(x) signif(x, 3)))))
cat("\n--- Alternative time / random-effect structures ---\n")
print(as.data.frame(al |> mutate(across(where(is.double), \(x) signif(x, 3)))))
dir.create("output/diagnostics", showWarnings = FALSE, recursive = TRUE); saveRDS(list(res = res, loto = lo, alt = al), "output/diagnostics/diag_lmm.rds")
