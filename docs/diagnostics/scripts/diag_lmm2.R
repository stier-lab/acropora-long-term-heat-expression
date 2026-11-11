suppressMessages(source(here::here("code", "00_setup.R")))
suppressMessages({library(lmerTest); library(nlme)})
ctrl <- lme4::lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4), optimizer = "bobyqa")
pam   <- readRDS(file.path(DATA_PROC, "pam_clean.rds")) |> mutate(thicket = factor(thicket)) |> filter(day >= 0)
color <- readRDS(file.path(DATA_PROC, "color_clean.rds")) |> mutate(thicket = factor(thicket)) |> filter(day >= 0)
sym   <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |> filter(is.finite(cells_per_cm2), cells_per_cm2 > 0) |>
  mutate(thicket = factor(thicket), biopsy_day_c = biopsy_day - 1, y = log(cells_per_cm2))
rep <- read_csv(file.path(TBL_DIR, "12_anova_summary.csv"), show_col_types = FALSE)
row <- function(model, structure, m, keys, extra = list()) {
  a <- as.data.frame(anova(m, type = 3))
  tibble(model, structure, term = keys, F = a[keys, "F value"], df_den = a[keys, "DenDF"], p = a[keys, "Pr(>F)"],
         singular = lme4::isSingular(m), conv = length(m@optinfo$conv$lme4$messages) > 0)
}
out <- list()
# Fv/Fm and color: slopes per tank and per fragment
for (nm in c("pam_fvfm", "color")) {
  d <- if (nm == "pam_fvfm") pam else color
  y <- if (nm == "pam_fvfm") "fv_fm" else "color_num"
  keys <- c("treatment:day", "treatment:wound:day", "treatment", "wound")
  f_both <- as.formula(paste(y, "~ treatment * wound * day * thicket + (1 + day | tank) + (1 + day | id)"))
  f_both_nc <- as.formula(paste(y, "~ treatment * wound * day * thicket + (1 + day || tank) + (1 + day || id)"))
  m <- lmerTest::lmer(f_both, data = d, control = ctrl); out[[length(out)+1]] <- row(nm, "day slopes per tank and fragment", m, keys)
  m2 <- lmerTest::lmer(f_both_nc, data = d, control = ctrl); out[[length(out)+1]] <- row(nm, "day slopes per tank and fragment (uncorrelated)", m2, keys)
  # Unequal variance by temperature + slopes (nlme; tank and fragment nested in tank)
  dd <- d |> mutate(id_in_tank = interaction(tank, id, drop = TRUE))
  fx <- as.formula(paste(y, "~ treatment * wound * day * thicket"))
  ml <- tryCatch(lme(fx, random = list(tank = ~ 1 + day, id_in_tank = ~ 1 + day), weights = varIdent(form = ~ 1 | treatment),
                     data = dd, method = "REML", control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)), error = function(e) {message(e$message); NULL})
  if (!is.null(ml)) {
    a <- anova(ml, type = "marginal")
    out[[length(out)+1]] <- tibble(model = nm, structure = "slopes + separate variance by temperature (nlme)", term = keys,
      F = a[keys, "F-value"], df_den = a[keys, "denDF"], p = a[keys, "p-value"], singular = NA, conv = FALSE)
    cat(nm, "nlme variance ratio heated/ambient (SD):", signif(coef(ml$modelStruct$varStruct, unconstrained = FALSE), 3), "\n")
  }
}
# Symbionts: slope per tank
keys <- c("treatment:biopsy_day_c", "treatment:wound:biopsy_day_c", "treatment")
m <- lmerTest::lmer(y ~ treatment * wound * biopsy_day_c * thicket + (1 + biopsy_day_c | tank), data = sym, control = ctrl)
out[[length(out)+1]] <- row("symbionts", "day slope per tank", m, keys)
m0 <- lmerTest::lmer(y ~ treatment * wound * biopsy_day_c * thicket + (1 | tank), data = sym, control = ctrl)
cat("symbionts AIC slope - intercept-only (ML):", AIC(refitML(m)) - AIC(refitML(m0)), "\n")
ms <- tryCatch(lme(y ~ treatment * wound * biopsy_day_c * thicket, random = ~ 1 + biopsy_day_c | tank, weights = varIdent(form = ~ 1 | treatment),
                   data = sym, control = lmeControl(opt = "optim")), error = function(e) NULL)
if (!is.null(ms)) { a <- anova(ms, type = "marginal"); out[[length(out)+1]] <- tibble(model = "symbionts", structure = "slope per tank + separate variance (nlme)", term = keys,
  F = a[keys, "F-value"], df_den = a[keys, "denDF"], p = a[keys, "p-value"], singular = NA, conv = FALSE) }
# Symbiont outliers: refit without DHARMa outliers
set.seed(42)
sim <- DHARMa::simulateResiduals(m0, n = 500)
oi <- DHARMa::outliers(sim)
cat("symbiont outlier rows:", length(oi), "\n"); print(as.data.frame(sym[oi, c("id", "treatment", "wound", "biopsy_day", "tank", "thicket", "cells_per_cm2", "count_source")]))
m3 <- lmerTest::lmer(y ~ treatment * wound * biopsy_day_c * thicket + (1 | tank), data = sym[-oi, ], control = ctrl)
out[[length(out)+1]] <- row("symbionts", "without DHARMa outliers", m3, keys)
base <- rep |> filter(term %in% c("treatment:day","treatment:wound:day","treatment","wound","treatment:biopsy_day_c","treatment:wound:biopsy_day_c")) |>
  transmute(model = recode(response_id, log_zoox_density = "symbionts", color_dscale = "color", growth_pct = "growth"), structure = "REPORTED (03_physiology.R)", term,
            F = `F value`, df_den = DenDF, p = `Pr(>F)`)
res <- bind_rows(base |> filter(model != "growth"), bind_rows(out)) |> arrange(model, term, structure)
options(width = 220); print(as.data.frame(res |> mutate(across(where(is.double), \(x) signif(x, 3)))))
dir.create("output/diagnostics", showWarnings = FALSE, recursive = TRUE); saveRDS(res, "output/diagnostics/diag_lmm2.rds")
