# =============================================================================
# Purpose: Penalized refit of the morphology GLMMs to handle complete /
#          quasi-complete separation flagged by the morphology diagnostics.
#
#          The raw `glmer(... family = binomial)` fits in script 12 produce
#          finite predictions (because the random-effect penalty already
#          regularizes them) and correct omnibus type-II Wald χ² from
#          `car::Anova`, but the *individual coefficient* Wald z's blow up
#          when any cell of treatment × day × thicket reaches 0 or 1
#          observed events. 7 of 8 traits hit this.
#
#          Fix: refit with `blme::bglmer` using weakly informative
#          Cauchy(0, 2.5) priors on the fixed-effect coefficients (Gelman
#          et al. 2008 default for logistic regression with separation).
#          Same data, same fixed and random structure — only the prior
#          is added. Coefficients become finite and Wald tests are
#          interpretable.
#
# Input:   data/processed/physio_clean.rds
# Output:  output/models/12c_morph_<trait>_blme.rds
#          output/tables/12c_morph_blme_fixed_effects.csv
#          output/tables/12c_morph_blme_anova.csv
# =============================================================================

source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages({ library(blme) })

ph <- readRDS(file.path(DATA_PROC, "physio_clean.rds")) |>
  filter(wound == "yes", !is.na(day), day >= 0) |>
  mutate(
    treatment = factor(treatment),
    thicket = factor(thicket)
  )

traits <- c("polyps_out", "hole_in_center", "polyp_in_hole",
            "wound_smoothed", "pigment_over_wound", "tip_exist",
            "tip_extension", "new_corallites_on_tip", "algae_on_wound")

fit_blme <- function(tr) {
  d <- ph |> mutate(y = .data[[tr]]) |> filter(!is.na(y))
  if (length(unique(d$y)) < 2 || nrow(d) < 30) return(NULL)
  contrasts(d$treatment) <- contr.treatment(nlevels(d$treatment))
  contrasts(d$thicket) <- contr.treatment(nlevels(d$thicket))
  # Cauchy(0, 2.5) prior on all fixed effects — Gelman 2008 default for
  # logistic regression with separation. Same structure as script 12.
  m <- tryCatch(
    suppressMessages(suppressWarnings(
      blme::bglmer(
        y ~ treatment * day * thicket + (1 | tank),
        family   = binomial,
        data     = d,
        fixef.prior = "t(scale = 2.5, df = 1)",   # Cauchy(0, 2.5) — Gelman 2008
        control  = lme4::glmerControl(optimizer = "bobyqa",
                                      optCtrl = list(maxfun = 1e5))
      )
    )),
    error = function(e) {
      message("bglmer failed for ", tr, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(m)) return(NULL)
  saveRDS(m, file.path(MOD_DIR, paste0("12c_morph_", tr, "_blme.rds")))
  av <- as.data.frame(car::Anova(m, type = 2)) |>
    tibble::rownames_to_column("term") |>
    mutate(trait = tr)
  tidy <- broom.mixed::tidy(m, effects = "fixed") |>
    mutate(trait = tr, max_se = max(std.error, na.rm = TRUE))
  list(anova = av, tidy = tidy)
}

cat("\n=== Penalized morphology GLMMs (blme::bglmer, Cauchy(0,2.5) priors) ===\n")
res <- map(traits, function(tr) {
  cat("Fitting:", tr, "... ")
  out <- fit_blme(tr)
  if (is.null(out)) cat("skipped\n")
  else cat("done (max fixef SE =",
           round(max(out$tidy$std.error, na.rm = TRUE), 2), ")\n")
  out
}) |> setNames(traits) |> compact()

# Aggregate
all_anova <- map_dfr(res, "anova")
all_tidy  <- map_dfr(res, "tidy")
write_csv(all_anova, file.path(TBL_DIR, "12c_morph_blme_anova.csv"))
write_csv(all_tidy,  file.path(TBL_DIR, "12c_morph_blme_fixed_effects.csv"))

# Compare: how many fixef SEs are still pathological?
cat("\n=== SE sanity (after Cauchy prior) ===\n")
ses <- all_tidy |> group_by(trait) |>
  summarise(n_fixed = n(),
            n_sep   = sum(std.error > 50, na.rm = TRUE),
            max_se  = max(std.error, na.rm = TRUE))
print(ses)

cat("\nWrote 12c_morph_blme_anova.csv (",
    nrow(all_anova), "rows ),\n",
    "      12c_morph_blme_fixed_effects.csv (",
    nrow(all_tidy), "rows),\n",
    "      output/models/12c_morph_<trait>_blme.rds (",
    length(res), "fits)\n", sep = "")
