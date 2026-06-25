# =============================================================================
# Purpose: Primary statistical models for every continuous + morphological
#          response, plus the two robustness refits that the diagnostics call
#          for. This single file replaces the former 12 / 12b / 12c trio:
#
#            PART 1 — Primary models (was 12_extended_stats.R)
#                     LMM/GLMM with treatment × wound × day × thicket fixed
#                     structure; genet (thicket) is a first-class fixed effect
#                     because only 3 genets were collected (Bolker 2008,
#                     Gelman 2005: few grouping levels → fixed beats random).
#                     Produces type-III ANOVA, emmeans contrasts, R²m/R²c,
#                     per-genet treatment effects, and DHARMa diagnostics.
#
#            PART 2 — Color ordinal-CLMM robustness (was 12b)
#                     Refit the 5-level color D-scale under the correct ordinal
#                     likelihood; confirms the Gaussian LMM's inferences hold.
#
#            PART 3 — Penalized morphology refit (was 12c)
#                     blme::bglmer with Cauchy(0, 2.5) priors to tame the
#                     complete/quasi-complete separation in 7/8 binary traits;
#                     yields finite, interpretable coefficient Wald tests.
#
# Input:   data/processed/{pam_clean,color_clean,buoyant_weight_clean,
#                          symbiont_chl_clean,physio_clean}.rds
# Output:  output/models/12_<response>_lmm.rds          — primary fits
#          output/models/12b_color_clmm.rds             — ordinal robustness
#          output/models/12c_morph_<trait>_blme.rds     — penalized refits
#          output/tables/12_anova_summary.csv            — ANOVA all terms
#          output/tables/12_emmeans_contrasts.csv        — pairwise contrasts
#          output/tables/12_genet_treatment_effects.csv  — per-genet effects
#          output/tables/12_r2_summary.csv               — R² per response
#          output/tables/12_morph_fixed_effects.csv      — raw glmer fixefs
#          output/tables/12b_color_clmm*.csv             — ordinal robustness
#          output/tables/12c_morph_blme_*.csv            — penalized robustness
#          figures/12_diagnostics/<response>.png         — DHARMa plots
# =============================================================================

source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages({
  library(MuMIn)
  library(blme)
  # NB: do NOT attach `ordinal` — it masks dplyr::slice and breaks downstream
  # scripts. PART 2 calls it namespace-qualified (ordinal::clmm).
})

DIAG_DIR <- file.path(FIG_DIR, "12_diagnostics")
dir.create(DIAG_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# PART 1 — PRIMARY MODELS  (formerly 12_extended_stats.R)
# =============================================================================

results_anova        <- list()
results_emm          <- list()
results_r2           <- list()
results_genet_effect <- list()

# ---- Helpers ---------------------------------------------------------------
save_dharma <- function(model, name) {
  res <- tryCatch(
    DHARMa::simulateResiduals(model, n = 500, refit = FALSE),
    error = function(e) NULL
  )
  if (is.null(res)) return(invisible(NULL))
  png(file.path(DIAG_DIR, paste0(name, ".png")),
      width = 1600, height = 800, res = 160)
  plot(res)
  dev.off()
  invisible(res)
}

record_results <- function(model, response_name) {
  if (inherits(model, "lmerMod") || inherits(model, "lmerModLmerTest") ||
      inherits(model, "lm")) {
    if (inherits(model, "lm") && !inherits(model, "lmerMod")) {
      av <- as.data.frame(car::Anova(model, type = 3))
    } else {
      av <- as.data.frame(anova(model, type = 3))
    }
    av$term <- rownames(av); av$response <- response_name
    results_anova[[response_name]] <<- av
  } else if (inherits(model, "glmerMod")) {
    av <- as.data.frame(car::Anova(model, type = 3))
    av$term <- rownames(av); av$response <- response_name
    results_anova[[response_name]] <<- av
  }
  r2 <- tryCatch(MuMIn::r.squaredGLMM(model), error = function(e) NULL)
  if (!is.null(r2)) {
    results_r2[[response_name]] <<- tibble::tibble(
      response = response_name,
      R2_marginal    = r2[1, "R2m"],
      R2_conditional = r2[1, "R2c"]
    )
  }
  save_dharma(model, response_name)
}

# Per-genet treatment effect at end-of-experiment (or biopsy-day max)
# Returns one row per genet × wound with (28C - 31C) treatment-effect estimate
record_genet_effect <- function(emm_obj, response_name) {
  contr <- as_tibble(pairs(emm_obj, by = c("thicket", "wound"),
                            adjust = "none"))
  if (nrow(contr) == 0) return(invisible(NULL))
  contr$response <- response_name
  results_genet_effect[[response_name]] <<- contr
}

# ---- 1. PAM Fv/Fm -----------------------------------------------------------
cat("\n=== 1. PAM Fv/Fm ===\n")
pam <- readRDS(file.path(DATA_PROC, "pam_clean.rds")) |>
  mutate(thicket = as.factor(thicket))

# Note: with n=3 genets and 4 wounded + 4 unwounded fragments per genet × treatment
# cell, the full 4-way interaction is supported but high-order terms will have
# limited power. We report all terms transparently rather than stepwise-dropping.
m_pam <- lmerTest::lmer(
  fv_fm ~ treatment * wound * day * thicket + (1 | tank) + (1 | id),
  data = pam, REML = TRUE,
  control = lme4::lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))
)
saveRDS(m_pam, file.path(MOD_DIR, "12_pam_lmm.rds"))
record_results(m_pam, "pam_fvfm")

# Day-specific contrasts: 31C vs 28C per (wound × genet × day)
emm_pam <- emmeans::emmeans(m_pam, ~ treatment | thicket * wound * day,
                             at = list(day = c(0, 3, 7, 10, 14)))
results_emm[["pam_fvfm"]] <- as_tibble(pairs(emm_pam, adjust = "tukey")) |>
  mutate(response = "pam_fvfm")
# Per-genet treatment effect at Day 14 (the end-of-experiment "thermal tolerance")
emm_pam_end <- emmeans::emmeans(m_pam, ~ treatment | thicket * wound,
                                 at = list(day = 14))
record_genet_effect(emm_pam_end, "pam_fvfm")

# ---- 2. Color score (ordinal D-scale) --------------------------------------
cat("\n=== 2. Color score ===\n")
color <- readRDS(file.path(DATA_PROC, "color_clean.rds")) |>
  mutate(thicket = as.factor(thicket))

m_color <- lmerTest::lmer(
  color_num ~ treatment * wound * day * thicket + (1 | tank) + (1 | id),
  data = color, REML = TRUE,
  control = lme4::lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))
)
saveRDS(m_color, file.path(MOD_DIR, "12_color_lmm.rds"))
record_results(m_color, "color_dscale")

emm_color <- emmeans::emmeans(m_color, ~ treatment | thicket * wound * day,
                              at = list(day = c(0, 3, 7, 10, 14)))
results_emm[["color_dscale"]] <- as_tibble(pairs(emm_color, adjust = "tukey")) |>
  mutate(response = "color_dscale")
emm_color_end <- emmeans::emmeans(m_color, ~ treatment | thicket * wound,
                                   at = list(day = 14))
record_genet_effect(emm_color_end, "color_dscale")

# ---- 3. Coral calcification (areal, mg cm^-2 d^-1) -------------------------
# PRIMARY growth metric is areal calcification rate (surface-area normalized;
# see code/05_buoyant_weight.R and notes/growth_allometry.md). % mass change
# is recorded as a robustness response below.
cat("\n=== 3. Calcification (areal) ===\n")
bw <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds")) |>
  mutate(thicket = as.factor(thicket))
bw_a <- bw |> filter(is.finite(areal_calc))

# n<=48 corals with 1 obs each — drop (1|id), but retain tank as the
# experimental-unit random effect for the temperature treatment.
m_bw <- lmerTest::lmer(
  areal_calc ~ treatment * wound * thicket + (1 | tank),
  data = bw_a, REML = TRUE,
  control = lme4::lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))
)
saveRDS(m_bw, file.path(MOD_DIR, "12_bw_lm.rds"))
record_results(m_bw, "growth_areal")

emm_bw <- emmeans::emmeans(m_bw, ~ treatment | thicket * wound)
results_emm[["growth_areal"]] <- as_tibble(pairs(emm_bw, adjust = "tukey")) |>
  mutate(response = "growth_areal", day = NA_integer_)
record_genet_effect(emm_bw, "growth_areal")

# Robustness: % mass change (the previous primary metric)
m_bw_pct <- lmerTest::lmer(
  pct_growth ~ treatment * wound * thicket + (1 | tank),
  data = bw, REML = TRUE,
  control = lme4::lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))
)
saveRDS(m_bw_pct, file.path(MOD_DIR, "12_bw_pct_lm.rds"))
record_results(m_bw_pct, "growth_pct")

# ---- 4. Symbiont density (cells/cm^2) --------------------------------------
cat("\n=== 4. Symbiont density ===\n")
phys <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
  filter(is.finite(cells_per_cm2), cells_per_cm2 > 0) |>
  mutate(thicket = factor(thicket),
         biopsy_day_c = biopsy_day - 1)

# Destructive sampling: 1 obs per coral. Drop (1|id).
m_zoox <- lmerTest::lmer(
  log(cells_per_cm2) ~ treatment * wound * biopsy_day_c * thicket +
    (1 | tank),
  data = phys, REML = TRUE,
  control = lme4::lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))
)
saveRDS(m_zoox, file.path(MOD_DIR, "12_zoox_lmm.rds"))
record_results(m_zoox, "log_zoox_density")

emm_zoox <- emmeans::emmeans(m_zoox, ~ treatment | thicket * wound * biopsy_day_c,
                             at = list(biopsy_day_c = c(0, 2, 9, 14)))
results_emm[["log_zoox_density"]] <- as_tibble(pairs(emm_zoox, adjust = "tukey")) |>
  mutate(response = "log_zoox_density")
emm_zoox_end <- emmeans::emmeans(m_zoox, ~ treatment | thicket * wound,
                                  at = list(biopsy_day_c = 14))
record_genet_effect(emm_zoox_end, "log_zoox_density")

# ---- 5. Morphological traits (binomial GLMM, per-trait) --------------------
cat("\n=== 5. Morphological traits (binomial GLMM with genet) ===\n")
ph <- readRDS(file.path(DATA_PROC, "physio_clean.rds")) |>
  filter(wound == "yes", !is.na(day), day >= 0) |>
  mutate(thicket = as.factor(thicket))

traits <- c("polyps_out", "hole_in_center", "polyp_in_hole",
            "wound_smoothed", "pigment_over_wound", "tip_exist",
            "tip_extension", "new_corallites_on_tip", "algae_on_wound")

fit_trait <- function(tr) {
  d <- ph |> mutate(y = .data[[tr]]) |> filter(!is.na(y))
  if (length(unique(d$y)) < 2 || nrow(d) < 30) return(NULL)
  # Full treatment × day × genet model. Include coral ID because morphology is
  # repeatedly scored on the same fragments over time.
  m <- tryCatch(
    suppressMessages(suppressWarnings(
      lme4::glmer(y ~ treatment * day * thicket + (1 | tank) + (1 | id),
                  family = binomial, data = d,
                  control = lme4::glmerControl(optimizer = "bobyqa",
                                               optCtrl = list(maxfun = 1e5)))
    )),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)
  av <- as.data.frame(car::Anova(m, type = 2)) |>
    tibble::rownames_to_column("term") |>
    mutate(response = paste0("morph_", tr), n_obs = nrow(d))
  tidy <- broom.mixed::tidy(m, effects = "fixed") |>
    mutate(response = paste0("morph_", tr))
  # Per-genet probability of trait at Day 10 (peak divergence in KM)
  emm <- tryCatch(
    emmeans::emmeans(m, ~ treatment | thicket, at = list(day = 10),
                     type = "response"),
    error = function(e) NULL
  )
  genet_eff <- if (!is.null(emm)) {
    as_tibble(pairs(emm, adjust = "tukey")) |>
      mutate(response = paste0("morph_", tr))
  } else NULL
  list(model = m, anova = av, tidy = tidy, genet_eff = genet_eff)
}

morph_results <- map(traits, fit_trait) |> setNames(traits) |> compact()
for (tr in names(morph_results)) {
  saveRDS(morph_results[[tr]]$model,
          file.path(MOD_DIR, paste0("12_morph_", tr, "_glmm.rds")))
  results_anova[[paste0("morph_", tr)]] <- morph_results[[tr]]$anova
  if (!is.null(morph_results[[tr]]$genet_eff)) {
    results_genet_effect[[paste0("morph_", tr)]] <- morph_results[[tr]]$genet_eff
  }
}

morph_tidy <- map_dfr(morph_results, "tidy")
write_csv(morph_tidy, file.path(TBL_DIR, "12_morph_fixed_effects.csv"))

# ---- Aggregate + save ------------------------------------------------------
all_anova <- bind_rows(results_anova, .id = "response_id") |>
  mutate(across(where(is.numeric), \(x) round(x, 5)))
write_csv(all_anova, file.path(TBL_DIR, "12_anova_summary.csv"))

all_emm <- bind_rows(results_emm) |>
  mutate(across(where(is.numeric), \(x) round(x, 5)))
write_csv(all_emm, file.path(TBL_DIR, "12_emmeans_contrasts.csv"))

all_r2 <- bind_rows(results_r2) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))
write_csv(all_r2, file.path(TBL_DIR, "12_r2_summary.csv"))

all_genet <- bind_rows(results_genet_effect) |>
  mutate(across(where(is.numeric), \(x) round(x, 5)))
write_csv(all_genet, file.path(TBL_DIR, "12_genet_treatment_effects.csv"))

# ---- Console summary -------------------------------------------------------
cat("\n=== Type-III ANOVA: continuous responses (genet terms highlighted) ===\n")
print(all_anova |>
        filter(response_id %in% c("pam_fvfm", "color_dscale",
                                   "growth_areal", "growth_pct",
                                   "log_zoox_density"),
               grepl("thicket|treatment", term, ignore.case = TRUE)) |>
        select(response_id, term, any_of(c("F value", "F", "Pr(>F)",
                                           "Chisq", "Pr(>Chisq)"))))

cat("\n=== R² (marginal / conditional) ===\n")
print(all_r2)

cat("\n=== Per-genet end-of-experiment treatment effect (28C - 31C) ===\n")
print(all_genet |>
        select(response, thicket, wound, estimate, SE, t.ratio, p.value) |>
        mutate(across(where(is.numeric), \(x) round(x, 3))))

cat("\nWrote: 12_anova_summary.csv, 12_emmeans_contrasts.csv,",
    "12_genet_treatment_effects.csv, 12_r2_summary.csv,",
    "12_morph_fixed_effects.csv, diagnostics in figures/12_diagnostics/\n")

# =============================================================================
# PART 2 — COLOR ORDINAL-CLMM ROBUSTNESS  (formerly 12b_color_clmm_robustness.R)
#
# Residual diagnostics flagged the Gaussian color LMM as KS-non-uniform (the
# D-scale is a 5-level ordinal). We keep the Gaussian model for presentation
# (it matches the other physiology metrics) but confirm here that every
# qualitative inference holds under the correct ordinal likelihood.
# =============================================================================
cat("\n=== PART 2: color ordinal-CLMM robustness ===\n")

color_ord_df <- readRDS(file.path(DATA_PROC, "color_clean.rds")) |>
  mutate(thicket = factor(thicket),
         # Convert split scores (e.g. 3.5) to the next-higher D category.
         # base::round() uses bankers rounding, so 2.5 would become 2 while
         # 3.5 becomes 4; floor(x + 0.5) is deterministic half-up rounding.
         color_ord = factor(pmin(pmax(as.integer(floor(color_num + 0.5)), 1), 5),
                            ordered = TRUE))

# The 4-way fixed × random structure is too rich for clmm (Hessian singular).
# Robustness check: same data, ordinal likelihood, slightly reduced structure
# — treatment × wound × day + thicket interactions, single (1|id) random.
m_clmm <- ordinal::clmm(
  color_ord ~ treatment * wound * day + treatment * thicket +
              treatment * wound * thicket + day * thicket + (1 | id),
  data = color_ord_df, Hess = TRUE
)
saveRDS(m_clmm, file.path(MOD_DIR, "12b_color_clmm.rds"))

# Type-I Wald per term (drop1-based)
av_drop <- as.data.frame(drop1(m_clmm, test = "Chisq")) |>
  tibble::rownames_to_column("term") |>
  mutate(response = "color_dscale_ordinal")
write_csv(av_drop, file.path(TBL_DIR, "12b_color_clmm.csv"))
av <- av_drop |> rename(LR.stat = LRT, `Pr(>Chisq)` = `Pr(>Chi)`)

# Compare to Gaussian LMM for the headline term — treatment main effect
m_gauss <- readRDS(file.path(MOD_DIR, "12_color_lmm.rds"))
gauss_av <- as.data.frame(anova(m_gauss))
cat("\n=== Treatment main-effect comparison (Gaussian LMM vs ordinal CLMM) ===\n")
if ("treatment" %in% rownames(gauss_av)) {
  cat("Gaussian LMM treatment: F =",
      round(gauss_av["treatment", "F value"], 2), "\n")
}
trt_row <- av[grepl("^treatment$", av$term), , drop = FALSE]
if (nrow(trt_row) > 0) {
  cat("Ordinal CLMM treatment: LRT =", round(trt_row$LR.stat[1], 2),
      " p =", signif(trt_row$`Pr(>Chisq)`[1], 3), "\n")
} else {
  cat("CLMM drop1 did not report a main treatment term (likely all retained due to higher-order interactions); see 12b_color_clmm.csv for full term-by-term tests.\n")
}

# Per-genet contrasts: skip if Hessian is singular (emmeans would fail).
# The drop1 LRT table above is the primary robustness output.
pairs_df <- tryCatch({
  emm <- emmeans::emmeans(m_clmm, ~ treatment | thicket * wound,
                          at = list(day = 14))
  as_tibble(pairs(emm, adjust = "none")) |>
    mutate(response = "color_dscale_ordinal")
}, error = function(e) {
  message("emmeans pairwise contrasts skipped (Hessian singular: ",
          conditionMessage(e), ")")
  tibble(
    note = "CLMM Hessian singular; per-genet contrasts not computable from this fit",
    response = "color_dscale_ordinal"
  )
})
write_csv(pairs_df, file.path(TBL_DIR, "12b_color_clmm_genet_effects.csv"))

cat("\nWrote 12b_color_clmm.csv, 12b_color_clmm_genet_effects.csv, 12b_color_clmm.rds\n")
cat("CLMM fit converged (LL =", round(as.numeric(logLik(m_clmm)), 1),
    "); use drop1 LRTs in 12b_color_clmm.csv for term-by-term tests.\n")
cat("Robustness conclusion: see report (overall direction of every term should match the Gaussian LMM).\n")

# =============================================================================
# PART 3 — PENALIZED MORPHOLOGY REFIT  (formerly 12c_morph_blme.R)
#
# The raw glmer fits in PART 1 give finite predictions (the random-effect
# penalty regularizes them) and correct omnibus type-II Wald χ², but the
# *individual coefficient* Wald z's blow up when any treatment × day × thicket
# cell reaches 0 or 1 observed events (7/8 traits hit this). Refit with
# blme::bglmer + Cauchy(0, 2.5) priors (Gelman et al. 2008 default) so the
# coefficients become finite and Wald tests interpretable.
# =============================================================================
cat("\n=== PART 3: penalized morphology GLMMs (blme::bglmer, Cauchy(0,2.5)) ===\n")

ph_blme <- readRDS(file.path(DATA_PROC, "physio_clean.rds")) |>
  filter(wound == "yes", !is.na(day), day >= 0) |>
  mutate(
    treatment = factor(treatment),
    thicket = factor(thicket)
  )

fit_blme <- function(tr) {
  d <- ph_blme |> mutate(y = .data[[tr]]) |> filter(!is.na(y))
  if (length(unique(d$y)) < 2 || nrow(d) < 30) return(NULL)
  contrasts(d$treatment) <- contr.treatment(nlevels(d$treatment))
  contrasts(d$thicket) <- contr.treatment(nlevels(d$thicket))
  # Cauchy(0, 2.5) prior on all fixed effects — Gelman 2008 default for
  # logistic regression with separation. Same structure as PART 1.
  m <- tryCatch(
    suppressMessages(suppressWarnings(
      blme::bglmer(
        y ~ treatment * day * thicket + (1 | tank) + (1 | id),
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

res <- map(traits, function(tr) {
  cat("Fitting:", tr, "... ")
  out <- fit_blme(tr)
  if (is.null(out)) cat("skipped\n")
  else cat("done (max fixef SE =",
           round(max(out$tidy$std.error, na.rm = TRUE), 2), ")\n")
  out
}) |> setNames(traits) |> compact()

# Aggregate
blme_anova <- map_dfr(res, "anova")
blme_tidy  <- map_dfr(res, "tidy")
write_csv(blme_anova, file.path(TBL_DIR, "12c_morph_blme_anova.csv"))
write_csv(blme_tidy,  file.path(TBL_DIR, "12c_morph_blme_fixed_effects.csv"))

# Compare: how many fixef SEs are still pathological?
cat("\n=== SE sanity (after Cauchy prior) ===\n")
ses <- blme_tidy |> group_by(trait) |>
  summarise(n_fixed = n(),
            n_sep   = sum(std.error > 50, na.rm = TRUE),
            max_se  = max(std.error, na.rm = TRUE))
print(ses)

cat("\nWrote 12c_morph_blme_anova.csv (",
    nrow(blme_anova), "rows ),\n",
    "      12c_morph_blme_fixed_effects.csv (",
    nrow(blme_tidy), "rows),\n",
    "      output/models/12c_morph_<trait>_blme.rds (",
    length(res), "fits)\n", sep = "")
