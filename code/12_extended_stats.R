# =============================================================================
# Purpose: Primary statistical analysis with genet (thicket) as a first-class
#          fixed effect. Models full treatment × wound × day × thicket
#          interactions where data permit; tank and individual coral retained
#          as random effects.
#
#          With only 3 field-collected genets (a, c, d), the design treats
#          thicket as fixed rather than random (Bolker 2008, Gelman 2005):
#          when the number of levels of a grouping factor is small, a fixed-
#          effects representation has better estimation properties and
#          surfaces per-genet effect sizes directly.
#
#          For each response, this script produces:
#            - LMM/GLMM with treatment × wound × day × thicket fixed structure
#            - Type-III ANOVA for all main and interaction terms
#            - emmeans estimates per genet × treatment × wound (key contrasts)
#            - DHARMa residual diagnostics
#            - R²m / R²c (MuMIn::r.squaredGLMM)
#            - Per-genet treatment-effect summary (forest-plot-ready)
#
# Input:   data/processed/{pam_clean,color_clean,buoyant_weight_clean,
#                          symbiont_chl_clean,physio_clean}.rds
# Output:  output/models/12_<response>_lmm.rds         — fitted model objects
#          output/tables/12_anova_summary.csv           — ANOVA all terms
#          output/tables/12_emmeans_contrasts.csv       — pairwise contrasts
#          output/tables/12_genet_treatment_effects.csv — per-genet treatment effect
#          output/tables/12_r2_summary.csv              — R² per response
#          figures/12_diagnostics/<response>.png        — DHARMa plots
# =============================================================================

source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages({
  library(MuMIn)
})

DIAG_DIR <- file.path(FIG_DIR, "12_diagnostics")
dir.create(DIAG_DIR, recursive = TRUE, showWarnings = FALSE)

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
