# =============================================================================
# Purpose: Formal LRT test of whether adding genet × treatment improves fit
#          over a baseline without genet structure, plus a reaction-norm
#          figure showing each genet's mean response shift under heating.
#
#          Note: the BROADER genet analysis (per-response interaction tests,
#          per-genet emmeans contrasts, full ANOVA) is integrated into the
#          primary statistics in code/12_models.R. This script is the
#          focused formal test that the Progress Notes asked for (2026-04-29:
#          "Add in genet effects and see if there is a difference"), kept as
#          a separate file so reviewers can replicate the comparison directly.
#
#          For each response we compare:
#            null:  response ~ treatment * wound * day + thicket
#                                + (1|tank) + (1|id)
#            genet: response ~ treatment * wound * day * thicket
#                                + (1|tank) + (1|id)         # genet as fixed
#          The likelihood-ratio test uses the same random structure in both
#          models; the additional fixed-effect terms test genet variation in
#          plasticity. Significant LRT
#          → adding genet × treatment significantly improves fit.
#
# What & why: code/12 already estimates the genet effects inside each full model.
#   This script asks one separable question that a reviewer can run on its own:
#   does genotype matter at all? We compare two nested models that differ ONLY in
#   whether genet interacts with the other predictors. The null lets each genet
#   have its own intercept (thicket as an additive main effect) but forces all
#   genets to respond to heat/wound/time the SAME way. The genet model lets
#   genet × everything vary (genet as a fully crossed fixed effect). A likelihood-
#   ratio test of null vs genet then tests whether genotypes differ in their
#   PLASTICITY (their reaction to treatment), i.e. G × E. Both models are fit with
#   REML = FALSE (ML) because an LRT comparing models with DIFFERENT fixed effects
#   is only valid under ML — REML likelihoods are not comparable across different
#   fixed structures. The companion figure draws the reaction norms (genet means
#   at 28 vs 31 °C): parallel lines = additive temperature effect; crossing/fanning
#   lines = genotype × environment.
#
# Input:   data/processed/{pam_clean,color_clean,buoyant_weight_clean,
#                          symbiont_chl_clean}.rds
# Output:  output/tables/13_genet_anova.csv             — LRT comparison
#          output/tables/13_genet_emmeans.csv           — per-genet end-of-exp means
#          figures/13_genet_response_panel.{pdf,png}    — 4-panel reaction norms
#          output/models/13_<response>_genet_lmm.rds   — alt-model objects
# =============================================================================

# 00_setup.R: packages, shared paths (DATA_PROC, MOD_DIR, TBL_DIR), palettes
# (PAL_GENO) and the theme_pub() / save_fig() figure helpers.
source(here::here("code", "00_setup.R"))

# Load the four continuous responses; thicket (genet) → factor so it enters the
# models as discrete fixed levels rather than a number.
pam   <- readRDS(file.path(DATA_PROC, "pam_clean.rds")) |>
  mutate(thicket = factor(thicket)) |>
  filter(day >= 0)   # exclude the pre-treatment baseline from the day-slope models
color <- readRDS(file.path(DATA_PROC, "color_clean.rds")) |>
  mutate(thicket = factor(thicket)) |>
  filter(day >= 0)   # (matches 12_models, so the genet LRTs use the same data)
bw    <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds")) |>
  mutate(thicket = factor(thicket))
phys  <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
  filter(is.finite(cells_per_cm2), cells_per_cm2 > 0) |>
  mutate(thicket = factor(thicket),
         biopsy_day_c = biopsy_day - 1)   # centred time axis (matches code/12)

# fit_and_report(): build the null and genet models for one response, fit both by
# ML, save the genet model, and return a one-row LRT summary.
#   fixed_extra  — swap the time variable name (e.g. "day" → "biopsy_day_c") for
#                  responses whose time column differs.
#   include_id   — drop the (1|id) random intercept for destructively-sampled
#                  responses (one observation per coral → id not identifiable).
fit_and_report <- function(data, response, name, fixed_extra = NULL,
                           include_id = TRUE) {
  id_term <- if (include_id) " + (1 | id)" else ""
  # NULL model: genet is additive only (thicket main effect) — all genets share
  # one treatment × wound × day response surface.
  rhs_null  <- paste0("treatment * wound * day + thicket + (1 | tank)",
                      id_term)
  # GENET model: genet fully crossed — each genet gets its own response surface.
  rhs_genet <- paste0("treatment * wound * day * thicket + (1 | tank)",
                      id_term)
  if (!is.null(fixed_extra)) {
    rhs_null  <- str_replace(rhs_null,  "day", fixed_extra)
    rhs_genet <- str_replace(rhs_genet, "day", fixed_extra)
  }
  f0 <- as.formula(paste(response, "~", rhs_null))
  f1 <- as.formula(paste(response, "~", rhs_genet))
  # REML = FALSE (ML): required so the two models' likelihoods are comparable —
  # they differ in fixed effects, and REML is invalid for that comparison.
  m0 <- suppressWarnings(suppressMessages(
    lme4::lmer(f0, data = data, REML = FALSE)))
  m1 <- suppressWarnings(suppressMessages(
    lme4::lmer(f1, data = data, REML = FALSE)))
  saveRDS(m1, file.path(MOD_DIR, paste0("13_", name, "_genet_lmm.rds")))
  # LRT comparison: anova(m0, m1) gives the χ², df (= number of extra genet
  # interaction terms), and p-value; we also report ΔAIC as a fit-vs-complexity
  # check that does not rely on a p-value threshold.
  lrt <- anova(m0, m1)
  tibble(
    response      = name,
    n_obs         = nrow(data),
    aic_null      = AIC(m0),
    aic_genet     = AIC(m1),
    delta_aic     = AIC(m1) - AIC(m0),
    lrt_chisq     = lrt$Chisq[2],
    lrt_df        = lrt$Df[2],
    lrt_p         = lrt$`Pr(>Chisq)`[2]
  )
}

# Run the LRT for the three time-series responses and stack the summary rows.
genet_tests <- bind_rows(
  fit_and_report(pam,   "fv_fm",              "pam_fvfm"),
  fit_and_report(color, "color_num",          "color_dscale"),
  # Symbionts: each id has 1 obs per biopsy day; drop (1|id) since each coral
  # was destructively sampled at a single timepoint.
  fit_and_report(phys,  "log(cells_per_cm2)", "log_zoox",
                 fixed_extra = "biopsy_day_c", include_id = FALSE)
)
# Growth is handled separately because it has NO time dimension (one calcification
# rate per coral), so the null/genet pair drops the `day` terms entirely. Retain
# (1|tank) as the experimental block for the treatment assignment; same ML fit so
# the LRT is valid.
bw_a       <- bw |> filter(is.finite(pct_growth))
m_bw_null  <- lme4::lmer(
  pct_growth ~ treatment * wound + thicket + (1 | tank),
  data = bw_a, REML = FALSE,
  control = lme4::lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))
)
m_bw_genet <- lme4::lmer(
  pct_growth ~ treatment * wound * thicket + (1 | tank),
  data = bw_a, REML = FALSE,
  control = lme4::lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))
)
bw_lrt <- anova(m_bw_null, m_bw_genet)
# Guard for lme4 version differences: the LRT df column is sometimes "Chi Df",
# sometimes "Df" — pick whichever is present.
bw_lrt_df <- if ("Chi Df" %in% names(bw_lrt)) bw_lrt$`Chi Df`[2] else bw_lrt$Df[2]
genet_tests <- bind_rows(genet_tests, tibble(
  response  = "growth_pct",
  n_obs     = nrow(bw_a),
  aic_null  = AIC(m_bw_null),
  aic_genet = AIC(m_bw_genet),
  delta_aic = AIC(m_bw_genet) - AIC(m_bw_null),
  lrt_chisq = bw_lrt$Chisq[2],
  lrt_df    = bw_lrt_df,
  lrt_p     = bw_lrt$`Pr(>Chisq)`[2]
))
saveRDS(m_bw_genet, file.path(MOD_DIR, "13_growth_genet_lmm.rds"))

write_csv(genet_tests, file.path(TBL_DIR, "13_genet_anova.csv"))

# ---- Per-genet means at end of experiment ---------------------------------
# Reaction norms: response ~ treatment, one slope per genet
# make_react_norm(): take the last timepoint per response (the end-of-experiment
# state, where treatment divergence is largest) and compute the genet × treatment
# mean ± SE. These raw cell means are what the reaction-norm figure plots; they
# visualize the same G × E that the LRT above tests formally. Raw summaries, not
# model-adjusted emmeans, so the figure shows the data directly.
make_react_norm <- function(data, response_col, label) {
  last_day <- if ("day" %in% names(data)) {
    data |> filter(day == max(day, na.rm = TRUE))
  } else if ("biopsy_day" %in% names(data)) {
    data |> filter(biopsy_day == max(biopsy_day, na.rm = TRUE))
  } else {
    data
  }
  last_day |>
    group_by(treatment, thicket) |>
    summarise(mean = mean(.data[[response_col]], na.rm = TRUE),
              se   = sd(.data[[response_col]], na.rm = TRUE) / sqrt(n()),
              n    = n(),
              .groups = "drop") |>
    mutate(response = label)
}

# Assemble the four-panel reaction-norm table. PAM and color use the helper;
# growth and symbionts are summarized inline because growth has no time axis and
# symbionts are shown on a log10 scale for readability.
emm_long <- bind_rows(
  make_react_norm(pam,   "fv_fm",                "PAM Fv/Fm"),
  make_react_norm(color, "color_num",            "Color (D-scale)"),
  bw_a |> group_by(treatment, thicket) |>
    summarise(mean = mean(pct_growth, na.rm = TRUE),
              se   = sd(pct_growth, na.rm = TRUE) / sqrt(n()),
              n    = n(),
              .groups = "drop") |>
    mutate(response = "Growth (% mass change)"),
  phys |> mutate(log10_cells = log10(cells_per_cm2)) |>
    group_by(treatment, thicket) |>
    summarise(mean = mean(log10_cells, na.rm = TRUE),
              se   = sd(log10_cells, na.rm = TRUE) / sqrt(n()),
              n    = n(),
              .groups = "drop") |>
    mutate(response = "log10 symbionts cm⁻² (raw summary)")
)
write_csv(emm_long, file.path(TBL_DIR, "13_genet_emmeans.csv"))

# ---- Reaction-norm figure --------------------------------------------------
# One panel per response; x = temperature, one coloured line per genet. The
# slope of each line is that genet's heat response; differences in slope between
# genets are the visual indication of the G × E that the LRT tests.
p_react <- ggplot(emm_long, aes(treatment, mean,
                                 colour = thicket, group = thicket)) +
  geom_line(linewidth = 0.6) +
  geom_pointrange(aes(ymin = mean - se, ymax = mean + se),
                  size = 0.35, fatten = 2.5) +
  facet_wrap(~ response, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = PAL_GENO, name = "Genotype",
                      labels = c("a", "c", "d")) +
  labs(x = NULL, y = "End-of-experiment mean (± 1 SE)",
       title = "Genet-level reaction norms: 28 °C vs 31 °C",
       subtitle = "Parallel lines = additive temperature effect; crossing = G × E") +
  theme_pub(9)

save_fig(p_react, "13_genet_response_panel", width = 200, height = 95)

cat("\n=== Genet × treatment interaction tests ===\n")
print(genet_tests |> mutate(across(where(is.numeric), \(x) round(x, 3))))
