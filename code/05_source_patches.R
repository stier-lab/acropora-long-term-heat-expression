# =============================================================================
# 05_source_patches.R: does the heat effect differ among source patches?
#
# A source patch (A, C, D) is the field patch a fragment's parent branch came
# from. Patches are field labels, not confirmed genotypes (see 06_genetics.R).
#
# For each physiology response, the 03_physiology.R model is compared with the
# same model minus every term containing source patch x temperature, using a
# likelihood-ratio test (LRT: does the bigger model fit the data significantly
# better than the smaller one?). LRTs need ML fits, not REML; 8 df, 4 for
# growth. Until 2026-11-11 the reduced model also dropped patch x wound and
# patch x day terms, so the test was broader than its name; see
# docs/diagnostics/. The per-patch heat effects (ambient minus heated, Day 14)
# come from the 03_physiology.R models.
#
# Sections
#   1. Likelihood-ratio tests
#   2. End-of-experiment means by patch (figure)
#   3. Numbers for the results text
#
# Reads:  data/processed/{pam,color,buoyant_weight,symbiont_chl}_clean.rds,
#         output/tables/12_genet_treatment_effects.csv (from 03_physiology.R)
# Writes: output/tables/13_genet_anova.csv, 13_genet_emmeans.csv
#         output/models/13_*_genet_lmm.rds
#         figures/13_genet_response_panel
#         output/paper/source_patches_numbers.csv
# Run:    sourced by code/run_all.R (or: Rscript code/run_all.R)
# =============================================================================

source(here::here("code", "00_setup.R"))


# ---- 1. Likelihood-ratio tests ---------------------------------------------
# Does adding source patch x temperature terms improve the physiology models?
# In the data `thicket` is the source patch column.

## 1a. Data, filtered exactly as in 03_physiology.R --------------------------
pam   <- readRDS(file.path(DATA_PROC, "pam_clean.rds")) |>
  mutate(thicket = factor(thicket)) |>
  filter(day >= 0)   # exclude the pre-treatment baseline from the day-slope models
color <- readRDS(file.path(DATA_PROC, "color_clean.rds")) |>
  mutate(thicket = factor(thicket)) |>
  filter(day >= 0)   # same filter as the 03_physiology.R models (formerly 12_models)
bw    <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds")) |>
  mutate(thicket = factor(thicket))
phys  <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
  filter(is.finite(cells_per_cm2), cells_per_cm2 > 0) |>
  mutate(thicket = factor(thicket),
         biopsy_day_c = biopsy_day - 1)   # centred on Day 1, as in 03_physiology.R

## 1b. The test --------------------------------------------------------------
# Full model = the 03_physiology.R model (same random effects, ML fit).
# Reduced model = full minus every term that contains source patch x
# temperature: patch:temp, patch:temp:wound, patch:temp:day and
# patch:temp:wound:day (8 df; 4 for growth, which has no day term).
# Patch x wound and patch x day terms stay in both models, so the test is only
# about whether the heat effect differs among patches.
# Singular-fit messages are silenced (a zero variance estimate is expected for
# some tank intercepts); bobyqa is the optimiser.
ctrl <- lme4::lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4),
                          optimizer = "bobyqa")
lrt_patch_temp <- function(data, response, name, time = "day", re) {
  drop <- c("treatment:thicket", "treatment:wound:thicket")

  # Build the full formula as text; time = NA means no day term (growth).
  if (is.na(time)) {
    f1 <- as.formula(paste(response, "~ treatment * wound * thicket +", re))
  } else {
    f1 <- as.formula(paste(response, "~ treatment * wound *", time, "* thicket +", re))
    drop <- c(drop, paste0("treatment:", time, ":thicket"),
              paste0("treatment:wound:", time, ":thicket"))
  }

  # update(f1, . ~ . - a - b) = the same formula with those terms removed.
  f0 <- update(f1, as.formula(paste(". ~ . -", paste(drop, collapse = " - "))))
  m1 <- lme4::lmer(f1, data = data, REML = FALSE, control = ctrl)
  m0 <- lme4::lmer(f0, data = data, REML = FALSE, control = ctrl)
  saveRDS(m1, file.path(MOD_DIR, paste0("13_", name, "_genet_lmm.rds")))

  # anova(m0, m1) on two nested models = the likelihood-ratio test; row 2
  # holds the comparison.
  lrt <- anova(m0, m1)
  tibble(
    response  = name,
    test      = paste("LRT, drop", paste(drop, collapse = ", ")),
    random    = re,
    n_obs     = nrow(data),
    aic_null  = AIC(m0),
    aic_genet = AIC(m1),
    delta_aic = AIC(m1) - AIC(m0),
    lrt_chisq = lrt$Chisq[2],
    lrt_df    = lrt$Df[2],
    lrt_p     = lrt$`Pr(>Chisq)`[2]
  )
}

## 1c. Run it for the four responses -----------------------------------------
# Random effects match 03_physiology.R: `(1 + day || tank)` gives each tank its
# own intercept and day slope (`||` = uncorrelated), because temperature was
# applied to tanks; `(1 | id)` handles repeated measures on each fragment.
bw_a <- bw |> filter(is.finite(pct_growth))
genet_tests <- bind_rows(
  lrt_patch_temp(pam,   "fv_fm",              "pam_fvfm",
                 re = "(1 + day || tank) + (1 | id)"),
  lrt_patch_temp(color, "color_num",          "color_dscale",
                 re = "(1 + day || tank) + (1 | id)"),
  lrt_patch_temp(phys,  "log(cells_per_cm2)", "log_zoox", time = "biopsy_day_c",
                 re = "(1 + biopsy_day_c || tank)"),
  lrt_patch_temp(bw_a,  "pct_growth",         "growth_pct", time = NA,
                 re = "(1 | tank)")
)
write_csv(genet_tests, file.path(TBL_DIR, "13_genet_anova.csv"))


# ---- 2. End-of-experiment means by patch (figure) --------------------------
# Raw means (not model estimates) at the last measurement of each response,
# with SE across fragments. Lines that are not parallel mean the heat effect
# differs among patches (a "reaction norm" plot).

# Mean and SE by temperature x patch on the last measurement day. Uses `day`
# if present, else `biopsy_day`, else all rows.
make_react_norm <- function(data, response_col, label) {
  last_day <- if ("day" %in% names(data)) {
    data |> filter(day == max(day, na.rm = TRUE))
  } else if ("biopsy_day" %in% names(data)) {
    data |> filter(biopsy_day == max(biopsy_day, na.rm = TRUE))
  } else {
    data
  }

  # .data[[response_col]] = the column whose name is stored in response_col.
  last_day |>
    group_by(treatment, thicket) |>
    summarise(mean = mean(.data[[response_col]], na.rm = TRUE),
              se   = sd(.data[[response_col]], na.rm = TRUE) / sqrt(n()),
              n    = n(),
              .groups = "drop") |>
    mutate(response = label)
}

# Panel labels, in panel order (the \u escapes are a superscript 2 and an en dash).
RESP_LEVELS <- c("Fv/Fm (Day 14)", "Color score (Day 14)",
                 "Symbionts, log10 cells per cm\u00b2 (Day 15)",
                 "Growth, % mass gained (Days 0\u201315)")

# Growth has a single end-point measurement, so it is summarised directly.
emm_long <- bind_rows(
  make_react_norm(pam,   "fv_fm",     RESP_LEVELS[1]),
  make_react_norm(color, "color_num", RESP_LEVELS[2]),
  phys |>
    mutate(log10_cells = log10(cells_per_cm2)) |>
    make_react_norm("log10_cells", RESP_LEVELS[3]),
  bw_a |>
    group_by(treatment, thicket) |>
    summarise(mean = mean(pct_growth, na.rm = TRUE),
              se   = sd(pct_growth, na.rm = TRUE) / sqrt(n()),
              n    = n(),
              .groups = "drop") |>
    mutate(response = RESP_LEVELS[4])
) |>
  mutate(response = factor(response, levels = RESP_LEVELS))
write_csv(emm_long, file.path(TBL_DIR, "13_genet_emmeans.csv"))

p_react <- ggplot(emm_long, aes(treatment, mean, colour = thicket, group = thicket)) +
  geom_line(linewidth = 0.6) +
  geom_linerange(aes(ymin = mean - se, ymax = mean + se), linewidth = 0.5) +
  geom_point(size = 2) +
  facet_wrap(~ response, scales = "free_y", ncol = 4,
             labeller = label_wrap_gen(width = 24)) +
  scale_x_discrete(labels = c("28C" = "28 \u00b0C", "31C" = "31 \u00b0C")) +
  scale_colour_manual(values = PAL_GENO, name = "Source patch",
                      labels = c("A", "C", "D")) +
  labs(x = NULL, y = "Mean at end of experiment (\u00b1 1 SE)") +
  theme_pub(9)
save_fig(p_react, "13_genet_response_panel", width = 200, height = 90)


# ---- 3. Numbers for the results text ---------------------------------------
# Re-read the written tables so the reported numbers come from the same files
# a reader can check.
lrt <- readr::read_csv(file.path(TBL_DIR, "13_genet_anova.csv"), show_col_types = FALSE)
eff <- readr::read_csv(file.path(TBL_DIR, "12_genet_treatment_effects.csv"),
                       show_col_types = FALSE)

# One row per response: the source patch x temperature LRT.
labels <- c(pam_fvfm = "Fv/Fm", color_dscale = "Color", log_zoox = "Symbionts",
            growth_pct = "Growth")
rows <- list()
for (r in names(labels)) {
  d <- lrt[lrt$response == r, ]
  if (nrow(d) != 1) next
  rows[[length(rows) + 1]] <- paper_row(
    paste0("patch_x_temp_", r),
    paste0(labels[[r]], ": does the heat effect differ among source patches?"),
    test = "likelihood-ratio test of all source patch x temperature terms",
    statistic = d$lrt_chisq, df1 = d$lrt_df, p_value = d$lrt_p,
    n = d$n_obs, source = "output/tables/13_genet_anova.csv")
}

# Heat effect (ambient minus heated) within each patch on Day 14, averaged over
# the unwounded and wounded estimates from 03_physiology.R.
for (r in c("pam_fvfm", "color_dscale")) {
  for (p in c("a", "c", "d")) {
    d <- eff[eff$response == r & eff$thicket == p & eff$contrast == "28C - 31C", ]
    if (nrow(d) == 0) next
    rows[[length(rows) + 1]] <- paper_row(
      paste0("heat_drop_", r, "_patch_", toupper(p)),
      paste0(labels[[r]], ": ambient minus heated on Day 14, patch ", toupper(p),
             " (mean over wound status)"),
      estimate = mean(d$estimate), test = "model estimate (emmeans)",
      source = "output/tables/12_genet_treatment_effects.csv")
  }
}
write_paper(rows, "source_patches")
