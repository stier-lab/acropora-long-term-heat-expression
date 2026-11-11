# =============================================================================
# 05_source_patches.R: does the heat effect differ among source patches?
#
# A source patch (A, C, D) is the field patch a fragment's parent branch came
# from. Patches are field labels, not confirmed genotypes (see 06_genetics.R).
#
# For each physiology response, a model with all source patch x temperature
# terms is compared with one without them (likelihood-ratio test, ML fits).
# Fixed effects otherwise match 03_physiology.R. The per-patch heat effects
# (ambient minus heated, Day 14) come from the 03_physiology.R models.
#
# Writes: output/tables/13_genet_anova.csv, 13_genet_emmeans.csv
#         figures/13_genet_response_panel
#         output/paper/source_patches_numbers.csv
# =============================================================================

source(here::here("code", "00_setup.R"))

# ---- 1. Likelihood-ratio tests ---------------------------------------------
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
fit_and_report <- function(data, response, name, fixed_extra = NULL,
                           include_id = TRUE) {
  id_term <- if (include_id) " + (1 | id)" else ""
  rhs_null  <- paste0("treatment * wound * day + thicket + (1 | tank)",
                      id_term)
  rhs_genet <- paste0("treatment * wound * day * thicket + (1 | tank)",
                      id_term)
  if (!is.null(fixed_extra)) {
    rhs_null  <- str_replace(rhs_null,  "day", fixed_extra)
    rhs_genet <- str_replace(rhs_genet, "day", fixed_extra)
  }
  f0 <- as.formula(paste(response, "~", rhs_null))
  f1 <- as.formula(paste(response, "~", rhs_genet))
  m0 <- suppressWarnings(suppressMessages(
    lme4::lmer(f0, data = data, REML = FALSE)))
  m1 <- suppressWarnings(suppressMessages(
    lme4::lmer(f1, data = data, REML = FALSE)))
  saveRDS(m1, file.path(MOD_DIR, paste0("13_", name, "_genet_lmm.rds")))
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
genet_tests <- bind_rows(
  fit_and_report(pam,   "fv_fm",              "pam_fvfm"),
  fit_and_report(color, "color_num",          "color_dscale"),
  fit_and_report(phys,  "log(cells_per_cm2)", "log_zoox",
                 fixed_extra = "biopsy_day_c", include_id = FALSE)
)
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

# ---- 2. End-of-experiment means by patch (figure) --------------------------
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

# ---- 3. Numbers for the results text ---------------------------------------
lrt <- readr::read_csv(file.path(TBL_DIR, "13_genet_anova.csv"), show_col_types = FALSE)
eff <- readr::read_csv(file.path(TBL_DIR, "12_genet_treatment_effects.csv"), show_col_types = FALSE)

labels <- c(pam_fvfm = "Fv/Fm", color_dscale = "Color", log_zoox = "Symbionts",
            growth_pct = "Growth")
rows <- list()
for (r in names(labels)) {
  d <- lrt[lrt$response == r, ]
  if (nrow(d) != 1) next
  rows[[length(rows) + 1]] <- paper_row(
    paste0("patch_x_temp_", r), paste0(labels[[r]], ": does the heat effect differ among source patches?"),
    test = "likelihood-ratio test", statistic = d$lrt_chisq, df1 = d$lrt_df, p_value = d$lrt_p,
    n = d$n_obs, source = "output/tables/13_genet_anova.csv")
}
for (r in c("pam_fvfm", "color_dscale")) {
  for (p in c("a", "c", "d")) {
    d <- eff[eff$response == r & eff$thicket == p & eff$contrast == "28C - 31C", ]
    if (nrow(d) == 0) next
    rows[[length(rows) + 1]] <- paper_row(
      paste0("heat_drop_", r, "_patch_", toupper(p)),
      paste0(labels[[r]], ": ambient minus heated on Day 14, patch ", toupper(p), " (mean over wound status)"),
      estimate = mean(d$estimate), test = "model estimate (emmeans)",
      source = "output/tables/12_genet_treatment_effects.csv")
  }
}
write_paper(rows, "source_patches")
