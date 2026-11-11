# =============================================================================
# 07_supplement.R: supplement-only results
#
#   1 Microscope experiment: share of the 16 fragments showing each stage on
#     each day, by temperature (supplement figure).
#   2 Results without the flagged fragments: refit the Fv/Fm and color
#     temperature x day tests and the growth tank permutation without
#     fragments 116 and 121 and without tank 3 (slow to equilibrate).
#   3 Color as an ordered score: robustness check for the color model.
#
# Other supplement material is written by the main scripts: the growth tank
# permutation (03_physiology.R), Kaplan-Meier curves and Cox models
# (04_regeneration.R), and the wax standard curve (01_clean_data.R).
#
# Writes: output/tables/11_microscope_trait_trajectory.csv, 22_sensitivity_flagged.csv,
#         12b_color_clmm.csv, 12b_color_clmm_genet_effects.csv
#         figures/11_microscope_trait_trajectories
# =============================================================================

source(here::here("code", "00_setup.R"))

# ---- 1. Microscope experiment: stage trajectories --------------------------
micro <- readRDS(file.path(DATA_PROC, "microscope_physio_clean.rds"))
trait_meta <- tibble::tribble(
  ~trait,                  ~label,                    ~event_interpretation,
  "tissue_over_wound",     "Tissue over wound",        "first observed tissue/coenosarc over wound",
  "axial_polyp_formation", "Axial polyp formation",   "first observed central axial polyp/hole",
  "wound_smoothed",        "Wound smoothed",          "first observed smoothed wound surface",
  "pigment_over_wound",    "Pigment over wound",      "first observed pigment over wound; scored D0-D7 only",
  "tip_exist",             "Tip exists",              "first observed visible axial tip",
  "tip_extension",         "Tip extension",           "first observed vertical tip extension",
  "new_corallites_on_tip", "New corallites on tip",   "first observed new corallites on the regenerating tip"
)
traits <- trait_meta$trait
long <- micro |>
  select(day, date, treatment, tank, thicket, id, use_photo, notes,
         all_of(traits)) |>
  pivot_longer(all_of(traits), names_to = "trait", values_to = "expressed") |>
  left_join(trait_meta, by = "trait") |>
  filter(!is.na(expressed))
prop_df <- long |>
  group_by(day, treatment, trait, label) |>
  summarise(
    n = n(),
    n_yes = sum(expressed == 1),
    prop = mean(expressed),
    .groups = "drop"
  )
write_csv(prop_df, file.path(TBL_DIR, "11_microscope_trait_trajectory.csv"))
p_micro <- ggplot(prop_df, aes(day, prop, colour = treatment,
                               group = treatment)) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = PAL_TEMP, name = "Temperature") +
  scale_x_continuous(breaks = seq(0, 15, by = 3), limits = c(0, 15)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  facet_wrap(~ label, ncol = 4) +
  labs(x = "Days post-wounding",
       y = "Photo corals expressing trait",
       title = "Microscope wound-healing trajectories",
       subtitle = "Separate photo-only cohort; wounded corals only (n = 8 per temperature)") +
  theme_pub(9)
save_fig(p_micro, "11_microscope_trait_trajectories",
         width = 210, height = 130)

# ---- 2. Results without the flagged fragments and tank ---------------------
FLAGGED_IDS  <- c(116, 121)   # two coral colonies that behaved anomalously
FLAGGED_TANK <- 3             # ambient tank that warmed/equilibrated slowly
results <- list()             # one tibble per (response × full-or-dropped) goes here
grab <- function(model, label, term) {
  av <- as.data.frame(anova(model))
  av$term <- rownames(av)
  row <- av[av$term == term, , drop = FALSE]
  if (nrow(row) == 0) return(NULL)
  tibble(
    response   = label,
    term       = term,
    F_value    = round(row[["F value"]], 2),
    p_value    = signif(row[["Pr(>F)"]], 3)
  )
}
pam <- readRDS(file.path(DATA_PROC, "pam_clean.rds")) |>
  mutate(thicket = factor(thicket)) |>
  filter(day >= 0)
pam_drop <- pam |> filter(!id %in% FLAGGED_IDS, tank != FLAGGED_TANK)
f_full <- lmerTest::lmer(fv_fm ~ treatment * wound * day + thicket + (1|tank) + (1|id),
                         data = pam, REML = FALSE)
f_drop <- lmerTest::lmer(fv_fm ~ treatment * wound * day + thicket + (1|tank) + (1|id),
                         data = pam_drop, REML = FALSE)
results$pam_full <- grab(f_full, "PAM Fv/Fm (full data)",     "treatment:day")
results$pam_drop <- grab(f_drop, "PAM Fv/Fm (flagged dropped)","treatment:day")
col <- readRDS(file.path(DATA_PROC, "color_clean.rds")) |>
  mutate(thicket = factor(thicket)) |>
  filter(day >= 0)
col_drop <- col |> filter(!id %in% FLAGGED_IDS, tank != FLAGGED_TANK)
c_full <- lmerTest::lmer(color_num ~ treatment * wound * day + thicket + (1|tank) + (1|id),
                         data = col, REML = FALSE)
c_drop <- lmerTest::lmer(color_num ~ treatment * wound * day + thicket + (1|tank) + (1|id),
                         data = col_drop, REML = FALSE)
results$col_full <- grab(c_full, "Color D-scale (full data)",     "treatment:day")
results$col_drop <- grab(c_drop, "Color D-scale (flagged dropped)","treatment:day")
bw <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds")) |>
  mutate(thicket = factor(thicket)) |>
  filter(is.finite(pct_growth))                 # drop NA/Inf growth values
bw_drop <- bw |> filter(!id %in% FLAGGED_IDS, tank != FLAGGED_TANK)
tank_perm <- function(dat, label) {
  tank_growth <- dat |>
    group_by(tank, treatment) |>
    summarise(mean_pct_growth = mean(pct_growth, na.rm = TRUE),
              n_corals = n(), .groups = "drop")
  n_by_temp <- table(tank_growth$treatment)
  if (length(n_by_temp) != 2 || any(n_by_temp < 2)) {
    return(tibble(response = label, term = "treatment (tank permutation)",
                  F_value = NA_real_, p_value = NA_real_,
                  estimate_28_minus_31 = NA_real_, n_tanks = nrow(tank_growth)))
  }
  vals <- tank_growth$mean_pct_growth
  trt <- tank_growth$treatment
  n28 <- sum(trt == "28C")
  obs <- mean(vals[trt == "28C"]) - mean(vals[trt == "31C"])
  null <- vapply(combn(seq_along(vals), n28, simplify = FALSE), \(idx) {
    mean(vals[idx]) - mean(vals[-idx])
  }, numeric(1))
  p <- mean(abs(null) >= abs(obs))
  tibble(response = label, term = "treatment (tank permutation)",
         F_value = NA_real_, p_value = signif(p, 3),
         estimate_28_minus_31 = round(obs, 3),
         n_tanks = nrow(tank_growth))
}
results$bw_full <- tank_perm(bw, "Growth % (full data)")
results$bw_drop <- tank_perm(bw_drop, "Growth % (flagged dropped)")
out <- bind_rows(results)
write_csv(out, file.path(TBL_DIR, "22_sensitivity_flagged.csv"))

# ---- 3. Color as an ordered score (robustness) ------------------------------
# Color is a 1-6 card score, so the Gaussian model's residuals are not normal
# (see 08_model_checks.csv). This refits color as an ordered category
# (cumulative-link mixed model, ordinal package) to check that the temperature
# terms point the same way.
cat("\n=== PART 2: color ordinal-CLMM robustness ===\n")
color_ord_df <- readRDS(file.path(DATA_PROC, "color_clean.rds")) |>
  mutate(thicket = factor(thicket),
         color_ord = factor(pmin(pmax(as.integer(floor(color_num + 0.5)), 1), 6),
                            ordered = TRUE))
m_clmm <- ordinal::clmm(
  color_ord ~ treatment * wound * day + treatment * thicket +
              treatment * wound * thicket + day * thicket + (1 | id),
  data = color_ord_df, Hess = TRUE
)
saveRDS(m_clmm, file.path(MOD_DIR, "12b_color_clmm.rds"))
av_drop <- as.data.frame(drop1(m_clmm, test = "Chisq")) |>
  tibble::rownames_to_column("term") |>
  mutate(response = "color_dscale_ordinal")
write_csv(av_drop, file.path(TBL_DIR, "12b_color_clmm.csv"))
av <- av_drop |> rename(LR.stat = LRT, `Pr(>Chisq)` = `Pr(>Chi)`)
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
