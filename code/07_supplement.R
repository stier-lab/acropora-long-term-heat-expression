# =============================================================================
# 07_supplement.R: supplement-only results
#
#   1 Microscope experiment: share of the 16 fragments showing each stage on
#     each day, by temperature (supplement figure).
#   2 Results without the flagged fragments: refit the Fv/Fm and color
#     temperature x day tests and the growth tank permutation without
#     fragments 116 and 121 and without tank 3 (slow to equilibrate).
#   3 Color: tank-level permutation check (no distribution assumed).
#
# Other supplement material is written by the main scripts: the growth tank
# permutation (03_physiology.R), Kaplan-Meier curves and Cox models
# (04_regeneration.R), and the wax standard curve (01_clean_data.R).
#
# Writes: output/tables/11_microscope_trait_trajectory.csv, 22_sensitivity_flagged.csv,
#         12b_color_tank_test.csv
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

# ---- 3. Color: tank-level check ------------------------------------------
# Color is a 1-6 card score and ambient fragments almost all stayed at 4, so
# the Gaussian model's residuals are not normal (18_model_checks.csv) and an
# ordered-category model cannot be estimated (near-complete separation). The
# check that needs no distribution: compare final mean color of the 4 ambient
# and 4 heated tanks against all 70 ways of splitting 8 tanks into two groups
# of 4, as for growth (03_physiology.R) and budding (04_regeneration.R).
color_final <- readRDS(file.path(DATA_PROC, "color_clean.rds")) |>
  filter(day == max(day, na.rm = TRUE))
tank_color <- color_final |>
  group_by(tank, treatment) |>
  summarise(mean_color = mean(color_num), n_corals = n(),
            n_paled = sum(color_num < 4), .groups = "drop")
obs_diff <- with(tank_color, mean(mean_color[treatment == "28C"]) - mean(mean_color[treatment == "31C"]))
vals <- tank_color$mean_color
assignments <- combn(seq_along(vals), sum(tank_color$treatment == "28C"), simplify = FALSE)
perm_diffs <- vapply(assignments, function(idx) mean(vals[idx]) - mean(vals[-idx]), numeric(1))
color_tank_test <- tibble(
  response = "color_dscale",
  day = unique(color_final$day),
  test = "tank-level exact permutation",
  n_tanks = nrow(tank_color),
  estimate_28_minus_31 = obs_diff,
  p_two_sided = mean(abs(perm_diffs) >= abs(obs_diff)),
  paled_per_tank_28C = paste(tank_color$n_paled[tank_color$treatment == "28C"], collapse = ", "),
  paled_per_tank_31C = paste(tank_color$n_paled[tank_color$treatment == "31C"], collapse = ", "),
  note = "Paled = card score below 4. 2/70 = 0.029 is the smallest possible two-sided p with 4 tanks per temperature."
) |>
  mutate(across(where(is.double), \(x) signif(x, 4)))
write_csv(color_tank_test, file.path(TBL_DIR, "12b_color_tank_test.csv"))
print(as.data.frame(color_tank_test))
