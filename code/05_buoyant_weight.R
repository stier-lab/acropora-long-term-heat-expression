# =============================================================================
# Purpose: Coral growth from buoyant weight.
#
#          PRIMARY metric: areal calcification rate (mg CaCO3 cm^-2 d^-1) —
#          dry-mass gain normalized to the coral's calcifying SURFACE AREA and
#          to time. Calcification is a surface-mediated process (new aragonite
#          is deposited under the living tissue veneer), so surface area, not
#          skeletal mass, is the mechanistically correct denominator, and
#          mg cm^-2 d^-1 is the field-standard unit (Jokiel et al. 1978).
#
#          Surface area comes from the day-15 wax-dipping standard-curve SA
#          (data/processed/wax_clean.rds). The 48 growth corals were biopsied
#          only terminally (all biopsy_day == 15), so day-15 SA is the SA at the
#          end of a clean 15-day growth window; SA is stable across biopsy days
#          and does not differ by treatment (28C 4.46 vs 31C 4.73 cm2, p=0.27),
#          so it is an unbiased denominator. See notes/growth_allometry.md.
#
#          Reported alongside for robustness: % mass change (delta_M / M_init)
#          and specific growth rate (SGR, % d^-1). The heat effect is the same
#          under all three; see output/tables/05b_growth_metric_comparison.csv.
#
#          Buoyant-weight algorithm: Davies 1989 / Jokiel et al. 1978; the raw
#          sheet stores Excel formulas, so we recompute the dry-mass conversion.
# Input:   data/raw/buoyant_weight/data.csv
#          data/processed/wax_clean.rds  (day-15 surface area)
# Output:  data/processed/buoyant_weight_clean.rds
#          figures/05_buoyant_weight_growth.{pdf,png}
#          output/tables/05_buoyant_weight_lm.csv
#          output/tables/05_buoyant_weight_tank_test.csv
#          output/tables/05b_growth_metric_comparison.csv
# =============================================================================

source(here::here("code", "00_setup.R"))

raw <- read_csv(file.path(DATA_RAW, "buoyant_weight", "data.csv"),
                show_col_types = FALSE) |>
  janitor::clean_names()

# Reference (calibration) measurements — these don't change between weighings
# so we use the row-level reference columns as supplied.
bw <- raw |>
  rename(thicket = matches("^thicket")) |>
  mutate(
    thicket    = str_to_lower(str_squish(thicket)),
    id         = as.integer(id),
    wound      = factor(wound, levels = c("no", "yes")),
    tank       = as.integer(tank),
    treatment  = factor(as.integer(treatment), levels = c(28, 31),
                        labels = c("28C", "31C")),
    biopsy_day = as.integer(biopsy_day),
    initial_w  = as.numeric(initial_coral_buoyant_weigh_uncorrected),
    final_w    = as.numeric(final_coral_weight_uncorrected),
    initial_temp = as.numeric(initial_temp_coral),
    final_temp   = as.numeric(final_temp_coral)
  ) |>
  filter(!is.na(initial_w), !is.na(final_w))

# Davies 1989 buoyant weight conversion:
# m_dry = m_buoy / (1 - rho_sw / rho_arag)
# Using rho_arag = 2.93 g/cm^3, rho_sw computed from temp + salinity reference.
# Here we approximate rho_sw at the recorded temperatures using a polynomial
# from the raw sheet's formula `-5e-6*T^2 + 7e-6*T + 1.0001` — same constants
# the field team used. Salinity correction omitted because S_field == S_ref.
rho_sw_temp <- function(t_c) -5e-6 * t_c^2 + 7e-6 * t_c + 1.0001
ARAG <- 2.93

bw <- bw |>
  mutate(
    initial_dry = initial_w / (1 - rho_sw_temp(initial_temp) / ARAG),
    final_dry   = final_w   / (1 - rho_sw_temp(final_temp)   / ARAG),
    delta_g     = final_dry - initial_dry,
    days        = pmax(biopsy_day, 1),
    pct_growth  = 100 * delta_g / initial_dry,                 # robustness metric
    sgr         = 100 * (log(final_dry) - log(initial_dry)) / days,  # % d^-1
    g_per_day   = delta_g / days
  )

# ---- Join surface area & compute areal calcification (PRIMARY) -------------
# Day-15 wax standard-curve SA, the calcifying-surface denominator.
sa <- readRDS(file.path(DATA_PROC, "wax_clean.rds")) |>
  filter(biopsy_day == 15, is.finite(sa_curve_cm2), sa_curve_cm2 > 0) |>
  distinct(id, .keep_all = TRUE) |>
  select(id, sa_cm2 = sa_curve_cm2)

bw <- bw |>
  left_join(sa, by = "id") |>
  mutate(
    # mg CaCO3 per cm^2 per day:  (delta_g * 1000 mg/g) / SA / days
    areal_calc = 1000 * delta_g / sa_cm2 / days
  )

n_missing_sa <- sum(is.na(bw$sa_cm2))
if (n_missing_sa > 0)
  message("NOTE: ", n_missing_sa, " growth corals lack a day-15 SA; ",
          "areal_calc is NA for these (pct_growth/sgr still available).")

saveRDS(bw, file.path(DATA_PROC, "buoyant_weight_clean.rds"))

# ---- Primary model: areal calcification ------------------------------------
bw_a <- bw |> filter(is.finite(areal_calc))
# Coral-level coefficients are useful effect-size summaries, but temperature is
# assigned at the tank level. The inferential heat-effect p-value below is
# therefore a tank-level permutation test.
m_growth <- lm(areal_calc ~ treatment * wound + thicket, data = bw_a)
res <- broom::tidy(m_growth, conf.int = TRUE) |>
  mutate(model_scope = "coral-level descriptive coefficients",
         across(where(is.numeric), \(x) round(x, 4))) |>
  relocate(model_scope, .after = term)
write_csv(res, file.path(TBL_DIR, "05_buoyant_weight_lm.csv"))

tank_growth <- bw_a |>
  group_by(tank, treatment) |>
  summarise(mean_areal_calc = mean(areal_calc, na.rm = TRUE),
            n_corals = n(), .groups = "drop")
obs_diff <- with(tank_growth,
                 mean(mean_areal_calc[treatment == "28C"]) -
                   mean(mean_areal_calc[treatment == "31C"]))
vals <- tank_growth$mean_areal_calc
assignments <- combn(seq_along(vals), sum(tank_growth$treatment == "28C"),
                     simplify = FALSE)
perm_diffs <- vapply(assignments, function(idx) {
  mean(vals[idx]) - mean(vals[-idx])
}, numeric(1))
tank_test <- tibble(
  response = "areal_calc",
  test = "tank-level exact permutation",
  n_tanks = nrow(tank_growth),
  n_28_tanks = sum(tank_growth$treatment == "28C"),
  n_31_tanks = sum(tank_growth$treatment == "31C"),
  estimate_28_minus_31 = obs_diff,
  p_two_sided = mean(abs(perm_diffs) >= abs(obs_diff)),
  note = "Temperature was randomized at the tank level; this p-value treats tanks, not corals, as exchangeable units."
) |>
  mutate(across(where(is.numeric), \(x) signif(x, 4)))
write_csv(tank_test, file.path(TBL_DIR, "05_buoyant_weight_tank_test.csv"))

# ---- Metric comparison (heat effect is robust to metric choice) ------------
metric_compare <- purrr::map_dfr(
  c(areal_calc = "areal_calc", pct_growth = "pct_growth", sgr = "sgr"),
  function(v) {
    d <- bw |> filter(is.finite(.data[[v]]))
    a <- car::Anova(lm(reformulate("treatment * wound + thicket", v), data = d),
                    type = 2)
    tibble(
      term      = rownames(a)[1:3],
      F_value   = round(a[1:3, "F value"], 2),
      p_value   = signif(a[1:3, "Pr(>F)"], 3)
    )
  }, .id = "metric"
)
write_csv(metric_compare, file.path(TBL_DIR, "05b_growth_metric_comparison.csv"))

# ---- Plot: areal calcification ---------------------------------------------
p_bw <- ggplot(bw_a, aes(interaction(treatment, wound, sep = " · "),
                         areal_calc, fill = treatment)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(colour = thicket), width = 0.18, height = 0,
              size = 1.6, alpha = 0.85) +
  scale_fill_manual(values  = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                    name = "Temperature") +
  scale_colour_manual(values = PAL_GENO, name = "Genotype") +
  labs(x = NULL,
       y = expression("Calcification ("*mg~cm^{-2}~d^{-1}*")"),
       title = "Coral calcification under heating × wounding",
       subtitle = "Buoyant-weight mass gain per unit surface area (Jokiel et al. 1978)") +
  theme_pub(10) +
  theme(axis.text.x = element_text(angle = 0))

save_fig(p_bw, "05_buoyant_weight_growth", width = 150, height = 95)

cat("\n=== Areal calcification LM (coral-level descriptive) ===\n")
print(res)
cat("\n=== Areal calcification tank-level temperature test ===\n")
print(tank_test)
cat("\n=== Heat effect is robust to metric choice ===\n")
print(metric_compare |> filter(term == "treatment"))
cat("\nWrote: buoyant_weight_clean.rds, 05_buoyant_weight_lm.csv,",
    "05_buoyant_weight_tank_test.csv,",
    "05b_growth_metric_comparison.csv, 05_buoyant_weight_growth.{pdf,png}\n")
