# =============================================================================
# Purpose: Wax-dipping surface-area calibration and per-coral SA estimates.
#          - Fit the standard curve: wax mass (g) -> SA (cm^2) from cylinders
#            of known diameter and height.
#          - Re-derive every coral's SA from its (wax_2 - wax_1) and the curve.
#          - Cross-check with caliper-based SA = 2π(d/2)·h.
#
# What & why: many coral metrics (calcification, symbiont density) must be
#   normalized to living surface area, but coral skeletons are too irregular to
#   measure with calipers alone. Wax dipping is the field-standard fix (Stimson &
#   Kinzie 1991): dip the object in molten wax, and the wax mass that sticks is
#   proportional to its surface area. We build that proportionality once, from
#   smooth CYLINDERS whose true SA we can compute geometrically (2π·r·h), by
#   regressing known SA on wax mass — this is the "standard curve". We then dip
#   every coral, take the wax it picked up (wax_2 - wax_1), and read its SA off
#   that fitted line. A caliper-based SA (treating each coral as a cylinder) is
#   computed in parallel as a sanity check. The resulting per-coral SA feeds the
#   calcification (code/05) and symbiont (code/06) denominators.
# Input:   data/raw/wax_dipping/data.csv
#          data/raw/wax_dipping/Standard_curve.csv
# Output:  data/processed/wax_clean.rds
#          figures/07_wax_standard_curve.{pdf,png}
#          output/tables/07_wax_curve_fit.csv
# =============================================================================

# 00_setup.R loads packages, shared paths (DATA_RAW, DATA_PROC, TBL_DIR),
# theme_pub(), and save_fig().
source(here::here("code", "00_setup.R"))

# ---- Load ------------------------------------------------------------------
# sc = calibration cylinders (known geometry); wax = the experimental corals.
# suppressWarnings + guess_max absorb the mixed/blank cells in these exports;
# clean_names() makes every header snake_case.
sc <- suppressWarnings(
  read_csv(file.path(DATA_RAW, "wax_dipping", "Standard_curve.csv"),
           show_col_types = FALSE, guess_max = 2000)
) |>
  janitor::clean_names()

wax <- suppressWarnings(
  read_csv(file.path(DATA_RAW, "wax_dipping", "data.csv"),
           show_col_types = FALSE, guess_max = 2000)
) |>
  janitor::clean_names() |>
  filter(!is.na(species) & !is.na(id))   # drop blank/trailing spreadsheet rows

# ---- Standard curve --------------------------------------------------------
# Sheet columns: id, height_mm, diameter_mm, surface_area_mm (formula),
# dry_weight_g, wax_weight_g, wax_minus_dry_g, ..., wax_dry, surface_area
# The trailing `wax_dry` / `surface_area` pair is the values already evaluated.
sc_clean <- sc |>
  mutate(
    height   = as.numeric(height_mm),
    diameter = as.numeric(diameter_mm),
    dry_g    = as.numeric(dry_weight_g),
    wax_g    = as.numeric(wax_weight_g),
    # `wax_dry_g` is a spreadsheet formula string; compute directly.
    # Wax that adhered = dipped weight minus the bare dry weight.
    wax_minus_dry = wax_g - dry_g,
    # Lateral surface area of a cylinder = circumference x height = 2π(d/2)·h.
    # d and h are in mm so the product is mm^2; /100 converts mm^2 -> cm^2.
    sa_cm2   = 2 * pi * (diameter / 2) * height / 100
  ) |>
  filter(!is.na(wax_minus_dry), !is.na(sa_cm2),
         wax_minus_dry > 0, sa_cm2 > 0)   # keep only physically valid calibration points

fit_df <- sc_clean |>
  transmute(wax_g = wax_minus_dry, sa = sa_cm2)   # predictor = wax mass, response = true SA

# The calibration line: true SA regressed on wax mass. Its slope/intercept are
# what convert a coral's wax pickup into an SA estimate below.
cal <- lm(sa ~ wax_g, data = fit_df)
fit_tbl <- broom::tidy(cal, conf.int = TRUE)   # tidy coefficient table with 95% CIs
write_csv(fit_tbl, file.path(TBL_DIR, "07_wax_curve_fit.csv"))

# ---- Apply to corals -------------------------------------------------------
wax_proc <- wax |>
  rename(thicket = matches("^thicket")) |>
  mutate(
    id           = as.integer(id),
    treatment    = factor(as.integer(treatment), levels = c(28, 31),
                          labels = c("28C", "31C")),    # baseline = 28C ambient
    biopsy_day   = as.integer(biopsy_day),
    diameter_mm  = as.numeric(diameter),
    height_mm    = as.numeric(height),
    dry_g        = as.numeric(dry_weight),
    wax1_g       = as.numeric(wax_weight_1),   # coral weight before dipping
    wax2_g       = as.numeric(wax_weight_2),   # coral weight after dipping
    # wax_2_minus_wax_1 column is a string formula in the export; recompute
    wax_g        = wax2_g - wax1_g,            # wax the coral picked up
    # Two SA estimates per coral: the caliper/cylinder geometry (cross-check)...
    sa_caliper_cm2 = 2 * pi * (diameter_mm / 2) * height_mm / 100,
    # ...and the wax-curve SA, read off the fitted line at this coral's wax_g.
    # This curve-based SA is the one used downstream (code/05, code/06).
    sa_curve_cm2   = predict(cal, newdata = tibble(wax_g = wax_g))
  ) |>
  select(id, treatment, biopsy_day, thicket, wound, dry_g, wax_g,
         sa_caliper_cm2, sa_curve_cm2)

# Per-coral SA table consumed by the calcification and symbiont scripts.
saveRDS(wax_proc, file.path(DATA_PROC, "wax_clean.rds"))

# ---- Plot the standard curve ----------------------------------------------
# Calibration points with the fitted line + CI band; subtitle reports the fitted
# equation, R^2 and n so the curve's quality is visible on the figure itself.
p_sc <- ggplot(fit_df, aes(wax_g, sa)) +
  geom_point(size = 2, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, colour = "#0072B2") +   # same lm as `cal`
  labs(x = "Wax mass (g)", y = expression(Surface~area~(cm^2)),
       title = "Wax-dipping standard curve",
       subtitle = sprintf("SA = %.1f + %.1f · wax  (R² = %.3f, n = %d)",
                          coef(cal)[1], coef(cal)[2],
                          summary(cal)$r.squared, nobs(cal))) +
  theme_pub(10)

save_fig(p_sc, "07_wax_standard_curve", width = 120, height = 95)

cat("\n=== Wax dipping ===\n")
print(fit_tbl)
