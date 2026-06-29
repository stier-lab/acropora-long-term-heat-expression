# =============================================================================
# Purpose: Coral growth from buoyant weight.
#
#          PRIMARY metric: % mass change over the 15-day window
#          (100 · (dry_final - dry_initial) / dry_initial). Reported alongside:
#          specific growth rate (SGR, % d^-1). The heat effect is
#          the same under both metrics; see output/tables/05b_growth_metric_comparison.csv.
#
#          NOTE on areal calcification: an areal rate (mg CaCO3 cm^-2 d^-1) is NOT
#          reported. It would require the calcifying SURFACE AREA of the whole
#          weighed fragment, but the wax-dipping SA we have is for the small
#          sub-fragment taken for the symbiont slurry — the rest of each fragment
#          was chopped for transcriptomics, so no whole-fragment SA exists (M.
#          Brzezinski, pers. comm., 2026). Dividing whole-fragment mass gain by a
#          sub-fragment's SA mixes two objects, so growth is expressed as % mass
#          change (and SGR), which need no surface-area normalization. (The same
#          wax SA IS valid for symbiont density in code/06, where the cell count
#          and the SA come from the same sub-fragment.)
#
#          Buoyant-weight algorithm: Davies 1989 / Jokiel et al. 1978; the raw
#          sheet stores Excel formulas, so we recompute the dry-mass conversion.
#
# What & why: buoyant weighing is the non-destructive way to track how much new
#   skeleton a coral lays down. A coral is weighed while submerged at the start
#   and end of the experiment; because seawater buoys the skeleton, the buoyant
#   weight can be converted to dry skeletal (aragonite) mass with Archimedes'
#   principle (Davies 1989). The mass gained over the 15-day window, expressed as
#   a percentage of the coral's starting mass, is the growth metric. This script
#   estimates and tests the main result: heating reduces growth by ~34%. The
#   design has one consequence: temperature was applied to
#   whole TANKS, not to individual corals, so the coral-level linear model gives
#   valid effect-size estimates but not a valid p-value (corals in a
#   tank are not independent). The inferential test for the heat effect is
#   therefore a tank-level permutation test that treats the 8 tanks (not the 48
#   corals) as the units that were randomized.
# Input:   data/raw/buoyant_weight/data.csv
# Output:  data/processed/buoyant_weight_clean.rds
#          figures/05_buoyant_weight_growth.{pdf,png}
#          output/tables/05_buoyant_weight_lm.csv
#          output/tables/05_buoyant_weight_tank_test.csv
#          output/tables/05b_growth_metric_comparison.csv
# =============================================================================

# 00_setup.R loads packages (tidyverse etc.), defines shared paths (DATA_RAW,
# DATA_PROC, TBL_DIR, the PAL_GENO genotype palette), theme_pub(), and save_fig().
source(here::here("code", "00_setup.R"))

# ---- Load ------------------------------------------------------------------
# clean_names() rewrites the spreadsheet headers as snake_case so the long raw
# column names below ("initial_coral_buoyant_weigh_uncorrected", etc.) are stable.
raw <- read_csv(file.path(DATA_RAW, "buoyant_weight", "data.csv"),
                show_col_types = FALSE) |>
  janitor::clean_names()

# ---- Clean + type columns --------------------------------------------------
# Coerce each field to its proper type and set the experimental factors. Factor
# LEVEL ORDER matters: it fixes the model's reference (baseline) group. wound's
# baseline is "no" and treatment's baseline is 28C, so model coefficients read as
# the effect of wounding and of heating relative to the unstressed control.
bw <- raw |>
  rename(thicket = matches("^thicket")) |>
  mutate(
    # "thicket" is the source colony / genet (A, C, D); lower-case + squish
    # whitespace so spelling variants collapse to one label per genet.
    thicket    = str_to_lower(str_squish(thicket)),
    id         = as.integer(id),
    wound      = factor(wound, levels = c("no", "yes")),     # baseline = unwounded
    tank       = as.integer(tank),
    treatment  = factor(as.integer(treatment), levels = c(28, 31),
                        labels = c("28C", "31C")),           # baseline = 28C ambient
    biopsy_day = as.integer(biopsy_day),
    initial_w  = as.numeric(initial_coral_buoyant_weigh_uncorrected),  # buoyant g, t0
    final_w    = as.numeric(final_coral_weight_uncorrected),           # buoyant g, t-end
    initial_temp = as.numeric(initial_temp_coral),   # water temp at each weighing,
    final_temp   = as.numeric(final_temp_coral)      # needed for the density correction
  ) |>
  filter(!is.na(initial_w), !is.na(final_w))          # need both weights to get a gain

# ---- Buoyant weight -> dry skeletal mass -----------------------------------
# Davies 1989 buoyant weight conversion:
# m_dry = m_buoy / (1 - rho_sw / rho_arag)
# Archimedes: a submerged skeleton weighs less by the mass of water it displaces,
# so dividing by (1 - rho_seawater/rho_aragonite) scales the buoyant weight back
# up to true dry mass. Using rho_arag = 2.93 g/cm^3, rho_sw computed from temp +
# salinity reference. Here we approximate rho_sw at the recorded temperatures
# using a polynomial from the raw sheet's formula `-5e-6*T^2 + 7e-6*T + 1.0001`
# — same constants the field team used. Salinity correction omitted because
# S_field == S_ref (seawater density depends on both T and S; here only T varied).
rho_sw_temp <- function(t_c) -5e-6 * t_c^2 + 7e-6 * t_c + 1.0001
ARAG <- 2.93   # density of aragonite, the CaCO3 polymorph corals build (g/cm^3)

# ---- Derive the growth metrics ---------------------------------------------
# delta_g (absolute mass gain) is the raw signal; the two rate metrics below
# normalize it by starting mass (no surface area required), so we can show the
# heat effect is not an artifact of how growth is expressed (see the metric-
# comparison table later).
bw <- bw |>
  mutate(
    initial_dry = initial_w / (1 - rho_sw_temp(initial_temp) / ARAG),  # g, t0
    final_dry   = final_w   / (1 - rho_sw_temp(final_temp)   / ARAG),  # g, t-end
    delta_g     = final_dry - initial_dry,    # skeletal mass gained over the window
    days        = pmax(biopsy_day, 1),        # growth-window length; floor at 1 to
                                              # avoid divide-by-zero for any day-0 row
    pct_growth  = 100 * delta_g / initial_dry,                 # PRIMARY metric (% mass change)
    sgr         = 100 * (log(final_dry) - log(initial_dry)) / days,  # % d^-1,
                                              # size-independent exponential growth rate
    g_per_day   = delta_g / days
  )

# Save the cleaned one-row-per-coral table for downstream use.
saveRDS(bw, file.path(DATA_PROC, "buoyant_weight_clean.rds"))

# ---- Primary model: % mass change ------------------------------------------
# Fit on corals with a finite % mass change.
bw_a <- bw |> filter(is.finite(pct_growth))
# Coral-level coefficients are effect-size summaries, but temperature is
# assigned at the tank level. The inferential heat-effect p-value below is
# therefore a tank-level permutation test.
# Model: treatment * wound interaction (does wounding change the heat effect?)
# plus thicket as an additive blocking term to absorb genet-to-genet differences.
m_growth <- lm(pct_growth ~ treatment * wound + thicket, data = bw_a)
res <- broom::tidy(m_growth, conf.int = TRUE) |>   # coefficient table with 95% CIs
  mutate(model_scope = "coral-level descriptive coefficients",
         across(where(is.numeric), \(x) round(x, 4))) |>
  relocate(model_scope, .after = term)
write_csv(res, file.path(TBL_DIR, "05_buoyant_weight_lm.csv"))

# ---- Tank-level permutation test for the heat effect ----------
# Collapse to one mean per tank so each randomized unit contributes one number.
tank_growth <- bw_a |>
  group_by(tank, treatment) |>
  summarise(mean_pct_growth = mean(pct_growth, na.rm = TRUE),
            n_corals = n(), .groups = "drop")
# Observed effect: mean tank % growth at 28C minus at 31C.
obs_diff <- with(tank_growth,
                 mean(mean_pct_growth[treatment == "28C"]) -
                   mean(mean_pct_growth[treatment == "31C"]))
# Exact permutation null: re-label which tanks are "28C" in every possible way
# (combn enumerates all choices of the 28C set), recompute the difference each
# time, and ask how often a re-labeling is as extreme as what we observed. This
# needs no distributional assumption — it just uses the randomization itself.
vals <- tank_growth$mean_pct_growth
assignments <- combn(seq_along(vals), sum(tank_growth$treatment == "28C"),
                     simplify = FALSE)
perm_diffs <- vapply(assignments, function(idx) {
  mean(vals[idx]) - mean(vals[-idx])
}, numeric(1))
tank_test <- tibble(
  response = "pct_growth",
  test = "tank-level exact permutation",
  n_tanks = nrow(tank_growth),
  n_28_tanks = sum(tank_growth$treatment == "28C"),
  n_31_tanks = sum(tank_growth$treatment == "31C"),
  estimate_28_minus_31 = obs_diff,
  # two-sided p = fraction of permuted differences at least as large in magnitude
  # as the observed one (|perm| >= |obs|).
  p_two_sided = mean(abs(perm_diffs) >= abs(obs_diff)),
  note = "Temperature was randomized at the tank level; this p-value treats tanks, not corals, as exchangeable units."
) |>
  mutate(across(where(is.numeric), \(x) signif(x, 4)))
write_csv(tank_test, file.path(TBL_DIR, "05_buoyant_weight_tank_test.csv"))

# ---- Metric comparison (heat effect is robust to metric choice) ------------
# Refit the same fixed-effects model for each surface-area-free growth metric and
# extract the treatment / wound / interaction tests, to show the heat effect does
# not depend on whether growth is expressed as % mass change or SGR.
metric_compare <- purrr::map_dfr(   # map over metrics, stack results into one df
  c(pct_growth = "pct_growth", sgr = "sgr"),
  function(v) {
    d <- bw |> filter(is.finite(.data[[v]]))   # .data[[v]] selects the metric by name
    # car::Anova type 2 = each main effect tested after the other main effect but
    # ignoring the interaction; order-invariant (unlike base
    # anova()'s type-1 sequential sums of squares).
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

# ---- Plot: % mass change ---------------------------------------------------
# x-axis is the 4 treatment x wound combinations (interaction() pastes them into
# one factor); boxes coloured by temperature, jittered points coloured by genet.
p_bw <- ggplot(bw_a, aes(interaction(treatment, wound, sep = " · "),
                         pct_growth, fill = treatment)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(colour = thicket), width = 0.18, height = 0,
              size = 1.6, alpha = 0.85) +
  scale_fill_manual(values  = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                    name = "Temperature") +
  scale_colour_manual(values = PAL_GENO, name = "Genotype") +
  labs(x = NULL,
       y = "Growth (% skeletal mass change over 15 d)",
       title = "Coral growth under heating × wounding",
       subtitle = "Buoyant-weight % mass change (Davies 1989; Jokiel et al. 1978)") +
  theme_pub(10) +
  theme(axis.text.x = element_text(angle = 0))

save_fig(p_bw, "05_buoyant_weight_growth", width = 150, height = 95)

# ---- Console summary -------------------------------------------------------
# Echo the tables so a run shows the main result without opening CSVs.
cat("\n=== % mass change LM (coral-level descriptive) ===\n")
print(res)
cat("\n=== % mass change tank-level temperature test ===\n")
print(tank_test)
cat("\n=== Heat effect is robust to metric choice ===\n")
print(metric_compare |> filter(term == "treatment"))
cat("\nWrote: buoyant_weight_clean.rds, 05_buoyant_weight_lm.csv,",
    "05_buoyant_weight_tank_test.csv,",
    "05b_growth_metric_comparison.csv, 05_buoyant_weight_growth.{pdf,png}\n")
