# =============================================================================
# Purpose: Daily YSI spot checks — TEMP, DO%, DO mg/L, SAL, pH.
#          - Convert temperature from °F (recorded) to °C
#          - Multi-panel daily-trace figure
#
# What & why: the YSI is a handheld probe dipped in each tank once a day (~1800h)
#   to confirm the tanks stayed within safe, comparable water-quality bounds.
#   This is a QA/QC check, not a hypothesis test — we want to show that the only
#   thing that differed between the 28 °C and 31 °C tanks was temperature, not
#   oxygen, salinity, or pH. (Continuous high-resolution temperature comes from
#   the Apex loggers in code/08; this is the daily manual cross-check.)
# Input:   data/raw/ysi/Sheet1.csv
# Output:  data/processed/ysi_clean.rds
#          figures/09_ysi_water_chem.{pdf,png}
# =============================================================================

# 00_setup.R loads packages (tidyverse etc.) and defines shared paths
# (DATA_RAW, DATA_PROC, ...), the theme_pub() plot theme, and save_fig().
source(here::here("code", "00_setup.R"))

# ---- Load ------------------------------------------------------------------
# clean_names() standardises the spreadsheet's column headers to snake_case
# (e.g. "DO mg/L" -> do_mg_l) so we can refer to them predictably below.
raw <- read_csv(file.path(DATA_RAW, "ysi", "Sheet1.csv"),
                show_col_types = FALSE) |>
  janitor::clean_names()

# ---- Clean + unit-convert --------------------------------------------------
# Coerce every column to its proper type (the CSV reads some as text), convert
# temperature to °C, and tag each row with its experimental treatment.
ysi <- raw |>
  filter(!is.na(date)) |>                       # drop blank/trailing rows
  mutate(
    date     = as_date(date),
    tank     = as.integer(tank),
    temp_f   = as.numeric(temp),                # YSI recorded temperature in °F
    temp_c   = (temp_f - 32) * 5 / 9,           # standard °F -> °C conversion
    do_pct   = as.numeric(do_percent),          # dissolved oxygen, % saturation
    do_mgl   = as.numeric(do_mg_l),             # dissolved oxygen, mg/L
    sal      = as.numeric(sal),                 # salinity (psu)
    ph       = as.numeric(ph),
    # Map tank number -> treatment via the shared plumbing layout (tank_treatment()
    # in 00_setup.R; tanks 3,6,9,12 = ambient 28 °C, 4,5,10,11 = heated 31 °C,
    # anything else -> NA). Same assignment used by 08 (Apex) and 16 (Fig 1).
    treatment = tank_treatment(tank)
  ) |>
  filter(!is.na(treatment))                     # keep only the experimental tanks

# Save the cleaned, one-row-per-reading table for any downstream use.
saveRDS(ysi, file.path(DATA_PROC, "ysi_clean.rds"))

# ---- Reshape for plotting --------------------------------------------------
# ggplot facets need "long" data: one row per (reading × variable). pivot_longer
# stacks the four measured variables into a var/value pair, and we relabel them
# with publication-ready facet titles (incl. units) in a fixed panel order.
ysi_long <- ysi |>
  select(date, tank, treatment, temp_c, do_mgl, sal, ph) |>
  pivot_longer(c(temp_c, do_mgl, sal, ph),
               names_to = "var", values_to = "value") |>
  mutate(var = factor(var,
                      levels = c("temp_c", "do_mgl", "sal", "ph"),
                      labels = c("Temperature (°C)", "DO (mg/L)",
                                 "Salinity (psu)", "pH")))

# ---- Figure ----------------------------------------------------------------
# One faceted panel per variable; one line per tank, coloured by treatment.
# scales = "free_y" lets each panel use its own y-range (the four variables are
# on completely different scales). Blue = ambient, orange/red = heated.
p_ysi <- ggplot(ysi_long, aes(date, value,
                              group = tank, colour = treatment)) +
  geom_line(linewidth = 0.3, alpha = 0.8) +
  geom_point(size = 0.9, alpha = 0.85) +
  scale_colour_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                      name = "Treatment") +
  facet_wrap(~ var, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = NULL,
       title = "Tank water-quality spot checks (YSI)",
       subtitle = "One line per tank, daily 1800 readings") +
  theme_pub(9)

# save_fig() writes both a .pdf (vector, for the manuscript) and a .png (preview)
# at the given size in mm. See 00_setup.R.
save_fig(p_ysi, "09_ysi_water_chem", width = 170, height = 130)

# ---- Done ------------------------------------------------------------------
cat("Wrote ysi_clean.rds and 09_ysi_water_chem.{pdf,png}\n")
