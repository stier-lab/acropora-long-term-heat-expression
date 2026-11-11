# =============================================================================
# Purpose: Numbers for the "tank temperature and water quality" results.
# Inputs:  output/tables/08_apex_treatment_summary.csv   (analysis/08_temperature.R)
#          output/tables/09_ysi_study_summary.csv        (analysis/09_water_quality.R)
#          output/tables/09_ysi_treatment_summary.csv
# Output:  output/paper/tanks_water_numbers.csv
# Methods: Apex logs every 10 min, averaged per tank per day, then across the
#          four tanks per temperature (Days 0-15). Handheld YSI readings once
#          per tank on 14 of 16 days (none on June 4 or 7).
# =============================================================================
source(here::here("code", "00_setup.R"))

apex <- readr::read_csv(file.path(TBL_DIR, "08_apex_treatment_summary.csv"), show_col_types = FALSE)
ysi  <- readr::read_csv(file.path(TBL_DIR, "09_ysi_study_summary.csv"), show_col_types = FALSE)
ysi_t <- readr::read_csv(file.path(TBL_DIR, "09_ysi_treatment_summary.csv"), show_col_types = FALSE)

a <- function(tr) apex[apex$treatment == tr, ]
y <- function(v) ysi[ysi$variable == v, ]
yt <- function(tr, v) ysi_t[ysi_t$treatment == tr & ysi_t$variable == v, ]
src_apex <- "output/tables/08_apex_treatment_summary.csv"
src_ysi  <- "output/tables/09_ysi_study_summary.csv"

rows <- list()
for (tr in c("28C", "31C")) {
  lab <- if (tr == "28C") "ambient" else "heated"
  rows[[length(rows) + 1]] <- paper_row(
    paste0("temp_mean_", lab), paste("Mean tank temperature (°C),", lab),
    estimate = a(tr)$mean_c, ci_low = a(tr)$min_daily_c, ci_high = a(tr)$max_daily_c,
    test = "mean; ci_low/ci_high = lowest and highest daily tank mean",
    n = a(tr)$n_tank_days, source = src_apex)
}
rows[[length(rows) + 1]] <- paper_row(
  "temp_difference", "Heated minus ambient mean temperature (°C)",
  estimate = a("31C")$mean_c - a("28C")$mean_c, n = a("28C")$n_tanks, source = src_apex)
for (v in c("sal", "do_mgl", "ph")) {
  rows[[length(rows) + 1]] <- paper_row(
    paste0("ysi_", v), paste("Meter reading, mean (statistic = SD):", v),
    estimate = y(v)$mean, statistic = y(v)$sd_readings, test = "mean and SD across readings",
    n = y(v)$n_readings, df1 = y(v)$n_dates, source = src_ysi)
}
for (tr in c("28C", "31C")) {
  lab <- if (tr == "28C") "ambient" else "heated"
  rows[[length(rows) + 1]] <- paper_row(
    paste0("ysi_do_", lab), paste("Dissolved oxygen (mg/L),", lab, "(statistic = SD)"),
    estimate = yt(tr, "do_mgl")$mean, statistic = yt(tr, "do_mgl")$sd_readings,
    test = "mean and SD across readings", n = yt(tr, "do_mgl")$n_readings,
    source = "output/tables/09_ysi_treatment_summary.csv")
}
write_paper(rows, "tanks_water")
