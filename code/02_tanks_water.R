# =============================================================================
# 02_tanks_water.R: tank temperature and water quality
#
# The experiment used eight tanks: four held at 28 °C (ambient) and four at
# 31 °C (heated). Because tanks, not fragments, received the temperature, the
# tank is the replicate when comparing treatments. Section 1 rebuilds tank
# temperatures from the Apex controller logs (one reading every 10 min) and
# checks them against Molly's compiled sheet and the tank assignments.
# Section 2 cleans the handheld YSI meter readings (salinity, dissolved
# oxygen, pH). Section 3 writes the numbers for the results text, including
# the oxygen drop expected from warmer water alone. Day 0 = 4 June 2025, the
# day tips were clipped.
#
# Sections
#   1. Tank temperature (Apex logs)
#   2. Water quality (YSI meter)
#   3. Numbers for the results text
#
# Reads:  data/raw/apex/**/*.xml, data/raw/apex/recovered_2026-09-10/
#           (compiled_temperature.csv, tank_assignments.csv)
#         data/raw/ysi/daily_log_2026-10-27.csv, daily_log_units_2026-10-27.csv,
#           Sheet1.csv
# Writes: data/processed/apex_temperature(_daily).rds, ysi_clean.rds
#         output/tables/08_apex_*.csv, 09_ysi_*.csv
#         figures/08_apex_temperature, 09_ysi_water_chem
#         output/paper/tanks_water_numbers.csv
# Run:    sourced by code/run_all.R (or: Rscript code/run_all.R)
# =============================================================================

source(here::here("code", "00_setup.R"))

suppressPackageStartupMessages(library(xml2))   # reads the Apex .xml logs


# ---- 1. Tank temperature (Apex logs) ---------------------------------------
# The Apex aquarium controller logs every temperature probe (one per tank)
# every 10 minutes. Readings above 60 are Fahrenheit and are converted.
# Readings outside 15-40 °C are probe errors and are dropped (listed in
# 08_apex_excluded_readings.csv). Every step leaves an audit table in
# output/tables/ so the cleaning can be checked without rerunning it.

## 1a. Parse every Apex XML log ----------------------------------------------
recovered <- file.path(DATA_RAW, "apex", "recovered_2026-09-10")
xml_files <- list.files(file.path(DATA_RAW, "apex"), pattern = "\\.xml$",
                        recursive = TRUE, full.names = TRUE)
stopifnot(length(xml_files) > 0)   # guarantee the logs were found

# A malformed file is skipped (and recorded) instead of stopping the run.
# `<<-` writes to parse_failures outside the function, so the failures
# collected inside tryCatch() survive after the function returns.
parse_failures <- list()
parse_apex <- function(path) {
  doc <- tryCatch(read_xml(path), error = function(e) {
    parse_failures[[path]] <<- tibble(source_file = basename(path),
                                    reason = conditionMessage(e))
    message("Excluded malformed XML: ", basename(path))
    NULL
  })
  if (is.null(doc)) return(tibble())

  # Each <record> holds a timestamp and one <probe> per sensor; keep the
  # temperature probes named Temp<tank number>.
  probes <- xml_find_all(doc, "//record/probe[type='Temp']")
  tibble(
    source_path = path,
    serial = xml_text(xml_find_first(doc, "//serial")),
    logger_timezone = xml_text(xml_find_first(doc, "//timezone")),
    datetime = mdy_hms(xml_text(xml_find_first(probes, "../date")), tz = "UTC"),
    probe = str_squish(xml_text(xml_find_first(probes, "./name"))),
    raw_value = as.numeric(xml_text(xml_find_first(probes, "./value")))
  ) |>
    filter(str_detect(probe, "^Temp[0-9]+$")) |>
    mutate(value = if_else(raw_value > 60, (raw_value - 32) * 5 / 9, raw_value))
}

# map_dfr() runs parse_apex() on every file and stacks the results.
raw <- map_dfr(xml_files, parse_apex)
write_csv(bind_rows(parse_failures), file.path(TBL_DIR, "08_apex_parse_failures.csv"))
stopifnot(!anyNA(raw$datetime), !anyNA(raw$value))   # every reading has a time and value

## 1b. Source manifest and duplicate audit -----------------------------------
# The manifest lists each log file with its date range, reading count and
# md5 checksum (a fingerprint that changes if the file is edited).
manifest <- raw |>
  group_by(source_path, serial, logger_timezone) |>
  summarise(first_record = min(datetime), last_record = max(datetime),
            n_readings = n(), .groups = "drop") |>
  mutate(md5 = unname(tools::md5sum(source_path)),
         source_path = sub(paste0(DATA_RAW, "/"), "", source_path, fixed = TRUE))
write_csv(manifest, file.path(TBL_DIR, "08_apex_source_manifest.csv"))

# Log files overlap in time, so the same reading can appear more than once.
# Copies must agree (one distinct value), and all files must come from a
# single controller (one serial number); then the copies are collapsed.
duplicates <- raw |>
  group_by(serial, datetime, probe) |>
  summarise(n_copies = n(), n_values = n_distinct(value), .groups = "drop") |>
  filter(n_copies > 1)
write_csv(duplicates, file.path(TBL_DIR, "08_apex_duplicate_audit.csv"))
stopifnot(all(duplicates$n_values == 1), n_distinct(raw$serial) == 1)

## 1c. Drop probe errors -----------------------------------------------------
readings <- raw |>
  distinct(serial, datetime, probe, value) |>
  mutate(tank = as.integer(str_extract(probe, "[0-9]+")),
         valid = is.finite(value) & value > 15 & value < 40)
write_csv(filter(readings, !valid), file.path(TBL_DIR, "08_apex_excluded_readings.csv"))
readings <- filter(readings, valid)

## 1d. Cross-check against Molly's compiled temperature sheet ----------------
# Every row of the compiled sheet must match an XML reading exactly.
# relationship = "many-to-one" makes the join stop with an error if one
# compiled row would match more than one XML reading.
compiled <- read_csv(file.path(recovered, "compiled_temperature.csv"),
                     show_col_types = FALSE) |>
  transmute(datetime = mdy_hms(.data[["/record/date"]], tz = "UTC"),
            probe = .data[["/record/probe/name"]],
            raw_value = .data[["/record/probe/value"]],
            compiled_c = if_else(raw_value > 60, (raw_value - 32) * 5 / 9, raw_value))

compiled_check <- compiled |>
  left_join(select(readings, datetime, probe, value),
            by = c("datetime", "probe"), relationship = "many-to-one") |>
  mutate(matches_xml = !is.na(value) & abs(compiled_c - value) < 1e-8)
write_csv(tibble(n_compiled = nrow(compiled_check),
                 n_duplicate_rows = sum(duplicated(compiled)),
                 n_matches_xml = sum(compiled_check$matches_xml),
                 n_unmatched = sum(!compiled_check$matches_xml)),
          file.path(TBL_DIR, "08_apex_compiled_crosscheck.csv"))
stopifnot(all(compiled_check$matches_xml))   # no compiled value disagrees with the XML

# The reverse check: every XML reading from the study tanks between 28 May
# (Day -7) and 19 June 2025 must also be in the compiled sheet.
# anti_join() keeps XML rows that have no match in the compiled sheet.
compiled_unique <- compiled |> distinct(datetime, probe, compiled_c)
study_xml <- readings |>
  filter(as.Date(datetime) >= as.Date("2025-05-28"),
         as.Date(datetime) <= as.Date("2025-06-19"),
         tank %in% c(TANK_28C, TANK_31C)) |>
  select(datetime, probe, value)
missing_compiled <- anti_join(study_xml, compiled_unique,
                              by = c("datetime", "probe", "value" = "compiled_c"))
write_csv(tibble(first_compiled = min(compiled_unique$datetime),
                 last_compiled = max(compiled_unique$datetime),
                 n_unique_compiled = nrow(compiled_unique),
                 n_study_xml = nrow(study_xml),
                 n_study_missing_from_compiled = nrow(missing_compiled)),
          file.path(TBL_DIR, "08_apex_compiled_standalone_check.csv"))
stopifnot(nrow(missing_compiled) == 0)

## 1e. Hourly and daily tank means -------------------------------------------
# floor_date(..., "hour") rounds each time down to the start of its hour.
hourly <- readings |>
  mutate(datetime = floor_date(datetime, "hour")) |>
  group_by(datetime, probe) |>
  summarise(value_mean = mean(value), value_sd = sd(value), n = n(), .groups = "drop")
daily <- readings |>
  mutate(date = as.Date(datetime)) |>
  group_by(date, probe) |>
  summarise(value_mean = mean(value), value_sd = sd(value), n = n(),
            n_hours = n_distinct(floor_date(datetime, "hour")), .groups = "drop")
saveRDS(hourly, file.path(DATA_PROC, "apex_temperature.rds"))
saveRDS(daily, file.path(DATA_PROC, "apex_temperature_daily.rds"))

## 1f. Match daily means to the tank assignment sheet ------------------------
# The assignment sheet says which treatment and phase (acclimation, ramp or
# hold) each tank was in on each date. Its treatments must agree with the
# fixed tank map in tank_treatment() (00_setup.R), with one row per tank per
# date.
assignments <- read_csv(file.path(recovered, "tank_assignments.csv"),
                       show_col_types = FALSE) |>
  filter(project == "17. LTH expression by temperature") |>
  transmute(date = mdy(date), tank = as.integer(tank),
            treatment = paste0(treatment, "C"), phase = str_trim(.data[["acc.ramp.hold"]]))
stopifnot(!anyDuplicated(select(assignments, date, tank)),
          all(assignments$treatment == tank_treatment(assignments$tank)))

# relationship = "one-to-one": error if any tank-day would match twice.
temp_daily <- daily |>
  mutate(tank = as.integer(str_extract(probe, "[0-9]+"))) |>
  inner_join(assignments, by = c("date", "tank"), relationship = "one-to-one") |>
  mutate(day = as.integer(date - as.Date("2025-06-04")))

## 1g. Coverage and treatment summary, Days 0-15 -----------------------------
# One row per study tank per day (Days 0-16), with 0 readings where the log
# has a gap. On Days 0-15 every tank must have all 144 ten-minute readings
# (6 per hour x 24 hours) on every day.
coverage <- expand_grid(date = as.Date("2025-06-04") + 0:16,
                        tank = sort(c(TANK_28C, TANK_31C))) |>
  left_join(daily |> mutate(tank = as.integer(str_extract(probe, "[0-9]+"))),
            by = c("date", "tank")) |>
  mutate(day = as.integer(date - as.Date("2025-06-04")),
         treatment = tank_treatment(tank),
         across(c(n, n_hours), ~ replace_na(.x, 0L)))
write_csv(coverage, file.path(TBL_DIR, "08_apex_daily_coverage.csv"))

# Mean, lowest and highest daily tank mean for each treatment.
apex_summary <- temp_daily |>
  filter(day >= 0, day <= 15) |>
  group_by(treatment) |>
  summarise(n_tanks = n_distinct(tank), n_tank_days = n(),
            mean_c = mean(value_mean), min_daily_c = min(value_mean),
            max_daily_c = max(value_mean), min_hours = min(n_hours), .groups = "drop")
write_csv(apex_summary, file.path(TBL_DIR, "08_apex_treatment_summary.csv"))
stopifnot(all(filter(coverage, day <= 15)$n == 144),
          all(filter(coverage, day <= 15)$n_hours == 24))

## 1h. Clock sensitivity -----------------------------------------------------
# Timestamps are parsed as UTC. To check that the treatment means do not
# depend on where the day boundary falls, they are recomputed with every
# time shifted back 3 hours. Each tank-day mean is taken first, then averaged per treatment.
clock_sensitivity <- map_dfr(c(0, -3), function(offset) {
  readings |>
    mutate(date = as.Date(datetime + hours(offset)),
           treatment = tank_treatment(tank)) |>
    filter(!is.na(treatment), date >= as.Date("2025-06-04"),
           date <= as.Date("2025-06-19")) |>
    group_by(date, tank, treatment) |>
    summarise(mean_c = mean(value), .groups = "drop") |>
    group_by(treatment) |>
    summarise(mean_c = mean(mean_c), .groups = "drop") |>
    mutate(clock_shift_hours = offset)
})
write_csv(clock_sensitivity, file.path(TBL_DIR, "08_apex_clock_sensitivity.csv"))

## 1i. Figure: daily mean temperature of each tank ---------------------------
# The grey band marks the pre-clipping period (Day -7 to 0); dashed lines
# mark the 28 and 31 °C targets.
p_apex <- ggplot(temp_daily, aes(day, value_mean, group = tank, colour = treatment)) +
  annotate("rect", xmin = -7, xmax = 0, ymin = -Inf, ymax = Inf,
           fill = "grey94", colour = NA) +
  geom_hline(yintercept = c(28, 31), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey40") +
  geom_line(linewidth = 0.45, alpha = 0.85) +
  geom_point(size = 1.1) +
  scale_colour_manual(values = c('28C' = "#56B4E9", '31C' = "#D55E00"),
                      labels = c('28C' = "28 °C", '31C' = "31 °C"), name = "Treatment") +
  scale_x_continuous(breaks = c(-7, -3, 0, 5, 10, 15)) +
  labs(x = "Day relative to tip clipping (Day 0 = 4 June 2025)",
       y = "Daily mean tank temperature (°C)") +
  theme_pub(10)
save_fig(p_apex, "08_apex_temperature", width = 180, height = 105)


# ---- 2. Water quality (YSI meter) ------------------------------------------
# One reading per tank per day on 14 of the 16 study days (none on June 4 or 7).
# The October 27 daily log is the source; the older Sheet1 export must agree
# wherever the two overlap. Note that `raw` and `coverage` are reused names:
# from here on they hold the YSI data, not the Apex data.

## 2a. Read the daily log and parse numbers ----------------------------------
# Every column is read as text first, so cells that will not convert to a
# number can be listed (09_ysi_parse_audit.csv) instead of silently becoming NA.
raw <- read_csv(file.path(DATA_RAW, "ysi", "daily_log_2026-10-27.csv"),
                col_types = cols(.default = col_character())) |>
  janitor::clean_names()

# Guarantee the units are what the code assumes (oxygen in mg/L, salinity
# unitless on the practical salinity scale).
units <- read_csv(file.path(DATA_RAW, "ysi", "daily_log_units_2026-10-27.csv"),
                  show_col_types = FALSE)
stopifnot(units$unit[units$column == "DO mg/L"] == "mg/L",
          units$unit[units$column == "SAL"] == "unitless")

# source_row = the row number in the spreadsheet (+1 for the header row).
raw <- raw |>
  mutate(source_row = row_number() + 1L, date = mdy(date),
         tank = as.integer(tank))
# Every date parsed, and exactly one reading per tank per day.
stopifnot(!anyNA(raw$date), !anyDuplicated(raw[c("date", "tank")]))

numeric_fields <- c("temp", "do_percent", "do_mg_l", "sal", "ph")
parse_audit <- raw |>
  select(source_row, date, tank, all_of(numeric_fields)) |>
  pivot_longer(all_of(numeric_fields), names_to = "variable", values_to = "recorded") |>
  mutate(parsed = suppressWarnings(as.numeric(recorded))) |>
  filter(!is.na(recorded), is.na(parsed))
write_csv(parse_audit, file.path(TBL_DIR, "09_ysi_parse_audit.csv"))

# The meter records temperature in °F; temp_c converts it. Treatment is taken
# from the tank number and must match the treatment written in the log.
ysi <- raw |>
  mutate(across(all_of(numeric_fields), ~ suppressWarnings(as.numeric(.x))),
         recorded_treatment = as.integer(treatment), treatment = tank_treatment(tank),
         temp_f = temp, temp_c = (temp - 32) * 5 / 9, do_pct = do_percent, do_mgl = do_mg_l)
stopifnot(!anyNA(ysi$treatment),
          all(as.character(ysi$treatment) == paste0(ysi$recorded_treatment, "C")))
saveRDS(ysi, file.path(DATA_PROC, "ysi_clean.rds"))

## 2b. Study-period coverage and range checks --------------------------------
# Study period = Days 0-15 (4-19 June 2025), eight tanks. Expect 112 readings
# (8 tanks x 14 days), no missing values and physically possible numbers.
study <- ysi |> filter(date >= as.Date("2025-06-04"), date <= as.Date("2025-06-19"))
coverage <- expand_grid(date = seq(as.Date("2025-06-04"), as.Date("2025-06-19"), by = "day"),
                        tank = sort(c(TANK_28C, TANK_31C))) |>
  left_join(study |> transmute(date, tank, observed = TRUE), by = c("date", "tank")) |>
  mutate(observed = replace_na(observed, FALSE))
write_csv(coverage, file.path(TBL_DIR, "09_ysi_study_coverage.csv"))
stopifnot(nrow(study) == 112L, n_distinct(study$date) == 14L,
          all(count(study, tank)$n == 14L),
          all(complete.cases(study[c("temp_c", "do_mgl", "do_pct", "sal", "ph")])),
          all(study$do_mgl > 0), all(study$sal > 0),
          all(study$ph > 0 & study$ph < 14))

## 2c. The older Sheet1 export must agree with the daily log -----------------
# Join the two on date and tank; every value present in both must be equal.
old <- read_csv(file.path(DATA_RAW, "ysi", "Sheet1.csv"), show_col_types = FALSE) |>
  janitor::clean_names() |>
  mutate(date = as_date(date), tank = as.integer(tank)) |>
  select(date, tank, all_of(numeric_fields))
overlap <- inner_join(old, raw |> select(date, tank, all_of(numeric_fields)),
                       by = c("date", "tank"), suffix = c("_old", "_new"))
write_csv(overlap, file.path(TBL_DIR, "09_ysi_legacy_overlap.csv"))
for (v in numeric_fields) {
  a <- suppressWarnings(as.numeric(overlap[[paste0(v, "_old")]]))
  b <- suppressWarnings(as.numeric(overlap[[paste0(v, "_new")]]))
  stopifnot(all(is.na(a) | is.na(b) | abs(a - b) < 1e-10))
}

## 2d. Summaries: pooled, by treatment and by tank ---------------------------
# Long format: one row per reading per variable.
long <- study |>
  select(date, tank, treatment, temp_c, do_mgl, do_pct, sal, ph) |>
  pivot_longer(c(temp_c, do_mgl, do_pct, sal, ph),
               names_to = "variable", values_to = "value")

# pool() summarises every reading of each variable (n, mean, SD, range).
pool <- function(d) d |>
  group_by(variable) |>
  summarise(n_readings = n(), n_dates = n_distinct(date), n_tanks = n_distinct(tank),
            mean = mean(value), sd_readings = sd(value), minimum = min(value),
            maximum = max(value), .groups = "drop")
overall <- pool(long)
# group_modify() runs pool() separately within each treatment.
by_treatment <- long |> group_by(treatment) |> group_modify(~ pool(.x)) |> ungroup()

# Tank means first, then the SD and SE of those means: tanks, not single
# readings, are the replicate, so this is the uncertainty to report when
# comparing treatments. Each tank mean appears twice: once under "All tanks"
# and once under its own treatment.
tank_means <- long |>
  group_by(variable, treatment, tank) |>
  summarise(mean = mean(value), n_dates = n(), .groups = "drop")
tank_summary <- bind_rows(
  tank_means |> mutate(scope = "All tanks"),
  tank_means |> mutate(scope = as.character(treatment))) |>
  group_by(scope, variable) |>
  summarise(n_tanks = n(), sd_tank_means = sd(mean), mean = mean(mean),
            se_tank_means = sd_tank_means / sqrt(n_tanks), .groups = "drop")
stopifnot(all(is.finite(tank_summary$sd_tank_means)),
          all(is.finite(tank_summary$se_tank_means)))

write_csv(overall, file.path(TBL_DIR, "09_ysi_study_summary.csv"))
write_csv(by_treatment, file.path(TBL_DIR, "09_ysi_treatment_summary.csv"))
write_csv(tank_means, file.path(TBL_DIR, "09_ysi_tank_means.csv"))
write_csv(tank_summary, file.path(TBL_DIR, "09_ysi_tank_mean_uncertainty.csv"))

## 2e. Figure: daily readings of each tank, one panel per variable -----------
# Percent oxygen saturation is left out of the figure (not in `labels`).
labels <- c(temp_c = "Temperature (degrees C)", do_mgl = "Dissolved oxygen (mg/L)",
            sal = "Salinity (unitless)", ph = "pH")
plot_data <- long |>
  filter(variable %in% names(labels)) |>
  mutate(variable = factor(variable, levels = names(labels), labels = labels))

p <- ggplot(plot_data, aes(date, value, colour = treatment, group = tank)) +
  geom_line(linewidth = 0.35, alpha = 0.7) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = PAL_TEMP, name = "Temperature treatment") +
  facet_wrap(~ variable, scales = "free_y", ncol = 2) +
  scale_x_date(limits = as.Date(c("2025-06-04", "2025-06-19")), date_breaks = "5 days",
               date_labels = "%b %d") +
  labs(x = NULL, y = NULL) +
  theme_pub(9)
save_fig(p, "09_ysi_water_chem", width = 185, height = 140)


# ---- 3. Numbers for the results text ---------------------------------------
# Reads back the summary tables written above and turns them into rows of
# output/paper/tanks_water_numbers.csv via paper_row() and write_paper()
# (both in 00_setup.R). "ambient" = 28 °C, "heated" = 31 °C. Note that `ysi`
# is reused here for the summary table, not the clean readings.
apex <- readr::read_csv(file.path(TBL_DIR, "08_apex_treatment_summary.csv"),
                        show_col_types = FALSE)
ysi  <- readr::read_csv(file.path(TBL_DIR, "09_ysi_study_summary.csv"),
                        show_col_types = FALSE)
ysi_t <- readr::read_csv(file.path(TBL_DIR, "09_ysi_treatment_summary.csv"),
                         show_col_types = FALSE)

# Small helpers that pick one row: a() by treatment, y() by variable, yt() by
# both.
a <- function(tr) apex[apex$treatment == tr, ]
y <- function(v) ysi[ysi$variable == v, ]
yt <- function(tr, v) ysi_t[ysi_t$treatment == tr & ysi_t$variable == v, ]
src_apex <- "output/tables/08_apex_treatment_summary.csv"
src_ysi  <- "output/tables/09_ysi_study_summary.csv"

## 3a. Temperature -----------------------------------------------------------
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

## 3b. Water quality, all tanks pooled, and oxygen by treatment --------------
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

## 3c. Oxygen in the heated tanks: how much is temperature alone? ------------
# Oxygen in the heated tanks (open question). Warm water holds less oxygen, so
# part of the lower concentration is expected at any saturation. Solubility
# from Garcia & Gordon (1992), their recommended fit to the Benson & Krause
# data (the coefficients below), converted from ml/L to mg/L (x 1.42905).
# t_c = temperature (°C), s = salinity; ts is their scaled temperature.
o2_sat_mgl <- function(t_c, s) {
  ts <- log((298.15 - t_c) / (273.15 + t_c))
  ln <- 2.00907 + 3.22014 * ts + 4.0501 * ts^2 + 4.94457 * ts^3 - 0.256847 * ts^4 +
    3.88767 * ts^5 + s * (-6.24523e-3 - 7.37614e-3 * ts - 1.0341e-2 * ts^2 -
    8.17083e-3 * ts^3) - 4.88682e-7 * s^2
  exp(ln) * 1.42905
}

for (tr in c("28C", "31C")) {
  lab <- if (tr == "28C") "ambient" else "heated"
  rows[[length(rows) + 1]] <- paper_row(
    paste0("ysi_do_pct_", lab), paste("Dissolved oxygen (% saturation),", lab,
                                      "(ci_low/ci_high = lowest and highest reading)"),
    estimate = yt(tr, "do_pct")$mean, ci_low = yt(tr, "do_pct")$minimum,
    ci_high = yt(tr, "do_pct")$maximum, n = yt(tr, "do_pct")$n_readings,
    source = "output/tables/09_ysi_treatment_summary.csv")
}

# Expected drop = saturation at the ambient tanks' mean temperature and
# salinity minus saturation at the heated tanks'. It is compared with the
# observed drop in mean mg/L (stored as `statistic`).
sol_drop <- o2_sat_mgl(yt("28C", "temp_c")$mean, yt("28C", "sal")$mean) -
  o2_sat_mgl(yt("31C", "temp_c")$mean, yt("31C", "sal")$mean)
rows[[length(rows) + 1]] <- paper_row(
  "do_solubility_drop",
  "Oxygen drop expected from warmer water alone at 100% saturation (mg/L); statistic = observed drop",
  estimate = sol_drop, statistic = yt("28C", "do_mgl")$mean - yt("31C", "do_mgl")$mean,
  test = "Garcia & Gordon (1992) solubility at mean tank temperature and salinity",
  source = "output/tables/09_ysi_treatment_summary.csv")

write_paper(rows, "tanks_water")
