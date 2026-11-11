# =============================================================================
# 02_tanks_water.R: tank temperature and water quality
#
# Section 1 rebuilds tank temperatures from the Apex controller logs (every
# 10 min) and checks them against Molly's compiled sheet and the tank
# assignments. Section 2 cleans the handheld YSI meter readings (salinity,
# dissolved oxygen, pH). Section 3 writes the numbers for the results text.
#
# Writes: data/processed/apex_temperature(_daily).rds, ysi_clean.rds
#         output/tables/08_apex_*.csv, 09_ysi_*.csv
#         figures/08_apex_temperature, 09_ysi_water_chem
#         output/paper/tanks_water_numbers.csv
# =============================================================================

source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages(library(xml2))

# ---- 1. Tank temperature (Apex logs) ---------------------------------------
# Readings above 60 are Fahrenheit and are converted. Readings outside 15-40 °C
# are probe errors and are dropped (listed in 08_apex_excluded_readings.csv).
# Day 0 = 4 June 2025, the day tips were clipped.
recovered <- file.path(DATA_RAW, "apex", "recovered_2026-09-10")
xml_files <- list.files(file.path(DATA_RAW, "apex"), pattern = "\\.xml$",
                        recursive = TRUE, full.names = TRUE)
stopifnot(length(xml_files) > 0)
parse_failures <- list()
parse_apex <- function(path) {
  doc <- tryCatch(read_xml(path), error = function(e) {
    parse_failures[[path]] <<- tibble(source_file = basename(path),
                                    reason = conditionMessage(e))
    message("Excluded malformed XML: ", basename(path))
    NULL
  })
  if (is.null(doc)) return(tibble())
  probes <- xml_find_all(doc, "//record/probe[type='Temp']")
  tibble(
    source_path = path,
    serial = xml_text(xml_find_first(doc, "//serial")),
    logger_timezone = xml_text(xml_find_first(doc, "//timezone")),
    datetime = mdy_hms(xml_text(xml_find_first(probes, "../date")), tz = "UTC"),
    probe = str_squish(xml_text(xml_find_first(probes, "./name"))),
    raw_value = as.numeric(xml_text(xml_find_first(probes, "./value")))
  ) |> filter(str_detect(probe, "^Temp[0-9]+$")) |>
    mutate(value = if_else(raw_value > 60, (raw_value - 32) * 5 / 9, raw_value))
}
raw <- map_dfr(xml_files, parse_apex)
write_csv(bind_rows(parse_failures), file.path(TBL_DIR, "08_apex_parse_failures.csv"))
stopifnot(!anyNA(raw$datetime), !anyNA(raw$value))
manifest <- raw |> group_by(source_path, serial, logger_timezone) |>
  summarise(first_record = min(datetime), last_record = max(datetime),
            n_readings = n(), .groups = "drop") |>
  mutate(md5 = unname(tools::md5sum(source_path)),
         source_path = sub(paste0(DATA_RAW, "/"), "", source_path, fixed = TRUE))
write_csv(manifest, file.path(TBL_DIR, "08_apex_source_manifest.csv"))
duplicates <- raw |> group_by(serial, datetime, probe) |>
  summarise(n_copies = n(), n_values = n_distinct(value), .groups = "drop") |>
  filter(n_copies > 1)
write_csv(duplicates, file.path(TBL_DIR, "08_apex_duplicate_audit.csv"))
stopifnot(all(duplicates$n_values == 1), n_distinct(raw$serial) == 1)
readings <- raw |> distinct(serial, datetime, probe, value) |>
  mutate(tank = as.integer(str_extract(probe, "[0-9]+")),
         valid = is.finite(value) & value > 15 & value < 40)
write_csv(filter(readings, !valid), file.path(TBL_DIR, "08_apex_excluded_readings.csv"))
readings <- filter(readings, valid)
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
stopifnot(all(compiled_check$matches_xml))
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
hourly <- readings |> mutate(datetime = floor_date(datetime, "hour")) |>
  group_by(datetime, probe) |>
  summarise(value_mean = mean(value), value_sd = sd(value), n = n(), .groups = "drop")
daily <- readings |> mutate(date = as.Date(datetime)) |>
  group_by(date, probe) |>
  summarise(value_mean = mean(value), value_sd = sd(value), n = n(),
            n_hours = n_distinct(floor_date(datetime, "hour")), .groups = "drop")
saveRDS(hourly, file.path(DATA_PROC, "apex_temperature.rds"))
saveRDS(daily, file.path(DATA_PROC, "apex_temperature_daily.rds"))
assignments <- read_csv(file.path(recovered, "tank_assignments.csv"),
                       show_col_types = FALSE) |>
  filter(project == "17. LTH expression by temperature") |>
  transmute(date = mdy(date), tank = as.integer(tank),
            treatment = paste0(treatment, "C"), phase = str_trim(.data[["acc.ramp.hold"]]))
stopifnot(!anyDuplicated(select(assignments, date, tank)),
          all(assignments$treatment == tank_treatment(assignments$tank)))
temp_daily <- daily |>
  mutate(tank = as.integer(str_extract(probe, "[0-9]+"))) |>
  inner_join(assignments, by = c("date", "tank"), relationship = "one-to-one") |>
  mutate(day = as.integer(date - as.Date("2025-06-04")))
coverage <- expand_grid(date = as.Date("2025-06-04") + 0:16,
                        tank = sort(c(TANK_28C, TANK_31C))) |>
  left_join(daily |> mutate(tank = as.integer(str_extract(probe, "[0-9]+"))),
            by = c("date", "tank")) |>
  mutate(day = as.integer(date - as.Date("2025-06-04")),
         treatment = tank_treatment(tank),
         across(c(n, n_hours), ~replace_na(.x, 0L)))
write_csv(coverage, file.path(TBL_DIR, "08_apex_daily_coverage.csv"))
summary <- temp_daily |> filter(day >= 0, day <= 15) |>
  group_by(treatment) |>
  summarise(n_tanks = n_distinct(tank), n_tank_days = n(),
            mean_c = mean(value_mean), min_daily_c = min(value_mean),
            max_daily_c = max(value_mean), min_hours = min(n_hours), .groups = "drop")
write_csv(summary, file.path(TBL_DIR, "08_apex_treatment_summary.csv"))
stopifnot(all(filter(coverage, day <= 15)$n == 144),
          all(filter(coverage, day <= 15)$n_hours == 24))
clock_sensitivity <- map_dfr(c(0, -3), function(offset) {
  readings |> mutate(date = as.Date(datetime + hours(offset)),
                     treatment = tank_treatment(tank)) |>
    filter(!is.na(treatment), date >= as.Date("2025-06-04"),
           date <= as.Date("2025-06-19")) |>
    group_by(date, tank, treatment) |> summarise(mean_c = mean(value), .groups = "drop") |>
    group_by(treatment) |> summarise(mean_c = mean(mean_c), .groups = "drop") |>
    mutate(clock_shift_hours = offset)
})
write_csv(clock_sensitivity, file.path(TBL_DIR, "08_apex_clock_sensitivity.csv"))
p_apex <- ggplot(temp_daily, aes(day, value_mean, group = tank, colour = treatment)) +
  annotate("rect", xmin = -7, xmax = 0, ymin = -Inf, ymax = Inf,
           fill = "grey94", colour = NA) +
  geom_hline(yintercept = c(28, 31), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey40") +
  geom_line(linewidth = 0.45, alpha = 0.85) + geom_point(size = 1.1) +
  scale_colour_manual(values = c('28C' = "#56B4E9", '31C' = "#D55E00"),
                      labels = c('28C' = "28 °C", '31C' = "31 °C"), name = "Treatment") +
  scale_x_continuous(breaks = c(-7, -3, 0, 5, 10, 15)) +
  labs(x = "Day relative to tip clipping (Day 0 = 4 June 2025)",
       y = "Daily mean tank temperature (°C)",
       title = "Tank temperatures through the final sampling day",
       subtitle = "One line per tank; shaded area precedes clipping; dashed lines mark targets") +
  theme_pub(10)
save_fig(p_apex, "08_apex_temperature", width = 180, height = 105)

# ---- 2. Water quality (YSI meter) ------------------------------------------
# One reading per tank per day on 14 of the 16 study days (none on June 4 or 7).
# The October 27 daily log is the source; the older Sheet1 export must agree
# wherever the two overlap.
raw <- read_csv(file.path(DATA_RAW, "ysi", "daily_log_2026-10-27.csv"),
                col_types = cols(.default = col_character())) |> janitor::clean_names()
units <- read_csv(file.path(DATA_RAW, "ysi", "daily_log_units_2026-10-27.csv"),
                  show_col_types = FALSE)
stopifnot(units$unit[units$column == "DO mg/L"] == "mg/L",
          units$unit[units$column == "SAL"] == "unitless")
raw <- raw |> mutate(source_row = row_number() + 1L, date = mdy(date),
                     tank = as.integer(tank))
stopifnot(!anyNA(raw$date), !anyDuplicated(raw[c("date", "tank")]))
numeric_fields <- c("temp", "do_percent", "do_mg_l", "sal", "ph")
parse_audit <- raw |> select(source_row, date, tank, all_of(numeric_fields)) |>
  pivot_longer(all_of(numeric_fields), names_to = "variable", values_to = "recorded") |>
  mutate(parsed = suppressWarnings(as.numeric(recorded))) |>
  filter(!is.na(recorded), is.na(parsed))
write_csv(parse_audit, file.path(TBL_DIR, "09_ysi_parse_audit.csv"))
ysi <- raw |> mutate(across(all_of(numeric_fields), ~suppressWarnings(as.numeric(.x))),
  recorded_treatment = as.integer(treatment), treatment = tank_treatment(tank),
  temp_f = temp, temp_c = (temp - 32) * 5 / 9, do_pct = do_percent, do_mgl = do_mg_l)
stopifnot(!anyNA(ysi$treatment),
          all(as.character(ysi$treatment) == paste0(ysi$recorded_treatment, "C")))
saveRDS(ysi, file.path(DATA_PROC, "ysi_clean.rds"))
study <- ysi |> filter(date >= as.Date("2025-06-04"), date <= as.Date("2025-06-19"))
coverage <- expand_grid(date = seq(as.Date("2025-06-04"), as.Date("2025-06-19"), by = "day"),
                        tank = c(3L,4L,5L,6L,9L,10L,11L,12L)) |>
  left_join(study |> transmute(date, tank, observed = TRUE), by = c("date", "tank")) |>
  mutate(observed = replace_na(observed, FALSE))
write_csv(coverage, file.path(TBL_DIR, "09_ysi_study_coverage.csv"))
stopifnot(nrow(study) == 112L, n_distinct(study$date) == 14L,
          all(count(study, tank)$n == 14L),
          all(complete.cases(study[c("temp_c","do_mgl","do_pct","sal","ph")])),
          all(study$do_mgl > 0), all(study$sal > 0),
          all(study$ph > 0 & study$ph < 14))
old <- read_csv(file.path(DATA_RAW, "ysi", "Sheet1.csv"), show_col_types = FALSE) |>
  janitor::clean_names() |> mutate(date = as_date(date), tank = as.integer(tank)) |>
  select(date, tank, all_of(numeric_fields))
overlap <- inner_join(old, raw |> select(date,tank,all_of(numeric_fields)),
                       by = c("date","tank"), suffix = c("_old","_new"))
write_csv(overlap, file.path(TBL_DIR, "09_ysi_legacy_overlap.csv"))
for (v in numeric_fields) {
  a <- suppressWarnings(as.numeric(overlap[[paste0(v,"_old")]]))
  b <- suppressWarnings(as.numeric(overlap[[paste0(v,"_new")]]))
  stopifnot(all(is.na(a) | is.na(b) | abs(a-b) < 1e-10))
}
long <- study |> select(date,tank,treatment,temp_c,do_mgl,do_pct,sal,ph) |>
  pivot_longer(c(temp_c,do_mgl,do_pct,sal,ph), names_to = "variable", values_to = "value")
pool <- function(d) d |> group_by(variable) |>
  summarise(n_readings = n(), n_dates = n_distinct(date), n_tanks = n_distinct(tank),
            mean = mean(value), sd_readings = sd(value), minimum = min(value),
            maximum = max(value), .groups = "drop")
overall <- pool(long)
by_treatment <- long |> group_by(treatment) |> group_modify(~pool(.x)) |> ungroup()
tank_means <- long |> group_by(variable,treatment,tank) |>
  summarise(mean = mean(value), n_dates = n(), .groups = "drop")
tank_summary <- bind_rows(
  tank_means |> mutate(scope = "All tanks"),
  tank_means |> mutate(scope = as.character(treatment))) |>
  group_by(scope, variable) |>
  summarise(n_tanks = n(), sd_tank_means = sd(mean), mean = mean(mean),
            se_tank_means = sd_tank_means / sqrt(n_tanks), .groups = "drop")
stopifnot(all(is.finite(tank_summary$sd_tank_means)),
          all(is.finite(tank_summary$se_tank_means)))
write_csv(overall, file.path(TBL_DIR,"09_ysi_study_summary.csv"))
write_csv(by_treatment, file.path(TBL_DIR,"09_ysi_treatment_summary.csv"))
write_csv(tank_means, file.path(TBL_DIR,"09_ysi_tank_means.csv"))
write_csv(tank_summary, file.path(TBL_DIR,"09_ysi_tank_mean_uncertainty.csv"))
labels <- c(temp_c="Temperature (degrees C)",do_mgl="Dissolved oxygen (mg/L)",
            sal="Salinity (unitless)",ph="pH")
plot_data <- long |> filter(variable %in% names(labels)) |>
  mutate(variable = factor(variable, levels = names(labels), labels = labels))
p <- ggplot(plot_data, aes(date,value,colour=treatment,group=tank)) +
  geom_line(linewidth=0.35,alpha=0.7) + geom_point(size=1.2) +
  scale_colour_manual(values=PAL_TEMP, name="Temperature treatment") +
  facet_wrap(~variable,scales="free_y",ncol=2) +
  scale_x_date(limits=as.Date(c("2025-06-04","2025-06-19")),date_breaks="5 days",date_labels="%b %d") +
  labs(x=NULL,y=NULL,title="Water-quality spot checks during the study",
       subtitle="Eight tanks; 14 sampled dates. June 4 and June 7 were not recorded.") +
  theme_pub(9)
save_fig(p,"09_ysi_water_chem",width=185,height=140)

# ---- 3. Numbers for the results text ---------------------------------------
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
