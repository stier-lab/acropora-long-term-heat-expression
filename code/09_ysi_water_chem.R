# YSI spot checks: observed conditions, not a test of treatment equivalence.
# Source snapshot and units are retained unchanged under data/raw/ysi.
source(here::here("code", "00_setup.R"))
raw <- read_csv(file.path(DATA_RAW, "ysi", "daily_log_2026-10-27.csv"),
                col_types = cols(.default = col_character())) |> janitor::clean_names()
units <- read_csv(file.path(DATA_RAW, "ysi", "daily_log_units_2026-10-27.csv"),
                  show_col_types = FALSE)
stopifnot(units$unit[units$column == "DO mg/L"] == "mg/L",
          units$unit[units$column == "SAL"] == "unitless")
raw <- raw |> mutate(source_row = row_number() + 1L, date = mdy(date),
                     tank = as.integer(tank))
stopifnot(!anyNA(raw$date), !anyDuplicated(raw[c("date", "tank")]))
# Invalid numeric strings remain missing and are reported, not guessed.
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
# Check the overlap without appending the old export (which would duplicate it).
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
fmt <- function(v) {
  r <- filter(overall,variable==v)
  sprintf("%.2f +/- %.2f",r$mean,r$sd_readings)
}
methods <- paste0("We measured salinity, dissolved oxygen, and pH using a YSI meter in all eight tanks on 14 dates within June 4-19, 2025 (112 tank-date observations; no readings were available for June 4 or June 7). Across the recorded observations, mean +/- SD was ",fmt("sal")," for salinity (unitless), ",fmt("do_mgl")," mg/L for dissolved oxygen, and ",fmt("ph")," for pH. Standard deviations describe variation across tanks and sampled dates.")
writeLines(methods,file.path(TBL_DIR,"09_ysi_methods_wording.txt"))
print(overall)
print(by_treatment)
print(tank_summary)
