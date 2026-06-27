# =============================================================================
# Purpose: Parse Neptune Apex XML datalogs into a per-tank temperature time
#          series and visualize the heat ramp + steady-state period.
#          Apex logs are large (one <record> per probe poll); we stream-parse
#          file-by-file, aggregate to hourly means inside each file, and
#          discard raw records before moving on. This keeps memory bounded.
#
# What & why: the Neptune Apex is the aquarium controller that ran the heating
#   system; it polls every tank's temperature probe continuously and writes the
#   readings to XML log files. This is the high-resolution, ground-truth record
#   of the thermal treatment — it confirms the heated tanks actually held ~+3 °C
#   (31 vs 28 °C) for the whole experiment, not just at the daily YSI spot check
#   (code/09). The catch: these logs are huge (a reading every minute or two,
#   for every probe, for weeks). Loading them all into memory at once would
#   blow up R, so we parse one file at a time, immediately collapse each file to
#   hourly means, throw away the raw records, then move to the next file. The
#   hourly series is rolled up again to daily means for the trend figures.
# Input:   data/raw/apex/datalog*.xml
# Output:  data/processed/apex_temperature.rds  (hourly per probe)
#          data/processed/apex_temperature_daily.rds (daily per probe)
#          figures/08_apex_temperature.{pdf,png}
# =============================================================================

# ---- Setup -----------------------------------------------------------------
# 00_setup.R loads packages and defines shared paths (DATA_RAW, DATA_PROC, ...),
# theme_pub(), and save_fig(). xml2 is the parser we use to walk the Apex XML.
source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages(library(xml2))

# Collect every Apex datalog in the raw folder. stopifnot() halts immediately
# with a clear error if none are found, rather than failing cryptically later.
xml_files <- list.files(file.path(DATA_RAW, "apex"),
                        pattern = "\\.xml$", full.names = TRUE)
stopifnot(length(xml_files) > 0)

# ---- Parser: one XML file -> hourly mean per probe -------------------------
# This is the memory-bounded workhorse. It reads a single Apex file, pulls out
# every probe reading, then immediately summarises to hourly means so the raw
# (minute-resolution) records can be discarded before the next file is read.
# Returns tibble(datetime [hour], probe, value_mean, value_sd, n).
parse_apex_hourly <- function(path) {
  message("Parsing ", basename(path), " (", round(file.size(path) / 1e6, 1), " MB)")
  # RECOVER/NOWARNING/NOERROR tell libxml2 to salvage what it can from a
  # malformed/truncated log instead of aborting; tryCatch returns NULL on a
  # hard failure so one bad file can't kill the whole run.
  doc <- tryCatch(
    read_xml(path, options = c("RECOVER", "NOWARNING", "NOERROR")),
    error = function(e) NULL
  )
  if (is.null(doc)) return(tibble())

  # Each <record> is one poll of all probes at one timestamp.
  records <- xml_find_all(doc, "//record")
  if (length(records) == 0) return(tibble())

  # Extract record dates en masse (vectorised over all records at once).
  # Apex timestamps are US-format mdy_hms; tz set to Tahiti (the field site) so
  # hours line up with local day/night, not UTC.
  dates_text <- xml_text(xml_find_first(records, "./date"))
  dates <- suppressWarnings(mdy_hms(dates_text, tz = "Pacific/Tahiti"))

  # For each record extract probe name+value vectors; collapse to one row per
  # (record, probe). Records with an unparseable date or no probes are skipped.
  out_list <- vector("list", length(records))
  for (i in seq_along(records)) {
    if (is.na(dates[i])) next
    probes <- xml_find_all(records[[i]], "./probe")
    if (length(probes) == 0) next
    nms <- xml_text(xml_find_first(probes, "./name"))
    # Probe values arrive as text; as.numeric coerces, with non-numeric -> NA.
    vals <- suppressWarnings(as.numeric(xml_text(xml_find_first(probes, "./value"))))
    # Convert °F -> °C at the RAW reading level, for TEMPERATURE probes only. The
    # Apex ran US firmware and logged °F on some units (switched mid-deployment);
    # a seawater temperature reading > 60 can only be Fahrenheit. Doing this BEFORE
    # hourly/daily aggregation is what matters: averaging first would mix °F and °C
    # on any firmware-switch day and corrupt that day's mean. We restrict to Temp
    # probes so non-temp channels (ORP/pH, whose values can legitimately exceed 60)
    # are never mis-converted.
    is_temp <- grepl("^Temp\\d+$", trimws(nms))
    vals <- if_else(is_temp & !is.na(vals) & vals > 60, (vals - 32) * 5 / 9, vals)
    out_list[[i]] <- tibble(datetime = dates[i], probe = nms, value = vals)
  }

  # Free the parsed document and force garbage collection before returning —
  # this is what keeps total memory flat across hundreds of large files.
  records <- NULL; doc <- NULL; gc()

  # Stack all records, drop missing values, bin to the hour, and summarise.
  bind_rows(out_list) |>
    filter(!is.na(value)) |>
    mutate(hour = lubridate::floor_date(datetime, "hour")) |>  # truncate to top of hour
    group_by(hour, probe) |>
    summarise(value_mean = mean(value, na.rm = TRUE),
              value_sd   = sd(value, na.rm = TRUE),
              n          = n(),                                 # readings per hour
              .groups = "drop") |>
    rename(datetime = hour)
}

# Run the parser over every file and row-bind the hourly results into one table.
hourly <- map_dfr(xml_files, parse_apex_hourly)
saveRDS(hourly, file.path(DATA_PROC, "apex_temperature.rds"))

# ---- Daily aggregation for plotting ---------------------------------------
# Roll the hourly series up to one mean per probe per day — smooth enough to
# read the heat ramp and steady state on a multi-week figure. str_squish()
# trims stray whitespace in probe names so e.g. " Temp4" and "Temp4" group as one.
daily <- hourly |>
  mutate(date = as_date(datetime),
         probe = str_squish(probe)) |>
  group_by(date, probe) |>
  summarise(value_mean = mean(value_mean, na.rm = TRUE),
            value_sd   = mean(value_sd, na.rm = TRUE),
            .groups = "drop")
saveRDS(daily, file.path(DATA_PROC, "apex_temperature_daily.rds"))

# ---- Filter to temperature-like probes ------------------------------------
# The controller logs many probe types (pH, ORP, etc.); keep only the tank
# temperature probes, which Apex names `Temp1`-`Temp12` (one per outlet/tank).
# The regex "^Temp\\d+$" matches exactly that pattern (Temp + digits, nothing else).
temp_daily <- daily |>
  filter(str_detect(probe, "^Temp\\d+$"),
         is.finite(value_mean), value_mean > 0) |>
  mutate(
    # Pull the trailing number out of the probe name to recover the tank ID.
    tank = as.integer(str_extract(probe, "\\d+")),
    # value_mean is ALREADY in °C — the °F->°C conversion happens at the raw-
    # reading level inside parse_apex_hourly(), so the daily means never mix units.
    value_c = value_mean,
    # Map tank -> treatment via the shared plumbing layout (tank_treatment() in
    # 00_setup.R; same assignment used by 09 (YSI) and 16 (Fig 1 panel A)).
    treatment = tank_treatment(tank)
  ) |>
  # Keep experimental tanks only, and drop physically impossible temperatures
  # (sensor dropouts/air exposure) by bounding to a plausible 15-40 °C window.
  filter(!is.na(treatment), value_c > 15, value_c < 40)

# ---- Plot: experimental window --------------------------------------------
# Main figure. Trim to the experiment dates so the steady-state contrast is the
# focus (the full logger record, with ramps/cooldowns, is the companion below).
exp_window <- as_date(c("2025-05-25", "2025-06-25"))
d_plot <- temp_daily |>
  filter(between(date, exp_window[1], exp_window[2]))

# One line per tank, coloured by treatment. Dashed reference lines mark the two
# target set-points (28 and 31 °C) so the reader can see how tightly tanks held.
p_apex <- ggplot(d_plot, aes(date, value_c, group = tank,
                              colour = treatment)) +
  geom_hline(yintercept = c(28, 31), linetype = "dashed",
             colour = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.4, alpha = 0.85) +
  geom_point(size = 1.0, alpha = 0.75) +
  # Direct-label the two target lines at the right edge instead of adding a
  # second legend (inherit.aes = FALSE so this layer ignores the main aes()).
  geom_text(data = data.frame(date = exp_window[2], y = c(28, 31),
                              lab = c("28 °C target", "31 °C target")),
            aes(x = date, y = y, label = lab),
            inherit.aes = FALSE, hjust = 1.1, vjust = -0.4,
            size = 2.6, colour = "grey40") +
  scale_colour_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),  # blue=ambient, orange=heated
                      name = "Treatment") +
  scale_x_date(date_breaks = "5 days", date_labels = "%b %d") +
  labs(x = NULL, y = "Daily mean tank T (°C)",
       title = "Apex tank temperatures across the experimental window",
       subtitle = "One line per tank; heat ramp begins late May, sustained 31 °C through experiment") +
  theme_pub(10)

save_fig(p_apex, "08_apex_temperature", width = 170, height = 100)

# ---- Plot: full datalog period (companion) --------------------------------
# Same data, no date filter — shows the whole record including the initial ramp
# up and the post-experiment cooldown, for the supplement / sanity checking.
p_full <- ggplot(temp_daily, aes(date, value_c, group = tank,
                                  colour = treatment)) +
  geom_hline(yintercept = c(28, 31), linetype = "dashed",
             colour = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.3, alpha = 0.75) +
  scale_colour_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                      name = "Treatment") +
  labs(x = NULL, y = "Daily mean tank T (°C)",
       title = "Apex tank temperatures — full datalog period",
       subtitle = "Tanks 3, 6, 9, 12 = 28 °C; tanks 4, 5, 10, 11 = 31 °C") +
  theme_pub(10)
save_fig(p_full, "08b_apex_temperature_full", width = 200, height = 100)

# ---- Console report --------------------------------------------------------
# Print the most-logged probes as a quick sanity check that the expected Temp
# probes dominate the file (and to spot any unexpected probe names).
cat("\nProbes detected (top 20 by record count):\n")
hourly |> count(probe, sort = TRUE) |> head(20) |> print()

cat("\nWrote apex_temperature.rds (", nrow(hourly), " hourly records),",
    " apex_temperature_daily.rds (", nrow(daily), "),",
    " 08_apex_temperature.{pdf,png}\n", sep = "")
