# =============================================================================
# 01_clean_data.R: turn every raw coral spreadsheet into a clean table
#
# Reads each raw spreadsheet in data/raw/ and writes one clean .rds file per
# data stream to data/processed/. There are no statistics here; the later
# scripts read these files. Cleaning means: standard column names, the right
# type for each column (numbers, dates, factors), lower-case source patch codes,
# and the same 28C/31C and no/yes factor levels in every table so they join
# cleanly. Fragment 119's wound label is corrected at import by
# correct_wound_labels() in 00_setup.R (confirmed by Molly, 2026-10-19).
# Tank temperature and water quality are cleaned in 02_tanks_water.R.
#
# Sections (raw input -> clean output, both relative to data/raw/ and
# data/processed/)
#   1. Fragment metadata       metadata/metadata.csv     -> coral_metadata.rds (+ .csv)
#   2. Photosynthesis (Fv/Fm)  pam/PAM_data.csv          -> pam_clean.rds
#   3. Color score             color_card/data.csv       -> color_clean.rds
#   4. Regeneration scores     physio_morphology/data.csv -> physio_clean.rds
#   5. Growth                  buoyant_weight/data.csv   -> buoyant_weight_clean.rds
#   6. Symbiont density        symbiont_counts/*.csv     -> symbiont_chl_clean.rds
#   7. Surface area (wax)      wax_dipping/*.csv         -> wax_clean.rds
#   8. Microscope experiment   microscope_physio/data.csv -> microscope_physio_clean.rds
#   9. Flatworm checks         worm_presence/Sheet1.csv  -> worm_clean.rds
#
# Reads:  data/raw/ (the files listed above)
# Writes: data/processed/*.rds (listed above), data/processed/coral_metadata.csv
#         output/tables/01_metadata_summary.csv, 07_wax_curve_fit.csv,
#           11_microscope_design_crosswalk.csv, 10_worm_summary.csv
#         figures/07_wax_standard_curve
# Run:    sourced by code/run_all.R (or: Rscript code/run_all.R)
# =============================================================================

source(here::here("code", "00_setup.R"))

# Make sure the output folders exist before anything is written.
dir.create(DATA_PROC, recursive = TRUE, showWarnings = FALSE)
dir.create(TBL_DIR, recursive = TRUE, showWarnings = FALSE)


# ---- 1. Fragment metadata --------------------------------------------------
# One row per fragment: id, source patch (`thicket` = a/c/d), temperature,
# wound status, tank and biopsy day. Every other data stream joins back to
# this table by fragment id, so it is cleaned first.

# suppressWarnings() hides readr's column-type guessing messages; guess_max
# makes readr look at more rows before deciding each column's type.
meta_raw <- suppressWarnings(
  read_csv(file.path(DATA_RAW, "metadata", "metadata.csv"),
           show_col_types = FALSE, guess_max = 5000)
) |>
  janitor::clean_names() |>
  filter(!is.na(species) & !is.na(id))          # drop blank trailing rows

# Temperature is stored as the number 28 or 31; it becomes the factor
# 28C/31C used everywhere downstream. Chlorophyll and symbiont columns lose
# the awkward "_cm_2" suffix that clean_names() produced.
meta <- meta_raw |>
  rename(thicket = matches("^thicket")) |>      # header is "thicket " with a space
  mutate(
    species        = str_squish(species),
    thicket        = str_to_lower(str_squish(thicket)),
    id             = as.integer(id),
    sample         = str_squish(sample),
    wound          = factor(correct_wound_labels(id, wound), levels = c("no", "yes")),
    tank           = as.integer(tank),
    treatment_c    = as.integer(treatment),
    treatment      = factor(treatment_c, levels = c(28, 31), labels = c("28C", "31C")),
    biopsy_day     = as.integer(biopsy_day),
    collection_date = as_date(collection_date),
    biopsy_date    = as_date(biopsy_date),
    coord_lat      = as.numeric(coord_lat),
    coord_long     = as.numeric(coord_long),
    chlorophyll_ug_cm2 = as.numeric(chlorophyll_ug_cm_2),
    zooxanthellae_cells_cm2 = as.numeric(zooxanthellae_cells_cm_2),
    percent_growth_bw = as.numeric(percent_growth_bw),
    calculated_sa  = as.numeric(calculated_sa)
  ) |>
  select(-treatment_c, -ends_with("_cm_2"))
message(glue::glue("Loaded metadata: {nrow(meta)} rows"))

# Save the clean table, plus a count of fragments in each design cell
# (temperature x wound x source patch) as a quick check of the design.
saveRDS(meta, file.path(DATA_PROC, "coral_metadata.rds"))
write_csv(meta, file.path(DATA_PROC, "coral_metadata.csv"))
write_csv(meta |> count(treatment, wound, thicket, name = "n") |>
            arrange(treatment, wound, thicket),
          file.path(TBL_DIR, "01_metadata_summary.csv"))


# ---- 2. Photosynthesis (Fv/Fm) ---------------------------------------------
# Fv/Fm is the maximum efficiency of photosystem II, measured with a PAM
# fluorometer: about 0.6-0.7 in a healthy fragment, falling under heat stress.
# Some Fv/Fm cells were spreadsheet formulas that export as text; for those
# rows Fv/Fm is recomputed from the raw Y reading (Y / 1000). Values outside
# (0, 1) are reading errors and are dropped.
pam <- read_csv(file.path(DATA_RAW, "pam", "PAM_data.csv"), show_col_types = FALSE) |>
  janitor::clean_names() |>
  rename(thicket = matches("^thicket"), fv_fm = matches("^fv_fm")) |>
  mutate(
    date      = as_date(date),
    day       = as.integer(day),
    treatment = factor(as.integer(treatment), levels = c(28, 31), labels = c("28C", "31C")),
    tank      = as.integer(tank),
    wound     = factor(correct_wound_labels(id, wound), levels = c("no", "yes")),
    thicket   = str_to_lower(str_squish(thicket)),
    id        = as.integer(id),
    location  = str_to_lower(location),
    f = as.numeric(f), m = as.numeric(m), y = as.numeric(y), e = as.numeric(e),
    # If the sheet's Fv/Fm will not convert to a number, fall back to Y / 1000.
    fv_fm     = if_else(is.na(suppressWarnings(as.numeric(fv_fm))),
                        y / 1000, suppressWarnings(as.numeric(fv_fm)))
  ) |>
  filter(!is.na(fv_fm), fv_fm > 0, fv_fm < 1)

# The two probe positions per fragment (top and base) are averaged to one
# value per fragment per day, so each fragment counts once in the models.
pam_avg <- pam |>
  group_by(date, day, treatment, tank, thicket, wound, id) |>
  summarise(fv_fm = mean(fv_fm, na.rm = TRUE), .groups = "drop")
saveRDS(pam_avg, file.path(DATA_PROC, "pam_clean.rds"))


# ---- 3. Color score --------------------------------------------------------
# Fragment color scored against the CoralWatch color card (1 = palest,
# 6 = darkest), a field proxy for bleaching. Card codes such as "D4" or
# "D3/D4" become numbers: the mean of the digits, so "D3/D4" = 3.5.

# str_extract_all() pulls out every run of digits; vapply() then averages
# them per cell and returns NA when a cell has no digits at all.
convert_color <- function(x) {
  nums <- stringr::str_extract_all(x, "\\d+")
  vapply(nums, function(v) if (length(v) == 0) NA_real_ else mean(as.numeric(v)),
         numeric(1))
}

cc <- read_csv(file.path(DATA_RAW, "color_card", "data.csv"), show_col_types = FALSE) |>
  janitor::clean_names() |>
  rename(thicket = matches("^thicket"), wound = wounded) |>
  mutate(
    date      = as_date(date),
    day       = as.integer(day),
    treatment = factor(as.integer(treatment), levels = c(28, 31), labels = c("28C", "31C")),
    tank      = as.integer(tank),
    thicket   = str_to_lower(str_squish(thicket)),
    id        = as.integer(id),
    wound     = factor(correct_wound_labels(id, wound), levels = c("no", "yes")),
    health_status  = str_to_lower(str_squish(health_status)),
    paling         = factor(str_to_lower(str_squish(paling)), levels = c("no", "yes")),
    hole_at_center = factor(str_to_lower(str_squish(hole_at_center)), levels = c("no", "yes")),
    color_num = convert_color(color)
  ) |>
  filter(!is.na(color_num))

# The color sheet records every fragment as unwounded at the Day -4 baseline
# (before clipping). Use each fragment's assigned wound group from the metadata
# so the baseline is grouped like every other day (as in the Fv/Fm data).
# stopifnot() stops the run if a check fails: here, every color-sheet
# fragment must exist in the metadata, so the join below cannot leave NAs.
stopifnot(all(cc$id %in% meta$id))
cc <- cc |>
  select(-wound) |>
  left_join(meta |> distinct(id, wound), by = "id") |>
  relocate(wound, .after = id)
saveRDS(cc, file.path(DATA_PROC, "color_clean.rds"))


# ---- 4. Daily regeneration scores (main experiment) ------------------------
# Each regeneration stage (polyps out, wound smoothed, pigment over the wound,
# new tip, ...) is scored yes/no per fragment per day; stored here as 1/0 so
# later scripts can sum and average them. "hole_in_center" and "polyp_in_hole"
# are the same score in the sheet; the check below stops the run if they ever
# differ.

# The date column mixes true dates, date-times and text, depending on how
# each row was entered. This helper turns any of them into a Date.
parse_lth_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(lubridate::as_date(x))
  lubridate::as_date(lubridate::parse_date_time(
    as.character(x), orders = c("mdy HMS", "mdy", "ymd HMS", "ymd"), quiet = TRUE))
}

# across(..., ~ ...) applies the same yes -> 1, anything else -> 0 rule to
# every stage column at once.
ph <- suppressWarnings(
  read_csv(file.path(DATA_RAW, "physio_morphology", "data.csv"),
           show_col_types = FALSE, guess_max = 2000)
) |>
  janitor::clean_names() |>
  filter(!is.na(species) & !is.na(id)) |>
  rename(thicket = matches("^thicket"), wound = wounded) |>
  mutate(
    date      = parse_lth_date(date),
    day       = as.integer(day),
    treatment = factor(as.integer(treatment), levels = c(28, 31), labels = c("28C", "31C")),
    tank      = as.integer(tank),
    thicket   = str_to_lower(str_squish(thicket)),
    id        = as.integer(id),
    wound     = factor(correct_wound_labels(id, wound), levels = c("no", "yes")),
    across(c(polyps_out, hole_in_center, polyp_in_hole, wound_smoothed,
             pigment_over_wound, tip_exist, tip_extension,
             new_corallites_on_tip, algae_on_wound),
           ~ as.integer(str_to_lower(str_squish(.x)) == "yes"))
  )

# Guarantee the two duplicate columns are identical, then keep one copy under
# a clearer name: the axial (tip) polyp forming in the cut.
stopifnot(identical(ph$hole_in_center, ph$polyp_in_hole))
ph <- ph |> mutate(axial_polyp_formation = hole_in_center)
saveRDS(ph, file.path(DATA_PROC, "physio_clean.rds"))


# ---- 5. Growth (buoyant weight) --------------------------------------------
# Buoyant weight is a fragment's weight while submerged in seawater. It is
# converted to dry skeletal mass with the buoyant-weight method (Davies 1989),
# following Molly's calculation sheet (`calcuations` in
# Buoyant_Weight_-_complete.xlsx) step by step:
#   1. Freshwater density at the freshwater temperature (the quadratic below).
#   2. Density of a glass reference object, from its weight in air and in
#      freshwater.
#   3. Seawater density, from the same object's weight in seawater. This is
#      measured at every weighing, so it tracks the real salinity and
#      temperature of the tank water.
#   4. Dry mass = buoyant weight / (1 - seawater density / aragonite density),
#      with aragonite density 2.93 g/cm3.
#   5. Each fragment is glued to an aragonite plug, weighed separately. The
#      plug's dry mass is subtracted from the Day 0 and Day 15 masses, so
#      growth is measured on the coral's own skeleton (plugs are ~60% of the
#      Day 0 total; see metadata codebook, percent_growth_bw).
# The reported growth measure is % change in skeletal mass, Day 0 to Day 15.
#
# Fixed 2026-11-13: earlier versions used the freshwater curve as if it were
# seawater density and did not subtract the plug. The plug omission shrank
# % growth about 2.6-fold (e.g. 6.1% vs 16.0% at 28 C). The calculation now
# reproduces Molly's sheet and the metadata column percent_growth_bw.

# Freshwater density (g/cm3) as a function of temperature (degrees C).
rho_fw_temp <- function(t_c) -5e-6 * t_c^2 + 7e-6 * t_c + 1.0001
ARAG <- 2.93

# Seawater density from the glass reference object: its density comes from
# its weights in air and freshwater; seawater density then follows from its
# weight in seawater.
rho_sw_ref <- function(ref_dry, ref_fw, ref_sw, t_fw) {
  rho_ref <- ref_dry * rho_fw_temp(t_fw) / (ref_dry - ref_fw)
  rho_ref * (ref_dry - ref_sw) / ref_dry
}

# Buoyant weight in seawater -> dry aragonite mass.
dry_mass <- function(buoyant_w, rho_sw) buoyant_w / (1 - rho_sw / ARAG)

bw <- read_csv(file.path(DATA_RAW, "buoyant_weight", "data.csv"), show_col_types = FALSE) |>
  janitor::clean_names() |>
  rename(thicket = matches("^thicket")) |>
  mutate(
    thicket    = str_to_lower(str_squish(thicket)),
    id         = as.integer(id),
    wound      = factor(correct_wound_labels(id, wound), levels = c("no", "yes")),
    tank       = as.integer(tank),
    treatment  = factor(as.integer(treatment), levels = c(28, 31), labels = c("28C", "31C")),
    biopsy_day = as.integer(biopsy_day),
    initial_w  = as.numeric(initial_coral_buoyant_weigh_uncorrected),
    final_w    = as.numeric(final_coral_weight_uncorrected),
    initial_temp = as.numeric(initial_temp_coral),
    final_temp   = as.numeric(final_temp_coral)
  ) |>
  filter(!is.na(initial_w), !is.na(final_w)) |>
  # Seawater density at each of the three weighings (plug, Day 0, Day 15).
  # The same glass reference object (plug_ref_dry, in air) was used for all.
  mutate(
    plug_rho_sw    = rho_sw_ref(plug_ref_dry, plug_ref_fw, plug_ref_sw, temp_fw),
    initial_rho_sw = rho_sw_ref(plug_ref_dry, initial_ref_fw, initial_ref_sw,
                                initial_temp_fw),
    final_rho_sw   = rho_sw_ref(plug_ref_dry, final_ref_fw, final_ref_sw,
                                final_temp_fw)
  ) |>
  # Dry masses with the plug removed. Three growth measures follow: % change
  # (the one reported), specific growth rate (sgr, % per day on the log scale)
  # and grams per day. pmax(biopsy_day, 1) guards against dividing by zero days.
  mutate(
    plug_dry    = dry_mass(plug_buoyant_weight_uncorrected, plug_rho_sw),
    initial_dry = dry_mass(initial_w, initial_rho_sw) - plug_dry,
    final_dry   = dry_mass(final_w, final_rho_sw) - plug_dry,
    delta_g     = final_dry - initial_dry,
    days        = pmax(biopsy_day, 1),
    pct_growth  = 100 * delta_g / initial_dry,
    sgr         = 100 * (log(final_dry) - log(initial_dry)) / days,
    g_per_day   = delta_g / days
  )
# Guarantee the calculation matches the growth Molly reported in the
# metadata sheet (percent_growth_bw, plug subtracted) for every fragment.
bw_check <- bw |>
  select(id, pct_growth) |>
  inner_join(select(meta, id, percent_growth_bw), by = "id")
stopifnot(nrow(bw_check) == nrow(bw),
          all(abs(bw_check$pct_growth - bw_check$percent_growth_bw) < 1e-3))

saveRDS(bw, file.path(DATA_PROC, "buoyant_weight_clean.rds"))


# ---- 6. Symbiont density ---------------------------------------------------
# Symbionts (zooxanthellae) were counted on a hemocytometer, a counting slide
# with a grid of four quadrants. Average the four quadrants within each count,
# then the replicate counts within each fragment. Where raw counts are missing,
# the spreadsheet's average is used (count_source says which).
# cells per cm2 = mean count x 1e4 x slurry volume (mL) / surface area (cm2).

## 6a. Mean hemocytometer count per fragment ---------------------------------
# rowwise() makes mean() work across the four quadrant columns within a row;
# ungroup() switches that off again.
reps <- suppressWarnings(
  read_csv(file.path(DATA_RAW, "symbiont_counts", "Raw_counts.csv"),
           show_col_types = FALSE, guess_max = 2000)
) |>
  janitor::clean_names() |>
  rename(coral_id = coral_id_3) |>
  filter(!is.na(coral_id)) |>
  mutate(coral_id = as.integer(coral_id), across(c(q1, q2, q3, q4), as.numeric)) |>
  rowwise() |>
  mutate(quad_mean = mean(c(q1, q2, q3, q4), na.rm = TRUE)) |>
  ungroup() |>
  group_by(coral_id) |>
  summarise(mean_quad_count = mean(quad_mean, na.rm = TRUE), n_reps = n(), .groups = "drop") |>
  filter(is.finite(mean_quad_count))

## 6b. Cells per cm2, preferring raw counts over the sheet average -----------
# coalesce() takes the first non-missing value: the raw-count mean if there
# is one, otherwise the average typed into the spreadsheet.
zoox <- suppressWarnings(
  read_csv(file.path(DATA_RAW, "symbiont_counts", "metadata_ordered_merge.csv"),
           show_col_types = FALSE, guess_max = 2000)
) |>
  janitor::clean_names() |>
  rename(thicket = matches("^thicket"), id = id_3) |>
  mutate(
    id         = as.integer(id),
    treatment  = factor(as.integer(treatment), levels = c(28, 31), labels = c("28C", "31C")),
    biopsy_day = as.integer(biopsy_day),
    thicket    = str_to_lower(str_squish(thicket)),
    wound      = factor(correct_wound_labels(id, wound), levels = c("no", "yes")),
    sa_cm2     = as.numeric(calculated_sa_standard_curve),
    slurry_ml  = as.numeric(slurry_volume_m_l),
    zoox_avg_sheet = as.numeric(average),
    cells_per_cm2_sheet = (zoox_avg_sheet * 1e4 * slurry_ml) / sa_cm2
  ) |>
  filter(!is.na(id)) |>
  left_join(reps, by = c("id" = "coral_id")) |>
  mutate(
    zoox_avg_hemo = coalesce(mean_quad_count, zoox_avg_sheet),
    count_source  = if_else(!is.na(mean_quad_count), "raw_counts", "metadata_average"),
    cells_per_cm2 = (zoox_avg_hemo * 1e4 * slurry_ml) / sa_cm2,
    tank = as.integer(tank)
  ) |>
  select(id, treatment, wound, biopsy_day, thicket, tank, sa_cm2,
         cells_per_cm2, count_source, n_reps)

## 6c. Add chlorophyll from the metadata -------------------------------------
sym <- zoox |>
  left_join(meta |> filter(!is.na(chlorophyll_ug_cm2)) |> select(id, chlorophyll_ug_cm2),
            by = "id")
saveRDS(sym, file.path(DATA_PROC, "symbiont_chl_clean.rds"))


# ---- 7. Surface area (wax dipping) -----------------------------------------
# Fragment surface area comes from wax dipping: the mass of wax a fragment
# picks up is proportional to its surface area. A standard curve, fitted on
# calibration cylinders of known surface area, converts wax mass to cm2.

## 7a. Standard curve: wax mass on cylinders vs their known surface area -----
# Cylinder side area = 2 * pi * radius * height; dividing by 100 turns mm2
# into cm2. Only positive wax gains and areas are kept.
sc_fit <- suppressWarnings(
  read_csv(file.path(DATA_RAW, "wax_dipping", "Standard_curve.csv"),
           show_col_types = FALSE, guess_max = 2000)
) |>
  janitor::clean_names() |>
  mutate(
    height = as.numeric(height_mm), diameter = as.numeric(diameter_mm),
    dry_g = as.numeric(dry_weight_g), wax_g = as.numeric(wax_weight_g),
    wax_minus_dry = wax_g - dry_g,
    sa_cm2 = 2 * pi * (diameter / 2) * height / 100
  ) |>
  filter(!is.na(wax_minus_dry), !is.na(sa_cm2), wax_minus_dry > 0, sa_cm2 > 0) |>
  transmute(wax_g = wax_minus_dry, sa = sa_cm2)

# Linear fit of surface area on wax mass; broom::tidy() saves the slope and
# intercept with 95% confidence intervals.
cal <- lm(sa ~ wax_g, data = sc_fit)
write_csv(broom::tidy(cal, conf.int = TRUE), file.path(TBL_DIR, "07_wax_curve_fit.csv"))

## 7b. Apply the curve to the experimental fragments -------------------------
# Wax gained = second weighing minus first. Two areas are kept: one from the
# standard curve (sa_curve_cm2) and one treating the fragment as a cylinder
# from caliper measurements (sa_caliper_cm2), as a cross-check.
wax <- suppressWarnings(
  read_csv(file.path(DATA_RAW, "wax_dipping", "data.csv"),
           show_col_types = FALSE, guess_max = 2000)
) |>
  janitor::clean_names() |>
  filter(!is.na(species) & !is.na(id)) |>
  rename(thicket = matches("^thicket")) |>
  mutate(
    id = as.integer(id),
    # Fragment 119 is wounded (see correct_wound_labels() in 00_setup.R).
    wound = factor(correct_wound_labels(id, wound), levels = c("no", "yes")),
    treatment = factor(as.integer(treatment), levels = c(28, 31), labels = c("28C", "31C")),
    biopsy_day = as.integer(biopsy_day),
    diameter_mm = as.numeric(diameter), height_mm = as.numeric(height),
    dry_g = as.numeric(dry_weight),
    wax1_g = as.numeric(wax_weight_1), wax2_g = as.numeric(wax_weight_2),
    wax_g = wax2_g - wax1_g,
    sa_caliper_cm2 = 2 * pi * (diameter_mm / 2) * height_mm / 100,
    sa_curve_cm2   = predict(cal, newdata = tibble(wax_g = wax_g))
  ) |>
  select(id, treatment, biopsy_day, thicket, wound, dry_g, wax_g,
         sa_caliper_cm2, sa_curve_cm2)
saveRDS(wax, file.path(DATA_PROC, "wax_clean.rds"))

## 7c. Figure: the standard curve --------------------------------------------
p_sc <- ggplot(sc_fit, aes(wax_g, sa)) +
  geom_point(size = 2, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, colour = "#0072B2") +
  labs(x = "Wax mass (g)", y = expression(Surface~area~(cm^2))) +
  theme_pub(10)
save_fig(p_sc, "07_wax_standard_curve", width = 120, height = 95)


# ---- 8. Microscope experiment ----------------------------------------------
# A separate set of 16 wounded fragments photographed daily under a dissecting
# microscope (8 per temperature, source patches A and C). Only this experiment
# scored living tissue over the cut. Kept separate from the main experiment.
# This import is stricter than the others: it stops with an error if a column
# is missing or a yes/no cell holds anything unexpected.

## 8a. Read the export and check its columns ---------------------------------
micro_raw <- readr::read_csv(file.path(DATA_RAW, "microscope_physio", "data.csv"),
                             show_col_types = FALSE, guess_max = 1000, na = c("", "NA")) |>
  janitor::clean_names()

micro_traits <- c("tissue_over_wound", "hole_in_center", "polyp_in_hole", "wound_smoothed",
                  "pigment_over_wound", "tip_exist", "tip_extension", "new_corallites_on_tip")

missing_cols <- setdiff(c("thicket", "wound", "tank", "sample", "treatment", "id", "date",
                          "day", micro_traits, "use_photo", "notes"), names(micro_raw))
if (length(missing_cols)) stop("Microscope export missing columns: ",
                               paste(missing_cols, collapse = ", "), call. = FALSE)

## 8b. Convert yes/no scores to 1/0, refusing unexpected codes ---------------
# Blank, "na" and "nan" become NA (not scored); "yes" = 1, "no" = 0; any
# other code stops the run and names the column it came from.
yn_to_int <- function(x, column) {
  y <- stringr::str_to_lower(stringr::str_squish(as.character(x)))
  y[y %in% c("", "na", "nan")] <- NA_character_
  unexpected <- setdiff(stats::na.omit(unique(y)), c("yes", "no"))
  if (length(unexpected)) stop(sprintf("Unexpected yes/no code in %s: %s", column,
                                       paste(unexpected, collapse = ", ")), call. = FALSE)
  ifelse(is.na(y), NA_integer_, as.integer(y == "yes"))
}

# dplyr::cur_column() passes the name of the column being converted, so an
# error message can say which column held the bad code.
micro <- micro_raw |>
  mutate(
    date = lubridate::mdy(date),
    day = as.integer(day),
    treatment = factor(paste0(as.integer(treatment), "C"), levels = c("28C", "31C")),
    tank = as.integer(tank),
    thicket = stringr::str_to_lower(stringr::str_squish(thicket)),
    id = as.integer(id),
    wound = factor(stringr::str_to_lower(stringr::str_squish(wound)), levels = c("no", "yes")),
    sample = stringr::str_to_lower(stringr::str_squish(sample)),
    notes = dplyr::na_if(as.character(notes), "NA"),
    across(all_of(c(micro_traits, "use_photo")), ~ yn_to_int(.x, dplyr::cur_column()))
  )

# As in Section 4, the two duplicate columns must be identical.
if (!identical(micro$hole_in_center, micro$polyp_in_hole)) {
  stop("hole_in_center and polyp_in_hole differ in the microscope data.", call. = FALSE)
}
micro <- micro |> mutate(axial_polyp_formation = hole_in_center)
saveRDS(micro, file.path(DATA_PROC, "microscope_physio_clean.rds"))

## 8c. Check every microscope fragment against the master metadata -----------
# One row per microscope fragment, with its tank, temperature, source patch,
# wound and sample type from both sources side by side. matches_metadata is
# TRUE only when all five agree. A mismatch gives a warning, not an error,
# so the rest of the pipeline still runs; the table shows which ids differ.
crosswalk <- micro |>
  distinct(id, tank, treatment, thicket, wound, sample) |>
  mutate(treatment = as.character(treatment), wound = as.character(wound)) |>
  left_join(meta |> filter(sample == "photo") |>
              transmute(id,
                        meta_tank = as.integer(tank),
                        meta_treatment = as.character(treatment),
                        meta_thicket = as.character(thicket),
                        meta_wound = as.character(wound),
                        meta_sample = as.character(sample)),
            by = "id") |>
  mutate(matches_metadata = !is.na(meta_tank) & tank == meta_tank &
           treatment == meta_treatment & thicket == meta_thicket &
           wound == meta_wound & sample == meta_sample) |>
  arrange(treatment, thicket, tank, id)
write_csv(crosswalk, file.path(TBL_DIR, "11_microscope_design_crosswalk.csv"))

if (any(!crosswalk$matches_metadata)) {
  warning("Microscope IDs do not all match the metadata; see 11_microscope_design_crosswalk.csv",
          call. = FALSE)
}


# ---- 9. Flatworm checks ----------------------------------------------------
# Acropora-eating flatworm counts from three inspections (methods detail).
# The sheet is wide (one worms_<date> column per inspection);
# pivot_longer() turns it into one row per fragment per inspection.
# Temperature comes from the tank number via tank_treatment() in 00_setup.R;
# rows from tanks outside the experiment get NA and are dropped.
worms <- read_csv(file.path(DATA_RAW, "worm_presence", "Sheet1.csv"),
                  show_col_types = FALSE) |>
  janitor::clean_names() |>
  rename(thicket = matches("^thicket")) |>
  mutate(thicket = str_to_lower(str_squish(thicket)), tank = as.integer(tank),
         id = as.integer(id), treatment = tank_treatment(tank)) |>
  pivot_longer(starts_with("worms_"), names_to = "date_raw", values_to = "n_worms") |>
  mutate(date = mdy(str_remove(date_raw, "^worms_")), present = as.integer(n_worms > 0)) |>
  filter(!is.na(treatment))
saveRDS(worms, file.path(DATA_PROC, "worm_clean.rds"))

# Share of fragments carrying at least one worm, per inspection and tank.
write_csv(worms |> group_by(date, treatment, tank) |>
            summarise(n_corals = n(), n_with_worms = sum(present, na.rm = TRUE),
                      pct = 100 * n_with_worms / n_corals, .groups = "drop"),
          file.path(TBL_DIR, "10_worm_summary.csv"))

cat("01_clean_data: wrote", length(list.files(DATA_PROC, "\\.rds$")), "clean .rds files\n")
