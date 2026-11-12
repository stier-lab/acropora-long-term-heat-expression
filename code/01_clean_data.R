# =============================================================================
# 01_clean_data.R: turn every raw coral spreadsheet into a clean table
#
# Reads data/raw/ and writes one clean .rds per data stream to data/processed/.
# No statistics here; later scripts read these files.
#
#   Section                    Raw input                                  Clean output
#   1 Fragment metadata        metadata/metadata.csv                      coral_metadata.rds (+ .csv)
#   2 Photosynthesis (Fv/Fm)   pam/PAM_data.csv                           pam_clean.rds
#   3 Color score              color_card/data.csv                        color_clean.rds
#   4 Regeneration scores      physio_morphology/data.csv                 physio_clean.rds
#   5 Growth                   buoyant_weight/data.csv                    buoyant_weight_clean.rds
#   6 Symbiont density         symbiont_counts/*.csv                      symbiont_chl_clean.rds
#   7 Surface area (wax)       wax_dipping/*.csv                          wax_clean.rds
#   8 Microscope experiment    microscope_physio/data.csv                 microscope_physio_clean.rds
#   9 Flatworm checks          worm_presence/Sheet1.csv                   worm_clean.rds
#
# Tank temperature and water quality are cleaned in 02_tanks_water.R.
# Also writes: output/tables/01_metadata_summary.csv, 07_wax_curve_fit.csv,
#   11_microscope_design_crosswalk.csv, 10_worm_summary.csv;
#   figures/07_wax_standard_curve.
# =============================================================================

source(here::here("code", "00_setup.R"))
dir.create(DATA_PROC, recursive = TRUE, showWarnings = FALSE)
dir.create(TBL_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Fragment metadata --------------------------------------------------
# One row per fragment: id, source patch (`thicket` = a/c/d), temperature,
# wound status, tank and biopsy day. Every other stream joins back to this.
meta_raw <- suppressWarnings(
  read_csv(file.path(DATA_RAW, "metadata", "metadata.csv"),
           show_col_types = FALSE, guess_max = 5000)
) |>
  janitor::clean_names() |>
  filter(!is.na(species) & !is.na(id))          # drop blank trailing rows

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
    collection_date= as_date(collection_date),
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

saveRDS(meta, file.path(DATA_PROC, "coral_metadata.rds"))
write_csv(meta, file.path(DATA_PROC, "coral_metadata.csv"))
write_csv(meta |> count(treatment, wound, thicket, name = "n") |>
            arrange(treatment, wound, thicket),
          file.path(TBL_DIR, "01_metadata_summary.csv"))

# ---- 2. Photosynthesis (Fv/Fm) ---------------------------------------------
# Some Fv/Fm cells were spreadsheet formulas that export as text; for those
# rows Fv/Fm is recomputed from the raw Y reading (Y / 1000). Values outside
# (0, 1) are reading errors and are dropped. The two probe positions per
# fragment (top and base) are averaged to one value per fragment per day.
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
    fv_fm     = if_else(is.na(suppressWarnings(as.numeric(fv_fm))),
                        y / 1000, suppressWarnings(as.numeric(fv_fm)))
  ) |>
  filter(!is.na(fv_fm), fv_fm > 0, fv_fm < 1)

pam_avg <- pam |>
  group_by(date, day, treatment, tank, thicket, wound, id) |>
  summarise(fv_fm = mean(fv_fm, na.rm = TRUE), .groups = "drop")
saveRDS(pam_avg, file.path(DATA_PROC, "pam_clean.rds"))

# ---- 3. Color score --------------------------------------------------------
# CoralWatch card codes such as "D4" or "D3/D4" become numbers (the mean of
# the digits, so "D3/D4" = 3.5).
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
stopifnot(all(cc$id %in% meta$id))
cc <- cc |>
  select(-wound) |>
  left_join(meta |> distinct(id, wound), by = "id") |>
  relocate(wound, .after = id)
saveRDS(cc, file.path(DATA_PROC, "color_clean.rds"))

# ---- 4. Daily regeneration scores (main experiment) ------------------------
# Each stage is scored yes/no per fragment per day; stored here as 1/0.
# "hole_in_center" and "polyp_in_hole" are the same score in the sheet; the
# check below stops the run if they ever differ.
parse_lth_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(lubridate::as_date(x))
  lubridate::as_date(lubridate::parse_date_time(
    as.character(x), orders = c("mdy HMS", "mdy", "ymd HMS", "ymd"), quiet = TRUE))
}
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
    wound     = factor(wound, levels = c("no", "yes")),
    across(c(polyps_out, hole_in_center, polyp_in_hole, wound_smoothed,
             pigment_over_wound, tip_exist, tip_extension,
             new_corallites_on_tip, algae_on_wound),
           ~ as.integer(str_to_lower(str_squish(.x)) == "yes"))
  )
stopifnot(identical(ph$hole_in_center, ph$polyp_in_hole))
ph <- ph |> mutate(axial_polyp_formation = hole_in_center)
saveRDS(ph, file.path(DATA_PROC, "physio_clean.rds"))

# ---- 5. Growth (buoyant weight) --------------------------------------------
# Buoyant weight is converted to dry skeletal mass with seawater density at the
# measured water temperature and aragonite density 2.93 g/cm3 (Davies 1989).
# The reported growth measure is % change in skeletal mass, Day 0 to Day 15.
rho_sw_temp <- function(t_c) -5e-6 * t_c^2 + 7e-6 * t_c + 1.0001
ARAG <- 2.93
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
  mutate(
    initial_dry = initial_w / (1 - rho_sw_temp(initial_temp) / ARAG),
    final_dry   = final_w   / (1 - rho_sw_temp(final_temp)   / ARAG),
    delta_g     = final_dry - initial_dry,
    days        = pmax(biopsy_day, 1),
    pct_growth  = 100 * delta_g / initial_dry,
    sgr         = 100 * (log(final_dry) - log(initial_dry)) / days,
    g_per_day   = delta_g / days
  )
saveRDS(bw, file.path(DATA_PROC, "buoyant_weight_clean.rds"))

# ---- 6. Symbiont density ---------------------------------------------------
# Hemocytometer counts: average the four quadrants within each count, then the
# replicate counts within each fragment. Where raw counts are missing, the
# spreadsheet's average is used (count_source says which).
# cells per cm2 = mean count x 1e4 x slurry volume (mL) / surface area (cm2).
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

sym <- zoox |>
  left_join(meta |> filter(!is.na(chlorophyll_ug_cm2)) |> select(id, chlorophyll_ug_cm2),
            by = "id")
saveRDS(sym, file.path(DATA_PROC, "symbiont_chl_clean.rds"))

# ---- 7. Surface area (wax dipping) -----------------------------------------
# Standard curve: wax mass on calibration cylinders vs their known surface area.
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
cal <- lm(sa ~ wax_g, data = sc_fit)
write_csv(broom::tidy(cal, conf.int = TRUE), file.path(TBL_DIR, "07_wax_curve_fit.csv"))

wax <- suppressWarnings(
  read_csv(file.path(DATA_RAW, "wax_dipping", "data.csv"),
           show_col_types = FALSE, guess_max = 2000)
) |>
  janitor::clean_names() |>
  filter(!is.na(species) & !is.na(id)) |>
  rename(thicket = matches("^thicket")) |>
  mutate(
    id = as.integer(id),
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
micro_raw <- readr::read_csv(file.path(DATA_RAW, "microscope_physio", "data.csv"),
                             show_col_types = FALSE, guess_max = 1000, na = c("", "NA")) |>
  janitor::clean_names()
micro_traits <- c("tissue_over_wound", "hole_in_center", "polyp_in_hole", "wound_smoothed",
                  "pigment_over_wound", "tip_exist", "tip_extension", "new_corallites_on_tip")
missing_cols <- setdiff(c("thicket", "wound", "tank", "sample", "treatment", "id", "date",
                          "day", micro_traits, "use_photo", "notes"), names(micro_raw))
if (length(missing_cols)) stop("Microscope export missing columns: ",
                               paste(missing_cols, collapse = ", "), call. = FALSE)
yn_to_int <- function(x, column) {
  y <- stringr::str_to_lower(stringr::str_squish(as.character(x)))
  y[y %in% c("", "na", "nan")] <- NA_character_
  unexpected <- setdiff(stats::na.omit(unique(y)), c("yes", "no"))
  if (length(unexpected)) stop(sprintf("Unexpected yes/no code in %s: %s", column,
                                       paste(unexpected, collapse = ", ")), call. = FALSE)
  ifelse(is.na(y), NA_integer_, as.integer(y == "yes"))
}
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
if (!identical(micro$hole_in_center, micro$polyp_in_hole)) {
  stop("hole_in_center and polyp_in_hole differ in the microscope data.", call. = FALSE)
}
micro <- micro |> mutate(axial_polyp_formation = hole_in_center)
saveRDS(micro, file.path(DATA_PROC, "microscope_physio_clean.rds"))

# Check every microscope fragment against the master metadata.
crosswalk <- micro |>
  distinct(id, tank, treatment, thicket, wound, sample) |>
  mutate(treatment = as.character(treatment), wound = as.character(wound)) |>
  left_join(meta |> filter(sample == "photo") |>
              transmute(id, meta_tank = as.integer(tank), meta_treatment = as.character(treatment),
                        meta_thicket = as.character(thicket), meta_wound = as.character(wound),
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
worms <- read_csv(file.path(DATA_RAW, "worm_presence", "Sheet1.csv"), show_col_types = FALSE) |>
  janitor::clean_names() |>
  rename(thicket = matches("^thicket")) |>
  mutate(thicket = str_to_lower(str_squish(thicket)), tank = as.integer(tank),
         id = as.integer(id), treatment = tank_treatment(tank)) |>
  pivot_longer(starts_with("worms_"), names_to = "date_raw", values_to = "n_worms") |>
  mutate(date = mdy(str_remove(date_raw, "^worms_")), present = as.integer(n_worms > 0)) |>
  filter(!is.na(treatment))
saveRDS(worms, file.path(DATA_PROC, "worm_clean.rds"))
write_csv(worms |> group_by(date, treatment, tank) |>
            summarise(n_corals = n(), n_with_worms = sum(present, na.rm = TRUE),
                      pct = 100 * n_with_worms / n_corals, .groups = "drop"),
          file.path(TBL_DIR, "10_worm_summary.csv"))

cat("01_clean_data: wrote", length(list.files(DATA_PROC, "\\.rds$")), "clean .rds files\n")
