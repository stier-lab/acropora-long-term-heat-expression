# =============================================================================
# 08_checks.R: data validation and model checks
#
# Checks that the cleaned data look the way the design says they should, and
# that the main models meet their assumptions well enough to trust. Nothing
# here changes a result; it writes tables that say which checks pass.
# Each validation row gets a status: PASS (as expected), FAIL (must be fixed
# before results are shared), WARN (worth a look), HANDLED (a known issue
# with a documented fix) or INFO (a value printed for the record).
#
# Sections
#   1. Data validation: sample sizes, tank-temperature assignments, value
#      ranges, IDs that match the metadata, and provenance files
#   2. Mixed-model checks for the four 03_physiology.R models: convergence,
#      singular fits, and simulated-residual tests (DHARMa; uniformity,
#      dispersion, outliers). Residual plots go to figures/diagnostics/ (not
#      tracked)
#   3. Proportional-hazards checks for the Cox models in 04_regeneration.R
#
# Reads:  data/processed/*.rds, data/metadata/source_inventory.csv and
#         codebooks, docs/provenance/drive_methods/, output/models/12_*.rds,
#         output/tables/14_event_intervals.csv
# Writes: output/tables/18_validation_summary.csv, 18_model_checks.csv,
#         14_cox_ph_tests.csv; figures/diagnostics/*.png
# Run:    sourced by code/run_all.R (or: Rscript code/run_all.R)
# =============================================================================

source(here::here("code", "00_setup.R"))

suppressPackageStartupMessages(library(survival))
DIAG_DIR <- file.path(FIG_DIR, "diagnostics")
dir.create(DIAG_DIR, recursive = TRUE, showWarnings = FALSE)


# ---- 1. Data validation ----------------------------------------------------
# Each check adds one row to a validation table. The checks guard against
# silent errors in 01_clean_data.R: a lost or duplicated fragment, a tank
# assigned to the wrong temperature, a value outside its possible range, or
# an ID that does not match the metadata (a typo that would drop a fragment
# from a join).

## 1a. Helpers ---------------------------------------------------------------
# add(): append one row to `results`. `<<-` writes to the `results` list
# defined here, outside the function, rather than to a local copy.
results <- list()
add <- function(stream, check, status, value, expected = "", notes = "") {
  results[[length(results) + 1]] <<- tibble(
    stream = stream, check = check, status = status,
    value = as.character(value), expected = as.character(expected),
    notes = notes
  )
}

# The design, typed out independently of 00_setup.R's TANK_28C / TANK_31C so
# a mistake there would show up here as a FAIL.
EXPECTED_28C_TANKS <- c(3L, 6L, 9L, 12L)
EXPECTED_31C_TANKS <- c(4L, 5L, 10L, 11L)
EXPECTED_GENETS    <- c("a", "c", "d")

# check_tanks(): do the tanks labelled 28C (and 31C) in a data set match the
# design exactly? A mismatch means temperature labels were joined wrongly.
# `if (length(obs_28))` skips a temperature that is absent from the data.
check_tanks <- function(stream, d) {
  obs_28 <- sort(unique(d$tank[d$treatment == "28C"]))
  obs_31 <- sort(unique(d$tank[d$treatment == "31C"]))
  if (length(obs_28)) {
    if (identical(obs_28, EXPECTED_28C_TANKS)) {
      add(stream, "28C tank IDs", "PASS",
          paste(obs_28, collapse = ","), paste(EXPECTED_28C_TANKS, collapse = ","))
    } else {
      add(stream, "28C tank IDs", "FAIL",
          paste(obs_28, collapse = ","),
          paste(EXPECTED_28C_TANKS, collapse = ","),
          "tank-treatment assignment mismatch")
    }
  }
  if (length(obs_31)) {
    if (identical(obs_31, EXPECTED_31C_TANKS)) {
      add(stream, "31C tank IDs", "PASS",
          paste(obs_31, collapse = ","), paste(EXPECTED_31C_TANKS, collapse = ","))
    } else {
      add(stream, "31C tank IDs", "FAIL",
          paste(obs_31, collapse = ","),
          paste(EXPECTED_31C_TANKS, collapse = ","),
          "tank-treatment assignment mismatch")
    }
  }
}

## 1b. Fragment metadata -----------------------------------------------------
# The master list of fragments: right total, right factor levels, right
# genets and tanks. Every other data stream is checked against it in 1f.
meta <- readRDS(file.path(DATA_PROC, "coral_metadata.rds"))
add("metadata", "n total fragments", if (nrow(meta) == 208) "PASS" else "FAIL",
    nrow(meta), 208, "192 destructive + 16 microscope = 208")
add("metadata", "treatments",
    if (setequal(levels(meta$treatment), c("28C", "31C"))) "PASS" else "FAIL",
    paste(levels(meta$treatment), collapse = ","), "28C, 31C")
add("metadata", "wound levels",
    if (setequal(levels(meta$wound), c("no", "yes"))) "PASS" else "FAIL",
    paste(levels(meta$wound), collapse = ","), "no, yes")
add("metadata", "genets",
    if (setequal(unique(meta$thicket), EXPECTED_GENETS)) "PASS" else "FAIL",
    paste(sort(unique(meta$thicket)), collapse = ","),
    paste(EXPECTED_GENETS, collapse = ","))
check_tanks("metadata", meta)

# Chlorophyll a was never assayed, so the column is empty by design.
n_chl_missing <- sum(is.na(meta$chlorophyll_ug_cm2))
add("metadata", "chl-a populated",
    if (n_chl_missing == 0) "PASS" else "HANDLED",
    sprintf("%d/%d missing", n_chl_missing, nrow(meta)),
    "all missing because assay was not run",
    "confirmed by Molly; analysis uses PAM, color, and symbiont density")

## 1c. Physiology streams: PAM, color, growth, symbionts, wax ----------------
# Row counts guard against lost or duplicated observations; ranges catch unit
# or data-entry errors (Fv/Fm must lie in 0-1, the color card runs 1-6).
pam <- readRDS(file.path(DATA_PROC, "pam_clean.rds"))
add("pam", "n observations (top+bottom averaged)",
    if (nrow(pam) > 300) "PASS" else "WARN",
    nrow(pam), "~336",
    "raw obs ≈ 672; cleaned data averages top/bottom probe locations per coral")
add("pam", "unique corals", if (length(unique(pam$id)) >= 48) "PASS" else "WARN",
    length(unique(pam$id)), 48)
add("pam", "Fv/Fm range", "INFO",
    sprintf("[%.2f, %.2f]", min(pam$fv_fm), max(pam$fv_fm)),
    "[0, 1]")
check_tanks("pam", pam)

color <- readRDS(file.path(DATA_PROC, "color_clean.rds"))
add("color", "n observations", if (nrow(color) > 300) "PASS" else "WARN",
    nrow(color), "~336")
add("color", "D-scale range", "INFO",
    sprintf("[%g, %g]", min(color$color_num), max(color$color_num)),
    "[1, 6]")
check_tanks("color", color)

bw <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds"))
add("buoyant_weight", "n corals", if (nrow(bw) == 48) "PASS" else "WARN",
    nrow(bw), 48, "physiology subset")
add("buoyant_weight", "% growth range", "INFO",
    sprintf("[%.1f, %.1f]", min(bw$pct_growth), max(bw$pct_growth)),
    "positive for growing corals")

# Symbiont densities: count finite values separately from the total, since
# NA or Inf densities are dropped from the models.
sym <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds"))
add("symbionts", "n biopsies", if (nrow(sym) >= 150) "PASS" else "WARN",
    nrow(sym), "~192")
sym_finite <- sum(is.finite(sym$cells_per_cm2))
add("symbionts", "finite densities",
    if (sym_finite >= 150) "PASS" else "WARN",
    sprintf("%d/%d finite", sym_finite, nrow(sym)), "192")
add("symbionts", "median (cells/cm²)", "INFO",
    sprintf("%.2e", median(sym$cells_per_cm2[is.finite(sym$cells_per_cm2)])),
    "~7e5–1e6")

# Wax dipping: surface area estimated two ways (calipers vs the wax standard
# curve) should agree; poor agreement is known and documented (HANDLED).
wax <- readRDS(file.path(DATA_PROC, "wax_clean.rds"))
add("wax_dipping", "n corals", if (nrow(wax) >= 150) "PASS" else "WARN",
    nrow(wax), "~192")
agreement <- with(wax, cor(sa_caliper_cm2, sa_curve_cm2,
                            use = "complete.obs"))
add("wax_dipping", "caliper-curve SA correlation",
    if (agreement > 0.7) "PASS" else "HANDLED",
    sprintf("r = %.3f", agreement), "> 0.7",
    "below threshold; standard curve retained and plotted in figure 07 for transparency")

## 1d. Daily morphology and the microscope experiment ------------------------
# Morphology: enough rows, every date parsed (a failed date parse would
# silently drop rows from any day-based analysis), and right tanks.
ph <- readRDS(file.path(DATA_PROC, "physio_clean.rds"))
add("morphology", "n observations", if (nrow(ph) > 700) "PASS" else "WARN",
    nrow(ph), "~768")
add("morphology", "date parsing",
    if (all(!is.na(ph$date))) "PASS" else "FAIL",
    if (all(!is.na(ph$date))) paste(range(ph$date), collapse = " to ")
    else paste0(sum(is.na(ph$date)), " missing"),
    "no missing dates", "all physio morphology rows retain calendar dates")
add("morphology", "wounded corals tracked", "INFO",
    length(unique(ph$id[ph$wound == "yes"])), "~24 (per treatment)")
check_tanks("morphology", ph)

# Microscope experiment: a separate cohort of 16 wounded "photo" fragments
# from genets a and c in two tanks per temperature, scored every day 0-15.
# These checks guard against it being mixed up with the main experiment.
micro <- readRDS(file.path(DATA_PROC, "microscope_physio_clean.rds"))
add("microscope", "n observations", if (nrow(micro) == 256) "PASS" else "WARN",
    nrow(micro), 256, "16 photo corals x 16 daily observations")
add("microscope", "n photo corals",
    if (length(unique(micro$id)) == 16) "PASS" else "FAIL",
    length(unique(micro$id)), 16)
add("microscope", "days scored",
    if (identical(sort(unique(micro$day)), 0:15)) "PASS" else "FAIL",
    paste(range(micro$day), collapse = "-"), "0-15")
add("microscope", "cohort",
    if (all(micro$sample == "photo") && all(micro$wound == "yes")) "PASS" else "FAIL",
    paste(sort(unique(micro$sample)), collapse = ","),
    "sample=photo; wound=yes")
add("microscope", "genets",
    if (setequal(unique(micro$thicket), c("a", "c"))) "PASS" else "FAIL",
    paste(sort(unique(micro$thicket)), collapse = ","), "a,c")
add("microscope", "28C tank IDs",
    if (identical(sort(unique(micro$tank[micro$treatment == "28C"])), c(6L, 9L))) "PASS"
    else "FAIL",
    paste(sort(unique(micro$tank[micro$treatment == "28C"])), collapse = ","),
    "6,9", "photo-only subset uses two of the four ambient tanks")
add("microscope", "31C tank IDs",
    if (identical(sort(unique(micro$tank[micro$treatment == "31C"])), c(4L, 5L))) "PASS"
    else "FAIL",
    paste(sort(unique(micro$tank[micro$treatment == "31C"])), collapse = ","),
    "4,5", "photo-only subset uses two of the four heated tanks")

# Pigment was not scored after Day 7 in the microscope set. Blanks there mean
# "not scored", not "no pigment", so they must stay NA (not 0).
add("microscope", "pigment D8-D15 blanks",
    if (all(is.na(micro$pigment_over_wound[micro$day > 7]))) "HANDLED" else "WARN",
    sprintf("%d/%d missing after D7",
            sum(is.na(micro$pigment_over_wound[micro$day > 7])),
            sum(micro$day > 7)),
    "blank = not scored, not no")

## 1e. Water quality, temperature logs and worms -----------------------------
# YSI handheld spot checks during the study window (June 4-19, 2025).
ysi <- readRDS(file.path(DATA_PROC, "ysi_clean.rds"))
ysi_study <- ysi |> filter(date >= as.Date("2025-06-04"), date <= as.Date("2025-06-19"))
add("ysi", "n study spot checks", if (nrow(ysi_study) == 112) "PASS" else "WARN",
    nrow(ysi_study), "112", "8 tanks × 14 dates; June 4 and 7 missing")
add("ysi", "temperature (°C) range", "INFO",
    sprintf("[%.1f, %.1f]", min(ysi_study$temp_c, na.rm = TRUE),
                            max(ysi_study$temp_c, na.rm = TRUE)),
    "[27, 33] reasonable")

# Apex controller probe logs, summarised to daily records.
apex <- readRDS(file.path(DATA_PROC, "apex_temperature_daily.rds"))
add("apex", "n daily probe-records", "INFO", nrow(apex),
    "thousands across 6 datalog files")

# Acropora-eating flatworm sightings.
worms <- readRDS(file.path(DATA_PROC, "worm_clean.rds"))
total_pos <- sum(worms$present, na.rm = TRUE)
add("worms", "total worm-positive obs", "INFO",
    total_pos, "0–few",
    "21 unique corals contaminated; concentrated in 31C tanks on 06/07 — treated with Worm Exit")

## 1f. Cross-stream ID matching ----------------------------------------------
# "Orphan" IDs appear in a data stream but not in the metadata. An orphan
# usually means a typo in a fragment ID, and that fragment would silently
# lose its treatment labels when the streams are joined.
ids_in_meta <- meta$id
ids_in_pam  <- unique(pam$id)
ids_in_bw   <- unique(bw$id)
ids_in_wax  <- unique(wax$id)
ids_in_micro <- unique(micro$id)
orphan_pam <- setdiff(ids_in_pam, ids_in_meta)
orphan_bw  <- setdiff(ids_in_bw,  ids_in_meta)
orphan_wax <- setdiff(ids_in_wax, ids_in_meta)
orphan_micro <- setdiff(ids_in_micro, ids_in_meta)
add("cross-stream", "PAM ids in metadata",
    if (length(orphan_pam) == 0) "PASS" else "FAIL",
    length(orphan_pam), 0,
    paste("Orphans:", paste(head(orphan_pam, 5), collapse = ",")))
add("cross-stream", "buoyant-weight ids in metadata",
    if (length(orphan_bw) == 0) "PASS" else "FAIL",
    length(orphan_bw), 0)
add("cross-stream", "wax-dipping ids in metadata",
    if (length(orphan_wax) == 0) "PASS" else "FAIL",
    length(orphan_wax), 0)
add("cross-stream", "microscope ids in metadata",
    if (length(orphan_micro) == 0) "PASS" else "FAIL",
    length(orphan_micro), 0)

## 1g. Provenance files ------------------------------------------------------
# Can a reader trace every input back to its source? Checks that the source
# inventory exists, that every raw-data folder has a codebook (column
# definitions), that the Drive method documents were copied in, and that the
# paths and scripts the inventory lists actually exist.
source_inventory_path <- file.path(DATA_META, "source_inventory.csv")
add("provenance", "source inventory exists",
    if (file.exists(source_inventory_path)) "PASS" else "FAIL",
    basename(source_inventory_path), "data/metadata/source_inventory.csv")

# Each folder in data/raw/ should have data/metadata/<folder>_codebook.csv.
raw_dirs <- list.dirs(DATA_RAW, full.names = FALSE, recursive = FALSE)
expected_codebooks <- file.path(DATA_META, paste0(raw_dirs, "_codebook.csv"))
missing_codebooks <- raw_dirs[!file.exists(expected_codebooks)]
add("provenance", "raw folders have codebooks",
    if (length(missing_codebooks) == 0) "PASS" else "FAIL",
    if (length(missing_codebooks) == 0) length(raw_dirs)
    else paste(missing_codebooks, collapse = ","),
    paste0(length(raw_dirs), " raw folders documented"))

method_dir <- here::here("docs", "provenance", "drive_methods")
method_files <- list.files(method_dir, all.files = FALSE, recursive = FALSE,
                           full.names = FALSE)
add("provenance", "Drive method docs imported",
    if (length(method_files) >= 8) "PASS" else "WARN",
    length(method_files), ">= 8")

if (file.exists(source_inventory_path)) {
  source_inventory <- read_csv(source_inventory_path, show_col_types = FALSE)

  # path_exists(): local_path can hold several paths separated by ";". TRUE
  # only if there is at least one path and every one exists on disk.
  path_exists <- function(path_string) {
    paths <- str_split(path_string, ";", simplify = FALSE)[[1]] |>
      str_squish()
    paths <- paths[nzchar(paths)]
    length(paths) > 0 && all(file.exists(paths))
  }

  # Sources marked "imported" must really be on disk.
  imported <- source_inventory |>
    filter(import_status == "imported")
  missing_imports <- imported$source_title[!map_lgl(imported$local_path, path_exists)]
  add("provenance", "imported source paths exist",
      if (length(missing_imports) == 0) "PASS" else "FAIL",
      if (length(missing_imports) == 0) nrow(imported)
      else paste(missing_imports, collapse = ","),
      paste0(nrow(imported), " imported sources"))

  # Enough analysis inputs are inventoried to cover every data stream.
  analysis_inputs <- source_inventory |>
    filter(pipeline_role == "analysis_input")
  add("provenance", "analysis inputs inventoried",
      if (nrow(analysis_inputs) >= 12) "PASS" else "WARN",
      nrow(analysis_inputs), ">= 12",
      "metadata, physiology, morphology, temperature, worms, and RNA-seq handoff inputs")

  # Every script the inventory names as reading a source must exist.
  scripted_sources <- source_inventory |>
    filter(!is.na(analysis_script), nzchar(analysis_script))
  missing_scripts <-
    scripted_sources$analysis_script[!file.exists(scripted_sources$analysis_script)]
  add("provenance", "inventoried analysis scripts exist",
      if (length(missing_scripts) == 0) "PASS" else "FAIL",
      if (length(missing_scripts) == 0) nrow(scripted_sources)
      else paste(unique(missing_scripts), collapse = ","),
      paste0(nrow(scripted_sources), " scripted sources"))
}

## 1h. Write and print the validation table ----------------------------------
out <- bind_rows(results)
write_csv(out, file.path(TBL_DIR, "18_validation_summary.csv"))

# Print one aligned line per check with a status symbol, so problems stand
# out in the run log. rowwise() lets sprintf() and nzchar() work row by row.
cat("Data validation summary:\n\n")
out |>
  mutate(marker = case_when(
    status == "PASS" ~ "✓",
    status == "FAIL" ~ "✘",
    status == "WARN" ~ "⚠",
    status == "HANDLED" ~ "↻",
    TRUE              ~ "ℹ"
  )) |>
  rowwise() |>
  mutate(line = sprintf("%s [%-15s] %-32s : %-30s | expected: %s%s",
                        marker, stream, check, value, expected,
                        if (nzchar(notes)) paste0("  — ", notes) else "")) |>
  pull(line) |> cat(sep = "\n")

n_fail <- sum(out$status == "FAIL")
n_handled <- sum(out$status == "HANDLED")
n_warn <- sum(out$status == "WARN")
n_pass <- sum(out$status == "PASS")
cat(sprintf("\n\nTotal: %d PASS, %d HANDLED, %d WARN, %d FAIL\n",
            n_pass, n_handled, n_warn, n_fail))
if (n_fail > 0) cat("⚠ Address FAILs before publishing results.\n")


# ---- 2. Mixed-model checks -------------------------------------------------
# For the four 03_physiology.R models (all lmer fits, read back from
# output/models/):
#   - convergence messages: the optimizer may not have found the best fit;
#   - singular fit: a random-effect variance estimated at (or near) zero, so
#     that random effect is not doing anything;
#   - DHARMa residuals: simulate 500 new data sets from the fitted model and
#     ask where each real observation falls among them. If the model is
#     right, these "scaled residuals" are uniform between 0 and 1.
#     Uniformity tests the overall shape, dispersion tests whether the data
#     spread more or less than the model expects, and outliers counts points
#     outside every simulated value.
# The seed is reset so the simulations, and hence the p-values, repeat.
model_files <- c(pam_fvfm = "12_pam_lmm.rds", color_dscale = "12_color_lmm.rds",
                 log_zoox_density = "12_zoox_lmm.rds", growth_pct = "12_bw_lmm.rds")
set.seed(42)

# imap_dfr(): loop over model_files with both the value (file) and the name,
# and stack the one-row results into one table.
model_checks <- purrr::imap_dfr(model_files, function(file, name) {
  m <- readRDS(file.path(MOD_DIR, file))
  conv <- m@optinfo$conv$lme4$messages   # `@` reads a slot of an S4 object

  # Save the standard DHARMa residual plots for a visual check.
  sim <- DHARMa::simulateResiduals(m, n = 500)
  png(file.path(DIAG_DIR, paste0("08_dharma_", name, ".png")),
      width = 1600, height = 800, res = 160)
  plot(sim)
  dev.off()

  tibble(
    model = name,
    n_obs = nobs(m),
    singular_fit = lme4::isSingular(m),
    convergence_messages = if (length(conv)) paste(conv, collapse = "; ") else "none",
    residual_uniformity_p = DHARMa::testUniformity(sim, plot = FALSE)$p.value,
    residual_dispersion_p = DHARMa::testDispersion(sim, plot = FALSE)$p.value,
    residual_outliers_p   = DHARMa::testOutliers(sim, plot = FALSE)$p.value
  )
}) |>
  mutate(across(where(is.double), \(x) signif(x, 4)),
         note = "p < 0.05 flags a departure worth a look; with these sample sizes small departures are common")
write_csv(model_checks, file.path(TBL_DIR, "18_model_checks.csv"))
print(as.data.frame(model_checks))


# ---- 3. Proportional-hazards checks (Cox models) ---------------------------
# A Cox model summarises the temperature effect on the timing of each stage
# as one hazard ratio, which assumes the ratio stays constant over the 15
# days ("proportional hazards"). cox.zph() tests this with Schoenfeld
# residuals: if the temperature effect drifts with time, the residuals trend
# and p < 0.05. Event intervals come from 04_regeneration.R.
events <- read_csv(file.path(TBL_DIR, "14_event_intervals.csv"), show_col_types = FALSE) |>
  mutate(treatment = factor(treatment, levels = c("28C", "31C")), thicket = factor(thicket))

# Override the project-wide sum contrasts for this factor so the model has
# a single "treatment" term (31C vs 28C) that zph$table can be indexed by.
contrasts(events$treatment) <- contr.treatment(nlevels(events$treatment))
traits <- c("axial_polyp_formation", "wound_smoothed", "pigment_over_wound",
            "tip_exist", "tip_extension", "new_corallites_on_tip")
ph_rows <- list()

# One Cox model per stage, stratified by genet (strata(): each genet gets its
# own baseline hazard). Stages with fewer than 5 events, or fits that fail,
# are skipped; tryCatch() turns an error into NULL so the loop continues.
for (tr in traits) {
  d <- events |> filter(trait == tr)
  if (sum(d$event) < 5) next
  fit <- tryCatch(coxph(Surv(event_day, event) ~ treatment + strata(thicket),
                        data = d), error = function(e) NULL)
  if (is.null(fit)) next
  zph <- tryCatch(cox.zph(fit), error = function(e) NULL)  # the Schoenfeld PH test
  if (is.null(zph)) next
  ph_rows[[tr]] <- tibble(
    trait      = tr,
    n_event    = fit$nevent,
    zph_chisq  = zph$table["treatment", "chisq"],
    zph_p      = zph$table["treatment", "p"],
    ph_ok      = zph$table["treatment", "p"] >= 0.05
  )

  # Residual-vs-time plot. plot(zph) needs more than 4 distinct event times;
  # otherwise (or if it errors) fall back to a plain scatter with a lowess
  # trend line. A flat cloud around 0 supports proportional hazards.
  ttl <- sprintf("PH check: %s (zph p=%.3f)", tr, zph$table["treatment", "p"])
  tt  <- as.numeric(zph$time); yy <- as.numeric(zph$y)
  png(file.path(DIAG_DIR, paste0("14_cox_ph_", tr, ".png")),
      width = 800, height = 500, res = 130)
  if (length(unique(tt)) > 4) {
    ok <- tryCatch({ plot(zph, resid = TRUE, main = ttl); TRUE },
                   error = function(e) FALSE)
  } else ok <- FALSE
  if (!ok) {
    plot(tt, yy, xlab = "Time (days)", ylab = "Scaled Schoenfeld residual",
         main = ttl, pch = 19, col = "grey40",
         ylim = range(c(0, yy), na.rm = TRUE))
    abline(h = 0, lty = 2, col = "grey60")
    if (length(unique(tt)) > 2)
      try(lines(lowess(tt, yy), col = "#D55E00", lwd = 2), silent = TRUE)
  }
  dev.off()
}

cox_ph <- bind_rows(ph_rows)
write_csv(cox_ph, file.path(TBL_DIR, "14_cox_ph_tests.csv"))
cat("\n=== Proportional-hazards (cox.zph) tests, all overall Cox models ===\n")
print(as.data.frame(cox_ph |> mutate(across(where(is.numeric), \(x) round(x, 4)))))
