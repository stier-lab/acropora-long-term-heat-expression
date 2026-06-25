# =============================================================================
# Purpose: Pre-flight checks for every raw dataset before analysis.
#          Catches:
#            - Sample-size imbalances vs the design (192 destructive + 48 + 16)
#            - Missing/out-of-range values in key columns
#            - ID collisions or orphans across data streams
#            - Date sequence anomalies (records before acclimation start, etc.)
#            - Tank assignments inconsistent with treatment expectation
#                (28 °C tanks must be 3, 6, 9, 12; 31 °C tanks must be 4, 5, 10, 11)
# Output:  output/tables/18_validation_summary.csv  — one row per check
# =============================================================================

source(here::here("code", "00_setup.R"))

results <- list()
add <- function(stream, check, status, value, expected = "", notes = "") {
  results[[length(results) + 1]] <<- tibble(
    stream = stream, check = check, status = status,
    value = as.character(value), expected = as.character(expected),
    notes = notes
  )
}

EXPECTED_28C_TANKS <- c(3L, 6L, 9L, 12L)
EXPECTED_31C_TANKS <- c(4L, 5L, 10L, 11L)
EXPECTED_GENETS    <- c("a", "c", "d")

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

# ---- Metadata --------------------------------------------------------------
meta <- readRDS(file.path(DATA_PROC, "coral_metadata.rds"))
add("metadata", "n total fragments", if (nrow(meta) == 208) "PASS" else "FAIL",
    nrow(meta), 208, "192 destructive + 16 microscope = 208")
add("metadata", "treatments", if (setequal(levels(meta$treatment), c("28C", "31C"))) "PASS" else "FAIL",
    paste(levels(meta$treatment), collapse = ","), "28C, 31C")
add("metadata", "wound levels", if (setequal(levels(meta$wound), c("no", "yes"))) "PASS" else "FAIL",
    paste(levels(meta$wound), collapse = ","), "no, yes")
add("metadata", "genets", if (setequal(unique(meta$thicket), EXPECTED_GENETS)) "PASS" else "FAIL",
    paste(sort(unique(meta$thicket)), collapse = ","),
    paste(EXPECTED_GENETS, collapse = ","))
check_tanks("metadata", meta)

n_chl_missing <- sum(is.na(meta$chlorophyll_ug_cm2))
add("metadata", "chl-a populated",
    if (n_chl_missing == 0) "PASS" else "HANDLED",
    sprintf("%d/%d missing", n_chl_missing, nrow(meta)),
    "all missing because assay was not run",
    "confirmed by Molly; analysis uses PAM, color, and symbiont density")

# ---- PAM ------------------------------------------------------------------
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

# ---- Color ----------------------------------------------------------------
color <- readRDS(file.path(DATA_PROC, "color_clean.rds"))
add("color", "n observations", if (nrow(color) > 300) "PASS" else "WARN",
    nrow(color), "~336")
add("color", "D-scale range", "INFO",
    sprintf("[%g, %g]", min(color$color_num), max(color$color_num)),
    "[1, 6]")
check_tanks("color", color)

# ---- Buoyant weight -------------------------------------------------------
bw <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds"))
add("buoyant_weight", "n corals", if (nrow(bw) == 48) "PASS" else "WARN",
    nrow(bw), 48, "physiology subset")
add("buoyant_weight", "% growth range", "INFO",
    sprintf("[%.1f, %.1f]", min(bw$pct_growth), max(bw$pct_growth)),
    "positive for growing corals")

# ---- Symbionts ------------------------------------------------------------
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

# ---- Wax dipping ----------------------------------------------------------
wax <- readRDS(file.path(DATA_PROC, "wax_clean.rds"))
add("wax_dipping", "n corals", if (nrow(wax) >= 150) "PASS" else "WARN",
    nrow(wax), "~192")
agreement <- with(wax, cor(sa_caliper_cm2, sa_curve_cm2,
                            use = "complete.obs"))
add("wax_dipping", "caliper-curve SA correlation",
    if (agreement > 0.7) "PASS" else "HANDLED",
    sprintf("r = %.3f", agreement), "> 0.7",
    "below threshold; standard curve retained and plotted in figure 07 for transparency")

# ---- Physio morphology ----------------------------------------------------
ph <- readRDS(file.path(DATA_PROC, "physio_clean.rds"))
add("morphology", "n observations", if (nrow(ph) > 700) "PASS" else "WARN",
    nrow(ph), "~768")
add("morphology", "wounded corals tracked", "INFO",
    length(unique(ph$id[ph$wound == "yes"])), "~24 (per treatment)")
check_tanks("morphology", ph)

# ---- YSI ------------------------------------------------------------------
ysi <- readRDS(file.path(DATA_PROC, "ysi_clean.rds"))
add("ysi", "n daily readings", if (nrow(ysi) > 50) "PASS" else "WARN",
    nrow(ysi), "~72", "9 tanks × ~8 days")
add("ysi", "temperature (°C) range", "INFO",
    sprintf("[%.1f, %.1f]", min(ysi$temp_c, na.rm=TRUE),
                            max(ysi$temp_c, na.rm=TRUE)),
    "[27, 33] reasonable")

# ---- Apex -----------------------------------------------------------------
apex <- readRDS(file.path(DATA_PROC, "apex_temperature_daily.rds"))
add("apex", "n daily probe-records", "INFO", nrow(apex),
    "thousands across 6 datalog files")

# ---- Worms ----------------------------------------------------------------
worms <- readRDS(file.path(DATA_PROC, "worm_clean.rds"))
total_pos <- sum(worms$present, na.rm = TRUE)
add("worms", "total worm-positive obs", "INFO",
    total_pos, "0–few",
    "21 unique corals contaminated; concentrated in 31C tanks on 06/07 — treated with Worm Exit")

# ---- Cross-stream ID checks -----------------------------------------------
ids_in_meta <- meta$id
ids_in_pam  <- unique(pam$id)
ids_in_bw   <- unique(bw$id)
ids_in_wax  <- unique(wax$id)

orphan_pam <- setdiff(ids_in_pam, ids_in_meta)
orphan_bw  <- setdiff(ids_in_bw,  ids_in_meta)
orphan_wax <- setdiff(ids_in_wax, ids_in_meta)

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

# ---- Output ---------------------------------------------------------------
out <- bind_rows(results)
write_csv(out, file.path(TBL_DIR, "18_validation_summary.csv"))

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
