# =============================================================================
# Purpose: Figure render / sanity audit. Confirms that every figure in
#          figures/ actually rendered to a non-broken PNG, and flags files
#          that look wrong by size:
#            - too small (< 5 KB)   => probably a broken / empty render
#            - too large (> 600 KB) => bloated, consider downsizing for the PDF
#
# What & why: a manuscript has many figures, each produced by a different
#   script, so it is easy for one to silently fail to render (an empty or
#   truncated PNG) without anyone noticing until submission. This script is the
#   catch-all build check: it walks the figures/ folder, records each PNG's
#   size and last-modified time, and assigns a PASS / WARN / FAIL status from
#   the file size so a single glance confirms the whole figure set is present
#   and healthy. It is NOT a statistical analysis and it does not inspect the
#   pixels (no palette/legend/clipping checks here — those live in the
#   pub-figure-pipeline). The working diagnostic PNGs in figures/12_diagnostics/
#   are excluded: they are model-check plots, not manuscript figures, and have
#   their own size expectations.
# Input:   figures/**/*.png   (all PNGs under FIG_DIR, searched recursively)
# Output:  output/tables/17_figure_audit.csv
# =============================================================================

# 00_setup.R loads packages (tidyverse etc.) and the shared paths used below
# (FIG_DIR = figures/, TBL_DIR = output/tables/).
source(here::here("code", "00_setup.R"))

# ---- Find candidate figures ------------------------------------------------
# List every PNG under figures/ (recursive), then remove the working diagnostic
# plots in figures/12_diagnostics/ via setdiff() so only manuscript figures are
# held to the size rules below.
fig_files <- list.files(FIG_DIR, pattern = "\\.png$", full.names = TRUE,
                        recursive = TRUE) |>
  setdiff(list.files(file.path(FIG_DIR, "12_diagnostics"),
                     pattern = "\\.png$", full.names = TRUE))

# ---- Score each figure by file size ---------------------------------------
# One row per figure: file name, size in KB, last-modified time, and a status
# flag. The thresholds are heuristics — a near-empty PNG (< 5 KB) almost always
# means the render failed, and a very large PNG (> 600 KB) should be downsized
# before it goes into the vector PDF.
audit <- tibble(
  file        = basename(fig_files),
  size_kb     = round(file.size(fig_files) / 1024, 1),
  mtime       = file.mtime(fig_files),
  status      = case_when(
    file.size(fig_files) < 5  * 1024 ~ "FAIL:tiny — probably broken",
    file.size(fig_files) > 600 * 1024 ~ "WARN:large — consider downsizing for PDF",
    TRUE ~ "PASS"
  )
) |>
  arrange(file)

# Persist the audit table (machine-readable record of the figure inventory).
write_csv(audit, file.path(TBL_DIR, "17_figure_audit.csv"))

# ---- Print a human-readable summary to the console ------------------------
# Re-emit each row as a marked line (check / warn / cross) so whoever runs the
# script sees the verdict for every figure immediately, without opening the CSV.
cat("Figure audit (✓ PASS / ⚠ WARN / ✘ FAIL):\n\n")
audit |>
  mutate(marker = case_when(
    startsWith(status, "PASS") ~ "✓",
    startsWith(status, "WARN") ~ "⚠",
    TRUE                        ~ "✘"
  )) |>
  rowwise() |>
  mutate(line = sprintf("%s  %s (%5.1f KB) -- %s",
                        marker, file, size_kb, status)) |>
  pull(line) |>
  cat(sep = "\n")
cat("\n\nWrote output/tables/17_figure_audit.csv\n")
