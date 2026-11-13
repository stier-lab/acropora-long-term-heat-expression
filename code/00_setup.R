# =============================================================================
# 00_setup.R: shared setup that every analysis script loads first
#
# Each of 01-08 starts with `source(here::here("code", "00_setup.R"))`. This
# file loads the shared packages, fixes the random seed so results repeat
# exactly, sets how the models code treatment groups, and defines the folder
# paths, colour palettes, tank-to-temperature map, figure style and helper
# functions, so every script and figure follows the same rules. Nothing here
# is specific to one measurement; it is shared setup for the LTH (heat ×
# wound × source patch) scripts.
#
# Sections
#   1. Packages
#   2. Random seed and group coding
#   3. Folder paths
#   4. Colour palettes
#   5. Tank to treatment map
#   6. Plot theme
#   7. Saving figures
#   8. Label corrections at import
#   9. Setup confirmation
#   10. Helpers for the results-text numbers
#
# Reads:  nothing
# Writes: nothing by itself; defines theme_pub(), save_fig(), PAL_*, paths and
#         the paper_row() / write_paper() helpers used by the other scripts
# Run:    sourced at the top of code/01_clean_data.R ... code/08_checks.R
# =============================================================================


# ---- 1. Packages -----------------------------------------------------------
# suppressPackageStartupMessages() hides the messages packages print when
# they load. Grouped by job: reading, tidying and plotting data (tidyverse,
# lubridate, janitor, readxl, patchwork, scales); turning model output into
# tables (broom, broom.mixed); mixed models, which allow for repeated
# measurements on the same fragment and tank (lme4, plus lmerTest, which
# adds p-values using the Satterthwaite method to estimate how much
# independent information the data hold); the model's average for each
# group and comparisons between groups (emmeans); and checks of model fit
# that compare the data with data simulated from the model (DHARMa).
# Package versions are pinned in renv.lock.
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(lubridate)
  library(janitor)
  library(readxl)
  library(patchwork)
  library(scales)
  library(broom)
  library(broom.mixed)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(DHARMa)
})


# ---- 2. Random seed and group coding ---------------------------------------
# Fix the random-number generator so anything that uses random draws (e.g.
# DHARMa's simulated data, the reshuffling tests in later scripts) gives the
# same numbers on every run.
set.seed(42)

# Set how R codes treatment groups in models. By default R compares every
# group with one baseline group. In a model with interactions, that makes the
# test of each main effect depend on which group is the baseline.
# "Sum-to-zero" coding (contr.sum) compares each group with the overall
# average instead, so the Type-III ANOVA tables in later scripts (each term
# tested after allowing for all the others) give the same answer whatever
# order the groups are listed in. Set once here for all scripts.
options(contrasts = c("contr.sum", "contr.poly"))


# ---- 3. Folder paths -------------------------------------------------------
# here() builds every path from the project folder (the one holding the
# .Rproj / .git), so scripts work wherever R was started. Naming the folders
# once here means later scripts never type out folder paths themselves.
DATA_RAW  <- here("data", "raw")
DATA_PROC <- here("data", "processed")
DATA_META <- here("data", "metadata")
FIG_DIR   <- here("figures")
OUT_DIR   <- here("output")
TBL_DIR   <- here("output", "tables")
MOD_DIR   <- here("output", "models")


# ---- 4. Colour palettes ----------------------------------------------------
# The experiment's colour code is set once here so every figure gives the
# same temperature, wound status or source patch the same colour. Later
# scripts call scale_*_manual(values = PAL_*) instead of re-typing colour
# codes, so figures cannot drift apart.

# Okabe-Ito: a set of colours that people with colour blindness can tell
# apart.
PAL_OKABE <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#000000")

# Temperature: ambient (28 °C) vs heated (31 °C), the 2-level design.
PAL_TEMP <- c(`28C` = "#56B4E9",   # blue: ambient (28 °C)
              `31C` = "#D55E00")   # red:  heated  (31 °C)

# Wound: unwounded vs wounded. Neutral grey/black, NOT blue/red, so wound is
# not confused with temperature (wounded previously shared #D55E00 with 31C).
PAL_WOUND <- c(no = "#9E9E9E", yes = "#000000")

# Source patch colours and point shapes (the field patch a fragment's parent
# branch came from; older column names say `thicket` or `genet`). Lowercase
# a, c, d matches the raw data.
PAL_GENO <- c(a = "#E69F00", c = "#009E73", d = "#0072B2")
SHP_GENO <- c(a = 16, c = 17, d = 15)


# ---- 5. Tank to treatment map ----------------------------------------------
# Temperature was piped to tanks in a fixed layout: tanks 3/6/9/12 ran
# ambient (28 °C) and 4/5/10/11 ran heated (31 °C). Any other tank was not
# part of the experiment and gets NA (dropped later). The tank, not the
# fragment, received the temperature, so every script takes temperature from
# this one map. It used to be typed out separately in the Apex (tank
# controller temperature log), YSI (handheld water-quality meter) and Fig. 1A
# scripts (old numbering 08, 09 and 16); defining it once keeps them in step.
TANK_28C <- c(3, 6, 9, 12)
TANK_31C <- c(4, 5, 10, 11)

# case_when() gives the label of the first rule that matches; the last rule
# (TRUE ~ NA) gives NA to every other tank.
tank_treatment <- function(tank) {
  dplyr::case_when(
    tank %in% TANK_28C ~ "28C",
    tank %in% TANK_31C ~ "31C",
    TRUE ~ NA_character_
  )
}


# ---- 6. Plot theme ---------------------------------------------------------
# theme_pub(): the one place that sets figure text sizes. All text sizes
# follow base_size (axis titles = base, tick labels = base - 2, panel strip
# labels = base - 1, panel letters = base + 1). The default is 10; pick
# base_size from the width of one PANEL (not the whole figure):
#   panel narrower than 50 mm -> base_size = 8
#   panel 50-80 mm            -> base_size = 9
#   panel 80 mm or wider      -> base_size = 10
# Artwork rule (2026-11-12): no embedded titles, subtitles, captions or
# legends, so those elements are blanked and legend.position = "none". Colour
# and shape keys, sample sizes and findings go in the separate caption text
# (figures/FIGURE_INDEX.md, Key section; manuscript legends). The legend.*
# text settings below only matter if a script turns a legend back on.
# `%+replace%` replaces each listed theme_bw() setting outright instead of
# blending the new values into the old ones.
theme_pub <- function(base_size = 10, margin_mm = 10) {
  theme_bw(base_size = base_size, base_family = "") %+replace%
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(linewidth = 0.25, colour = "grey88"),
      panel.border      = element_rect(fill = NA, colour = "black", linewidth = 0.4),
      strip.background  = element_blank(),
      strip.text        = element_text(size = base_size - 1, face = "bold"),
      axis.text         = element_text(size = base_size - 2, colour = "black"),
      axis.title        = element_text(size = base_size, face = "plain"),
      axis.ticks        = element_line(colour = "black", linewidth = 0.3),
      plot.title        = element_blank(),
      plot.subtitle     = element_blank(),
      plot.tag          = element_text(size = base_size + 1, face = "bold"),
      plot.caption      = element_blank(),
      legend.position   = "none",
      legend.title      = element_text(size = base_size - 1),
      legend.text       = element_text(size = base_size - 2),
      legend.key        = element_blank(),
      legend.background = element_blank(),
      # 10 mm outer margin on every side (lab figure rule) so nothing is
      # clipped at the page edge. Multi-panel figures use patch_margins().
      plot.margin       = margin(margin_mm, margin_mm, margin_mm, margin_mm, "mm")
    )
}

# patch_margins(): spacing for patchwork (multi-panel) figures. Each panel
# gets a small margin (3 mm, so facing panels sit ~6 mm apart) and the whole
# assembled figure gets the 10 mm outer margin. Without this, every panel
# would carry the full 10 mm and the gaps between panels would be 20 mm.
# Use as: save_fig(patch_margins(p1 + p2 + plot_annotation(tag_levels = "A")), ...)
patch_margins <- function(patch, panel_mm = 3, outer_mm = 10) {
  (patch & theme(plot.margin = margin(panel_mm, panel_mm, panel_mm, panel_mm, "mm"))) +
    plot_annotation(theme = theme(plot.margin = margin(outer_mm, outer_mm, outer_mm,
                                                       outer_mm, "mm")))
}


# ---- 7. Saving figures -----------------------------------------------------
# save_fig() writes each figure twice: a .pdf that stays sharp at any zoom
# (for the manuscript) and a .png image (for previews and slides), both at
# the same physical size and 300 dpi. Scripts call save_fig() once and get
# both files. cairo_pdf stores the fonts inside the PDF and keeps symbols
# such as ° and × intact.
save_fig <- function(plot, name, width = 170, height = 110, units = "mm",
                     dpi = 300, dir = FIG_DIR) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (ext in c("pdf", "png")) {
    ggsave(
      filename = file.path(dir, paste0(name, ".", ext)),
      plot     = plot,
      width    = width,
      height   = height,
      units    = units,
      dpi      = dpi,
      device   = if (ext == "pdf") cairo_pdf else "png"
    )
  }
  invisible(plot)
}


# ---- 8. Label corrections at import ----------------------------------------
# Fragment 119 is wounded but was recorded as unwounded in the raw export.
# Confirmed by Molly on 2026-10-19 and checked against Drive metadata on
# 2026-10-27. Keep the archived raw exports unchanged and correct the label
# when the data are imported (01_clean_data.R).
correct_wound_labels <- function(id, wound) {
  ifelse(!is.na(as.integer(id)) & as.integer(id) == 119L,
         "yes", as.character(wound))
}


# ---- 9. Setup confirmation -------------------------------------------------
# (The R and package versions are recorded by run_all.R in
# output/session_info.txt, and pinned in renv.lock.)
# Print a one-line, timestamped confirmation so the calling script's log shows
# that setup ran. search() lists the loaded packages (plus a few built-in R
# workspaces); the message prints how many there are.
message("Setup complete. ", length(search()), " env entries; ",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"))


# ---- 10. Helpers for the results-text numbers ------------------------------
# The "Numbers for the results text" section of scripts 02-06 writes one CSV
# per results section to output/paper/. Each row is one number quoted in the
# text. paper_row() builds one row (unused fields stay NA); write_paper()
# stacks the rows, rounds values to 6 significant digits so tiny rounding
# differences in the last decimal places between runs do not show up as
# changes in git, and writes <section>_numbers.csv.
PAPER_DIR <- here("output", "paper")

paper_row <- function(key, description, estimate = NA_real_, ci_low = NA_real_,
                      ci_high = NA_real_, test = NA_character_, statistic = NA_real_,
                      df1 = NA_real_, df2 = NA_real_, p_value = NA_real_, n = NA_real_,
                      source = NA_character_) {
  tibble::tibble(key = key, description = description,
                 estimate = estimate, ci_low = ci_low, ci_high = ci_high,
                 test = test, statistic = statistic, df1 = df1, df2 = df2,
                 p_value = p_value, n = n, source = source)
}

write_paper <- function(rows, section) {
  dir.create(PAPER_DIR, recursive = TRUE, showWarnings = FALSE)
  out <- dplyr::bind_rows(rows) |>
    dplyr::mutate(dplyr::across(c(estimate, ci_low, ci_high, statistic, df1, df2, p_value, n),
                                ~ signif(as.numeric(.x), 6)))

  # Every key must be unique, so each number in the text has one source row.
  stopifnot(!anyDuplicated(out$key))

  path <- file.path(PAPER_DIR, paste0(section, "_numbers.csv"))
  readr::write_csv(out, path, na = "")
  message("Wrote ", path, " (", nrow(out), " numbers)")
  invisible(out)
}
