# =============================================================================
# Purpose: Project-wide setup — packages, paths, theme, palettes
# Input:   none
# Output:  attaches packages, defines theme_pub() and PAL_*, sets here() root
# Author:  Stier Lab
# =============================================================================

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

set.seed(42)

# ---- Paths -----------------------------------------------------------------
DATA_RAW  <- here("data", "raw")
DATA_PROC <- here("data", "processed")
DATA_META <- here("data", "metadata")
FIG_DIR   <- here("figures")
OUT_DIR   <- here("output")
TBL_DIR   <- here("output", "tables")
MOD_DIR   <- here("output", "models")

# ---- Palettes --------------------------------------------------------------
# Okabe-Ito (colorblind-friendly qualitative)
PAL_OKABE <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#000000")

# Treatment palette: ambient (28C) vs heated (31C) — the actual 2-level design.
# These exact hex codes are what every figure already draws, so this is the
# single source of truth: prefer `scale_*_manual(values = PAL_TEMP)` over
# re-typing the hexes in each script.
PAL_TEMP <- c(`28C` = "#56B4E9",   # blue — ambient (28 °C)
              `31C` = "#D55E00")   # red  — heated  (31 °C)

# Wound palette: unwounded vs wounded. Deliberately neutral grey/black — NOT
# blue/red — so wound is never visually confused with the temperature axis
# (28C blue / 31C red), which previously shared #D55E00 with "wounded".
PAL_WOUND <- c(no = "#9E9E9E", yes = "#000000")

# Genotype/thicket palette + point shapes: lowercase a, c, d (raw-data convention).
PAL_GENO <- c(a = "#E69F00", c = "#009E73", d = "#0072B2")
SHP_GENO <- c(a = 16, c = 17, d = 15)

# ---- Theme -----------------------------------------------------------------
# theme_pub(): single source of truth for figure text sizing.
# base_size is width-dependent — pass appropriate value per figure:
#   single-col (55-89 mm)  → base_size = 8
#   1.5-col   (115-170 mm) → base_size = 9-10
#   double-col(170-183 mm) → base_size = 10-11
theme_pub <- function(base_size = 10) {
  theme_bw(base_size = base_size, base_family = "") %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      panel.border     = element_rect(fill = NA, colour = "black", linewidth = 0.4),
      strip.background = element_blank(),
      strip.text       = element_text(size = base_size - 1, face = "bold"),
      axis.text        = element_text(size = base_size - 2, colour = "black"),
      axis.title       = element_text(size = base_size, face = "plain"),
      axis.ticks       = element_line(colour = "black", linewidth = 0.3),
      plot.title       = element_text(size = base_size, face = "bold", hjust = 0,
                                      margin = margin(b = 4)),
      plot.subtitle    = element_text(size = base_size - 1, hjust = 0,
                                      margin = margin(b = 6), colour = "grey30"),
      plot.tag         = element_text(size = base_size + 1, face = "bold"),
      plot.caption     = element_text(size = base_size - 2, hjust = 0, colour = "grey40"),
      legend.position  = "bottom",
      legend.title     = element_text(size = base_size - 1),
      legend.text      = element_text(size = base_size - 2),
      legend.key       = element_blank(),
      legend.background= element_blank(),
      plot.margin      = margin(8, 8, 8, 8, "pt")
    )
}

# ---- Save helper -----------------------------------------------------------
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

# ---- Session log -----------------------------------------------------------
log_session <- function(path = file.path(OUT_DIR, "session_info.txt")) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(sessionInfo()), path)
}

message("Setup complete. ", length(search()), " env entries; ",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
