# =============================================================================
# Purpose: Project-wide setup — packages, paths, theme, palettes
#
# What & why: every numbered analysis script starts with
#   `source(here::here("code", "00_setup.R"))`. It loads the shared packages,
#   fixes the random seed so results are reproducible, sets the project-wide
#   statistical contrast convention, and defines the folder paths, colour
#   palettes, and plotting theme so figures are consistent across the
#   manuscript. Nothing here is specific to one response variable; it is shared
#   infrastructure for the downstream LTH (heat × wound × genet) scripts.
# Input:   none
# Output:  attaches packages, defines theme_pub() and PAL_*, sets here() root
# Author:  Stier Lab
# =============================================================================

# ---- Packages --------------------------------------------------------------
# suppressPackageStartupMessages() hides the "attaching" banners so the
# console log stays readable. Grouped roughly by job: data wrangling/plotting
# (tidyverse, lubridate, janitor, readxl, patchwork, scales), tidy model output
# (broom, broom.mixed), mixed models (lme4 + lmerTest for p-values via
# Satterthwaite df), estimated marginal means / contrasts (emmeans), and
# simulation-based residual diagnostics for GLMMs (DHARMa).
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

# Fix the RNG so anything stochastic (e.g. DHARMa's simulated residuals,
# bootstrap/permutation steps in later scripts) returns the same numbers on
# every run — essential for a reproducible analysis.
set.seed(42)

# Sum-to-zero contrasts keep Type-III tests independent of reference levels.
# By default R uses treatment ("dummy") contrasts, where each main effect is
# tested against an arbitrary reference level — which makes the main-effect
# tests in a model WITH interactions depend on which level is
# the baseline. contr.sum codes factors so each level is compared to the grand
# mean instead, so the Type-III ANOVA tables in the downstream models are
# interpretable regardless of factor ordering. Set once here for all scripts.
options(contrasts = c("contr.sum", "contr.poly"))

# ---- Paths -----------------------------------------------------------------
# here() anchors every path to the project root (the folder with the .Rproj /
# .git), so scripts work no matter the working directory. Defining the folders
# once as constants means downstream scripts never hardcode relative paths.
DATA_RAW  <- here("data", "raw")
DATA_PROC <- here("data", "processed")
DATA_META <- here("data", "metadata")
FIG_DIR   <- here("figures")
OUT_DIR   <- here("output")
TBL_DIR   <- here("output", "tables")
MOD_DIR   <- here("output", "models")

# ---- Palettes --------------------------------------------------------------
# Define the experiment's colour code in one place so every figure maps the
# same treatment/wound/genet to the same colour. Centralising the hex codes
# here is what lets later scripts call scale_*_manual(values = PAL_*) instead of
# re-typing colours (and risking drift between figures).
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

# ---- Tank -> treatment map (single source of truth) ------------------------
# Temperature was plumbed to tanks by a fixed layout: tanks 3/6/9/12 ran ambient
# (28 °C) and 4/5/10/11 ran heated (31 °C); any other tank number is non-
# experimental and maps to NA (dropped downstream). This assignment was
# previously hardcoded identically in 08 (Apex), 09 (YSI), and 16 (Fig 1 panel
# A); defining it once here keeps those three in lock-step.
TANK_28C <- c(3, 6, 9, 12)
TANK_31C <- c(4, 5, 10, 11)
tank_treatment <- function(tank) {
  dplyr::case_when(
    tank %in% TANK_28C ~ "28C",
    tank %in% TANK_31C ~ "31C",
    TRUE ~ NA_character_
  )
}

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
# save_fig() writes each figure twice: a vector .pdf (what goes in the
# manuscript) and a raster .png (quick preview / slides), both at the same
# physical size and 300 dpi. Looping over the two extensions means downstream
# scripts call save_fig() once and get both files, always sized consistently.
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
# log_session() dumps sessionInfo() (R version + every package version) to a
# text file. Recording the exact software environment is part of making the
# analysis reproducible — a reader can recreate the package versions used.
log_session <- function(path = file.path(OUT_DIR, "session_info.txt")) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(sessionInfo()), path)
}

# Print a one-line confirmation (with a timestamp) so the sourcing script's log
# shows the setup ran. search() lists attached packages/environments.
message("Setup complete. ", length(search()), " env entries; ",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
