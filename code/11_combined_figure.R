# =============================================================================
# Purpose: Multi-panel publication figure stitching together the four main
#          response variables (PAM, color, growth, symbionts) into one Figure 1
#          for the manuscript.
#
# What & why: this script doesn't analyse anything new — it assembles the study's
#   four headline coral responses into a single Figure 1 so a reader can see the
#   whole story at a glance. The panels are: (A) photosynthetic efficiency
#   (PAM Fv/Fm), (B) colony colour score (a visual bleaching proxy),
#   (C) calcification rate from buoyant weighing (growth), and (D) symbiont
#   cell density. Each panel reads its own cleaned dataset (produced by earlier
#   scripts), builds a small standalone plot, and patchwork tiles them into a
#   2×2 grid with one shared legend. The recurring design choice — colour by
#   WOUND in A-C, faceted by treatment, but colour by TREATMENT in D — reflects
#   what varies within each dataset. Edit the upstream scripts to change the
#   data; edit here only to change the figure layout.
# Input:   data/processed/{pam_clean,color_clean,buoyant_weight_clean,
#                          symbiont_chl_clean}.rds
# Output:  figures/11_main_response_panel.{pdf,png}
# =============================================================================

# 00_setup.R loads packages, shared paths (DATA_PROC, ...), theme_pub(),
# save_fig(), and the shared palettes (PAL_WOUND) used across all panels.
source(here::here("code", "00_setup.R"))

# Load the four cleaned datasets, each written by its own upstream script.
pam   <- readRDS(file.path(DATA_PROC, "pam_clean.rds"))
color <- readRDS(file.path(DATA_PROC, "color_clean.rds"))
bw    <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds"))
phys  <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds"))

# ---- A: PAM trajectory -----------------------------------------------------
# Photosynthetic efficiency (Fv/Fm) over time. Collapse to a mean ± standard
# error per day × treatment × wound; the ribbon is the ±1 SE band around the mean.
A_df <- pam |>
  group_by(day, treatment, wound) |>
  summarise(m = mean(fv_fm, na.rm = TRUE),
            se = sd(fv_fm, na.rm = TRUE) / sqrt(n()),   # SE = SD / sqrt(n)
            .groups = "drop")
pA <- ggplot(A_df, aes(day, m, colour = wound, fill = wound)) +
  geom_ribbon(aes(ymin = m - se, ymax = m + se), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +  # day 0 = treatment start
  facet_wrap(~ treatment, ncol = 2) +                  # one sub-panel per temperature
  scale_colour_manual(values = PAL_WOUND, name = "Wound") +
  scale_fill_manual(values   = PAL_WOUND, guide = "none") +  # reuse palette for ribbon, no 2nd legend
  labs(x = "Day", y = expression(F[v]/F[m]), tag = "A") +     # expression() renders the Fv/Fm subscripts
  theme_pub(9)

# ---- B: Color trajectory ---------------------------------------------------
# Colony colour score over time (same mean ± SE treatment as panel A). Colour is
# a 1-6 ordinal bleaching proxy where LOWER = paler/more bleached, so the y-axis
# is reversed below to keep "healthy/darker" at the top, matching intuition.
B_df <- color |>
  group_by(day, treatment, wound) |>
  summarise(m = mean(color_num, na.rm = TRUE),
            se = sd(color_num, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")
pB <- ggplot(B_df, aes(day, m, colour = wound, fill = wound)) +
  geom_ribbon(aes(ymin = m - se, ymax = m + se), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  facet_wrap(~ treatment, ncol = 2) +
  scale_y_reverse(breaks = 1:6) +                      # reverse so darker/healthier is up
  scale_colour_manual(values = PAL_WOUND, name = "Wound") +
  scale_fill_manual(values   = PAL_WOUND, guide = "none") +
  labs(x = "Day", y = "Color (D)", tag = "B") +
  theme_pub(9)

# ---- C: Calcification ------------------------------------------------------
# End-of-experiment calcification (growth) rate by treatment × wound. Unlike the
# time series above, this is a single value per colony, so a box-and-jitter plot
# shows the distribution; is.finite() drops colonies with no valid growth rate.
pC <- ggplot(filter(bw, is.finite(areal_calc)),
             aes(treatment, areal_calc, fill = wound)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7,   # hide outlier pts; raw pts added next
               position = position_dodge(0.7)) +
  geom_point(aes(colour = wound),
             # jitterdodge spreads points so wound groups don't overplot within a treatment
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
             size = 1.1, alpha = 0.7) +
  scale_fill_manual(values   = PAL_WOUND, name = "Wound") +
  scale_colour_manual(values = PAL_WOUND, guide = "none") +
  labs(x = NULL, y = expression(Calcification~(mg~cm^{-2}~d^{-1})), tag = "C") +  # units via plotmath
  theme_pub(9)

# ---- D: Symbionts ----------------------------------------------------------
# Symbiont (zooxanthellae) cell density at each biopsy day, by treatment. Divide
# by 1e6 so the axis reads in millions of cells (raw counts are unwieldy). Note
# the deliberate switch: this panel colours by TREATMENT (temperature is the
# contrast here), whereas A-C coloured by wound.
pD <- ggplot(phys, aes(factor(biopsy_day), cells_per_cm2 / 1e6,   # factor() = discrete biopsy days
                       fill = treatment)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, alpha = 0.4, size = 0.7) +
  scale_fill_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),  # blue=ambient, orange=heated
                    name = "Treatment") +
  labs(x = "Biopsy day",
       y = expression(Symbionts~(10^6~cells~cm^{-2})),
       tag = "D") +
  theme_pub(9)

# ---- Compose ---------------------------------------------------------------
# patchwork arithmetic lays the panels out: (A | B) on top, (C | D) below, equal
# row heights. The `&` applies the bottom-legend theme to every panel, and
# patchwork collects the duplicate legends into one shared legend at the bottom.
fig1 <- (pA + pB) / (pC + pD) +
  patchwork::plot_layout(heights = c(1, 1)) &
  theme(legend.position = "bottom")

save_fig(fig1, "11_main_response_panel", width = 200, height = 165)

cat("Wrote 11_main_response_panel.{pdf,png}\n")
