# =============================================================================
# Purpose: Multi-panel publication figure stitching together the four main
#          response variables (PAM, color, growth, symbionts) into one Figure 1
#          for the manuscript.
# Input:   data/processed/{pam_clean,color_clean,buoyant_weight_clean,
#                          symbiont_chl_clean}.rds
# Output:  figures/11_main_response_panel.{pdf,png}
# =============================================================================

source(here::here("code", "00_setup.R"))

pam   <- readRDS(file.path(DATA_PROC, "pam_clean.rds"))
color <- readRDS(file.path(DATA_PROC, "color_clean.rds"))
bw    <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds"))
phys  <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds"))

# ---- A: PAM trajectory -----------------------------------------------------
A_df <- pam |>
  group_by(day, treatment, wound) |>
  summarise(m = mean(fv_fm, na.rm = TRUE),
            se = sd(fv_fm, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")
pA <- ggplot(A_df, aes(day, m, colour = wound, fill = wound)) +
  geom_ribbon(aes(ymin = m - se, ymax = m + se), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  facet_wrap(~ treatment, ncol = 2) +
  scale_colour_manual(values = PAL_WOUND, name = "Wound") +
  scale_fill_manual(values   = PAL_WOUND, guide = "none") +
  labs(x = "Day", y = expression(F[v]/F[m]), tag = "A") +
  theme_pub(9)

# ---- B: Color trajectory ---------------------------------------------------
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
  scale_y_reverse(breaks = 1:6) +
  scale_colour_manual(values = PAL_WOUND, name = "Wound") +
  scale_fill_manual(values   = PAL_WOUND, guide = "none") +
  labs(x = "Day", y = "Color (D)", tag = "B") +
  theme_pub(9)

# ---- C: Calcification ------------------------------------------------------
pC <- ggplot(filter(bw, is.finite(areal_calc)),
             aes(treatment, areal_calc, fill = wound)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7,
               position = position_dodge(0.7)) +
  geom_point(aes(colour = wound),
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
             size = 1.1, alpha = 0.7) +
  scale_fill_manual(values   = PAL_WOUND, name = "Wound") +
  scale_colour_manual(values = PAL_WOUND, guide = "none") +
  labs(x = NULL, y = expression(Calcification~(mg~cm^{-2}~d^{-1})), tag = "C") +
  theme_pub(9)

# ---- D: Symbionts ----------------------------------------------------------
pD <- ggplot(phys, aes(factor(biopsy_day), cells_per_cm2 / 1e6,
                       fill = treatment)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, alpha = 0.4, size = 0.7) +
  scale_fill_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                    name = "Treatment") +
  labs(x = "Biopsy day",
       y = expression(Symbionts~(10^6~cells~cm^{-2})),
       tag = "D") +
  theme_pub(9)

# ---- Compose ---------------------------------------------------------------
fig1 <- (pA + pB) / (pC + pD) +
  patchwork::plot_layout(heights = c(1, 1)) &
  theme(legend.position = "bottom")

save_fig(fig1, "11_main_response_panel", width = 200, height = 165)

cat("Wrote 11_main_response_panel.{pdf,png}\n")
