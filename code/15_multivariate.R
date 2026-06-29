# =============================================================================
# Purpose: Multivariate summary of the "heat stress syndrome" — collapse the
#          four key physiological responses (PAM, color, growth, symbionts)
#          to one biplot to visualize how thermal + wound treatments separate
#          in physiological state space.
#          Restricted to corals with measurements for all four responses.
#
# What & why: code/12 tested each response one at a time. But "heat stress" is a
#   syndrome — paling, falling Fv/Fm, slower growth and symbiont loss tend to move
#   together. PCA rotates the 4 correlated responses into uncorrelated axes ordered
#   by variance; PC1 typically becomes a single "overall heat-stress" axis. The
#   biplot then shows, in one picture, (a) whether 28 °C and 31 °C corals occupy
#   different regions of physiology space (treatment separation = a multivariate
#   heat effect) and (b) which original variables load on which axis (the arrows).
#   This is descriptive/visual, not a hypothesis test — it complements, not
#   replaces, the per-response models. We standardize before PCA (see below)
#   because the four variables are on wildly different units and scales.
# Input:   data/processed/{pam_clean,color_clean,buoyant_weight_clean,
#                          symbiont_chl_clean}.rds
# Output:  data/processed/coral_physio_wide.rds
#          figures/15_physio_PCA_biplot.{pdf,png}
#          output/tables/15_pca_loadings.csv
# =============================================================================

# 00_setup.R: packages, shared paths (DATA_PROC, TBL_DIR), palettes, theme_pub()
# and save_fig().
source(here::here("code", "00_setup.R"))

pam   <- readRDS(file.path(DATA_PROC, "pam_clean.rds"))
color <- readRDS(file.path(DATA_PROC, "color_clean.rds"))
bw    <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds"))
phys  <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
  filter(is.finite(cells_per_cm2))

# PCA needs one row per coral, so collapse each repeated-measures response to a
# single end-of-experiment value (the last day measured). Grouping by the design
# columns carries treatment/wound/genet through for colouring the biplot.
# Last-observation values per coral
pam_last <- pam |>
  group_by(id, treatment, wound, thicket) |>
  filter(day == max(day, na.rm = TRUE)) |>
  summarise(pam_end = mean(fv_fm), .groups = "drop")

color_last <- color |>
  group_by(id, treatment, wound, thicket) |>
  filter(day == max(day, na.rm = TRUE)) |>
  summarise(color_end = mean(color_num), .groups = "drop")

zoox_last <- phys |>
  group_by(id, treatment, wound, thicket) |>
  filter(biopsy_day == max(biopsy_day, na.rm = TRUE)) |>
  summarise(zoox_end = log10(mean(cells_per_cm2)), .groups = "drop")  # log10 to tame the skew

bw_summary <- bw |>
  group_by(id, treatment, wound, thicket) |>
  summarise(growth_pct = pct_growth[1], .groups = "drop")  # one growth rate per coral

# Join the four responses into one wide table (one row per coral). inner_join on
# PAM/color/growth keeps only corals with that always-present trio; zoox is
# left_join'd because symbiont density is missing for some corals.
# Wide-format physiology per coral
wide <- pam_last |>
  inner_join(color_last,  by = c("id", "treatment", "wound", "thicket")) |>
  inner_join(bw_summary,  by = c("id", "treatment", "wound", "thicket")) |>
  left_join(zoox_last,    by = c("id", "treatment", "wound", "thicket"))

saveRDS(wide, file.path(DATA_PROC, "coral_physio_wide.rds"))
cat("PCA-eligible corals:", nrow(wide),
    "(with all 4 responses:", sum(complete.cases(wide)), ")\n")

# PCA on the always-present trio (PAM, color, growth) plus symbiont density when
# available — so 4 responses here, since zoox_end is populated in this dataset.
pca_vars <- c("pam_end", "color_end", "growth_pct")
if (any(!is.na(wide$zoox_end))) pca_vars <- c(pca_vars, "zoox_end")

# PCA cannot handle missing values, so keep only complete-case corals. keep_rows
# is the matching logical index into `wide` so the design labels (groups) line up
# row-for-row with the PCA scores.
pca_input <- wide[, pca_vars, drop = FALSE] |> drop_na()
keep_rows <- complete.cases(wide[, pca_vars])
groups <- wide[keep_rows, c("treatment", "wound", "thicket")]

# center + scale. = TRUE → standardize each variable to mean 0, SD 1 before the
# rotation. Essential here: the variables are on different units (yield, D-score,
# mg cm⁻² d⁻¹, log10 cells), and without scaling the largest-variance variable
# would dominate PC1 purely because of its units. pca$rotation = the loadings
# (how each original variable maps onto each PC).
pca <- prcomp(pca_input, center = TRUE, scale. = TRUE)
loadings <- as.data.frame(pca$rotation) |>
  tibble::rownames_to_column("variable") |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))
write_csv(loadings, file.path(TBL_DIR, "15_pca_loadings.csv"))

# Proportion of total variance each PC captures (row 2 of the importance table),
# as a percentage — used for the axis labels below.
var_exp <- summary(pca)$importance[2, ] * 100

# pca$x = each coral's coordinates (scores) on the new axes; re-attach the design
# labels so we can colour/shape points by treatment, wound and genet.
scores <- as.data.frame(pca$x) |>
  bind_cols(groups)

# ---- Biplot (overall) -----------------------------------------------------
# A biplot overlays the loading vectors (arrows) on the score cloud. The ×2.8 is
# a purely cosmetic scaling so the unit-length loading arrows are visible against
# the spread of the scores — it does not change any inference.
loading_arrows <- as.data.frame(pca$rotation) |>
  tibble::rownames_to_column("variable") |>
  mutate(across(c(PC1, PC2), \(x) x * 2.8))

p_pca <- ggplot(scores, aes(PC1, PC2,
                            colour = treatment, shape = thicket)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.25) +
  # 68% (≈ 1 SD) data ellipse per temperature group — the visual test of whether
  # the two treatments occupy distinct regions of physiology space.
  stat_ellipse(aes(group = treatment, fill = treatment),
               geom = "polygon", level = 0.68,
               alpha = 0.1, colour = NA) +
  geom_point(size = 2.4, alpha = 0.9, stroke = 0.5) +
  geom_segment(data = loading_arrows,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(2, "mm")),
               colour = "grey30", linewidth = 0.35) +
  geom_text(data = loading_arrows,
            aes(PC1 * 1.15, PC2 * 1.15, label = variable),
            inherit.aes = FALSE,
            size = 3.2, colour = "grey15", fontface = "bold") +
  scale_colour_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                      name = "Temperature") +
  scale_fill_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                    guide = "none") +
  scale_shape_manual(values = c(a = 16, c = 17, d = 15), name = "Genet") +
  labs(x = sprintf("PC1 (%.0f%%)", var_exp[1]),
       y = sprintf("PC2 (%.0f%%)", var_exp[2]),
       title = "PCA of end-of-experiment physiology",
       subtitle = "Coloured ellipses: 68% confidence by treatment; shape = genet") +
  theme_pub(10)

save_fig(p_pca, "15_physio_PCA_biplot", width = 150, height = 130)

# ---- Per-genet faceted biplot — does each genet move the same way? --------
# Same biplot split into one panel per genet. If the 28C→31C shift is larger in
# one genet's panel than another's, that genet is more thermally plastic/sensitive
# — the multivariate echo of the genet effects in code/12 and code/13.
p_pca_facet <- ggplot(scores, aes(PC1, PC2,
                                   colour = treatment, shape = wound)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.25) +
  stat_ellipse(aes(group = treatment, fill = treatment),
               geom = "polygon", level = 0.68,
               alpha = 0.12, colour = NA) +
  geom_point(size = 2.2, alpha = 0.9, stroke = 0.4) +
  facet_wrap(~ thicket, labeller = labeller(thicket = function(x) paste("Genet", x))) +
  scale_colour_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                      name = "Temperature") +
  scale_fill_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                    guide = "none") +
  scale_shape_manual(values = c(no = 16, yes = 17), name = "Wound") +
  labs(x = sprintf("PC1 (%.0f%%)", var_exp[1]),
       y = sprintf("PC2 (%.0f%%)", var_exp[2]),
       title = "How each genet responds to heating in physiology space",
       subtitle = "Greater treatment separation = greater plasticity / thermal sensitivity") +
  theme_pub(10)

save_fig(p_pca_facet, "15b_physio_PCA_by_genet", width = 200, height = 100)

# ---- Per-genet centroid displacement under heat ---------------------------
# Compute centroid (PC1, PC2) per genet × treatment, then the vector length
# from 28C centroid to 31C centroid — this is the "thermal displacement" of
# each genet's physiology in PCA space.
# Centroid = mean PC1/PC2 position of each genet × treatment group.
centroids <- scores |>
  group_by(thicket, treatment) |>
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")

# Reshape to one row per genet with both treatment centroids side by side, then
# compute the Euclidean distance the genet's centroid moves from 28C to 31C.
# This `displacement` is a single-number index of how far heating pushes each
# genet through physiology space — bigger = more thermally responsive.
displ <- centroids |>
  pivot_wider(names_from = treatment, values_from = c(PC1, PC2),
              names_glue = "{.value}_{treatment}") |>
  mutate(
    dPC1 = `PC1_31C` - `PC1_28C`,
    dPC2 = `PC2_31C` - `PC2_28C`,
    displacement = sqrt(dPC1^2 + dPC2^2)   # vector length of the 28C→31C shift
  )
write_csv(displ, file.path(TBL_DIR, "15_genet_pca_displacement.csv"))

cat("\n=== PCA loadings ===\n")
print(loadings)
cat("\nPC variance explained: ",
    paste(sprintf("PC%d=%.0f%%", seq_along(var_exp), var_exp), collapse = ", "),
    "\n", sep = "")
cat("\n=== Per-genet thermal displacement in PCA space ===\n")
print(displ |> mutate(across(where(is.numeric), \(x) round(x, 2))))
