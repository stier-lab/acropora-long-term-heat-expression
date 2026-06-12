# =============================================================================
# Purpose: Multivariate summary of the "heat stress syndrome" — collapse the
#          four key physiological responses (PAM, color, growth, symbionts)
#          to one biplot to visualize how thermal + wound treatments separate
#          in physiological state space.
#          Restricted to corals with measurements for all four responses.
# Input:   data/processed/{pam_clean,color_clean,buoyant_weight_clean,
#                          symbiont_chl_clean}.rds
# Output:  data/processed/coral_physio_wide.rds
#          figures/15_physio_PCA_biplot.{pdf,png}
#          output/tables/15_pca_loadings.csv
# =============================================================================

source(here::here("code", "00_setup.R"))

pam   <- readRDS(file.path(DATA_PROC, "pam_clean.rds"))
color <- readRDS(file.path(DATA_PROC, "color_clean.rds"))
bw    <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds"))
phys  <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
  filter(is.finite(cells_per_cm2))

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
  summarise(zoox_end = log10(mean(cells_per_cm2)), .groups = "drop")

bw_summary <- bw |>
  group_by(id, treatment, wound, thicket) |>
  summarise(growth_areal = areal_calc[1], .groups = "drop")

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
pca_vars <- c("pam_end", "color_end", "growth_areal")
if (any(!is.na(wide$zoox_end))) pca_vars <- c(pca_vars, "zoox_end")

pca_input <- wide[, pca_vars, drop = FALSE] |> drop_na()
keep_rows <- complete.cases(wide[, pca_vars])
groups <- wide[keep_rows, c("treatment", "wound", "thicket")]

pca <- prcomp(pca_input, center = TRUE, scale. = TRUE)
loadings <- as.data.frame(pca$rotation) |>
  tibble::rownames_to_column("variable") |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))
write_csv(loadings, file.path(TBL_DIR, "15_pca_loadings.csv"))

# Variance explained
var_exp <- summary(pca)$importance[2, ] * 100

scores <- as.data.frame(pca$x) |>
  bind_cols(groups)

# ---- Biplot (overall) -----------------------------------------------------
loading_arrows <- as.data.frame(pca$rotation) |>
  tibble::rownames_to_column("variable") |>
  mutate(across(c(PC1, PC2), \(x) x * 2.8))

p_pca <- ggplot(scores, aes(PC1, PC2,
                            colour = treatment, shape = thicket)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.25) +
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
centroids <- scores |>
  group_by(thicket, treatment) |>
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")

displ <- centroids |>
  pivot_wider(names_from = treatment, values_from = c(PC1, PC2),
              names_glue = "{.value}_{treatment}") |>
  mutate(
    dPC1 = `PC1_31C` - `PC1_28C`,
    dPC2 = `PC2_31C` - `PC2_28C`,
    displacement = sqrt(dPC1^2 + dPC2^2)
  )
write_csv(displ, file.path(TBL_DIR, "15_genet_pca_displacement.csv"))

cat("\n=== PCA loadings ===\n")
print(loadings)
cat("\nPC variance explained: ",
    paste(sprintf("PC%d=%.0f%%", seq_along(var_exp), var_exp), collapse = ", "),
    "\n", sep = "")
cat("\n=== Per-genet thermal displacement in PCA space ===\n")
print(displ |> mutate(across(where(is.numeric), \(x) round(x, 2))))
