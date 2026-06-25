# =============================================================================
# Purpose: Symbiont density (cells/cm^2); chlorophyll-a was planned but not run.
#          Endpoint physiological response variables sampled destructively at
#          each biopsy timepoint.
#          Joins symbiont counts with metadata. The chl-a column is retained only
#          as provenance for the planned-but-not-collected assay.
# Input:   data/raw/symbiont_counts/Raw_counts.csv  (4 hemocytometer reps per coral)
#          data/raw/symbiont_counts/metadata_ordered_merge.csv
#          data/processed/coral_metadata.rds
# Output:  data/processed/symbiont_chl_clean.rds
#          figures/06_symbiont_chl_by_day.{pdf,png}
#          output/tables/06_symbiont_chl_summary.csv
# =============================================================================

source(here::here("code", "00_setup.R"))

raw_counts <- suppressWarnings(
  read_csv(file.path(DATA_RAW, "symbiont_counts", "Raw_counts.csv"),
           show_col_types = FALSE, guess_max = 2000)
) |>
  janitor::clean_names()
zoox_meta  <- suppressWarnings(
  read_csv(file.path(DATA_RAW, "symbiont_counts",
                     "metadata_ordered_merge.csv"),
           show_col_types = FALSE, guess_max = 2000)
) |>
  janitor::clean_names()
meta <- readRDS(file.path(DATA_PROC, "coral_metadata.rds"))

# ---- Reps per coral -------------------------------------------------------
# Raw_counts has 4 hemocytometer quadrant counts (Q1, Q2, Q3, Q4) per
# replicate count per coral, in a "long" tidy layout.
# Columns: day, (blank), coral_id, count, q1, q2, q3, q4, average, coral_id_2
# Some rows have Q values, some don't (partial cells). Aggregate per coral_id:
reps <- raw_counts |>
  rename(coral_id = coral_id_3) |>
  filter(!is.na(coral_id)) |>
  mutate(coral_id = as.integer(coral_id)) |>
  mutate(across(c(q1, q2, q3, q4), as.numeric)) |>
  rowwise() |>
  mutate(quad_mean = mean(c(q1, q2, q3, q4), na.rm = TRUE)) |>
  ungroup() |>
  group_by(coral_id) |>
  summarise(
    mean_quad_count = mean(quad_mean, na.rm = TRUE),
    n_reps          = n(),
    .groups = "drop"
  ) |>
  filter(is.finite(mean_quad_count))

# Join with per-coral metadata (treatment, day, SA, etc.)
# metadata_ordered_merge has columns: species, thicket, id_3, sample, wound,
# tank, treatment, biopsy_day, biopsy_date, average, coral_id, id_12,
# calculated_sa_standard_curve, slurry_volume_m_l
zoox <- zoox_meta |>
  rename(thicket = matches("^thicket"),
         id      = id_3) |>
  mutate(
    id              = as.integer(id),
    treatment       = factor(as.integer(treatment), levels = c(28, 31),
                              labels = c("28C", "31C")),
    biopsy_day      = as.integer(biopsy_day),
    thicket         = str_to_lower(str_squish(thicket)),
    wound           = factor(wound, levels = c("no", "yes")),
    sa_cm2          = as.numeric(calculated_sa_standard_curve),
    slurry_ml       = as.numeric(slurry_volume_m_l),
    zoox_avg_sheet  = as.numeric(average),
    # Hemocytometer formula: density = (cells/quadrant) * 10^4 (cells/mL)
    # Total in slurry = density * slurry_mL ; then per cm^2 = / SA
    cells_per_cm2_sheet = (zoox_avg_sheet * 1e4 * slurry_ml) / sa_cm2
  ) |>
  filter(!is.na(id)) |>
  left_join(reps, by = c("id" = "coral_id")) |>
  mutate(
    zoox_avg_hemo = coalesce(mean_quad_count, zoox_avg_sheet),
    count_source = if_else(!is.na(mean_quad_count), "raw_counts", "metadata_average"),
    # Hemocytometer formula: density = (cells/quadrant) * 10^4 (cells/mL)
    # Total in slurry = density * slurry_mL ; then per cm^2 = / SA
    cells_per_cm2 = (zoox_avg_hemo * 1e4 * slurry_ml) / sa_cm2
  ) |>
  mutate(tank = as.integer(tank)) |>
  select(id, treatment, wound, biopsy_day, thicket, tank, sa_cm2,
         cells_per_cm2, count_source, n_reps)

# Pull chlorophyll from master metadata
chl <- meta |>
  filter(!is.na(chlorophyll_ug_cm2)) |>
  select(id, chlorophyll_ug_cm2)

phys <- zoox |>
  left_join(chl, by = "id")

saveRDS(phys, file.path(DATA_PROC, "symbiont_chl_clean.rds"))

# ---- Summary --------------------------------------------------------------
summary_tbl <- phys |>
  group_by(treatment, wound, biopsy_day) |>
  summarise(
    n             = n(),
    mean_cells    = mean(cells_per_cm2, na.rm = TRUE),
    se_cells      = sd(cells_per_cm2, na.rm = TRUE) / sqrt(n()),
    mean_chl      = mean(chlorophyll_ug_cm2, na.rm = TRUE),
    se_chl        = sd(chlorophyll_ug_cm2, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
write_csv(summary_tbl, file.path(TBL_DIR, "06_symbiont_chl_summary.csv"))

# ---- Figure ----------------------------------------------------------------
# Express counts in millions for readable axes
p_zoox <- ggplot(phys, aes(factor(biopsy_day), cells_per_cm2 / 1e6,
                            fill = treatment)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.65) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.4, size = 0.8) +
  scale_fill_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                    name = "Temperature") +
  labs(x = "Biopsy day",
       y = expression(Symbionts~(10^6~cells~cm^{-2})),
       title = "Symbiont density loss under heating",
       subtitle = expression(italic(A.~pulchra)~"biopsies, n = 192 corals across 5 timepoints")) +
  theme_pub(10)

has_chl <- any(is.finite(phys$chlorophyll_ug_cm2))
if (has_chl) {
  p_chl <- ggplot(phys, aes(factor(biopsy_day), chlorophyll_ug_cm2,
                             fill = treatment)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.65) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.4, size = 0.8) +
    scale_fill_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                      name = "Temperature") +
    labs(x = "Biopsy day", y = expression(Chl-a~(mu*g~cm^{-2})),
         tag = "B") +
    theme_pub(10)
  p_combo <- (p_zoox + labs(tag = "A")) + p_chl +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  save_fig(p_combo, "06_symbiont_chl_by_day", width = 170, height = 95)
} else {
  message("Chlorophyll-a assay was not run for this project — ",
          "saving single-panel symbiont figure.")
  save_fig(p_zoox, "06_symbiont_density_by_day", width = 140, height = 100)
  # Also write the two-panel layout with an explicit note for chl-a
  p_chl_placeholder <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "Chlorophyll-a assay\nwas not run\nfor this project",
             size = 4, colour = "grey30", lineheight = 1.1) +
    theme_void() +
    theme(panel.background = element_rect(fill = "grey96", colour = NA),
          plot.tag = element_text(size = 11, face = "bold")) +
    labs(tag = "B")
  p_combo <- (p_zoox + labs(tag = "A")) + p_chl_placeholder +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  save_fig(p_combo, "06_symbiont_chl_by_day", width = 170, height = 95)
}

cat("Wrote symbiont_chl_clean.rds, 06_symbiont_chl_summary.csv,",
    "06_symbiont_chl_by_day.{pdf,png}\n")
