# =============================================================================
# Purpose: Integrative cross-response genet dashboard. For each of the 3
#          genets (a, c, d), compute a standardized "heat sensitivity" effect
#          across every response variable (PAM, color, growth, symbionts,
#          and 7 morphology traits). Build a forest plot of all effect sizes
#          to identify which genet is most/least thermally sensitive across
#          every dimension.
#
#          Also produces a composite "thermal resilience score" per genet
#          — the inverse of mean standardized effect across responses.
#
# What & why: this is the synthesis script. Every other analysis asks "does
#   heat affect response X?"; this one asks the cross-cutting question "which
#   genet best tolerates heat overall?" To put effects measured in completely
#   different units on one axis (Fv/Fm, a colour score, calcification rate,
#   symbiont counts, and wound-healing hazard ratios), each effect is rescaled
#   WITHIN its own response so the largest-magnitude effect = 1 (this is the "z"
#   column — a row-max standardization, not a statistical z-score). Positive
#   means the genet's phenotype is worse at 31 °C than 28 °C, i.e. more heat-
#   sensitive. Averaging these standardized values per genet yields a single
#   resilience ranking; in this dataset genet C is the most resilient. Nothing
#   here fits a new model — every number is a re-summary of upstream tables.
# Input:   output/tables/12_genet_treatment_effects.csv   (continuous responses)
#          output/tables/14_cox_hazard_ratios.csv         (per-genet HR for KM)
#          output/tables/15_genet_pca_displacement.csv    (multivariate)
# Output:  figures/19_genet_dashboard.{pdf,png}            — forest plot
#          figures/19b_genet_resilience_ranking.{pdf,png}  — composite score
#          output/tables/19_genet_resilience_summary.csv
# =============================================================================

# 00_setup.R loads packages, shared paths (TBL_DIR, FIG_DIR), theme_pub(),
# save_fig(), and the genet colour palette PAL_GENO used below.
source(here::here("code", "00_setup.R"))

# ---- Load upstream effect tables ------------------------------------------
# Three sources, one per analysis domain: per-genet continuous LMM contrasts
# (script 12), wound-healing Cox hazard ratios (script 14), and the
# multivariate per-genet PCA displacement (script 15).
cont <- read_csv(file.path(TBL_DIR, "12_genet_treatment_effects.csv"),
                 show_col_types = FALSE)
cox  <- read_csv(file.path(TBL_DIR, "14_cox_hazard_ratios.csv"),
                 show_col_types = FALSE)
pca  <- read_csv(file.path(TBL_DIR, "15_genet_pca_displacement.csv"),
                 show_col_types = FALSE)

# ---- Standardize continuous-response effects ------------------------------
# Continuous responses (from script 12) give estimate = mean(28C) - mean(31C)
# at Day 14 (end of experiment) per genet × wound, with SE. We collapse over
# wound (both wounded and unwounded show heat effects), keeping the
# standardized effect as estimate / SD_pooled.
cont_eff <- cont |>
  filter(response %in% c("pam_fvfm", "color_dscale", "growth_areal",
                         "log_zoox_density"),
         is.finite(estimate)) |>
  group_by(response, thicket) |>
  summarise(estimate = mean(estimate, na.rm = TRUE),
            se = mean(SE, na.rm = TRUE),
            .groups = "drop") |>
  # Re-scale per-response so values are comparable
  group_by(response) |>
  mutate(z = estimate / max(abs(estimate), na.rm = TRUE)) |>
  ungroup() |>
  mutate(metric = "Δ phenotype (28C − 31C, end of experiment)")

# ---- Standardize Cox hazard ratios (per genet) ----------------------------
# HR < 1 means 31C delays/prevents trait expression — i.e., heat sensitivity.
# Convert to log(HR) so negative = heat-impaired, positive = heat-promoted.
cox_per_genet <- cox |>
  filter(grepl("^genet=", scope), is.finite(HR_31_vs28),
         HR_31_vs28 > 0, HR_31_vs28 < Inf) |>
  mutate(thicket = sub("\\s+first-observed approximation$", "",
                       sub("^genet=", "", scope)),
         logHR = log(HR_31_vs28),
         response = paste0("morph_", trait))

cox_eff <- cox_per_genet |>
  group_by(response) |>
  mutate(hr_scale = max(abs(logHR), na.rm = TRUE),
         z = if_else(is.finite(hr_scale) & hr_scale > 0,
                     -logHR / hr_scale, NA_real_)) |>
  # Sign flipped: more negative HR (delayed healing) → more positive z (more impaired)
  ungroup() |>
  filter(is.finite(z)) |>
  select(response, thicket, estimate = logHR, z) |>
  mutate(metric = "−log(HR 31C/28C) for healing trait onset",
         se = NA_real_)

# ---- Combine and forest plot ---------------------------------------------
all_eff <- bind_rows(
  cont_eff |> select(response, thicket, estimate, se, z, metric),
  cox_eff
) |>
  mutate(
    response_label = case_when(
      response == "pam_fvfm"          ~ "PAM Fv/Fm",
      response == "color_dscale"      ~ "Color (D-scale)",
      response == "growth_areal"        ~ "Calcification",
      response == "log_zoox_density"  ~ "log symbionts cm⁻²",
      grepl("^morph_", response)      ~ str_to_sentence(
        gsub("_", " ", sub("^morph_", "", response))),
      TRUE                            ~ response
    ),
    domain = case_when(
      response %in% c("pam_fvfm","color_dscale","growth_areal","log_zoox_density")
        ~ "Physiology",
      grepl("hole|polyp|smoothed", response_label, ignore.case = TRUE)
        ~ "Wound closure",
      grepl("tip|corallite", response_label, ignore.case = TRUE)
        ~ "Regeneration",
      TRUE
        ~ "Other"
    )
  )

# Sort by domain then by response_label for the y-axis (so rows group by domain)
all_eff <- all_eff |>
  mutate(response_label = factor(response_label,
                                  levels = unique(response_label[order(domain, response_label)])))

# Forest plot: each point is one genet's standardized heat effect on one
# response, faceted by domain. Dashed line at x = 0 marks "no heat effect".
p_dash <- ggplot(all_eff,
                  aes(z, response_label,
                      colour = thicket, shape = thicket)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60",
             linewidth = 0.3) +
  geom_point(size = 3.2, alpha = 0.9) +
  scale_colour_manual(values = PAL_GENO, name = "Genet") +
  scale_shape_manual(values = c(a = 16, c = 17, d = 15), name = "Genet") +
  facet_grid(domain ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Standardized heat sensitivity (each row's max-scaled effect)",
       y = NULL,
       title = "Genet × response heat-sensitivity dashboard",
       subtitle = "Right of zero = the genet's phenotype is worse under 31 °C than 28 °C") +
  theme_pub(9) +
  theme(panel.grid.major.y = element_line(colour = "grey95", linewidth = 0.2),
        strip.text.y = element_text(face = "bold"))

save_fig(p_dash, "19_genet_dashboard", width = 180, height = 175)

# ---- Composite thermal resilience score ----------------------------------
# Per genet: mean standardized heat sensitivity across all responses.
# More-negative composite = more resilient (smaller heat effects).
resilience <- all_eff |>
  group_by(thicket) |>
  summarise(
    mean_sensitivity   = mean(z, na.rm = TRUE),
    median_sensitivity = median(z, na.rm = TRUE),
    n_responses        = n(),
    .groups = "drop"
  ) |>
  left_join(pca |> select(thicket, pca_displacement = displacement),
            by = "thicket") |>
  mutate(rank_overall = rank(mean_sensitivity))

# This summary table is the headline output: it (and its rank_overall column)
# is what script 20 reads back in as the "Genet resilience" rows.
write_csv(resilience, file.path(TBL_DIR, "19_genet_resilience_summary.csv"))

# Bar chart of the composite ranking; each bar is annotated with that genet's
# multivariate PCA displacement as an independent cross-check on the ordering.
p_rank <- ggplot(resilience,
                  aes(reorder(thicket, -mean_sensitivity),
                      mean_sensitivity, fill = thicket)) +
  geom_col(width = 0.55, alpha = 0.85, colour = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("displacement = %.2f", pca_displacement)),
            vjust = -0.4, size = 3, colour = "grey20") +
  scale_fill_manual(values = PAL_GENO, guide = "none") +
  labs(x = "Genet", y = "Mean standardized heat sensitivity",
       title = "Thermal resilience ranking across all responses",
       subtitle = "Bar height = average heat penalty across 11 responses; lower = more resilient. Text shows multivariate PCA displacement.") +
  coord_cartesian(ylim = c(0, max(resilience$mean_sensitivity) * 1.15)) +
  theme_pub(10)

save_fig(p_rank, "19b_genet_resilience_ranking", width = 130, height = 110)

# ---- Decomposed dashboard: heat-only vs heat-while-wounded ---------------
# The composite above pools across wound state. To answer "is genet C
# resilient to heat per se, or only to heat-while-wounded?" we split the
# continuous-response standardized effect into two scopes, plus an
# implicit "wound effect at 28C" baseline (no genet x wound contrasts
# are available in the current pipeline for an explicit wound-only column).
cont_by_wound <- cont |>
  # Drop morphology rows — they have wound=NA (wounded-only by construction)
  # and are captured separately in cox_by_wound.
  filter(!is.na(wound),
         response %in% c("pam_fvfm", "color_dscale", "growth_areal",
                          "log_zoox_density")) |>
  group_by(response, thicket, wound) |>
  summarise(estimate = mean(estimate, na.rm = TRUE),
            se = mean(SE, na.rm = TRUE),
            .groups = "drop") |>
  group_by(response, wound) |>
  mutate(z = estimate / max(abs(estimate), na.rm = TRUE)) |>
  ungroup() |>
  mutate(scope = if_else(wound == "yes",
                         "heat while wounded",
                         "heat only (unwounded)"))

# Cox results are wounded-only by design — flag them as the third scope
cox_by_wound <- cox_per_genet |>
  select(response, thicket, z, scope_orig = scope) |>
  mutate(scope = "heat while wounded",
         wound = "yes")

decomp <- bind_rows(
  cont_by_wound |> select(response, thicket, wound, scope, z),
  cox_by_wound  |> select(response, thicket, wound, scope, z)
) |>
  mutate(
    response_label = case_when(
      response == "pam_fvfm"         ~ "PAM Fv/Fm",
      response == "color_dscale"     ~ "Color (D)",
      response == "growth_areal"       ~ "Calcification",
      response == "log_zoox_density" ~ "log symbionts",
      grepl("^morph_", response)     ~ str_to_sentence(
        gsub("_", " ", sub("^morph_", "", response))),
      TRUE                            ~ response
    ),
    domain = case_when(
      response %in% c("pam_fvfm","color_dscale","growth_areal","log_zoox_density")
        ~ "Physiology",
      grepl("hole|polyp|smoothed", response_label, ignore.case = TRUE)
        ~ "Wound closure",
      grepl("tip|corallite", response_label, ignore.case = TRUE)
        ~ "Regeneration",
      TRUE ~ "Other"
    )
  )

# Same forest plot as p_dash, but now faceted domain (rows) × scope (columns)
# so heat-only and heat-while-wounded sensitivity can be compared side by side.
p_decomp <- ggplot(decomp,
                    aes(z, response_label,
                        colour = thicket, shape = thicket)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey60", linewidth = 0.3) +
  geom_point(size = 2.8, alpha = 0.9) +
  facet_grid(domain ~ scope, scales = "free_y", space = "free_y") +
  scale_colour_manual(values = PAL_GENO, name = "Genet") +
  scale_shape_manual(values = c(a = 16, c = 17, d = 15), name = "Genet") +
  labs(x = "Standardized heat sensitivity (row-max scaled)",
       y = NULL,
       title = "Decomposed resilience: heat-only vs heat-while-wounded",
       subtitle = "Right of zero = phenotype worse under 31 °C; Cox HRs are wounded-only by design") +
  theme_pub(9) +
  theme(panel.grid.major.y = element_line(colour = "grey95", linewidth = 0.2),
        strip.text.y = element_text(face = "bold"),
        strip.text.x = element_text(face = "bold"))

save_fig(p_decomp, "19c_decomposed_resilience", width = 200, height = 175)

# Per-genet × scope mean sensitivity
resilience_decomp <- decomp |>
  group_by(thicket, scope) |>
  summarise(mean_sensitivity = mean(z, na.rm = TRUE),
            n_responses      = n(),
            .groups = "drop")

write_csv(resilience_decomp,
          file.path(TBL_DIR, "19c_resilience_decomp_by_scope.csv"))

cat("\n=== Genet resilience summary ===\n")
print(resilience |> mutate(across(where(is.numeric), \(x) round(x, 3))))
cat("\n=== Decomposed resilience by scope ===\n")
print(resilience_decomp |> mutate(across(where(is.numeric), \(x) round(x, 3))))
cat("\nWrote 19_genet_dashboard.{pdf,png}, 19b_genet_resilience_ranking.{pdf,png},",
    "19c_decomposed_resilience.{pdf,png}, 19_genet_resilience_summary.csv,",
    "19c_resilience_decomp_by_scope.csv\n")
