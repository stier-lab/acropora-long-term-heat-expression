# =============================================================================
# Purpose: PAM Fv/Fm × treatment × wound × day
#          - Mixed model: Fv/Fm ~ treatment * wound * day + (1|tank) + (1|thicket)
#          - Within-coral repeated measures via random effect of id
#          - Figure: trajectories with 95% CI, faceted by treatment, color by wound
# Input:   data/raw/pam/PAM_data.csv
#          data/processed/coral_metadata.rds
# Output:  data/processed/pam_clean.rds
#          figures/02_pam_fvfm_trajectory.{pdf,png}
#          output/tables/02_pam_treatment_contrasts.csv
#          output/models/02_pam_lmer.rds
# =============================================================================

source(here::here("code", "00_setup.R"))

# ---- Load ------------------------------------------------------------------
pam_raw <- read_csv(file.path(DATA_RAW, "pam", "PAM_data.csv"),
                    show_col_types = FALSE) |>
  janitor::clean_names()

meta <- readRDS(file.path(DATA_PROC, "coral_metadata.rds"))

# ---- Clean -----------------------------------------------------------------
# The Fv/Fm column was a spreadsheet formula in some rows ("=L2/100"); recompute
# from Y and the divisor convention used (Y / 1000 in PAM convention).
pam <- pam_raw |>
  rename(thicket = matches("^thicket"),
         fv_fm   = matches("^fv_fm")) |>
  mutate(
    date       = as_date(date),
    day        = as.integer(day),
    treatment  = factor(as.integer(treatment), levels = c(28, 31),
                        labels = c("28C", "31C")),
    tank       = as.integer(tank),
    wound      = factor(wound, levels = c("no", "yes")),
    thicket    = str_to_lower(str_squish(thicket)),
    id         = as.integer(id),
    location   = str_to_lower(location),
    f          = as.numeric(f),
    m          = as.numeric(m),
    y          = as.numeric(y),
    e          = as.numeric(e),
    # PAM convention: Fv/Fm = Y / 1000 when reported as raw "Y"
    fv_fm      = if_else(is.na(suppressWarnings(as.numeric(fv_fm))),
                         y / 1000,
                         suppressWarnings(as.numeric(fv_fm)))
  ) |>
  filter(!is.na(fv_fm), fv_fm > 0, fv_fm < 1)

# Average top/bottom replicate measurements per coral-day
pam_avg <- pam |>
  group_by(date, day, treatment, tank, thicket, wound, id) |>
  summarise(fv_fm = mean(fv_fm, na.rm = TRUE), .groups = "drop")

saveRDS(pam_avg, file.path(DATA_PROC, "pam_clean.rds"))

# ---- Location (top vs bottom) sensitivity check ----------------------------
# Molly's original analysis (code/archive/molly_original/LTH_PAM.R) noted that
# top vs bottom probe placement appeared to differ. The primary pipeline
# averages the two as technical replicates (above). Here we test that decision
# directly: fit the model on the UN-averaged data with location as a fixed
# factor and report the location terms. If location and its interactions are
# non-significant, averaging is justified.
if ("location" %in% names(pam) &&
    length(unique(na.omit(pam$location))) > 1) {
  pam_loc <- pam |> mutate(location = factor(location))
  # lmerTest::lmer (not lme4::lmer) so anova() returns Satterthwaite p-values.
  mod_loc <- lmerTest::lmer(
    fv_fm ~ treatment * wound * day * location + (1 | tank) + (1 | id),
    data = pam_loc, REML = FALSE,
    control = lme4::lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))
  )
  loc_anova <- as.data.frame(anova(mod_loc)) |>
    tibble::rownames_to_column("term") |>
    filter(grepl("location", term))
  write_csv(loc_anova, file.path(TBL_DIR, "02b_pam_location_sensitivity.csv"))
  cat("\n=== PAM location (top/bottom) sensitivity — terms involving location ===\n")
  print(loc_anova)
  # The decision that matters: does location INTERACT with the experimental
  # factors? A pure main-effect offset is averaged out harmlessly; an
  # interaction would mean averaging distorts the treatment/wound/day effects.
  pcol <- intersect(c("Pr(>F)"), names(loc_anova))
  interaction_terms <- loc_anova[grepl(":", loc_anova$term), , drop = FALSE]
  inter_sig <- length(pcol) == 1 &&
    any(interaction_terms[[pcol]] < 0.05, na.rm = TRUE)
  main_sig <- length(pcol) == 1 &&
    any(loc_anova[loc_anova$term == "location", pcol] < 0.05, na.rm = TRUE)
  cat(sprintf(
    "  location main effect %s; location x (experimental factor) interactions %s.\n",
    if (main_sig) "SIGNIFICANT (real top/bottom offset)" else "n.s.",
    if (inter_sig) "SIGNIFICANT — averaging may distort effects; see 02b table"
    else "n.s. — averaging top/bottom is justified for treatment/wound/day effects"))
}

# ---- Mixed model -----------------------------------------------------------
# Continuous day, factor treatment×wound. Random: tank + thicket + id (within).
mod <- lme4::lmer(
  fv_fm ~ treatment * wound * day + (1 | tank) + (1 | thicket) + (1 | id),
  data = pam_avg, REML = TRUE
)
saveRDS(mod, file.path(MOD_DIR, "02_pam_lmer.rds"))

# Treatment contrasts at each integer day
emm <- emmeans::emmeans(mod, ~ treatment * wound | day,
                        at = list(day = sort(unique(pam_avg$day))))
contrasts_tbl <- as_tibble(pairs(emm, adjust = "tukey"))
write_csv(contrasts_tbl, file.path(TBL_DIR, "02_pam_treatment_contrasts.csv"))

# ---- Figure ----------------------------------------------------------------
plot_df <- pam_avg |>
  group_by(day, treatment, wound) |>
  summarise(
    mean = mean(fv_fm, na.rm = TRUE),
    se   = sd(fv_fm, na.rm = TRUE) / sqrt(n()),
    n    = n(),
    .groups = "drop"
  )

p_pam <- ggplot(plot_df, aes(day, mean,
                             colour = wound, fill = wound, group = wound)) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se),
              alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  facet_wrap(~ treatment, ncol = 2) +
  scale_colour_manual(values = PAL_WOUND,
                      name = "Wound") +
  scale_fill_manual(values   = PAL_WOUND,
                    guide = "none") +
  labs(x = "Day relative to wounding (D0)",
       y = expression(F[v]/F[m]),
       title = "Photochemical efficiency",
       subtitle = "Mean ± 1 SE across genets and tanks") +
  theme_pub(10)

save_fig(p_pam, "02_pam_fvfm_trajectory", width = 170, height = 90)

# Console summary
cat("\n=== PAM mixed model ===\n")
print(summary(mod)$coefficients)
cat("\nWrote: pam_clean.rds, 02_pam_lmer.rds, 02_pam_treatment_contrasts.csv,",
    "02_pam_fvfm_trajectory.{pdf,png}\n")
