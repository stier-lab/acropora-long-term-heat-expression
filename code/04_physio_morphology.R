# =============================================================================
# Purpose: 9 binary morphological characterizations of wound healing over time.
#          - For wounded corals only, plot cumulative proportion expressing each
#            trait by day × treatment, with genet (thicket) as a third
#            dimension (linetype) so genotype-specific responses are visible.
#          - Logistic mixed models per trait with treatment × thicket
#            interaction so genet-level differences in heat response are
#            tested directly (rather than absorbed in a random intercept).
#
# What & why: each wounded coral was photographed/scored over time for nine
#   yes/no signs of how the cut surface was healing. These split into two
#   biological phases: WOUND CLOSURE (the tissue sealing the lesion — e.g.
#   "wound smoothed", "pigment over wound", "hole in center" disappearing) and
#   REGENERATION (rebuilding lost skeletal structure — e.g. "tip exists", "tip
#   extension", "new corallites on tip"). This distinction is the heart of the
#   paper's headline: heat (31 °C) blocks REGENERATION but not basic wound
#   CLOSURE. Because the three genets (A/C/D) differ in heat resilience
#   (C > D > A), genet is modelled as a FIXED treatment × genet interaction —
#   we WANT to estimate each genotype's heat response, not average it away into
#   a random intercept. Each trait is a separate binary outcome, so each gets
#   its own logistic mixed model.
# Input:   data/raw/physio_morphology/data.csv
# Output:  data/processed/physio_clean.rds
#          figures/04_morphology_trajectories.{pdf,png}        — pooled-by-genet
#          figures/04b_morphology_trajectories_by_genet.{pdf,png} — split by genet
#          output/tables/04_morphology_trait_glmm_summaries.csv
# =============================================================================

# 00_setup.R loads packages, shared paths (DATA_RAW, DATA_PROC, TBL_DIR), the
# theme_pub() theme, and save_fig().
source(here::here("code", "00_setup.R"))

# ---- Load ------------------------------------------------------------------
# guess_max = 2000 makes readr scan many rows before guessing column types
# (some trait columns are blank early on, which would otherwise mis-type them).
# suppressWarnings hides the expected parse messages from that coercion.
# clean_names() -> snake_case; drop rows missing the essentials (species, id).
raw <- suppressWarnings(
  read_csv(file.path(DATA_RAW, "physio_morphology", "data.csv"),
           show_col_types = FALSE, guess_max = 2000)
) |>
  janitor::clean_names() |>
  filter(!is.na(species) & !is.na(id))

# ---- Clean -----------------------------------------------------------------
# Type-coerce the design columns and turn the nine yes/no trait columns into
# 0/1 integers so they can be modelled and averaged as proportions.
ph <- raw |>
  # Raw header "wounded" -> project-wide "wound"; match thicket by prefix.
  rename(thicket = matches("^thicket"),
         wound   = wounded) |>
  mutate(
    date      = as_date(date),
    day       = as.integer(day),                 # day relative to wounding (D0)
    # Treatment as set-point (28/31 °C); 28C is the reference level.
    treatment = factor(as.integer(treatment), levels = c(28, 31),
                       labels = c("28C", "31C")),
    tank      = as.integer(tank),
    thicket   = str_to_lower(str_squish(thicket)),  # genet ID; tidy case/whitespace
    id        = as.integer(id),                  # unique coral fragment ID
    wound     = factor(wound, levels = c("no", "yes")),  # "no" = reference level
    # Encode all nine trait columns at once: lower-case, trim, then test == "yes"
    # to get a clean 0/1 integer. Anything that isn't exactly "yes" -> 0 (NA
    # inputs stay NA via the comparison and are filtered later).
    across(c(polyps_out, hole_in_center, polyp_in_hole, wound_smoothed,
             pigment_over_wound, tip_exist, tip_extension,
             new_corallites_on_tip, algae_on_wound),
           ~ as.integer(str_to_lower(str_squish(.x)) == "yes"))
  )

# Save the cleaned, one-row-per-observation table for downstream scripts.
saveRDS(ph, file.path(DATA_PROC, "physio_clean.rds"))

# ---- Long-form for plotting ------------------------------------------------
# The nine trait columns, in the order they should appear in the facets below.
traits <- c("polyps_out", "hole_in_center", "polyp_in_hole",
            "wound_smoothed", "pigment_over_wound", "tip_exist",
            "tip_extension", "new_corallites_on_tip", "algae_on_wound")

# Restrict to WOUNDED corals (unwounded controls have no wound to heal), then
# stack the nine traits into one long trait/expressed column for faceting.
long <- ph |>
  filter(wound == "yes") |>
  select(day, treatment, tank, thicket, id, all_of(traits)) |>
  pivot_longer(all_of(traits), names_to = "trait", values_to = "expressed") |>
  filter(!is.na(expressed))                      # drop unscored coral-day-traits

# Proportion of wounded corals expressing each trait, by day × treatment. Because
# "expressed" is 0/1, its mean IS the proportion showing the trait.
prop_df <- long |>
  group_by(day, treatment, trait) |>
  summarise(prop = mean(expressed), n = n(), .groups = "drop") |>
  mutate(trait = factor(trait, levels = traits,
                        labels = c("Polyps out", "Hole in center",
                                   "Polyp in hole", "Wound smoothed",
                                   "Pigment over wound", "Tip exists",
                                   "Tip extension", "New corallites on tip",
                                   "Algae on wound")))

# ---- Pooled trajectory figure (one panel per trait) ------------------------
# Lines coloured by temperature; one panel per healing trait. Blue = 28 °C,
# orange/red = 31 °C. Dotted line at day 0 = wounding.
p_traits <- ggplot(prop_df, aes(day, prop,
                                colour = treatment, group = treatment)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.4) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  facet_wrap(~ trait, ncol = 4) +
  scale_colour_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                      name = "Temperature") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +          # y is a proportion, fixed 0-100%
  labs(x = "Day relative to wounding (D0)",
       y = "Wounded corals expressing trait",
       title = "Morphological wound-healing characteristics",
       subtitle = "Wounded corals only (n = 24 per treatment)") +
  theme_pub(9)

save_fig(p_traits, "04_morphology_trajectories", width = 200, height = 130)

# ---- Per-genet trajectory plot --------------------------------------------
# Same proportions, but now also split by genet (thicket) so genotype-specific
# heat responses (the C > D > A resilience ordering) are visible by eye.
prop_genet_df <- long |>
  group_by(day, treatment, thicket, trait) |>
  summarise(prop = mean(expressed), n = n(), .groups = "drop") |>
  mutate(trait = factor(trait, levels = traits,
                        labels = c("Polyps out", "Hole in center",
                                   "Polyp in hole", "Wound smoothed",
                                   "Pigment over wound", "Tip exists",
                                   "Tip extension", "New corallites on tip",
                                   "Algae on wound")))

# Temperature -> colour, genet -> linetype; group by their interaction so each
# temperature×genet combination is drawn as its own line.
p_traits_genet <- ggplot(prop_genet_df,
                          aes(day, prop, colour = treatment,
                              linetype = thicket,
                              group = interaction(treatment, thicket))) +
  geom_line(linewidth = 0.55, alpha = 0.85) +
  geom_point(size = 1.2, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  facet_wrap(~ trait, ncol = 4) +
  scale_colour_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                      name = "Temperature") +
  # Genets A/C/D get distinct dash patterns ("22"/"44" are on/off dash specs).
  scale_linetype_manual(values = c(a = "solid", c = "22", d = "44"),
                        name = "Genet") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(x = "Day relative to wounding (D0)",
       y = "Wounded corals expressing trait",
       title = "Morphological wound-healing trajectories by genet",
       subtitle = "Wounded corals only (n ≈ 8 per genet × treatment cell)") +
  theme_pub(9)

save_fig(p_traits_genet, "04b_morphology_trajectories_by_genet",
         width = 210, height = 135)

# ---- Per-trait GLMM with genet × treatment ---------------------------------
# Question: does the temperature effect on trait expression depend on genet?
# Each trait is a binary outcome, so we fit one logistic mixed model per trait.
# Fixed effects: treatment × day × thicket (genet is FIXED so we can read off
# each genotype's heat response, not absorb it into a random intercept).
# Random intercepts (1|tank) + (1|id) handle shared-tank and repeated-measures
# non-independence.
fit_one <- function(trait_name) {
  d <- long |>
    filter(trait == trait_name, day >= 0) |>    # post-wounding days only (healing window)
    mutate(thicket = factor(thicket))
  # Skip traits that are all-0 or all-1: a binary model needs both outcomes to
  # estimate anything (avoids complete separation / a degenerate fit).
  if (length(unique(d$expressed)) < 2) return(NULL)
  fit <- tryCatch(
    suppressMessages(suppressWarnings(
      # bobyqa optimizer + a high iteration cap (2e5) help these sparse binary
      # models converge; the whole fit is wrapped in tryCatch so one trait that
      # fails to converge returns NULL instead of stopping the loop.
      lme4::glmer(expressed ~ treatment * day * thicket + (1 | tank) + (1 | id),
                  family = binomial, data = d,
                  control = lme4::glmerControl(optimizer = "bobyqa",
                                               optCtrl = list(maxfun = 2e5)))
    )),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  # broom.mixed::tidy pulls the fixed-effect coefficient table into a tidy frame.
  tidy_fit <- broom.mixed::tidy(fit, effects = "fixed") |>
    mutate(trait = trait_name, n = nrow(d))
  # Type-II Wald ANOVA — reports significance of treatment, day, thicket and
  # their interactions. Type II respects marginality (tests each main effect
  # after the others but before its own higher-order interactions).
  av <- as.data.frame(car::Anova(fit, type = 2)) |>
    tibble::rownames_to_column("term") |>
    mutate(trait = trait_name)
  list(tidy = tidy_fit, anova = av)
}
# Fit every trait, name the list, and drop the NULLs (skipped/failed traits).
trait_results <- map(traits, fit_one) |> setNames(traits) |> compact()
# Stack the per-trait coefficient and ANOVA tables into two long data frames.
trait_tidy <- map_dfr(trait_results, "tidy")
trait_anova <- map_dfr(trait_results, "anova")
write_csv(trait_tidy,  file.path(TBL_DIR, "04_morphology_trait_glmm_summaries.csv"))
write_csv(trait_anova, file.path(TBL_DIR, "04_morphology_trait_anova_genet.csv"))

cat("Wrote physio_clean.rds, 04_morphology_trajectories.{pdf,png},",
    "04b_morphology_trajectories_by_genet.{pdf,png},",
    "04_morphology_trait_glmm_summaries.csv,",
    "04_morphology_trait_anova_genet.csv\n")
