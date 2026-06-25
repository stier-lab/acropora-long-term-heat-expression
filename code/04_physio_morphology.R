# =============================================================================
# Purpose: 9 binary morphological characterizations of wound healing over time.
#          - For wounded corals only, plot cumulative proportion expressing each
#            trait by day × treatment, with genet (thicket) as a third
#            dimension (linetype) so genotype-specific responses are visible.
#          - Logistic mixed models per trait with treatment × thicket
#            interaction so genet-level differences in heat response are
#            tested directly (rather than absorbed in a random intercept).
# Input:   data/raw/physio_morphology/data.csv
# Output:  data/processed/physio_clean.rds
#          figures/04_morphology_trajectories.{pdf,png}        — pooled-by-genet
#          figures/04b_morphology_trajectories_by_genet.{pdf,png} — split by genet
#          output/tables/04_morphology_trait_glmm_summaries.csv
# =============================================================================

source(here::here("code", "00_setup.R"))

raw <- suppressWarnings(
  read_csv(file.path(DATA_RAW, "physio_morphology", "data.csv"),
           show_col_types = FALSE, guess_max = 2000)
) |>
  janitor::clean_names() |>
  filter(!is.na(species) & !is.na(id))

ph <- raw |>
  rename(thicket = matches("^thicket"),
         wound   = wounded) |>
  mutate(
    date      = as_date(date),
    day       = as.integer(day),
    treatment = factor(as.integer(treatment), levels = c(28, 31),
                       labels = c("28C", "31C")),
    tank      = as.integer(tank),
    thicket   = str_to_lower(str_squish(thicket)),
    id        = as.integer(id),
    wound     = factor(wound, levels = c("no", "yes")),
    across(c(polyps_out, hole_in_center, polyp_in_hole, wound_smoothed,
             pigment_over_wound, tip_exist, tip_extension,
             new_corallites_on_tip, algae_on_wound),
           ~ as.integer(str_to_lower(str_squish(.x)) == "yes"))
  )

saveRDS(ph, file.path(DATA_PROC, "physio_clean.rds"))

# ---- Long-form for plotting ------------------------------------------------
traits <- c("polyps_out", "hole_in_center", "polyp_in_hole",
            "wound_smoothed", "pigment_over_wound", "tip_exist",
            "tip_extension", "new_corallites_on_tip", "algae_on_wound")

long <- ph |>
  filter(wound == "yes") |>
  select(day, treatment, tank, thicket, id, all_of(traits)) |>
  pivot_longer(all_of(traits), names_to = "trait", values_to = "expressed") |>
  filter(!is.na(expressed))

# proportion of wounded corals expressing trait by day × treatment
prop_df <- long |>
  group_by(day, treatment, trait) |>
  summarise(prop = mean(expressed), n = n(), .groups = "drop") |>
  mutate(trait = factor(trait, levels = traits,
                        labels = c("Polyps out", "Hole in center",
                                   "Polyp in hole", "Wound smoothed",
                                   "Pigment over wound", "Tip exists",
                                   "Tip extension", "New corallites on tip",
                                   "Algae on wound")))

p_traits <- ggplot(prop_df, aes(day, prop,
                                colour = treatment, group = treatment)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.4) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  facet_wrap(~ trait, ncol = 4) +
  scale_colour_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                      name = "Temperature") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(x = "Day relative to wounding (D0)",
       y = "Wounded corals expressing trait",
       title = "Morphological wound-healing characteristics",
       subtitle = "Wounded corals only (n = 24 per treatment)") +
  theme_pub(9)

save_fig(p_traits, "04_morphology_trajectories", width = 200, height = 130)

# ---- Per-genet trajectory plot --------------------------------------------
prop_genet_df <- long |>
  group_by(day, treatment, thicket, trait) |>
  summarise(prop = mean(expressed), n = n(), .groups = "drop") |>
  mutate(trait = factor(trait, levels = traits,
                        labels = c("Polyps out", "Hole in center",
                                   "Polyp in hole", "Wound smoothed",
                                   "Pigment over wound", "Tip exists",
                                   "Tip extension", "New corallites on tip",
                                   "Algae on wound")))

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
# Uses tank and coral ID random intercepts; thicket becomes a fixed factor.
fit_one <- function(trait_name) {
  d <- long |>
    filter(trait == trait_name, day >= 0) |>
    mutate(thicket = factor(thicket))
  if (length(unique(d$expressed)) < 2) return(NULL)
  fit <- tryCatch(
    suppressMessages(suppressWarnings(
      lme4::glmer(expressed ~ treatment * day * thicket + (1 | tank) + (1 | id),
                  family = binomial, data = d,
                  control = lme4::glmerControl(optimizer = "bobyqa",
                                               optCtrl = list(maxfun = 2e5)))
    )),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  tidy_fit <- broom.mixed::tidy(fit, effects = "fixed") |>
    mutate(trait = trait_name, n = nrow(d))
  # Type-II Wald ANOVA — reports significance of treatment, day, thicket and
  # their interactions
  av <- as.data.frame(car::Anova(fit, type = 2)) |>
    tibble::rownames_to_column("term") |>
    mutate(trait = trait_name)
  list(tidy = tidy_fit, anova = av)
}
trait_results <- map(traits, fit_one) |> setNames(traits) |> compact()
trait_tidy <- map_dfr(trait_results, "tidy")
trait_anova <- map_dfr(trait_results, "anova")
write_csv(trait_tidy,  file.path(TBL_DIR, "04_morphology_trait_glmm_summaries.csv"))
write_csv(trait_anova, file.path(TBL_DIR, "04_morphology_trait_anova_genet.csv"))

cat("Wrote physio_clean.rds, 04_morphology_trajectories.{pdf,png},",
    "04b_morphology_trajectories_by_genet.{pdf,png},",
    "04_morphology_trait_glmm_summaries.csv,",
    "04_morphology_trait_anova_genet.csv\n")
