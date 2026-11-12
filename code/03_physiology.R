# =============================================================================
# 03_physiology.R: photosynthesis, color, symbionts and growth (Figure 2)
#
# Models (lme4 + lmerTest, REML, Satterthwaite F tests, type III):
#   Fv/Fm, color: temperature x wound x day x source patch, day as a linear
#     trend (Day 0 onward). Random effects: an intercept and a day slope for
#     each tank (uncorrelated), and an intercept for each fragment.
#   Symbionts:    log(cells per cm2), same fixed effects with biopsy day;
#                 an intercept and a day slope for each tank (each fragment
#                 was sampled once).
#   Temperature was applied to tanks, so the temperature x day test has to be
#   judged against how much day slopes vary among tanks (about 6 denominator
#   df). Intercept-only models overstated it (F 95-288 on 160-230 df; see
#   docs/diagnostics/MODEL_DIAGNOSTICS_2026-11-11.md). Some fits put the tank
#   intercept variance at zero (singular); the slope variance is what matters.
#   Growth:       % skeletal mass change Day 0-15, temperature x wound x
#                 source patch, random intercept for tank. This is the
#                 reported temperature test (decided 2026-11-10); the exact
#                 tank permutation in section 4 is the supplement check. The
#                 result depends on ambient tank 3 (see 07_supplement.R).
#   Per-patch heat effects (emmeans, Day 14) carry a Holm-adjusted p across
#   the 6 patch x wound contrasts of each response.
#
# Writes: output/tables/12_anova_summary.csv, 12_genet_treatment_effects.csv,
#         03_color_end_proportions.csv, 05_buoyant_weight_tank_test.csv,
#         06_symbiont_chl_summary.csv, 45_figure2_summary.csv, 45_figure2_tank_means.csv
#         figures/45_physiology_figure2 (Figure 2; panels 45a-c), 03_color_trajectory,
#         05_buoyant_weight_growth, 06_symbiont_density_by_day
#         output/models/12_*.rds (not tracked; read by 08_checks.R)
#         output/paper/physiology_numbers.csv
# =============================================================================

source(here::here("code", "00_setup.R"))
dir.create(MOD_DIR, recursive = TRUE, showWarnings = FALSE)

pam   <- readRDS(file.path(DATA_PROC, "pam_clean.rds"))
color <- readRDS(file.path(DATA_PROC, "color_clean.rds"))
bw    <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds"))
sym   <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds"))

ctrl <- lme4::lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))
anova_rows <- list()   # type III F tests, one table per response
patch_rows <- list()   # heat effect (28C - 31C) within each source patch x wound

record_anova <- function(model, name) {
  av <- as.data.frame(anova(model, type = 3))
  av$term <- rownames(av); av$response <- name
  anova_rows[[name]] <<- av
}
record_patch_effect <- function(emm, name) {
  contr <- as_tibble(pairs(emm, by = c("thicket", "wound"), adjust = "none"))
  contr$response <- name
  contr$p_holm <- p.adjust(contr$p.value, "holm")   # 6 patch x wound contrasts
  patch_rows[[name]] <<- contr
}

# ---- 1. Models -------------------------------------------------------------
# Day -1 (Fv/Fm) and Day -4 (color) are pre-treatment baselines: shown in the
# figures but left out of the models, so the day slope is post-clipping only.
pam_m <- pam |> mutate(thicket = as.factor(thicket)) |> filter(day >= 0)
m_pam <- lmerTest::lmer(fv_fm ~ treatment * wound * day * thicket + (1 + day || tank) + (1 | id),
                        data = pam_m, REML = TRUE, control = ctrl)
saveRDS(m_pam, file.path(MOD_DIR, "12_pam_lmm.rds"))
record_anova(m_pam, "pam_fvfm")
record_patch_effect(emmeans::emmeans(m_pam, ~ treatment | thicket * wound, at = list(day = 14)),
                    "pam_fvfm")

color_m <- color |> mutate(thicket = as.factor(thicket)) |> filter(day >= 0)
m_color <- lmerTest::lmer(color_num ~ treatment * wound * day * thicket + (1 + day || tank) + (1 | id),
                          data = color_m, REML = TRUE, control = ctrl)
saveRDS(m_color, file.path(MOD_DIR, "12_color_lmm.rds"))
record_anova(m_color, "color_dscale")
record_patch_effect(emmeans::emmeans(m_color, ~ treatment | thicket * wound, at = list(day = 14)),
                    "color_dscale")

bw_m <- bw |> mutate(thicket = as.factor(thicket)) |> filter(is.finite(pct_growth))
m_bw <- lmerTest::lmer(pct_growth ~ treatment * wound * thicket + (1 | tank),
                       data = bw_m, REML = TRUE, control = ctrl)
saveRDS(m_bw, file.path(MOD_DIR, "12_bw_lm.rds"))
record_anova(m_bw, "growth_pct")
record_patch_effect(emmeans::emmeans(m_bw, ~ treatment | thicket * wound), "growth_pct")

# Biopsy day is centred on the first biopsy (Day 1).
sym_m <- sym |>
  filter(is.finite(cells_per_cm2), cells_per_cm2 > 0) |>
  mutate(thicket = factor(thicket), biopsy_day_c = biopsy_day - 1)
m_zoox <- lmerTest::lmer(log(cells_per_cm2) ~ treatment * wound * biopsy_day_c * thicket + (1 + biopsy_day_c || tank),
                         data = sym_m, REML = TRUE, control = ctrl)
saveRDS(m_zoox, file.path(MOD_DIR, "12_zoox_lmm.rds"))
record_anova(m_zoox, "log_zoox_density")
record_patch_effect(emmeans::emmeans(m_zoox, ~ treatment | thicket * wound,
                                     at = list(biopsy_day_c = 14)), "log_zoox_density")

write_csv(bind_rows(anova_rows, .id = "response_id") |>
            mutate(across(where(is.numeric), \(x) round(x, 5))),
          file.path(TBL_DIR, "12_anova_summary.csv"))
write_csv(bind_rows(patch_rows) |> mutate(across(where(is.numeric), \(x) round(x, 5))),
          file.path(TBL_DIR, "12_genet_treatment_effects.csv"))

# ---- 2. Color: end-of-experiment summary and trajectory figure -------------
end_day <- max(color$day, na.rm = TRUE)             # last day with color data
end_props <- color |>
  filter(day == end_day) |>
  group_by(treatment, wound) |>
  summarise(
    n_total       = n(),
    n_alive       = sum(health_status == "alive", na.rm = TRUE),
    n_paled       = sum(paling == "yes", na.rm = TRUE),
    mean_color    = mean(color_num, na.rm = TRUE),
    sd_color      = sd(color_num,   na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(prop_paled = n_paled / n_total)
write_csv(end_props, file.path(TBL_DIR, "03_color_end_proportions.csv"))
traj <- color |>
  group_by(day, treatment, wound) |>
  summarise(
    mean_color = mean(color_num, na.rm = TRUE),
    se = sd(color_num, na.rm = TRUE) / sqrt(n()),   # SE = SD / sqrt(n)
    n = n(),
    .groups = "drop"
  )
p_color <- ggplot(traj, aes(day, mean_color,
                            colour = wound, fill = wound, group = wound)) +
  geom_ribbon(aes(ymin = mean_color - se, ymax = mean_color + se),
              alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  facet_wrap(~ treatment, ncol = 2) +
  scale_y_reverse(breaks = 1:6, limits = c(6.2, 0.8)) +  # darker = higher D, plot 1 on top
  scale_colour_manual(values = PAL_WOUND,
                      name = "Wound") +
  scale_fill_manual(values   = PAL_WOUND,
                    guide = "none") +
  labs(x = "Day relative to wounding (D0)",
       y = "Color score (Siebeck D-scale; lower = paler)",
       title = "Pigmentation trajectory",
       subtitle = "Mean ± 1 SE; paler scores appear higher on the reversed axis") +
  theme_pub(10)
save_fig(p_color, "03_color_trajectory", width = 170, height = 90)

# ---- 3. Symbionts: summary by group and day --------------------------------
summary_tbl <- sym |>
  group_by(treatment, wound, biopsy_day) |>
  summarise(
    n             = n(),
    mean_cells    = mean(cells_per_cm2, na.rm = TRUE),
    se_cells      = sd(cells_per_cm2, na.rm = TRUE) / sqrt(n()),
    mean_chl      = mean(chlorophyll_ug_cm2, na.rm = TRUE),   # NaN when chl absent
    se_chl        = sd(chlorophyll_ug_cm2, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
write_csv(summary_tbl, file.path(TBL_DIR, "06_symbiont_chl_summary.csv"))
p_zoox <- ggplot(sym, aes(factor(biopsy_day), cells_per_cm2 / 1e6,
                            fill = treatment)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.65) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.4, size = 0.8) +
  scale_fill_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                    name = "Temperature") +
  labs(x = "Biopsy day",
       y = expression(Symbionts~(10^6~cells~cm^{-2})),
       title = "Symbiont density loss under heating",
       subtitle = expression(italic(A.~pulchra)~"biopsies, n = 192 across 4 timepoints")) +
  theme_pub(10)
set.seed(42)   # geom_jitter: keep point positions stable between runs
save_fig(p_zoox, "06_symbiont_density_by_day", width = 140, height = 100)

# ---- 4. Growth: tank-level permutation check and figure --------------------
# Temperature was set per tank, so this compares the 4 ambient and 4 heated
# tank means against all 70 ways of splitting 8 tanks into two groups of 4.
bw_a <- bw |> filter(is.finite(pct_growth))
tank_growth <- bw_a |>
  group_by(tank, treatment) |>
  summarise(mean_pct_growth = mean(pct_growth, na.rm = TRUE),
            n_corals = n(), .groups = "drop")
obs_diff <- with(tank_growth,
                 mean(mean_pct_growth[treatment == "28C"]) -
                   mean(mean_pct_growth[treatment == "31C"]))
vals <- tank_growth$mean_pct_growth
assignments <- combn(seq_along(vals), sum(tank_growth$treatment == "28C"),
                     simplify = FALSE)
perm_diffs <- vapply(assignments, function(idx) {
  mean(vals[idx]) - mean(vals[-idx])
}, numeric(1))
tank_test <- tibble(
  response = "pct_growth",
  test = "tank-level exact permutation",
  n_tanks = nrow(tank_growth),
  n_28_tanks = sum(tank_growth$treatment == "28C"),
  n_31_tanks = sum(tank_growth$treatment == "31C"),
  estimate_28_minus_31 = obs_diff,
  p_two_sided = mean(abs(perm_diffs) >= abs(obs_diff)),
  note = "Temperature was randomized at the tank level; this p-value treats tanks, not corals, as exchangeable units."
) |>
  mutate(across(where(is.numeric), \(x) signif(x, 4)))
write_csv(tank_test, file.path(TBL_DIR, "05_buoyant_weight_tank_test.csv"))
p_bw <- ggplot(bw_a, aes(interaction(treatment, wound, sep = " · "),
                         pct_growth, fill = treatment)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(colour = thicket), width = 0.18, height = 0,
              size = 1.6, alpha = 0.85) +
  scale_fill_manual(values  = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                    name = "Temperature") +
  scale_colour_manual(values = PAL_GENO, name = "Genotype") +
  labs(x = NULL,
       y = "Growth (% skeletal mass change over 15 d)",
       title = "Coral growth under heating × wounding",
       subtitle = "Buoyant-weight % mass change (Davies 1989; Jokiel et al. 1978)") +
  theme_pub(10) +
  theme(axis.text.x = element_text(angle = 0))
set.seed(42)   # geom_jitter
save_fig(p_bw, "05_buoyant_weight_growth", width = 150, height = 95)

# ---- 5. Figure 2 -----------------------------------------------------------
# A: Fv/Fm, B: symbiont density, C: growth. Error bars are SE across the four
# tank means per temperature (tanks are the replicates for temperature).
pam_f2 <- readRDS(file.path(DATA_PROC, "pam_clean.rds")) |>
  transmute(id, tank, treatment, wound, day, value = fv_fm)
sym_f2 <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
  transmute(id, tank, treatment, wound, day = biopsy_day, value = cells_per_cm2 / 1e6)
growth <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds")) |>
  filter(biopsy_day == 15, is.finite(pct_growth)) |>
  transmute(id, tank, treatment, wound, day = 15L, value = pct_growth)
stopifnot(nrow(pam_f2) == 336, nrow(sym_f2) == 192, nrow(growth) == 48,
          !anyDuplicated(pam_f2[c("id", "day")]), !anyDuplicated(sym_f2$id),
          !anyDuplicated(growth$id))
raw <- bind_rows(PAM = pam_f2, Symbionts = sym_f2, Growth = growth, .id = "response")
stopifnot(all(is.finite(raw$value)), all(!is.na(raw$wound)))
tanks <- raw |>
  group_by(response, treatment, wound, day, tank) |>
  summarise(value = mean(value), n_fragments = n(), .groups = "drop")
means <- tanks |>
  group_by(response, treatment, wound, day) |>
  summarise(mean = mean(value), se = sd(value) / sqrt(n()), n_tanks = n(),
            n_fragments = sum(n_fragments), .groups = "drop")
stopifnot(all(means$n_tanks == 4), all(tanks$n_fragments >= 2),
          all(tanks$n_fragments <= 4))
write_csv(tanks, file.path(TBL_DIR, "45_figure2_tank_means.csv"))
write_csv(means, file.path(TBL_DIR, "45_figure2_summary.csv"))
common <- function(lines = TRUE) list(
  facet_wrap(~treatment, ncol = 2,
             labeller = as_labeller(c(`28C` = "28 \u00b0C", `31C` = "31 \u00b0C"))),
  scale_colour_manual(values = PAL_TEMP, guide = "none"),
  scale_shape_manual(values = c(no = 1, yes = 16), name = NULL,
                     labels = c(no = "Unwounded", yes = "Wounded")),
  if (lines) scale_linetype_manual(values = c(no = "dashed", yes = "solid"), name = NULL,
                                  labels = c(no = "Unwounded", yes = "Wounded")),
  theme_pub(11),
  theme(legend.position = "bottom", plot.title = element_text(size = 13, face = "bold"),
        panel.grid.minor = element_blank(), strip.text = element_text(face = "bold")))
trajectory <- function(response_name, title, ylab, breaks) {
  d <- filter(means, response == response_name)
  ggplot(d, aes(day, mean, colour = treatment, shape = wound,
                linetype = wound, group = wound)) +
    geom_vline(xintercept = 0, colour = "grey70", linetype = "dotted") +
    geom_line(linewidth = 0.7) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  width = 0.25, position = position_dodge(width = 0.3),
                  linewidth = 0.45, show.legend = FALSE) +
    geom_point(size = 2.6, stroke = 0.9, position = position_dodge(width = 0.3)) +
    scale_x_continuous(breaks = breaks, limits = c(-1.5, 15.5)) +
    labs(title = title, x = "Day relative to clipping (Day 0)", y = ylab) + common()
}
a <- trajectory("PAM", "A  Photosynthesis efficiency", expression(F[v]/F[m]),
                c(-1, 0, 3, 6, 9, 12, 14))
box_style <- list(
  scale_fill_manual(values = c(no = "white", yes = "grey65"), name = NULL,
                    labels = c(no = "Unwounded", yes = "Wounded")),
  common(lines = FALSE))
b <- ggplot(sym_f2, aes(factor(day), value, colour = treatment,
                     fill = wound, group = interaction(day, wound))) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.65,
               outlier.shape = NA, linewidth = 0.5) +
  geom_point(aes(shape = wound), colour = "grey35", size = 1.2,
             position = position_jitterdodge(jitter.width = 0.12,
                                           dodge.width = 0.8, seed = 42)) +
  labs(title = "B  Algal symbiont density", x = "Sampling day",
       y = expression(10^6~cells~cm^{-2})) + box_style
c <- ggplot(growth, aes(wound, value, colour = treatment, fill = wound)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.5) +
  geom_point(aes(shape = wound), colour = "grey35", size = 1.7,
             position = position_jitter(width = 0.12, height = 0, seed = 42)) +
  scale_x_discrete(labels = c(no = "Unwounded", yes = "Wounded")) +
  labs(title = "C  Skeletal growth over 15 days", x = NULL,
       y = "Skeletal mass gain (%)") + box_style + theme(legend.position = "none")
figure2 <- (a / (b + guides(fill = "none", shape = "none")) /
              (c + guides(fill = "none", shape = "none"))) +
  plot_layout(guides = "collect", heights = c(1, 1, 0.95)) &
  theme(legend.position = "bottom")
save_fig(figure2, "45_physiology_figure2", width = 185, height = 235)
save_fig(a, "45a_photosynthesis", width = 185, height = 85)
save_fig(b, "45b_symbiont_boxplots", width = 185, height = 85)
save_fig(c, "45c_growth_boxplots", width = 185, height = 85)
stopifnot(all(raw$wound[raw$id == 119L] == "yes"),
          all(count(growth, treatment, wound)$n == 12L))
cat("Figure 2 verified: 48 repeated PAM fragments, 192 destructive symbiont samples,",
    "48 endpoint growth fragments; four tanks per temperature.\n")

# ---- 6. Numbers for the results text ---------------------------------------
anova <- readr::read_csv(file.path(TBL_DIR, "12_anova_summary.csv"), show_col_types = FALSE)
fig2  <- readr::read_csv(file.path(TBL_DIR, "45_figure2_summary.csv"), show_col_types = FALSE)
colr  <- readr::read_csv(file.path(TBL_DIR, "03_color_end_proportions.csv"), show_col_types = FALSE)
zoox  <- readr::read_csv(file.path(TBL_DIR, "06_symbiont_chl_summary.csv"), show_col_types = FALSE)
tank  <- readr::read_csv(file.path(TBL_DIR, "05_buoyant_weight_tank_test.csv"), show_col_types = FALSE)

f_row <- function(key, description, response, term) {
  d <- anova[anova$response_id == response & anova$term == term, ]
  stopifnot(nrow(d) == 1)
  paper_row(key, description, test = "F test (Satterthwaite)",
            statistic = d$`F value`, df1 = d$NumDF, df2 = d$DenDF, p_value = d$`Pr(>F)`,
            source = "output/tables/12_anova_summary.csv")
}
f2 <- function(resp, tr, day = NULL) {
  d <- fig2[fig2$response == resp & fig2$treatment == tr, ]
  if (!is.null(day)) d <- d[d$day == day, ]
  mean(d$mean)
}

rows <- list(
  f_row("pam_temp_x_day", "Fv/Fm: temperature x day", "pam_fvfm", "treatment:day"),
  f_row("pam_temp_x_wound_x_day", "Fv/Fm: temperature x wound x day", "pam_fvfm", "treatment:wound:day"),
  f_row("pam_wound", "Fv/Fm: wound", "pam_fvfm", "wound"),
  paper_row("pam_day14_ambient", "Fv/Fm on Day 14, ambient (mean of wound groups)",
            estimate = f2("PAM", "28C", 14), source = "output/tables/45_figure2_summary.csv"),
  paper_row("pam_day14_heated", "Fv/Fm on Day 14, heated (mean of wound groups)",
            estimate = f2("PAM", "31C", 14), source = "output/tables/45_figure2_summary.csv"),
  f_row("color_temp_x_day", "Color: temperature x day", "color_dscale", "treatment:day"),
  f_row("color_temp_x_wound_x_day", "Color: temperature x wound x day", "color_dscale", "treatment:wound:day"),
  paper_row("color_final_ambient", "Final color score (1-6), ambient",
            estimate = mean(colr$mean_color[colr$treatment == "28C"]),
            source = "output/tables/03_color_end_proportions.csv"),
  paper_row("color_final_heated", "Final color score (1-6), heated",
            estimate = mean(colr$mean_color[colr$treatment == "31C"]),
            source = "output/tables/03_color_end_proportions.csv"),
  f_row("zoox_temp_x_day", "Symbionts: temperature x day", "log_zoox_density", "treatment:biopsy_day_c"),
  f_row("zoox_temp_x_wound_x_day", "Symbionts: temperature x wound x day", "log_zoox_density", "treatment:wound:biopsy_day_c"),
  paper_row("zoox_day15_ambient", "Symbionts on Day 15, ambient (million cells per cm2)",
            estimate = mean(zoox$mean_cells[zoox$treatment == "28C" & zoox$biopsy_day == 15]) / 1e6,
            source = "output/tables/06_symbiont_chl_summary.csv"),
  paper_row("zoox_day15_heated", "Symbionts on Day 15, heated (million cells per cm2)",
            estimate = mean(zoox$mean_cells[zoox$treatment == "31C" & zoox$biopsy_day == 15]) / 1e6,
            source = "output/tables/06_symbiont_chl_summary.csv"),
  f_row("growth_temp", "Growth: temperature (reported test)", "growth_pct", "treatment"),
  f_row("growth_wound", "Growth: wound", "growth_pct", "wound"),
  f_row("growth_temp_x_wound", "Growth: temperature x wound", "growth_pct", "treatment:wound"),
  paper_row("growth_ambient", "Growth over 15 days (% mass change), ambient",
            estimate = f2("Growth", "28C"), source = "output/tables/45_figure2_summary.csv"),
  paper_row("growth_heated", "Growth over 15 days (% mass change), heated",
            estimate = f2("Growth", "31C"), source = "output/tables/45_figure2_summary.csv"),
  paper_row("growth_pct_lower", "Growth under heat, % lower than ambient",
            estimate = 100 * (1 - f2("Growth", "31C") / f2("Growth", "28C")),
            source = "output/tables/45_figure2_summary.csv"),
  paper_row("growth_tank_permutation", "Growth: tank permutation test (supplement check)",
            test = "exact permutation of 4 vs 4 tank means", p_value = tank$p_two_sided,
            source = "output/tables/05_buoyant_weight_tank_test.csv")
)
# Direction of the clipping difference in heated symbionts (open question):
# wounded / unwounded mean density on Days 10 and 15, and the number of heated
# tanks where wounded fragments kept more symbionts on both days.
wound_ratio <- sym |>
  filter(treatment == "31C", biopsy_day %in% c(10, 15)) |>
  group_by(tank, biopsy_day, wound) |>
  summarise(m = mean(cells_per_cm2), .groups = "drop") |>
  tidyr::pivot_wider(names_from = wound, values_from = m) |>
  group_by(tank) |>
  summarise(both_higher = all(yes > no), .groups = "drop")
for (d in c(10, 15)) {
  z <- sym |> filter(treatment == "31C", biopsy_day == d) |>
    group_by(wound) |> summarise(m = mean(cells_per_cm2))
  rows[[length(rows) + 1]] <- paper_row(
    paste0("zoox_heated_wound_ratio_day", d),
    paste0("Symbionts at 31 C, Day ", d, ": wounded / unwounded mean density"),
    estimate = z$m[z$wound == "yes"] / z$m[z$wound == "no"],
    source = "data/processed/symbiont_chl_clean.rds")
}
rows[[length(rows) + 1]] <- paper_row(
  "zoox_heated_tanks_wounded_higher", "Heated tanks where wounded fragments had more symbionts on both Day 10 and Day 15",
  estimate = sum(wound_ratio$both_higher), n = nrow(wound_ratio),
  source = "data/processed/symbiont_chl_clean.rds")
write_paper(rows, "physiology")
