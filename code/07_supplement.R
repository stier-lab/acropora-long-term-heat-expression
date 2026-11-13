# =============================================================================
# 07_supplement.R: supplement-only results
#
# Analyses that appear only in the supplement. They describe the microscope
# experiment and the other daily observations, and test how much the main
# results depend on a few fragments or tanks ("sensitivity analyses": refit
# the same model on a smaller data set and see whether the answer changes).
# Tank is the replicate for temperature: the 4 ambient and 4 heated tanks, not
# the fragments, received the temperature treatment.
#
# Sections
#   1. Microscope experiment: share of the 16 fragments showing each stage on
#      each day, by temperature (Supp. fig. 3)
#   2. Other daily observations: polyps, pigment, algae, flatworms
#      (Supp. fig. 2)
#   3. Results without the flagged fragments: refit the Fv/Fm and color
#      temperature x day tests and the growth tank-shuffling test without
#      fragments 116 and 121 and without tank 3 (slow to reach temperature)
#   4. How much the symbiont wound result and the growth result depend on
#      single fragments or tanks
#   5. Color: tank-level shuffling check (no bell curve assumed)
#
# Other supplement material is written by the main scripts: the growth
# tank-shuffling test (03_physiology.R), the time-to-stage curves and models
# (Kaplan-Meier and Cox; 04_regeneration.R), and the wax standard curve
# (01_clean_data.R).
#
# Reads:  data/processed/microscope_physio_clean.rds, physio_clean.rds,
#         pam_clean.rds, color_clean.rds, buoyant_weight_clean.rds,
#         symbiont_chl_clean.rds
# Writes: output/tables/11_microscope_trait_trajectory.csv,
#         47_other_daily_observations.csv, 22_sensitivity_flagged.csv,
#         22b_sensitivity_wound_growth.csv, 12b_color_tank_test.csv
#         figures/11_microscope_trait_trajectories (Supp. fig. 3),
#         47_other_daily_observations (Supp. fig. 2)
# Run:    sourced by code/run_all.R (or: Rscript code/run_all.R)
# =============================================================================

source(here::here("code", "00_setup.R"))


# ---- 1. Microscope experiment: stage trajectories --------------------------
# Supplementary figure 3 (Molly's numbering): healing plus the five
# regeneration stages, in manuscript order and terms. Plots the share of the
# 8 fragments per temperature scored "yes" on each day. In the labels,
# coenosarc = the living tissue sheet between polyps; corallite = the skeletal
# cup a polyp sits in (axial = the tip cup, radial = the side cups).
micro <- readRDS(file.path(DATA_PROC, "microscope_physio_clean.rds"))

# Lookup table: raw column name -> stage label. Setting the label order to
# this row order fixes the panel order in the figure.
trait_meta <- tibble::tribble(
  ~trait,                  ~label,
  "tissue_over_wound",     "Healing: coenosarc covers wound",
  "axial_polyp_formation", "1. Axial polyp forms",
  "wound_smoothed",        "2. Wound surface remodels",
  "tip_exist",             "3. Axial corallite forms",
  "tip_extension",         "4. Axial corallite extends",
  "new_corallites_on_tip", "5. Radial corallites form"
) |>
  mutate(label = factor(label, levels = label))
traits <- trait_meta$trait

# One row per fragment x day x stage (0/1); unscored cells are dropped.
long <- micro |>
  select(day, date, treatment, tank, thicket, id, all_of(traits)) |>
  pivot_longer(all_of(traits), names_to = "trait", values_to = "expressed") |>
  left_join(trait_meta, by = "trait") |>
  filter(!is.na(expressed))

# Share scored "yes" per day x temperature x stage (mean of a 0/1 column).
prop_df <- long |>
  group_by(day, treatment, trait, label) |>
  summarise(
    n = n(),
    n_yes = sum(expressed == 1),
    prop = mean(expressed),
    .groups = "drop"
  )
write_csv(prop_df, file.path(TBL_DIR, "11_microscope_trait_trajectory.csv"))

p_micro <- ggplot(prop_df, aes(day, prop, colour = treatment,
                               group = treatment)) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = PAL_TEMP, name = "Temperature",
                      labels = c(`28C` = "28 \u00b0C", `31C` = "31 \u00b0C")) +
  scale_x_continuous(breaks = seq(0, 15, by = 3), limits = c(0, 15)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  facet_wrap(~ label, ncol = 3, labeller = label_wrap_gen(width = 22)) +
  labs(x = "Days after wounding", y = "Fragments showing the stage (%)") +
  theme_pub(9)
save_fig(p_micro, "11_microscope_trait_trajectories",
         width = 190, height = 130)


# ---- 2. Other daily observations, main experiment (Supp. fig. 2) -----------
# Polyp extension (all 48 fragments; unwounded fragments were not scored on
# Days 6 and 10-15), pigment over the wound (wounded fragments; scored on
# Days 3, 4, 9, 12 and 15), algae on the wound (wounded; scored daily) and
# Acropora-eating flatworm (AEF) signs, which were recorded only when seen.
# On Day 4 every fragment was dipped in the AEF chemical treatment and all
# polyps were retracted (0/48); that dip is a handling effect, not heat.
ph_obs <- readRDS(file.path(DATA_PROC, "physio_clean.rds"))
obs_panels <- c("Polyps extended", "Pigment over wound (wounded only)",
                "Algae on wound (wounded only)", "Flatworm (AEF) signs recorded")

# For each observation: k of n fragments scored "yes" per day x temperature x
# wound, stacked into one table with a panel label. AEF has no NA filter
# because it was recorded only when seen (no record = not seen).
daily_obs <- bind_rows(
  ph_obs |> filter(!is.na(polyps_out)) |>
    group_by(day, treatment, wound) |>
    summarise(n = n(), k = sum(polyps_out), .groups = "drop") |>
    mutate(panel = obs_panels[1]),
  ph_obs |> filter(wound == "yes", !is.na(pigment_over_wound)) |>
    group_by(day, treatment, wound) |>
    summarise(n = n(), k = sum(pigment_over_wound), .groups = "drop") |>
    mutate(panel = obs_panels[2]),
  ph_obs |> filter(wound == "yes", !is.na(algae_on_wound)) |>
    group_by(day, treatment, wound) |>
    summarise(n = n(), k = sum(algae_on_wound), .groups = "drop") |>
    mutate(panel = obs_panels[3]),
  ph_obs |> group_by(day, treatment, wound) |>
    summarise(n = n(), k = sum(disease_status %in% "AEF"), .groups = "drop") |>
    mutate(panel = obs_panels[4])
) |>
  mutate(pct = 100 * k / n, panel = factor(panel, levels = obs_panels))
write_csv(daily_obs, file.path(TBL_DIR, "47_other_daily_observations.csv"))

p_obs <- ggplot(daily_obs, aes(day, pct, colour = treatment, linetype = wound,
                               group = interaction(treatment, wound))) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.4) +
  facet_wrap(~ panel, ncol = 2) +
  scale_colour_manual(values = PAL_TEMP, name = "Temperature",
                      labels = c(`28C` = "28 \u00b0C", `31C` = "31 \u00b0C")) +
  scale_linetype_manual(values = c(no = "dashed", yes = "solid"), name = "Wounded",
                        labels = c(no = "No", yes = "Yes")) +
  scale_x_continuous(breaks = seq(0, 15, by = 3), limits = c(0, 15)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Days after wounding", y = "Fragments (%)") +
  theme_pub(9)
save_fig(p_obs, "47_other_daily_observations", width = 180, height = 130)


# ---- 3. Results without the flagged fragments and tank ---------------------
# Refit each main test on the full data and on the data without fragments 116
# and 121 and tank 3, and put the two side by side. Fv/Fm (how efficiently the
# algae's photosynthesis works, 0-1, measured with the PAM fluorometer) and
# color use the 03_physiology.R models (each tank has its own change over
# days); growth uses the tank-shuffling test.

# Fitting settings shared by every lmer() fit in this script. lme4's
# "singular fit" warning (it estimated no difference among tanks) is switched
# off because the 03 models already report it.
ctrl <- lme4::lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4),
                          optimizer = "bobyqa")
FLAGGED_IDS  <- c(116, 121)   # two fragments flagged as behaving anomalously
FLAGGED_TANK <- 3             # ambient tank that was slow to reach temperature
results <- list()             # one table per (response × full-or-dropped) goes here

# grab(): pull the F and p for one term from the model's test table. type = 3
# tests each term after allowing for all the others; lmerTest supplies the
# p-values (using the Satterthwaite method to estimate how much information
# each test has). Returns NULL if the term is absent.
grab <- function(model, label, term) {
  av <- as.data.frame(anova(model, type = 3))
  av$term <- rownames(av)
  row <- av[av$term == term, , drop = FALSE]
  if (nrow(row) == 0) return(NULL)
  tibble(
    response   = label,
    term       = term,
    F_value    = round(row[["F value"]], 2),
    p_value    = signif(row[["Pr(>F)"]], 3)
  )
}

## 3a. Fv/Fm and color: temperature x day ------------------------------------
# Day -1 / Day -4 baselines are dropped (day >= 0), as in 03_physiology.R.
# In the formulas, `(1 + day || tank)` lets each tank have its own starting
# level and its own change over days (`||` = the two are estimated
# separately); `(1 | id)` accounts for measuring the same fragment many
# times. REML = TRUE is lme4's standard (restricted-likelihood) fit.
pam <- readRDS(file.path(DATA_PROC, "pam_clean.rds")) |>
  mutate(thicket = factor(thicket)) |>
  filter(day >= 0)
pam_drop <- pam |> filter(!id %in% FLAGGED_IDS, tank != FLAGGED_TANK)
f_full <- lmerTest::lmer(fv_fm ~ treatment * wound * day * thicket +
                           (1 + day || tank) + (1 | id),
                         data = pam, REML = TRUE, control = ctrl)
f_drop <- lmerTest::lmer(fv_fm ~ treatment * wound * day * thicket +
                           (1 + day || tank) + (1 | id),
                         data = pam_drop, REML = TRUE, control = ctrl)
results$pam_full <- grab(f_full, "PAM Fv/Fm (full data)",     "treatment:day")
results$pam_drop <- grab(f_drop, "PAM Fv/Fm (flagged dropped)","treatment:day")

col <- readRDS(file.path(DATA_PROC, "color_clean.rds")) |>
  mutate(thicket = factor(thicket)) |>
  filter(day >= 0)
col_drop <- col |> filter(!id %in% FLAGGED_IDS, tank != FLAGGED_TANK)
c_full <- lmerTest::lmer(color_num ~ treatment * wound * day * thicket +
                           (1 + day || tank) + (1 | id),
                         data = col, REML = TRUE, control = ctrl)
c_drop <- lmerTest::lmer(color_num ~ treatment * wound * day * thicket +
                           (1 + day || tank) + (1 | id),
                         data = col_drop, REML = TRUE, control = ctrl)
results$col_full <- grab(c_full, "Color D-scale (full data)",     "treatment:day")
results$col_drop <- grab(c_drop, "Color D-scale (flagged dropped)","treatment:day")

## 3b. Growth: tank-shuffling test -------------------------------------------
bw <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds")) |>
  mutate(thicket = factor(thicket)) |>
  filter(is.finite(pct_growth))                 # drop NA/Inf growth values
bw_drop <- bw |> filter(!id %in% FLAGGED_IDS, tank != FLAGGED_TANK)

# tank_perm(): a permutation (shuffling) test on tank means. Average growth
# within each tank, take the ambient-minus-heated difference, and compare it
# with the difference from every possible way of relabelling the tanks as
# ambient or heated (combn() lists them all). p = share of relabellings with a
# difference at least as large. Returns NA if either temperature has fewer
# than 2 tanks.
tank_perm <- function(dat, label) {
  tank_growth <- dat |>
    group_by(tank, treatment) |>
    summarise(mean_pct_growth = mean(pct_growth, na.rm = TRUE),
              n_corals = n(), .groups = "drop")
  n_by_temp <- table(tank_growth$treatment)
  if (length(n_by_temp) != 2 || any(n_by_temp < 2)) {
    return(tibble(response = label, term = "treatment (tank permutation)",
                  F_value = NA_real_, p_value = NA_real_,
                  estimate_28_minus_31 = NA_real_, n_tanks = nrow(tank_growth)))
  }
  vals <- tank_growth$mean_pct_growth
  trt <- tank_growth$treatment
  n28 <- sum(trt == "28C")
  obs <- mean(vals[trt == "28C"]) - mean(vals[trt == "31C"])
  null <- vapply(combn(seq_along(vals), n28, simplify = FALSE), \(idx) {
    mean(vals[idx]) - mean(vals[-idx])
  }, numeric(1))
  p <- mean(abs(null) >= abs(obs))
  tibble(response = label, term = "treatment (tank permutation)",
         F_value = NA_real_, p_value = signif(p, 3),
         estimate_28_minus_31 = round(obs, 3),
         n_tanks = nrow(tank_growth))
}
results$bw_full <- tank_perm(bw, "Growth % (full data)")
results$bw_drop <- tank_perm(bw_drop, "Growth % (flagged dropped)")

out <- bind_rows(results)
write_csv(out, file.path(TBL_DIR, "22_sensitivity_flagged.csv"))


# ---- 4. How much do two results depend on single fragments or tanks? -------
# (a) Symbionts (the algae in the coral tissue), temperature x wound x day:
#     the fragments the outlier check flags (two or three, depending on the
#     DHARMa package version; versions are pinned in renv.lock) are all
#     heated: 163 (patch A, Day 15, low), 212 (patch A, Day 3, high) and,
#     under DHARMa 0.5.0, 215 (patch D, Day 10). Nothing marks them as errors,
#     so they stay in the reported model; this shows how the wound result
#     moves without them.
# (b) Growth, temperature: ambient tank 3 (slow to reach temperature) grew
#     like the heated tanks.

## 4a. Setup: symbiont data, model formulas and helpers ----------------------
# Same data and model as 03_physiology.R: log density, biopsy day counted from
# the first biopsy (Day 1 = 0), and each tank with its own change over days.
sym <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
  filter(is.finite(cells_per_cm2), cells_per_cm2 > 0) |>
  mutate(thicket = factor(thicket), biopsy_day_c = biopsy_day - 1)
f_sym <- log(cells_per_cm2) ~ treatment * wound * biopsy_day_c * thicket +
  (1 + biopsy_day_c || tank)

# term_p(): F and p for one term from the same type = 3 test table.
term_p <- function(m, term) {
  a <- as.data.frame(anova(m, type = 3)); c(F = a[term, "F value"], p = a[term, "Pr(>F)"])
}

# add_sens(): append one row to the `sens` list. `<<-` writes to the `sens`
# defined outside the function rather than to a local copy.
sens <- list()
add_sens <- function(response, term, version, Fp, note = NA_character_)
  sens[[length(sens) + 1]] <<- tibble(response, term, version,
                                      F_value = unname(Fp[1]), p_value = unname(Fp[2]),
                                      note)
k3 <- "treatment:wound:biopsy_day_c"

## 4b. Symbionts: temperature x wound x day ----------------------------------
m_sym <- lmerTest::lmer(f_sym, data = sym, control = ctrl)
add_sens("symbionts", k3, "as reported (03_physiology.R)", term_p(m_sym, k3))

# DHARMa simulates 500 new data sets from the fitted model and flags any
# observation that falls outside every simulated value. The random seed is
# reset so the same rows are flagged on every run.
set.seed(42)
out_rows <- DHARMa::outliers(DHARMa::simulateResiduals(m_sym, n = 500))
# If nothing is flagged, sym[-integer(0), ] would return zero rows, so keep
# all rows explicitly in that case.
sym_no_out <- if (length(out_rows) > 0) sym[-out_rows, ] else sym
add_sens("symbionts", k3, "without residual outliers",
         term_p(lmerTest::lmer(f_sym, data = sym_no_out, control = ctrl), k3),
         if (length(out_rows) > 0) {
           paste("dropped fragments", paste(sort(sym$id[out_rows]), collapse = ", "))
         } else {
           "no fragments flagged; same as reported model"
         })

# The same model refitted with the nlme package, letting the scatter around
# the model differ between temperatures (varIdent), in case heated fragments
# are simply more variable. pdDiag gives each tank its own starting level and
# change over days, estimated separately, like `||` in lme4.
m_var <- nlme::lme(log(cells_per_cm2) ~ treatment * wound * biopsy_day_c * thicket,
                   random = list(tank = nlme::pdDiag(~ 1 + biopsy_day_c)),
                   weights = nlme::varIdent(form = ~ 1 | treatment), data = sym,
                   control = nlme::lmeControl(opt = "optim"))
a_var <- anova(m_var, type = "marginal")
add_sens("symbionts", k3, "separate residual variance by temperature (nlme)",
         c(a_var[k3, "F-value"], a_var[k3, "p-value"]))

# Leave one tank out: refit 8 times, each without one tank; report the
# lowest and highest p and which tank was dropped.
loto <- sapply(sort(unique(sym$tank)), \(tk)
  term_p(lmerTest::lmer(f_sym, data = sym[sym$tank != tk, ], control = ctrl), k3)["p"])
add_sens("symbionts", k3, "leave one tank out: lowest p", c(NA, min(loto)),
         paste("tank", sort(unique(sym$tank))[which.min(loto)], "dropped"))
add_sens("symbionts", k3, "leave one tank out: highest p", c(NA, max(loto)),
         paste("tank", sort(unique(sym$tank))[which.max(loto)], "dropped"))

# Direction check: on Day 15 in heated tanks, are wounded fragments higher or
# lower? Geometric mean = exp(mean(log x)), matching the log-scale model.
gm <- sym |> filter(biopsy_day == 15, treatment == "31C") |> group_by(wound) |>
  summarise(gm = exp(mean(log(cells_per_cm2))), .groups = "drop")
add_sens("symbionts", "Day 15 heated geometric mean, wounded / unwounded",
         "direction of the wound difference", c(NA, NA),
         sprintf("wounded %.0f vs unwounded %.0f cells per cm2",
                 gm$gm[gm$wound == "yes"], gm$gm[gm$wound == "no"]))

## 4c. Growth: temperature ---------------------------------------------------
f_bw <- pct_growth ~ treatment * wound * thicket + (1 | tank)
add_sens("growth", "treatment", "as reported (03_physiology.R)",
         term_p(lmerTest::lmer(f_bw, data = bw, control = ctrl), "treatment"))

# Leave one tank out, then drop tank 3 specifically.
tks <- sort(unique(bw$tank))
loto_bw <- sapply(tks, \(tk)
  term_p(lmerTest::lmer(f_bw, data = bw[bw$tank != tk, ], control = ctrl), "treatment")["p"])
add_sens("growth", "treatment", "leave one tank out: lowest p", c(NA, min(loto_bw)),
         paste("tank", tks[which.min(loto_bw)], "dropped"))
add_sens("growth", "treatment", "leave one tank out: highest p", c(NA, max(loto_bw)),
         paste("tank", tks[which.max(loto_bw)], "dropped"))
add_sens("growth", "treatment", "without ambient tank 3",
         term_p(lmerTest::lmer(f_bw, data = bw[bw$tank != 3, ], control = ctrl), "treatment"))

# Leave one fragment out: refit once per row (one row per fragment) and
# summarise the range of p and how many drops push p above 0.05.
lofo_bw <- sapply(seq_len(nrow(bw)), \(i)
  term_p(lmerTest::lmer(f_bw, data = bw[-i, ], control = ctrl), "treatment")["p"])
add_sens("growth", "treatment", "leave one fragment out: range of p", c(NA, NA),
         sprintf("p from %.3f to %.3f; %d of %d single-fragment drops give p > 0.05",
                 min(lofo_bw), max(lofo_bw), sum(lofo_bw > 0.05), length(lofo_bw)))

# Tank means side by side, to show tank 3 sits with the heated tanks.
tank_means <- bw |>
  group_by(treatment, tank) |>
  summarise(m = mean(pct_growth), .groups = "drop")
add_sens("growth", "tank means (% mass gained)", "ambient tanks vs heated tanks", c(NA, NA),
         paste0("ambient ",
                paste(sprintf("%d: %.1f",
                              tank_means$tank[tank_means$treatment == "28C"],
                              tank_means$m[tank_means$treatment == "28C"]),
                      collapse = ", "),
                "; heated ",
                paste(sprintf("%d: %.1f",
                              tank_means$tank[tank_means$treatment == "31C"],
                              tank_means$m[tank_means$treatment == "31C"]),
                      collapse = ", ")))

write_csv(bind_rows(sens) |> mutate(across(where(is.double), \(x) signif(x, 4))),
          file.path(TBL_DIR, "22b_sensitivity_wound_growth.csv"))


# ---- 5. Color: tank-level check --------------------------------------------
# Color is a 1-6 card score (the D-scale) and ambient fragments almost all
# stayed at 4. So the main model's leftovers are not bell-shaped, as that
# model assumes (18_model_checks.csv), and a model for ordered scores cannot
# be fitted (the ambient group sits almost entirely in one score). This check
# needs no bell curve: compare the final mean color of the 4 ambient and 4
# heated tanks against all 70 ways of splitting 8 tanks into two groups of 4,
# as for growth (03_physiology.R) and budding (04_regeneration.R).
color_final <- readRDS(file.path(DATA_PROC, "color_clean.rds")) |>
  filter(day == max(day, na.rm = TRUE))
tank_color <- color_final |>
  group_by(tank, treatment) |>
  summarise(mean_color = mean(color_num), n_corals = n(),
            n_paled = sum(color_num < 4), .groups = "drop")

# Observed ambient-minus-heated difference, then the same difference for
# every possible split of the 8 tank means (what chance alone would give).
obs_diff <- with(tank_color,
                 mean(mean_color[treatment == "28C"]) - mean(mean_color[treatment == "31C"]))
vals <- tank_color$mean_color
assignments <- combn(seq_along(vals), sum(tank_color$treatment == "28C"), simplify = FALSE)
perm_diffs <- vapply(assignments, function(idx) mean(vals[idx]) - mean(vals[-idx]), numeric(1))

color_tank_test <- tibble(
  response = "color_dscale",
  day = unique(color_final$day),
  test = "tank-level exact permutation",
  n_tanks = nrow(tank_color),
  estimate_28_minus_31 = obs_diff,
  p_two_sided = mean(abs(perm_diffs) >= abs(obs_diff)),
  paled_per_tank_28C = paste(tank_color$n_paled[tank_color$treatment == "28C"],
                             collapse = ", "),
  paled_per_tank_31C = paste(tank_color$n_paled[tank_color$treatment == "31C"],
                             collapse = ", "),
  note = "Paled = card score below 4. 2/70 = 0.029 is the smallest possible two-sided p with 4 tanks per temperature."
) |>
  mutate(across(where(is.double), \(x) signif(x, 4)))
write_csv(color_tank_test, file.path(TBL_DIR, "12b_color_tank_test.csv"))
print(as.data.frame(color_tank_test))
