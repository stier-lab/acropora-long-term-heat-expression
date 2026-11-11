# =============================================================================
# Purpose: Numbers for the healing and regeneration results (Figure 1).
# Inputs:  output/tables/11f_morphology_dataset_endpoint_tests.csv (analysis/11c_stage_counts.R)
#          output/tables/11f_morphology_dataset_logrank_tests.csv
#          output/tables/14_interval_survreg.csv, 14_cox_hazard_ratios.csv,
#            14_cox_genet_LRT.csv                    (analysis/14_stage_timing_models.R)
#          output/tables/36_stage_timing_mean_iqr.csv (analysis/36_stage_timing_tables.R)
#          data/processed/physio_clean.rds, microscope_physio_clean.rds
# Output:  output/paper/regeneration_numbers.csv
# Methods: 24 wounded fragments (12 per temperature) scored daily, Days 0-15,
#   for five stages. Stage counts compared with Fisher's exact test. Time to
#   each stage: interval-censored Weibull models (onset between the last day
#   without and the first day with the stage; not-yet-reached fragments kept),
#   with source patch; Cox models and log-rank tests as checks. Because
#   temperature was applied to tanks (3 wounded fragments per tank), budding is
#   also summarised tank by tank, with an exact permutation test of 4 vs 4 tanks.
#   Tissue cover comes from the separate 16-fragment microscope experiment.
# Do not report: p-values from the daily yes/no models in analysis/12_models.R
#   (they change with computer and package versions; see
#   docs/provenance/rerun_check_2026-11-10.md).
# =============================================================================
source(here::here("code", "00_setup.R"))

endpoint <- readr::read_csv(file.path(TBL_DIR, "11f_morphology_dataset_endpoint_tests.csv"), show_col_types = FALSE)
logrank  <- readr::read_csv(file.path(TBL_DIR, "11f_morphology_dataset_logrank_tests.csv"), show_col_types = FALSE)
aft      <- readr::read_csv(file.path(TBL_DIR, "14_interval_survreg.csv"), show_col_types = FALSE)
cox      <- readr::read_csv(file.path(TBL_DIR, "14_cox_hazard_ratios.csv"), show_col_types = FALSE)
cox_gx   <- readr::read_csv(file.path(TBL_DIR, "14_cox_genet_LRT.csv"), show_col_types = FALSE)
stages   <- readr::read_csv(file.path(TBL_DIR, "36_stage_timing_mean_iqr.csv"), show_col_types = FALSE)

traits <- c(axial_polyp_formation = "Axial polyp forms",
            wound_smoothed = "Cut surface smooths over",
            tip_exist = "New tip forms",
            tip_extension = "Tip extends",
            new_corallites_on_tip = "New radial corallites bud")
rows <- list()

# Stage present on the last scored day (Fisher's exact test), main experiment
for (tr in names(traits)) {
  d <- endpoint[endpoint$dataset == "physio_morphology" & endpoint$trait == tr, ]
  stopifnot(nrow(d) == 1)
  rows[[length(rows) + 1]] <- paper_row(
    paste0(tr, "_present_day15"),
    paste0(traits[[tr]], ": fragments with the stage on Day 15 (estimate = ambient, statistic = heated; n = per temperature)"),
    estimate = d$yes_28, statistic = d$yes_31, n = d$n_28,
    test = "Fisher's exact test", p_value = d$fisher_p,
    source = "output/tables/11f_morphology_dataset_endpoint_tests.csv")
}

# Stage reached at any point by Day 15, main experiment
stage_trait <- c(`2` = "axial_polyp_formation", `3` = "wound_smoothed", `4` = "tip_exist",
                 `5` = "tip_extension", `6` = "new_corallites_on_tip")
for (so in names(stage_trait)) {
  d28 <- stages[stages$stage_order == as.integer(so) & grepl("^28", stages$treatment), ]
  d31 <- stages[stages$stage_order == as.integer(so) & grepl("^31", stages$treatment), ]
  if (nrow(d28) != 1 || nrow(d31) != 1) next
  tr <- stage_trait[[so]]
  rows[[length(rows) + 1]] <- paper_row(
    paste0(tr, "_reached_by_day15"),
    paste0(traits[[tr]], ": fragments that reached the stage by Day 15 (estimate = ambient, statistic = heated)"),
    estimate = d28$n_reached, statistic = d31$n_reached, n = d28$n_scored,
    source = "output/tables/36_stage_timing_mean_iqr.csv")
}

d <- logrank[logrank$dataset == "physio_morphology" & logrank$trait == "new_corallites_on_tip", ]
rows[[length(rows) + 1]] <- paper_row(
  "budding_logrank", "Budding: log-rank test", test = "log-rank",
  statistic = d$logrank_chisq, df1 = d$logrank_df, p_value = d$logrank_p,
  source = "output/tables/11f_morphology_dataset_logrank_tests.csv")

# Time to each stage (interval-censored Weibull; time ratio > 1 = slower under heat)
for (tr in names(traits)) {
  d <- aft[aft$trait == tr, ]
  if (nrow(d) != 1) next
  rows[[length(rows) + 1]] <- paper_row(
    paste0(tr, "_time_ratio"), paste0(traits[[tr]], ": time ratio, heated / ambient"),
    estimate = d$time_ratio_31_vs28, ci_low = d$ratio_lo, ci_high = d$ratio_hi,
    test = "interval-censored Weibull model", statistic = d$z, p_value = d$p, n = d$n,
    source = "output/tables/14_interval_survreg.csv")
}
d <- cox[cox$trait == "new_corallites_on_tip" & grepl("^overall", cox$scope), ]
stopifnot(nrow(d) == 1)
rows[[length(rows) + 1]] <- paper_row(
  "budding_hazard_ratio", "Budding: hazard ratio, heated / ambient (Cox)",
  estimate = d$HR_31_vs28, ci_low = d$HR_lo, ci_high = d$HR_hi,
  test = "Cox model, stratified by source patch", statistic = d$z, p_value = d$p, n = d$n,
  source = "output/tables/14_cox_hazard_ratios.csv")
d <- cox_gx[cox_gx$trait == "new_corallites_on_tip", ]
rows[[length(rows) + 1]] <- paper_row(
  "budding_patch_x_temp", "Budding: source patch x temperature (exploratory; 4 fragments per patch and temperature)",
  test = "likelihood-ratio test", statistic = d$chisq_genet_trt_int, df1 = d$df_genet_trt_int,
  p_value = d$p_genet_trt_int, source = "output/tables/14_cox_genet_LRT.csv")

# Budding tank by tank (temperature was applied to tanks)
ph <- readRDS(file.path(DATA_PROC, "physio_clean.rds"))
by_tank <- ph |>
  dplyr::filter(wound == "yes") |>
  dplyr::group_by(treatment, tank, id) |>
  dplyr::summarise(budded = any(new_corallites_on_tip %in% c("yes", 1, TRUE)), .groups = "drop") |>
  dplyr::group_by(treatment, tank) |>
  dplyr::summarise(n = dplyr::n(), k = sum(budded), prop = k / n, .groups = "drop")
props <- setNames(by_tank$prop, by_tank$tank)
amb <- as.character(by_tank$tank[by_tank$treatment == "28C"])
diff_for <- function(set) mean(props[set]) - mean(props[setdiff(names(props), set)])
all_sets <- utils::combn(names(props), length(amb), simplify = FALSE)
obs <- diff_for(amb)
perm_p <- mean(vapply(all_sets, function(s) abs(diff_for(s)) >= abs(obs) - 1e-12, logical(1)))
for (tr in c("28C", "31C")) {
  lab <- if (tr == "28C") "ambient" else "heated"
  d <- by_tank[by_tank$treatment == tr, ]
  rows[[length(rows) + 1]] <- paper_row(
    paste0("budding_per_tank_", lab),
    paste0("Budding by Day 15, ", lab, " tanks: wounded fragments that budded per tank (estimate = fewest, statistic = most)"),
    estimate = min(d$k), statistic = max(d$k), n = d$n[1], df1 = nrow(d),
    source = "data/processed/physio_clean.rds")
}
rows[[length(rows) + 1]] <- paper_row(
  "budding_tank_permutation", "Budding: exact permutation test of tank proportions, 4 vs 4 tanks",
  estimate = obs, test = "exact permutation (70 assignments); smallest possible p = 2/70",
  p_value = perm_p, n = length(props), source = "data/processed/physio_clean.rds")

# Tissue cover, separate microscope experiment
mic <- readRDS(file.path(DATA_PROC, "microscope_physio_clean.rds"))
first_cover <- mic |>
  dplyr::filter(tissue_over_wound %in% c("yes", 1, TRUE)) |>
  dplyr::group_by(id) |>
  dplyr::summarise(first_day = min(day), .groups = "drop")
rows[[length(rows) + 1]] <- paper_row(
  "tissue_cover_by_day1", "Microscope experiment: fragments with tissue covering the cut by Day 1",
  estimate = sum(first_cover$first_day <= 1), n = dplyr::n_distinct(mic$id),
  source = "data/processed/microscope_physio_clean.rds")

write_paper(rows, "regeneration")
