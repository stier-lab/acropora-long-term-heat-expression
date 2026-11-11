# =============================================================================
# Purpose: Numbers for the physiology results (Figure 2): photosynthetic
#          efficiency (Fv/Fm), color, symbiont density and growth.
# Inputs:  output/tables/12_anova_summary.csv       (analysis/12_models.R)
#          output/tables/45_figure2_summary.csv     (analysis/45_figure2_physiology.R)
#          output/tables/03_color_end_proportions.csv (analysis/03_color.R)
#          output/tables/06_symbiont_chl_summary.csv  (analysis/06_symbionts.R)
#          output/tables/05_buoyant_weight_tank_test.csv (analysis/05_growth.R)
# Output:  output/paper/physiology_numbers.csv
# Methods (models fitted in analysis/12_models.R):
#   Fv/Fm and color: temperature x wound x day x source patch, day as a linear
#     trend, random intercepts for tank and fragment (lme4/lmerTest,
#     Satterthwaite df).
#   Symbionts: log(cells per cm2), same fixed effects, random intercept for
#     tank only (each fragment sampled once).
#   Growth: % change in buoyant mass Day 0-15, temperature x wound x source
#     patch, random intercept for tank. Reported test = this model
#     (decision 2026-11-10); the tank permutation test is a supplement check.
# =============================================================================
source(here::here("code", "00_setup.R"))

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
write_paper(rows, "physiology")
