# Verify the field-confirmed wound correction and compare physiology tests.
source(here::here("code", "00_setup.R"))

cohorts <- c("coral_metadata", "pam_clean", "color_clean",
             "buoyant_weight_clean", "symbiont_chl_clean", "physio_clean")
checks <- purrr::map_dfr(cohorts, function(name) {
  d <- readRDS(file.path(DATA_PROC, paste0(name, ".rds")))
  selected <- d |> filter(id == 119L)
  stopifnot(nrow(selected) > 0, all(as.character(selected$wound) == "yes"))
  tibble(dataset = name, n_rows_119 = nrow(selected), wound = "yes", status = "pass")
})
write_csv(checks, file.path(TBL_DIR, "46_wound_label_checks.csv"))

growth <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds")) |>
  filter(biopsy_day == 15, is.finite(pct_growth))
stopifnot(nrow(growth) == 48L, !anyDuplicated(growth$id),
          all(count(growth, treatment, wound)$n == 12L))
stopifnot(identical(correct_wound_labels(c(118, 119, 120), c("no", "no", "yes")),
                    c("no", "yes", "yes")))

old <- read_csv(here::here("data", "metadata", "119_before_correction_anova_2026-10-16.csv"),
                show_col_types = FALSE)
current <- read_csv(file.path(TBL_DIR, "12_anova_summary.csv"), show_col_types = FALSE)
select_tests <- function(d) d |>
  filter(response_id %in% c("pam_fvfm", "color_dscale", "growth_pct", "log_zoox_density"),
         term %in% c("treatment:wound", "treatment:wound:day",
                     "treatment:wound:biopsy_day_c")) |>
  select(response_id, term, p_value = `Pr(>F)`)
comparison <- inner_join(select_tests(old), select_tests(current),
                         by = c("response_id", "term"), suffix = c("_before", "_after"))
stopifnot(nrow(comparison) == 7L, all(is.finite(comparison$p_value_after)))
write_csv(comparison, file.path(TBL_DIR, "46_wound_label_model_comparison.csv"))
print(comparison)
