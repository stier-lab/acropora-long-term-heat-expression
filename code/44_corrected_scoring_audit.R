# Compare the October 16 Drive snapshots with the preserved earlier exports.
source(here::here("code", "00_setup.R"))
cohorts <- c("physio_morphology", "microscope_physio")
changes <- map_dfr(cohorts, function(cohort) {
  read_snapshot <- function(name) read_csv(file.path(DATA_RAW, cohort, name),
    col_types = cols(.default = col_character()), na = character()) |>
    janitor::clean_names() |> filter(id != "")
  old <- read_snapshot("data_before_2026-10-16.csv")
  new <- read_snapshot("data.csv")
  stopifnot(!anyDuplicated(old[c("id", "day")]),
            !anyDuplicated(new[c("id", "day")]),
            setequal(paste(old$id, old$day), paste(new$id, new$day)))
  a <- old |> pivot_longer(-c(id, day), names_to = "column", values_to = "previous")
  b <- new |> pivot_longer(-c(id, day), names_to = "column", values_to = "corrected")
  full_join(a, b, by = c("id", "day", "column")) |>
    mutate(previous = replace_na(previous, ""), corrected = replace_na(corrected, "")) |>
    filter(previous != corrected) |> mutate(cohort = cohort, .before = 1)
})
write_csv(changes, file.path(TBL_DIR, "44_scoring_changes.csv"))
stopifnot(!any(changes$column %in% c("id", "day", "tank", "treatment", "wounded", "wound", "thicket")))
main <- readRDS(file.path(DATA_PROC, "physio_clean.rds")) |> filter(wound == "yes")
stopifnot(nrow(main) == 384, n_distinct(main$id) == 24,
          all(is.na(main$pigment_over_wound[!main$day %in% c(3, 4, 9, 12, 15)])),
          all(!is.na(main$pigment_over_wound[main$day %in% c(3, 4, 9, 12, 15)])))
traits <- c("tip_exist", "tip_extension", "new_corallites_on_tip", "pigment_over_wound")
coverage <- main |> select(id, day, treatment, all_of(traits)) |>
  pivot_longer(all_of(traits), names_to = "trait", values_to = "value") |>
  group_by(trait, treatment, day) |>
  summarise(n = n(), n_scored = sum(!is.na(value)), n_yes = sum(value == 1, na.rm = TRUE),
            .groups = "drop")
write_csv(coverage, file.path(TBL_DIR, "44_corrected_score_coverage.csv"))
events <- main |> select(id, day, treatment, all_of(traits)) |>
  pivot_longer(all_of(traits), names_to = "trait", values_to = "value") |>
  filter(!is.na(value)) |> group_by(trait, treatment, id) |>
  summarise(reached = any(value == 1),
            first_day = if (any(value == 1)) min(day[value == 1]) else NA_real_,
            last_scored = max(day), .groups = "drop")
write_csv(events, file.path(TBL_DIR, "44_corrected_first_observations.csv"))
print(count(changes, cohort, column))
print(events |> group_by(trait, treatment) |> summarise(n = n(), reached = sum(reached), .groups = "drop"))

# Test the production censoring function without refitting its models.
definitions <- parse(file.path("code", "14_morphology_kaplan.R"))
is_event_function <- vapply(definitions, function(x) is.call(x) &&
  identical(x[[1]], as.name("<-")) && identical(x[[2]], as.name("compute_events")), logical(1))
stopifnot(sum(is_event_function) == 1)
test_env <- new.env(parent = environment())
eval(definitions[[which(is_event_function)]], envir = test_env)
fixture <- tibble(id = c(1L, 1L, 1L, 2L, 2L, 2L), treatment = "28C",
                  tank = 1L, thicket = "a", wound = "yes",
                  day = c(3, 8, 9, 3, 12, 15), score = c(0, NA, 1, 0, 0, NA))
check <- test_env$compute_events(fixture, "score")
stopifnot(check$event_lower[check$id == 1] == 3,
          check$event_upper[check$id == 1] == 9,
          check$event_lower[check$id == 2] == 12,
          is.infinite(check$event_upper[check$id == 2]))
cat("Onset regression checks passed: missing visits neither narrow intervals nor extend follow-up.\n")
