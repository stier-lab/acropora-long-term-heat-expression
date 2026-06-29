# =============================================================================
# Purpose: Confirmatory vs exploratory testing. A subset of the tests are
#          a-priori DIRECTED hypotheses grounded in the coral thermal-stress
#          literature; these are confirmatory and reported UNADJUSTED. The
#          remaining tests are exploratory (no strong prior prediction) and are
#          corrected for multiple comparisons with Benjamini-Hochberg (FDR)
#          within the exploratory family.
#
# A-priori directed hypotheses (reported raw; confirmatory):
#   * Sustained heat lowers photochemical efficiency over time
#     (Fv/Fm treatment x day) — Warner, Fitt & Schmidt 1999.
#   * Heat drives loss of pigmentation and symbionts over time
#     (color & symbiont treatment x time) — Hoegh-Guldberg 1999; Jokiel & Coles 1990.
#   * Heat reduces calcification (treatment) — Jokiel & Coles 1977.
#   * Heat impairs skeletal REGENERATION at the wound tip
#     (tip_exist, tip_extension, new_corallites_on_tip) — the study's central
#     directed hypothesis (energetic trade-off: heat diverts energy from
#     non-essential regrowth).
#
# Exploratory (BH-corrected): whether heat alters basic wound CLOSURE or other
# incidental traits (axial_polyp_formation, wound_smoothed, polyps_out,
# pigment_over_wound, algae_on_wound) — no strong a-priori prediction. (Note:
# hole_in_center + polyp_in_hole are combined upstream in code/04 into the single
# trait axial_polyp_formation, so each enters this BH family only once.)
#
# What & why: running many tests inflates the false-positive rate, but blanket
#   multiple-testing correction is too blunt — it penalizes a directed hypothesis
#   as harshly as an exploratory test with no prior prediction. Standard practice
#   (and what reviewers expect) is to split the tests into two families.
#   CONFIRMATORY tests are directed, pre-registered-in-spirit predictions backed
#   by specific prior literature (e.g. "heat lowers Fv/Fm", Warner et al. 1999);
#   each is a single planned comparison, so we report its RAW p-value. EXPLORATORY
#   tests have no strong prior; we treat them as a
#   family and apply Benjamini-Hochberg (BH) false-discovery-rate correction so
#   that, across that family, the expected proportion of false positives is
#   controlled. This script gathers p-values from the physiology, morphology, and
#   survival analyses, tags each as confirmatory or exploratory, BH-adjusts the
#   exploratory family only, and reports the appropriate p for each.
# Input:   output/tables/12_anova_summary.csv,
#          output/tables/12c_morph_blme_anova.csv,
#          output/tables/14_interval_survreg.csv
# Output:  output/tables/28_multiple_testing.csv
# =============================================================================

# 00_setup.R loads packages and shared paths (TBL_DIR, ...).
source(here::here("code", "00_setup.R"))

# The three skeletal-regeneration traits — the study's central directed
# hypothesis, so these are CONFIRMATORY wherever they appear below.
APRIORI_REGEN <- c("tip_exist", "tip_extension", "new_corallites_on_tip")

# Collect each test into a growing list; `add()` appends one tidy row. The `<<-`
# writes to the `rows` list defined in this (parent) scope, not a local copy.
rows <- list()
add <- function(test, p_value, hypothesis, rationale) {
  rows[[length(rows) + 1]] <<- tibble(test, p_value, hypothesis, rationale)
}

# ---- Physiology: all four are a-priori directed hypotheses ----------------
# Each physiology test maps to a specific prior prediction (the `rat` rationale),
# so all four are confirmatory and reported with raw p-values.
anova12 <- read_csv(file.path(TBL_DIR, "12_anova_summary.csv"),
                    show_col_types = FALSE) |>
  filter(!grepl("^morph_", response_id))
growth_tank <- read_csv(file.path(TBL_DIR, "05_buoyant_weight_tank_test.csv"),
                        show_col_types = FALSE)
phys <- list(
  c(resp = "pam_fvfm",         term = "treatment:day",
    rat = "Warner et al. 1999 (heat lowers Fv/Fm)"),
  c(resp = "color_dscale",     term = "treatment:day",
    rat = "Hoegh-Guldberg 1999 (heat-driven paling)"),
  c(resp = "log_zoox_density", term = "treatment:biopsy_day_c",
    rat = "Hoegh-Guldberg 1999; Jokiel & Coles 1990 (symbiont loss)"),
  c(resp = "growth_pct",     term = "treatment",
    rat = "Jokiel & Coles 1977 (heat reduces calcification); tank-level permutation p")
)
for (h in phys) {
  # Growth uses a tank-level permutation p (the only tank-randomized response),
  # so pull it from its own file rather than the model ANOVA table.
  if (h["resp"] == "growth_pct" && nrow(growth_tank) > 0) {
    add(paste0(h["resp"], " : ", h["term"], " (tank permutation)"),
        growth_tank$p_two_sided[1], "a priori (confirmatory)", h["rat"])
    next
  }
  # All other physiology responses: grab the row matching this response/term.
  r <- anova12 |> filter(response_id == h["resp"], term == h["term"])
  if (nrow(r) == 0) next                       # skip if the term wasn't fit
  # coalesce(): LMMs give an F-test p, GLMMs a Chisq p — take whichever exists.
  add(paste0(h["resp"], " : ", h["term"]),
      coalesce(r$`Pr(>F)`[1], r$`Pr(>Chisq)`[1]),
      "a priori (confirmatory)", h["rat"])
}

# ---- Morphology trait treatment effects -----------------------------------
# Each morphology trait's treatment effect: confirmatory if it's a regeneration
# trait (central hypothesis), otherwise exploratory (incidental closure trait).
morph <- read_csv(file.path(TBL_DIR, "12c_morph_blme_anova.csv"),
                  show_col_types = FALSE) |>
  filter(term == "treatment")
for (i in seq_len(nrow(morph))) {
  tr <- morph$trait[i]
  is_ap <- tr %in% APRIORI_REGEN               # is this a directed-hypothesis trait?
  add(paste0("morph:", tr, " : treatment"), morph$`Pr(>Chisq)`[i],
      if (is_ap) "a priori (confirmatory)" else "exploratory",
      if (is_ap) "central hypothesis: heat impairs tip regeneration"
      else "no a-priori prediction for this closure/incidental trait")
}

# ---- Survival: interval-censored timing per milestone ---------------------
# Time-to-event (when each milestone is reached) compared between 31C and 28C;
# same confirmatory/exploratory split as the morphology traits above.
surv <- read_csv(file.path(TBL_DIR, "14_interval_survreg.csv"),
                 show_col_types = FALSE) |>
  filter(!is.na(p))                            # keep only milestones with a fitted p
for (i in seq_len(nrow(surv))) {
  tr <- surv$trait[i]
  is_ap <- tr %in% APRIORI_REGEN
  add(paste0("survival_interval:", tr, " : 31C vs 28C"), surv$p[i],
      if (is_ap) "a priori (confirmatory)" else "exploratory",
      if (is_ap) "central hypothesis: heat delays/prevents tip regeneration"
      else "no a-priori prediction for this closure milestone")
}

# ---- BH correction within the EXPLORATORY family only ---------------------
# Design choice: BH is applied only across the exploratory tests, so the
# confirmatory tests keep their raw p and are not penalized for the exploratory ones.
out <- bind_rows(rows)                          # stack all collected test rows
out$p_BH <- NA_real_                            # confirmatory rows stay NA here
expl <- out$hypothesis == "exploratory"         # logical mask for the exploratory family
out$p_BH[expl] <- p.adjust(out$p_value[expl], method = "BH")   # FDR within that family
out <- out |>
  mutate(
    # a-priori tests are confirmatory -> use raw p; exploratory -> BH-adjusted
    p_reported = ifelse(hypothesis == "exploratory", p_BH, p_value),
    significant = p_reported < 0.05,
    across(c(p_value, p_BH, p_reported), \(x) signif(x, 4))
  ) |>
  select(hypothesis, test, rationale, p_value, p_BH, p_reported, significant)
write_csv(out, file.path(TBL_DIR, "28_multiple_testing.csv"))

n_ap  <- sum(out$hypothesis != "exploratory")
n_exp <- sum(expl)
cat("\n=== Confirmatory (a-priori, raw p) vs exploratory (BH-corrected) ===\n")
print(as.data.frame(out))
cat(sprintf("\nA-priori: %d tests (raw p), %d significant.\n",
            n_ap, sum(out$significant[out$hypothesis != "exploratory"], na.rm = TRUE)))
cat(sprintf("Exploratory: %d tests (BH-corrected), %d significant.\n",
            n_exp, sum(out$significant[expl], na.rm = TRUE)))
cat("Wrote output/tables/28_multiple_testing.csv\n")
