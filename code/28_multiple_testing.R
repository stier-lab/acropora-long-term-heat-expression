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
# incidental traits (hole_in_center, polyp_in_hole, wound_smoothed, polyps_out,
# pigment_over_wound, algae_on_wound) — no strong a-priori prediction.
#
# Input:   output/tables/12_anova_summary.csv,
#          output/tables/04_morphology_trait_anova_genet.csv,
#          output/tables/14_cox_hazard_ratios.csv
# Output:  output/tables/28_multiple_testing.csv
# =============================================================================

source(here::here("code", "00_setup.R"))

# Traits whose heat response was predicted a priori (skeletal regeneration)
APRIORI_REGEN <- c("tip_exist", "tip_extension", "new_corallites_on_tip")

rows <- list()
add <- function(test, p_value, hypothesis, rationale) {
  rows[[length(rows) + 1]] <<- tibble(test, p_value, hypothesis, rationale)
}

# ---- Physiology: all four are a-priori directed hypotheses ----------------
anova12 <- read_csv(file.path(TBL_DIR, "12_anova_summary.csv"),
                    show_col_types = FALSE) |>
  filter(!grepl("^morph_", response_id))
phys <- list(
  c(resp = "pam_fvfm",         term = "treatment:day",
    rat = "Warner et al. 1999 (heat lowers Fv/Fm)"),
  c(resp = "color_dscale",     term = "treatment:day",
    rat = "Hoegh-Guldberg 1999 (heat-driven paling)"),
  c(resp = "log_zoox_density", term = "treatment:biopsy_day_c",
    rat = "Hoegh-Guldberg 1999; Jokiel & Coles 1990 (symbiont loss)"),
  c(resp = "growth_areal",     term = "treatment",
    rat = "Jokiel & Coles 1977 (heat reduces calcification)")
)
for (h in phys) {
  r <- anova12 |> filter(response_id == h["resp"], term == h["term"])
  if (nrow(r) == 0) next
  add(paste0(h["resp"], " : ", h["term"]),
      coalesce(r$`Pr(>F)`[1], r$`Pr(>Chisq)`[1]),
      "a priori (confirmatory)", h["rat"])
}

# ---- Morphology trait treatment effects -----------------------------------
morph <- read_csv(file.path(TBL_DIR, "04_morphology_trait_anova_genet.csv"),
                  show_col_types = FALSE) |>
  filter(term == "treatment")
for (i in seq_len(nrow(morph))) {
  tr <- morph$trait[i]
  is_ap <- tr %in% APRIORI_REGEN
  add(paste0("morph:", tr, " : treatment"), morph$`Pr(>Chisq)`[i],
      if (is_ap) "a priori (confirmatory)" else "exploratory",
      if (is_ap) "central hypothesis: heat impairs tip regeneration"
      else "no a-priori prediction for this closure/incidental trait")
}

# ---- Survival: overall Cox HR per milestone -------------------------------
cox <- read_csv(file.path(TBL_DIR, "14_cox_hazard_ratios.csv"),
                show_col_types = FALSE) |>
  filter(grepl("overall", scope), !is.na(p))
for (i in seq_len(nrow(cox))) {
  tr <- cox$trait[i]
  is_ap <- tr %in% APRIORI_REGEN
  add(paste0("cox:", tr, " : 31C vs 28C"), cox$p[i],
      if (is_ap) "a priori (confirmatory)" else "exploratory",
      if (is_ap) "central hypothesis: heat delays/prevents tip regeneration"
      else "no a-priori prediction for this closure milestone")
}

# ---- BH correction within the EXPLORATORY family only ---------------------
out <- bind_rows(rows)
out$p_BH <- NA_real_
expl <- out$hypothesis == "exploratory"
out$p_BH[expl] <- p.adjust(out$p_value[expl], method = "BH")
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
