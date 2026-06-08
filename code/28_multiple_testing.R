# =============================================================================
# Purpose: Multiple-testing sensitivity. The study reports many tests across
#          several response families; primary inference is unadjusted (a-priori
#          directed hypotheses), but a reviewer will ask about family-wise
#          error. This applies Benjamini-Hochberg (FDR) and Bonferroni
#          corrections WITHIN coherent outcome families and reports whether each
#          primary result survives. (Defensive sensitivity, as in the
#          wound-type analyses — not the primary inference.)
# Families:
#   1. Physiology treatment effects (PAM, color, calcification, symbionts):
#      the key treatment / treatment×time interaction per response (script 12).
#   2. Morphology trait treatment effects: treatment term per binary trait
#      (script 04 / 12c).
#   3. Survival: overall Cox HR (31C vs 28C) per healing milestone (script 14).
# Input:   output/tables/12_anova_summary.csv,
#          output/tables/04_morphology_trait_anova_genet.csv,
#          output/tables/14_cox_hazard_ratios.csv
# Output:  output/tables/28_multiple_testing.csv
# =============================================================================

source(here::here("code", "00_setup.R"))

adjust_family <- function(df, family) {
  df |>
    mutate(family = family,
           p_BH = p.adjust(p_value, method = "BH"),
           p_bonferroni = p.adjust(p_value, method = "bonferroni"),
           sig_raw = p_value < 0.05,
           sig_BH = p_BH < 0.05,
           sig_bonferroni = p_bonferroni < 0.05)
}

# ---- Family 1: physiology (treatment×time interaction, or treatment) -------
anova12 <- read_csv(file.path(TBL_DIR, "12_anova_summary.csv"),
                    show_col_types = FALSE) |>
  filter(!grepl("^morph_", response_id))
phys_terms <- c(pam_fvfm = "treatment:day", color_dscale = "treatment:day",
                log_zoox_density = "treatment:biopsy_day_c",
                growth_areal = "treatment")
fam_phys <- imap_dfr(phys_terms, function(term, resp) {
  r <- anova12 |> filter(response_id == resp, term == !!term)
  if (nrow(r) == 0) return(NULL)
  tibble(test = paste0(resp, " : ", term),
         p_value = coalesce(r$`Pr(>F)`[1], r$`Pr(>Chisq)`[1]))
}) |> adjust_family("Physiology (treatment effect)")

# ---- Family 2: morphology trait treatment effects --------------------------
morph <- read_csv(file.path(TBL_DIR, "04_morphology_trait_anova_genet.csv"),
                  show_col_types = FALSE) |>
  filter(term == "treatment")
fam_morph <- tibble(test = paste0("morph:", morph$trait, " : treatment"),
                    p_value = morph$`Pr(>Chisq)`) |>
  adjust_family("Morphology (treatment effect)")

# ---- Family 3: survival (overall Cox HR per milestone) ---------------------
cox <- read_csv(file.path(TBL_DIR, "14_cox_hazard_ratios.csv"),
                show_col_types = FALSE) |>
  filter(grepl("overall", scope), !is.na(p))
fam_cox <- tibble(test = paste0("cox:", cox$trait, " : 31C vs 28C"),
                  p_value = cox$p) |>
  adjust_family("Survival (Cox HR)")

out <- bind_rows(fam_phys, fam_morph, fam_cox) |>
  mutate(across(c(p_value, p_BH, p_bonferroni), \(x) signif(x, 4))) |>
  select(family, test, p_value, p_BH, p_bonferroni,
         sig_raw, sig_BH, sig_bonferroni)
write_csv(out, file.path(TBL_DIR, "28_multiple_testing.csv"))

cat("\n=== Multiple-testing sensitivity (BH-FDR + Bonferroni within family) ===\n")
print(as.data.frame(out))
cat(sprintf("\nSurvive raw: %d | survive BH: %d | survive Bonferroni: %d (of %d tests)\n",
            sum(out$sig_raw, na.rm = TRUE), sum(out$sig_BH, na.rm = TRUE),
            sum(out$sig_bonferroni, na.rm = TRUE), nrow(out)))
cat("Wrote output/tables/28_multiple_testing.csv\n")
