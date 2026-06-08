# =============================================================================
# Purpose: Report morphology (binary trait) treatment contrasts on the
#          INTERPRETABLE probability scale (Δ probability, 31C vs 28C) in
#          addition to the odds-ratio scale, at a representative timepoint.
#          The primary morphology models report log-odds coefficients; readers
#          want "how much lower is the probability of regeneration under heat."
# Input:   output/models/12c_morph_*_blme.rds  (penalized binomial GLMMs)
# Output:  output/tables/29_morphology_prob_contrasts.csv
# =============================================================================

source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages(library(emmeans))

DAY <- 10   # peak divergence in the KM curves

morph_files <- list.files(MOD_DIR, pattern = "^12c_morph_.*_blme\\.rds$",
                          full.names = TRUE)

contrast_one <- function(f) {
  trait <- sub("^12c_morph_(.*)_blme\\.rds$", "\\1", basename(f))
  m <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(m)) return(NULL)

  # Probability scale: Δ probability (28C - 31C), averaged over genet
  prob <- tryCatch({
    emm <- emmeans(m, ~ treatment, at = list(day = DAY), type = "response")
    pr  <- as.data.frame(pairs(emm, reverse = FALSE))   # 28C / 31C
    em  <- as.data.frame(emm)
    p28 <- em$prob[em$treatment == "28C"]
    p31 <- em$prob[em$treatment == "31C"]
    tibble(trait = trait, day = DAY,
           prob_28 = p28, prob_31 = p31,
           delta_prob = p28 - p31)
  }, error = function(e) tibble(trait = trait, day = DAY,
                                prob_28 = NA, prob_31 = NA, delta_prob = NA))

  # Odds-ratio scale (31C vs 28C) with Wald CI + p
  orr <- tryCatch({
    emm_lo <- emmeans(m, ~ treatment, at = list(day = DAY))
    pc <- as.data.frame(pairs(emm_lo, adjust = "tukey"))   # log-odds diff
    # contrast "28C - 31C" -> OR(31 vs 28) = exp(-estimate)
    est <- pc$estimate[1]; se <- pc$SE[1]
    tibble(or_31_vs_28 = exp(-est),
           or_lo = exp(-(est + 1.96 * se)),
           or_hi = exp(-(est - 1.96 * se)),
           p_value = pc$p.value[1])
  }, error = function(e) tibble(or_31_vs_28 = NA, or_lo = NA, or_hi = NA,
                                p_value = NA))

  bind_cols(prob, orr)
}

out <- map_dfr(morph_files, contrast_one) |>
  mutate(across(where(is.numeric), \(x) round(x, 4)))
write_csv(out, file.path(TBL_DIR, "29_morphology_prob_contrasts.csv"))

cat("\n=== Morphology contrasts at Day", DAY,
    "(probability + odds-ratio scale, 31C vs 28C) ===\n")
print(as.data.frame(out))
cat("\ndelta_prob = P(trait | 28C) - P(trait | 31C); positive = heat lowers the",
    "probability. or_31_vs_28 < 1 = heat reduces the odds.\n")
cat("Wrote output/tables/29_morphology_prob_contrasts.csv\n")
