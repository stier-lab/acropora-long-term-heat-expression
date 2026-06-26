# =============================================================================
# Purpose: Report morphology (binary trait) treatment contrasts on the
#          INTERPRETABLE probability scale (Δ probability, 31C vs 28C) in
#          addition to the odds-ratio scale, at a representative timepoint.
#          The primary morphology models report log-odds coefficients; readers
#          want "how much lower is the probability of regeneration under heat."
#
# What & why: the morphology GLMMs are fit on the logit (log-odds) scale, where a
#   coefficient like -1.4 is mathematically clean but biologically opaque — no
#   reader thinks in log-odds. This script translates each treatment effect into
#   two scales a reader can actually picture, evaluated at day 10 (where the
#   Kaplan-Meier curves diverge most): (1) the Δ-probability — literally "the
#   probability of the trait under ambient minus under heat" (e.g. 0.30 = heat
#   lowers the chance by 30 percentage points); and (2) the odds ratio (31C vs
#   28C), where <1 means heat reduces the odds. We use emmeans to get model-based
#   marginal means averaged over genet. One numerical hazard: with binary data a
#   trait can be perfectly "separated" (probability pinned at 0 or 1 in a group),
#   which makes the log-odds — and thus the OR and its CI — blow up to absurd
#   values (~1e19). The bounded Δ-probability is still meaningful in that case,
#   so we keep it but null out the garbage OR/CI/p (see the separation guard).
# Input:   output/models/12c_morph_*_blme.rds  (penalized binomial GLMMs)
# Output:  output/tables/29_morphology_prob_contrasts.csv
# =============================================================================

# 00_setup.R loads packages and shared paths (MOD_DIR, TBL_DIR, ...).
source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages(library(emmeans))   # marginal means / contrasts

DAY <- 10   # peak divergence in the KM curves — the representative timepoint

# Glob every saved penalized-binomial morphology model (one per trait).
morph_files <- list.files(MOD_DIR, pattern = "^12c_morph_.*_blme\\.rds$",
                          full.names = TRUE)

# Helper: produce one row of probability + odds-ratio contrasts for one model.
contrast_one <- function(f) {
  trait <- sub("^12c_morph_(.*)_blme\\.rds$", "\\1", basename(f))   # trait from filename
  m <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(m)) return(NULL)                     # skip if the model is missing

  # Probability scale: Δ probability (28C - 31C), averaged over genet.
  # type = "response" back-transforms the marginal means from logit to probability.
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

  # Odds-ratio scale (31C vs 28C) with Wald CI + p. Here emmeans is left on the
  # default LINK (log-odds) scale, so the pairwise difference is a log-odds diff.
  orr <- tryCatch({
    emm_lo <- emmeans(m, ~ treatment, at = list(day = DAY))
    pc <- as.data.frame(pairs(emm_lo, adjust = "tukey"))   # log-odds diff (28C - 31C)
    # emmeans reports "28C - 31C"; flip the sign so the OR is heat-vs-ambient:
    # OR(31 vs 28) = exp(-estimate). CI from the Wald ±1.96*SE on the log scale.
    est <- pc$estimate[1]; se <- pc$SE[1]
    tibble(or_31_vs_28 = exp(-est),
           or_lo = exp(-(est + 1.96 * se)),
           or_hi = exp(-(est - 1.96 * se)),
           p_value = pc$p.value[1])
  }, error = function(e) tibble(or_31_vs_28 = NA, or_lo = NA, or_hi = NA,
                                p_value = NA))

  bind_cols(prob, orr) |>
    # Separation guard: when either probability sits at the {0,1} boundary, the
    # log-odds contrast (and therefore the odds ratio and its CI) is undefined
    # and explodes to ~1e19. Keep the bounded Δ-probability; null out the
    # meaningless OR/CI/p so these artifacts never reach a table or the paper.
    mutate(
      separated = coalesce(prob_28 <= 1e-6 | prob_28 >= 1 - 1e-6 |
                           prob_31 <= 1e-6 | prob_31 >= 1 - 1e-6, FALSE),
      across(c(or_31_vs_28, or_lo, or_hi, p_value),
             \(x) if_else(separated, NA_real_, x))
    ) |>
    select(-separated)
}

# Run the helper over every trait model and stack the one-row results.
out <- map_dfr(morph_files, contrast_one) |>
  mutate(across(where(is.numeric), \(x) round(x, 4)))
write_csv(out, file.path(TBL_DIR, "29_morphology_prob_contrasts.csv"))

cat("\n=== Morphology contrasts at Day", DAY,
    "(probability + odds-ratio scale, 31C vs 28C) ===\n")
print(as.data.frame(out))
cat("\ndelta_prob = P(trait | 28C) - P(trait | 31C); positive = heat lowers the",
    "probability. or_31_vs_28 < 1 = heat reduces the odds.\n")
cat("Wrote output/tables/29_morphology_prob_contrasts.csv\n")
