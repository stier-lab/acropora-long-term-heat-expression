# =============================================================================
# Purpose: Sensitivity analysis for the samples/tanks Molly flagged as odd
#          (see notes/QAQC_flagged_samples.md). Re-runs the key treatment-effect
#          tests with the flagged colonies (116, 121) and the slow ambient tank
#          (tank 3) removed, and compares the treatment effect to the full-data
#          model. If the conclusions hold, the flags are documented but
#          inconsequential.
#
# What & why: during the experiment Molly noticed a few suspect data points —
#   two coral colonies (IDs 116 and 121) that behaved oddly, and one ambient
#   tank (tank 3) that ran noticeably colder/slower than its siblings. A fair
#   worry is that the headline results are driven by these few oddballs rather
#   than by the real heat treatment. This script answers that worry the honest
#   way: re-fit each main model TWICE — once on the full data and once with the
#   flagged samples removed — and put the two treatment effects side by side. If
#   the effect (its size, sign, and significance) barely moves when the flagged
#   points are dropped, the conclusions are robust and the flags are just
#   footnotes. This is a sensitivity analysis, not a new hypothesis test.
# Input:   data/processed/{pam_clean,color_clean,buoyant_weight_clean}.rds
# Output:  output/tables/22_sensitivity_flagged.csv
# =============================================================================

# 00_setup.R loads packages and defines shared paths (DATA_PROC, TBL_DIR, ...).
source(here::here("code", "00_setup.R"))

# ---- The flagged samples ---------------------------------------------------
# The exact colonies and tank Molly flagged in QA/QC. Defined once here so the
# same removal rule (drop these IDs, drop this tank) is applied to every response.
FLAGGED_IDS  <- c(116, 121)   # two coral colonies that behaved anomalously
FLAGGED_TANK <- 3             # ambient tank that warmed/equilibrated slowly

results <- list()             # one tibble per (response × full-or-dropped) goes here

# ---- Helper: pull one ANOVA row out of a fitted model ----------------------
# For the mixed models below, lmerTest's anova() gives a Type-III F-test table.
# This grabs the single row we care about (e.g. the "treatment:day" interaction)
# and returns its F statistic and p-value as a tidy one-row tibble. Returning
# NULL when the term is absent lets bind_rows() simply skip it.
grab <- function(model, label, term) {
  av <- as.data.frame(anova(model))
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

# ---- PAM Fv/Fm: treatment × day --------------------------------------------
# Photosynthetic efficiency, measured repeatedly on each coral over time. The
# key effect is the treatment×day interaction: does the heated group's Fv/Fm
# decline differ from the ambient group's over the experiment? This MIRRORS the
# primary model (02/12): genet (thicket) is a FIXED blocking term (only 3 genets),
# random intercepts for tank and id absorb the nested repeated-measures structure,
# and the pre-treatment baseline (day < 0) is excluded so the comparison is
# apples-to-apples with the reported model. REML = FALSE (ML) is used so the same
# model fits the full vs dropped data comparably. pam_drop is the data minus flags.
pam <- readRDS(file.path(DATA_PROC, "pam_clean.rds")) |>
  mutate(thicket = factor(thicket)) |>
  filter(day >= 0)
pam_drop <- pam |> filter(!id %in% FLAGGED_IDS, tank != FLAGGED_TANK)

f_full <- lmerTest::lmer(fv_fm ~ treatment * wound * day + thicket + (1|tank) + (1|id),
                         data = pam, REML = FALSE)
f_drop <- lmerTest::lmer(fv_fm ~ treatment * wound * day + thicket + (1|tank) + (1|id),
                         data = pam_drop, REML = FALSE)
# Same model, two datasets -> compare the treatment:day F/p side by side.
results$pam_full <- grab(f_full, "PAM Fv/Fm (full data)",     "treatment:day")
results$pam_drop <- grab(f_drop, "PAM Fv/Fm (flagged dropped)","treatment:day")

# ---- Color D-scale: treatment × day ----------------------------------------
# Coral colour scored against the standard Coral Watch D-scale (a bleaching/
# paling proxy), again repeated over time. Identical model structure and
# full-vs-dropped logic as PAM above.
col <- readRDS(file.path(DATA_PROC, "color_clean.rds")) |>
  mutate(thicket = factor(thicket)) |>
  filter(day >= 0)
col_drop <- col |> filter(!id %in% FLAGGED_IDS, tank != FLAGGED_TANK)

c_full <- lmerTest::lmer(color_num ~ treatment * wound * day + thicket + (1|tank) + (1|id),
                         data = col, REML = FALSE)
c_drop <- lmerTest::lmer(color_num ~ treatment * wound * day + thicket + (1|tank) + (1|id),
                         data = col_drop, REML = FALSE)
results$col_full <- grab(c_full, "Color D-scale (full data)",     "treatment:day")
results$col_drop <- grab(c_drop, "Color D-scale (flagged dropped)","treatment:day")

# ---- Areal calcification: tank-level treatment permutation -----------------
# Growth (areal calcification, mg CaCO3 per cm2 per day, from buoyant weight).
# Unlike PAM/colour this is a single end-of-experiment value per coral, and the
# temperature treatment was applied at the TANK level (not the coral level), so
# the tank — not the coral — is the true unit of replication. With only ~4 tanks
# per treatment, a mixed model is overkill and fragile. Instead we average each
# tank to one number and run an exact permutation test on the tank means.
bw <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds")) |>
  mutate(thicket = factor(thicket)) |>
  filter(is.finite(areal_calc))                 # drop NA/Inf growth values
bw_drop <- bw |> filter(!id %in% FLAGGED_IDS, tank != FLAGGED_TANK)

# Exact (combinatorial) permutation test on tank-mean growth.
tank_perm <- function(dat, label) {
  # Collapse corals to one mean per tank -> tank is the unit of analysis.
  tank_growth <- dat |>
    group_by(tank, treatment) |>
    summarise(mean_areal_calc = mean(areal_calc, na.rm = TRUE),
              n_corals = n(), .groups = "drop")

  # Guard: need both treatments present with >=2 tanks each, else the test is
  # undefined (e.g. after dropping tank 3) -> return an all-NA placeholder row.
  n_by_temp <- table(tank_growth$treatment)
  if (length(n_by_temp) != 2 || any(n_by_temp < 2)) {
    return(tibble(response = label, term = "treatment (tank permutation)",
                  F_value = NA_real_, p_value = NA_real_,
                  estimate_28_minus_31 = NA_real_, n_tanks = nrow(tank_growth)))
  }

  vals <- tank_growth$mean_areal_calc
  trt <- tank_growth$treatment
  n28 <- sum(trt == "28C")
  # Observed effect: difference in mean growth, ambient minus heated.
  obs <- mean(vals[trt == "28C"]) - mean(vals[trt == "31C"])
  # Null distribution: re-label which tanks are "28C" in EVERY possible way
  # (all combn() choices of n28 tanks) and recompute the difference each time.
  # With so few tanks this enumerates the entire permutation null exactly — no
  # random sampling needed.
  null <- vapply(combn(seq_along(vals), n28, simplify = FALSE), \(idx) {
    mean(vals[idx]) - mean(vals[-idx])
  }, numeric(1))
  # Two-sided exact p: fraction of relabelings with an effect at least as
  # extreme (in absolute value) as the one actually observed.
  p <- mean(abs(null) >= abs(obs))

  tibble(response = label, term = "treatment (tank permutation)",
         F_value = NA_real_, p_value = signif(p, 3),
         estimate_28_minus_31 = round(obs, 3),
         n_tanks = nrow(tank_growth))
}
results$bw_full <- tank_perm(bw, "Areal calcification (full data)")
results$bw_drop <- tank_perm(bw_drop, "Areal calcification (flagged dropped)")

# ---- Collect + write -------------------------------------------------------
# Stack the full and dropped rows for every response into one table. Read it in
# full/dropped PAIRS: if F, p, and the growth estimate barely change, the
# conclusions are robust to the flagged samples.
out <- bind_rows(results)
write_csv(out, file.path(TBL_DIR, "22_sensitivity_flagged.csv"))

cat("\n=== Sensitivity to flagged samples (corals 116/121 + tank 3) ===\n")
print(as.data.frame(out))
cat("\nInterpret the PAM/color rows by their treatment:day F-tests and growth by",
    "the tank-level exact permutation p-value.\n")
cat("Wrote output/tables/22_sensitivity_flagged.csv\n")
