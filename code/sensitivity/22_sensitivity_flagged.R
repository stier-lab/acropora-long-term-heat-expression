# =============================================================================
# Purpose: Sensitivity analysis for the samples/tanks Molly flagged as odd
#          (see notes/QAQC_flagged_samples.md). Re-runs the key treatment-effect
#          tests with the flagged colonies (116, 121) and the slow ambient tank
#          (tank 3) removed, and compares the treatment effect to the full-data
#          model. If the conclusions hold, the flags are documented but
#          inconsequential.
# Input:   data/processed/{pam_clean,color_clean,buoyant_weight_clean}.rds
# Output:  output/tables/22_sensitivity_flagged.csv
# =============================================================================

source(here::here("code", "00_setup.R"))

FLAGGED_IDS  <- c(116, 121)
FLAGGED_TANK <- 3

results <- list()

# Helper: extract the treatment×day (or treatment) F and p from an lmerTest fit
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

# ---- PAM Fv/Fm: treatment × day ----
pam <- readRDS(file.path(DATA_PROC, "pam_clean.rds")) |>
  mutate(thicket = factor(thicket))
pam_drop <- pam |> filter(!id %in% FLAGGED_IDS, tank != FLAGGED_TANK)

f_full <- lmerTest::lmer(fv_fm ~ treatment * wound * day + (1|tank) + (1|thicket) + (1|id),
                         data = pam, REML = FALSE)
f_drop <- lmerTest::lmer(fv_fm ~ treatment * wound * day + (1|tank) + (1|thicket) + (1|id),
                         data = pam_drop, REML = FALSE)
results$pam_full <- grab(f_full, "PAM Fv/Fm (full data)",     "treatment:day")
results$pam_drop <- grab(f_drop, "PAM Fv/Fm (flagged dropped)","treatment:day")

# ---- Color D-scale: treatment × day ----
col <- readRDS(file.path(DATA_PROC, "color_clean.rds")) |>
  mutate(thicket = factor(thicket))
col_drop <- col |> filter(!id %in% FLAGGED_IDS, tank != FLAGGED_TANK)

c_full <- lmerTest::lmer(color_num ~ treatment * wound * day + (1|tank) + (1|thicket) + (1|id),
                         data = col, REML = FALSE)
c_drop <- lmerTest::lmer(color_num ~ treatment * wound * day + (1|tank) + (1|thicket) + (1|id),
                         data = col_drop, REML = FALSE)
results$col_full <- grab(c_full, "Color D-scale (full data)",     "treatment:day")
results$col_drop <- grab(c_drop, "Color D-scale (flagged dropped)","treatment:day")

# ---- Areal calcification: tank-level treatment permutation -----------------
bw <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds")) |>
  mutate(thicket = factor(thicket)) |>
  filter(is.finite(areal_calc))
bw_drop <- bw |> filter(!id %in% FLAGGED_IDS, tank != FLAGGED_TANK)

tank_perm <- function(dat, label) {
  tank_growth <- dat |>
    group_by(tank, treatment) |>
    summarise(mean_areal_calc = mean(areal_calc, na.rm = TRUE),
              n_corals = n(), .groups = "drop")

  n_by_temp <- table(tank_growth$treatment)
  if (length(n_by_temp) != 2 || any(n_by_temp < 2)) {
    return(tibble(response = label, term = "treatment (tank permutation)",
                  F_value = NA_real_, p_value = NA_real_,
                  estimate_28_minus_31 = NA_real_, n_tanks = nrow(tank_growth)))
  }

  vals <- tank_growth$mean_areal_calc
  trt <- tank_growth$treatment
  n28 <- sum(trt == "28C")
  obs <- mean(vals[trt == "28C"]) - mean(vals[trt == "31C"])
  null <- vapply(combn(seq_along(vals), n28, simplify = FALSE), \(idx) {
    mean(vals[idx]) - mean(vals[-idx])
  }, numeric(1))
  p <- mean(abs(null) >= abs(obs))

  tibble(response = label, term = "treatment (tank permutation)",
         F_value = NA_real_, p_value = signif(p, 3),
         estimate_28_minus_31 = round(obs, 3),
         n_tanks = nrow(tank_growth))
}
results$bw_full <- tank_perm(bw, "Areal calcification (full data)")
results$bw_drop <- tank_perm(bw_drop, "Areal calcification (flagged dropped)")

out <- bind_rows(results)
write_csv(out, file.path(TBL_DIR, "22_sensitivity_flagged.csv"))

cat("\n=== Sensitivity to flagged samples (corals 116/121 + tank 3) ===\n")
print(as.data.frame(out))
cat("\nInterpret the PAM/color rows by their treatment:day F-tests and growth by",
    "the tank-level exact permutation p-value.\n")
cat("Wrote output/tables/22_sensitivity_flagged.csv\n")
