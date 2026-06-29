# =============================================================================
# Purpose: REPRODUCIBILITY CHECK (ADVISORY) — verify that every reported
#          number in the manuscript still matches the freshly-regenerated
#          analysis outputs. Each check recomputes a canonical value from the
#          source tables (themselves regenerated upstream) and asserts the
#          manuscript text contains it. If a number has drifted, this script
#          emits a WARNING and writes a report — it does NOT stop the pipeline
#          (no build-fail), so the run still completes and regenerates.
#
#          SCOPE: this only covers the Stier-lab-authored **phenotype Methods
#          and Results** (physiology, morphology, growth, genet variation,
#          thermal context). It does NOT police the Introduction, Discussion,
#          Abstract, or the transcriptomics — those are the lead author's
#          (S. Banerjee) narrative and are not gated here.
#
#          To keep a new phenotype number checkable: add a check below. The
#          rule is one-directional — PASS = the manuscript contains the current
#          value; FLAG = it does not (stale, or formatted off).
#
# What & why: numbers in a manuscript drift. An analysis is rerun, an estimate
#   shifts in the third decimal, but the prose still quotes the old value.
#   This script checks for that, and it runs LAST in the pipeline (after every
#   table has been regenerated). For each headline phenotype number it recomputes
#   the value straight from the fresh outputs, then searches the manuscript text
#   for that value. If it is there -> PASS; if it is not -> FLAG (and a warning).
#   This is ADVISORY: it never errors out the build. A flag means "check this
#   sentence," not a failure — so a number that is mid-edit cannot block the
#   pipeline. It only covers the Stier-lab phenotype sections;
#   the lead author's narrative (Intro/Discussion/Abstract/transcriptomics) is
#   out of scope. Currently 15/15 checks pass.
#
# Input:   manuscript/Manuscript_LTH.md
#          output/tables/{12_anova_summary,14_interval_survreg,
#                         14_cox_hazard_ratios,14_cox_genet_LRT,
#                         14_milestone_lag_summary,15_genet_pca_displacement,
#                         26_thermal_context}.csv
#          data/processed/{buoyant_weight_clean,coral_physio_wide}.rds
# Output:  output/tables/30_manuscript_audit.csv   — one row per checked number
#          warns (does not stop) if any phenotype number has drifted
# =============================================================================

# 00_setup.R loads packages and defines shared paths (TBL_DIR for output/tables,
# DATA_PROC, ...).
source(here::here("code", "00_setup.R"))

# Read the whole manuscript into ONE string (lines collapsed with "\n") so the
# token searches below are simple substring tests against the full text.
ms_path <- here::here("manuscript", "Manuscript_LTH.md")
ms <- paste(readLines(ms_path, warn = FALSE), collapse = "\n")
# Normalize the typographic minus (U+2212, used in the prose) to ASCII "-" so a
# value like -0.03 matches formatC()'s "-0.03". (Numbers only; harmless to text.)
ms <- gsub("−", "-", ms)

# ---- Matching helpers ------------------------------------------------------
# A claim PASSES when the manuscript contains EVERY one of its identifying
# tokens (so a coincidental lone "0.22" elsewhere can't make it pass). Each
# token is matched against a few rounding candidates to tolerate the
# manuscript's display precision (e.g. 106.57 written as "106.6").
fmt1 <- function(v, d) formatC(round(v, d), format = "f", digits = d)  # scalar
fmt_cands <- function(x, digits) {
  # Match at the registered precision and one decimal finer. We do
  # NOT go coarser: "3.7" would substring-match inside a stale "3.74" and hide
  # real drift. Register `digits` = the precision the manuscript actually uses.
  ds <- unique(c(digits, digits + 1L))
  unique(vapply(ds, \(d) fmt1(x, d), character(1)))
}
token_found <- function(x, digits) {
  any(vapply(fmt_cands(x, digits), \(s) grepl(s, ms, fixed = TRUE), logical(1)))
}

# checks accumulates one tibble row per registered claim; add_check() appends to
# it via <<- (assign in the enclosing scope). Each check PASSES only when ALL of
# its identifying tokens are found, so a multi-part claim (e.g. estimate + both CI
# bounds) can't pass on a partial coincidental match.
checks <- list()
add_check <- function(label, values, digits, source) {
  # values: numeric vector of identifying tokens; digits: matching precision(s)
  digits <- rep_len(digits, length(values))  # recycle one digit spec across all values
  oks <- mapply(token_found, values, digits)
  checks[[length(checks) + 1]] <<- tibble(
    check    = label,
    expected = paste(mapply(fmt1, values, digits), collapse = " | "),
    found    = all(oks),
    missing  = paste(mapply(fmt1, values[!oks], digits[!oks]), collapse = ", "),
    source   = source
  )
}

# ---- 1. New-corallite regeneration timing (headline + Cox summary) ----------
aft <- read_csv(file.path(TBL_DIR, "14_interval_survreg.csv"), show_col_types = FALSE)
nc_aft <- aft |> filter(trait == "new_corallites_on_tip")
add_check("Interval AFT new corallites (time ratio + 95% CI bounds)",
          c(nc_aft$time_ratio_31_vs28[1], nc_aft$ratio_lo[1], nc_aft$ratio_hi[1]), 2,
          "14_interval_survreg.csv")

cox <- read_csv(file.path(TBL_DIR, "14_cox_hazard_ratios.csv"), show_col_types = FALSE)
nc  <- cox |> filter(trait == "new_corallites_on_tip", grepl("overall", scope))
add_check("Cox HR new corallites (HR + 95% CI bounds)",
          c(nc$HR_31_vs28[1], nc$HR_lo[1], nc$HR_hi[1]), 2,
          "14_cox_hazard_ratios.csv")
add_check("Cox p new corallites", nc$p[1], 3, "14_cox_hazard_ratios.csv")

# ---- 2. Growth (% skeletal mass change) reduction --------------------------
bw <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds"))
ac <- bw |> filter(is.finite(pct_growth)) |>
  group_by(treatment) |> summarise(m = mean(pct_growth), .groups = "drop")
m28 <- ac$m[ac$treatment == "28C"]; m31 <- ac$m[ac$treatment == "31C"]
add_check("Growth (% mass change) reduction % and means",
          c((1 - m31 / m28) * 100, m28, m31), c(0L, 2L, 2L),
          "buoyant_weight_clean.rds")

# ---- 3. Treatment x time F-stats (PAM, color, symbionts) -------------------
an <- read_csv(file.path(TBL_DIR, "12_anova_summary.csv"), show_col_types = FALSE)
getF <- function(resp, trm) {
  r <- an |> filter(response_id == resp, term == trm)
  coalesce(r$`F value`[1], r$Chisq[1])
}
add_check("PAM treatment x day F",          getF("pam_fvfm", "treatment:day"), 1, "12_anova_summary.csv")
add_check("Color treatment x day F",        getF("color_dscale", "treatment:day"), 1, "12_anova_summary.csv")
add_check("Symbiont treatment x day F",     getF("log_zoox_density", "treatment:biopsy_day_c"), 1, "12_anova_summary.csv")

# ---- 4. Per-genet PCA centroid displacement --------------------------------
disp <- read_csv(file.path(TBL_DIR, "15_genet_pca_displacement.csv"), show_col_types = FALSE)
gd <- \(g) disp$displacement[disp$thicket == g]
add_check("PCA displacement genets a / d / c",
          c(gd("a"), gd("d"), gd("c")), 2, "15_genet_pca_displacement.csv")

# ---- 5. PCA variance on PC1 (recomputed live) ------------------------------
w  <- as.data.frame(readRDS(file.path(DATA_PROC, "coral_physio_wide.rds")))
pv <- c("pam_end", "color_end", "growth_pct", "zoox_end")
pcm <- w[complete.cases(w[, pv]), pv]
ve  <- prcomp(pcm, center = TRUE, scale. = TRUE)$sdev^2
add_check("PCA variance PC1 (%)", 100 * ve[1] / sum(ve), 0, "coral_physio_wide.rds (prcomp)")

# ---- 6. Cunning ED50 thermal context ---------------------------------------
tc <- read_csv(file.path(TBL_DIR, "26_thermal_context.csv"), show_col_types = FALSE)
tv <- \(m) tc$value[tc$metric == m]
add_check("Acute ED50 mean (°C)",            tv("acute ED50 mean (°C)"), 1, "26_thermal_context.csv")
add_check("31 °C below mean ED50 (°C)",      tv("31 °C below mean ED50 (°C)"), 1, "26_thermal_context.csv")

# ---- 7. Healing-to-regeneration censored fraction --------------------------
lag <- read_csv(file.path(TBL_DIR, "14_milestone_lag_summary.csv"), show_col_types = FALSE)
pc31 <- lag$pct_closed_no_regen[lag$treatment == "31C"][1]
pc28 <- lag$pct_closed_no_regen[lag$treatment == "28C"][1]
add_check("Closed-but-never-regenerated % (31C / 28C)", c(pc31, pc28), 0,
          "14_milestone_lag_summary.csv")

# ---- 8. Genet x treatment LRTs the manuscript reports (continuous responses) -
g13 <- read_csv(file.path(TBL_DIR, "13_genet_anova.csv"), show_col_types = FALSE)
gchi <- \(resp) g13$lrt_chisq[g13$response == resp]
add_check("Genet x treatment LRT chi-sq (PAM / color / symbionts)",
          c(gchi("pam_fvfm"), gchi("color_dscale"), gchi("log_zoox")), 1,
          "13_genet_anova.csv")

# ---- 9. Color % paled at end of experiment ---------------------------------
cp <- read_csv(file.path(TBL_DIR, "03_color_end_proportions.csv"), show_col_types = FALSE)
pp <- \(tt, wd) 100 * cp$prop_paled[cp$treatment == tt & cp$wound == wd]
add_check("Color % paled (31C wounded / 31C unwounded / 28C unwounded)",
          c(pp("31C", "yes"), pp("31C", "no"), pp("28C", "no")), 0,
          "03_color_end_proportions.csv")

# ---- 10. Per-genet composite standardized heat sensitivity -----------------
rs <- read_csv(file.path(TBL_DIR, "19_genet_resilience_summary.csv"), show_col_types = FALSE)
sens <- \(g) rs$mean_sensitivity[rs$thicket == g]
add_check("Per-genet composite heat sensitivity (a / c / d)",
          c(sens("a"), sens("c"), sens("d")), 2,
          "19_genet_resilience_summary.csv")

# ===========================================================================
# Report (ADVISORY — warns, never stops)
# ===========================================================================
audit <- bind_rows(checks)
write_csv(audit, file.path(TBL_DIR, "30_manuscript_audit.csv"))

cat("\n=== Phenotype reproducibility check (manuscript vs regenerated tables) ===\n")
cat("Scope: Stier-lab phenotype Methods/Results only (not Intro/Discussion/Abstract/transcriptomics).\n")
print(as.data.frame(audit[, c("check", "expected", "found", "missing")]),
      row.names = FALSE)
n_fail <- sum(!audit$found)
cat(sprintf("\n%d/%d checks PASS.\n", sum(audit$found), nrow(audit)))

if (n_fail > 0) {
  bad <- audit |> filter(!found)
  warning(sprintf(
    paste0("PHENOTYPE NUMBERS MAY BE STALE (advisory — pipeline NOT failed): ",
           "%d phenotype number(s) in manuscript/Manuscript_LTH.md no longer ",
           "match the regenerated analysis outputs.\n%s\n",
           "Fix when convenient: update the flagged value(s) in the manuscript ",
           "from output/tables/20_master_results.csv. (This check covers only the ",
           "Stier-lab phenotype Methods/Results, not the lead author's narrative.)"),
    n_fail,
    paste(sprintf("  - %s: expected %s (source %s); manuscript missing: %s",
                  bad$check, bad$expected, bad$source, bad$missing),
          collapse = "\n")),
    call. = FALSE)
  cat(sprintf("\n[ADVISORY] %d phenotype number(s) flagged as possibly stale — see warning above and 30_manuscript_audit.csv.\n", n_fail))
} else {
  cat("All checked phenotype numbers are in sync.\n")
}
cat("Wrote output/tables/30_manuscript_audit.csv\n")
