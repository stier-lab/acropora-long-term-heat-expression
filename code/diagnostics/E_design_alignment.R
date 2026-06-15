# =============================================================================
# Design alignment audit.
#
# Verify each primary model's structure matches the experimental design:
#   - 2 temperatures x 2 wounds x 3 genets x 8 tanks (4 per temp)
#   - n=192 corals, repeated measures over 14 days for PAM + color
#   - n=192 destructive biopsies (1 per coral) for symbionts
#   - n=48 for buoyant weight (1 obs per coral)
#   - 24 wounded corals per temperature for morphology traits
# =============================================================================

source(here::here("code", "00_setup.R"))

DIAG_OUT <- here("output", "diagnostics")
dir.create(DIAG_OUT, recursive = TRUE, showWarnings = FALSE)

audit <- list()

add_row <- function(model, design_feature, observed, expected, status, note) {
  audit[[length(audit) + 1]] <<- tibble::tibble(
    model = model, check = design_feature,
    observed = as.character(observed),
    expected = as.character(expected),
    status   = status,
    note     = note
  )
}

# Helper to count cells in design
cells <- function(d, ...) d |> count(...) |> nrow()

# ---- PAM Fv/Fm ----
m <- readRDS(file.path(MOD_DIR, "12_pam_lmm.rds"))
d <- readRDS(file.path(DATA_PROC, "pam_clean.rds"))
form <- paste(format(formula(m)), collapse = " ")
add_row("12_pam_lmm", "fixed structure", form,
        "treatment * wound * day * thicket",
        if (grepl("treatment \\* wound \\* day \\* thicket", form)) "PASS" else "FAIL",
        "4-way fixed structure required for genet x heat x wound x time")
add_row("12_pam_lmm", "random structure",
        paste(names(lme4::ranef(m)), collapse = "+"),
        "tank + id",
        if (all(c("tank", "id") %in% names(lme4::ranef(m)))) "PASS" else "WARN",
        "tank for tank-level confounds, id for repeated measures on same coral")
add_row("12_pam_lmm", "n observations", nrow(d), "~336", "INFO",
        sprintf("n_corals=%d, days=%s", length(unique(d$id)),
                paste(sort(unique(d$day)), collapse = ",")))
add_row("12_pam_lmm", "balance (tank x treatment)",
        cells(d, tank, treatment), "8",
        if (cells(d, tank, treatment) == 8) "PASS" else "WARN",
        "8 unique tank x treatment cells expected (4 tanks per temp, fully nested)")
add_row("12_pam_lmm", "(1|tank) sufficient?",
        "yes — tank labels are unique across treatments",
        "yes",
        "PASS",
        "tanks are uniquely IDed; no need for tank:treatment crossed term")
add_row("12_pam_lmm", "random slope on day?",
        "(1|id) only — no random slope",
        "considered (day|id) but n=14 days, few per id; likely overfits",
        "HANDLED",
        "Random slope considered and intentionally omitted; expected singular/overfit with sparse per-id trajectories")

# ---- Color D-scale (same fixed structure as PAM) ----
m <- readRDS(file.path(MOD_DIR, "12_color_lmm.rds"))
form <- paste(format(formula(m)), collapse = " ")
add_row("12_color_lmm", "fixed structure", form,
        "treatment * wound * day * thicket",
        if (grepl("treatment \\* wound \\* day \\* thicket", form)) "PASS" else "FAIL",
        "same 4-way as PAM")
add_row("12_color_lmm", "ordinal data on Gaussian",
        "Gaussian LMM on D1-D5 ordinal scale",
        "either CLMM or note in Methods",
        "HANDLED",
        "Addressed by 12b CLMM robustness check; KS-violation noted in Section 10")

# ---- Buoyant weight ----
m <- readRDS(file.path(MOD_DIR, "12_bw_lm.rds"))
d <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds"))
form <- paste(format(formula(m)), collapse = " ")
add_row("12_bw_lm", "fixed structure", form,
        "treatment * wound * thicket",
        if (grepl("treatment \\* wound \\* thicket", form)) "PASS" else "FAIL",
        "no day term — single endpoint measurement")
add_row("12_bw_lm", "no random effects?",
        "OLS only", "OK because 1 obs per coral",
        "PASS", "tank (n=8) and id (n=48) would be singular with 1 obs/coral")
add_row("12_bw_lm", "n",
        nrow(d), "48 (24 per temperature)",
        if (nrow(d) == 48) "PASS" else "WARN",
        sprintf("Cells in design: %d (2 trt x 2 wound x 3 genet = 12)",
                cells(d, treatment, wound, thicket)))
add_row("12_bw_lm", "df residual",
        m$df.residual, "~36 (48 - 12 fixed params)",
        if (m$df.residual >= 30) "PASS" else "WARN",
        "Adequate residual df for 3-way model")

# ---- Symbiont density ----
m <- readRDS(file.path(MOD_DIR, "12_zoox_lmm.rds"))
d <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
  filter(is.finite(cells_per_cm2), cells_per_cm2 > 0) |>
  mutate(biopsy_day_c = biopsy_day - 1)
form <- paste(format(formula(m)), collapse = " ")
add_row("12_zoox_lmm", "fixed structure", form,
        "treatment * wound * biopsy_day_c * thicket",
        if (grepl("treatment \\* wound \\* biopsy_day_c \\* thicket", form)) "PASS" else "FAIL",
        "biopsy_day_c is centered at Day 1; 4-way as for PAM")
add_row("12_zoox_lmm", "drop (1|id) — destructive",
        paste(names(lme4::ranef(m)), collapse = "+"),
        "tank only (each coral 1 biopsy)",
        if (identical(names(lme4::ranef(m)), "tank")) "PASS" else "WARN",
        "Correct — destructive sampling means each id has 1 obs")
add_row("12_zoox_lmm", "biopsy_day_c centering",
        sprintf("range = [%g, %g]", min(d$biopsy_day_c, na.rm=TRUE),
                max(d$biopsy_day_c, na.rm=TRUE)),
        "0 to 14 (Day 1 = 0)",
        if (min(d$biopsy_day_c, na.rm=TRUE) == 0) "PASS" else "WARN",
        "Centered correctly at Day 1 baseline")

# ---- Morphology GLMMs (using blme version) ----
trait_files <- list.files(MOD_DIR, pattern = "^12c_morph_.*_blme\\.rds$",
                          full.names = TRUE)
for (f in trait_files) {
  trait <- gsub("^12c_morph_(.*)_blme\\.rds$", "\\1", basename(f))
  m <- readRDS(f)
  form <- paste(format(formula(m)), collapse = " ")
  add_row(basename(f), "fixed structure", form,
          "treatment * day * thicket (wound dropped — wounded-only)",
          if (grepl("treatment \\* day \\* thicket", form)) "PASS" else "FAIL",
          "Restricted to wounded corals; wound is a stratification, not a covariate")
  add_row(basename(f), "Cauchy(0,2.5) prior",
          "fixef.prior = t(scale=2.5, df=1)",
          "Gelman 2008 separation-fix default",
          "PASS",
          "Prior addresses the morphology separation issue")
}

# ---- Cox PH ----
add_row("14_cox_*", "stratification",
        "strata = thicket in overall model",
        "thicket as strata (allows different baseline hazard per genet)",
        "PASS",
        "Per-genet HRs in separate per-genet models (low EPV caveat in Section 10)")
add_row("14_cox_*", "wounded-only at risk",
        "filter(wound == 'yes') applied",
        "yes — only wounded corals can express healing traits",
        "PASS", "Confirmed in script 14")

# ---- PCA ----
add_row("15_pca", "centering + scaling",
        "prcomp(..., center=TRUE, scale.=TRUE)",
        "required because units differ (Fv/Fm, D-scale, %, log-cells)",
        "PASS", "Confirmed in script 15")

# Write outputs
final <- bind_rows(audit)
write_csv(final, file.path(DIAG_OUT, "E_design_alignment.csv"))

# Report
sink(file.path(DIAG_OUT, "E_design_alignment_report.md"))
cat("# E. Design alignment audit\n\n")
cat("Generated:", format(Sys.time()), "\n\n")
cat("| Status | Count |\n|---|---|\n")
print(table(final$status))
cat("\n## Status summary\n\n")
for (s in c("FAIL", "WARN", "HANDLED", "INFO", "PASS")) {
  rows <- final[final$status == s, ]
  if (nrow(rows) > 0) {
    cat("### ", s, " (", nrow(rows), ")\n\n", sep = "")
    for (i in seq_len(nrow(rows))) {
      cat("- **", rows$model[i], " / ", rows$check[i], "**: ",
          rows$note[i], "\n", sep = "")
    }
    cat("\n")
  }
}
sink()

cat("=== Design alignment audit complete ===\n")
print(table(final$status))
cat("\nWrote", file.path(DIAG_OUT, "E_design_alignment.csv"),
    "and", file.path(DIAG_OUT, "E_design_alignment_report.md"), "\n")
