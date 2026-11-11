# =============================================================================
# Purpose: Master spreadsheet coverage check (diagnostic suite H) — confirm that
#          every statistic from the per-script result CSVs made it into
#          the single master results table.
#
# What & why: each analysis script writes its own results CSV, and script 20
#   consolidates them into output/tables/20_master_results.csv — the one table the
#   manuscript draws from. A statistic that silently fails to copy across would
#   never appear in the paper. This script audits that hand-off: for every source
#   CSV it checks how many of its rows are represented in the master (MISSING /
#   PARTIAL / COVERED), then runs integrity checks on the master itself — duplicate
#   rows, and domain mislabeling (e.g. a morphology result tagged Physiology, or a
#   GLMM result tagged LMM). It's an accounting audit, not a statistical test.
# Input:   output/tables/20_master_results.csv (the consolidated table)
#          output/tables/<per-script>.csv      (the sources it should contain)
# Output:  output/diagnostics/H_coverage_by_source.csv
#          output/diagnostics/H_duplicates.csv
#          output/diagnostics/H_coverage_report.md
# =============================================================================

# 00_setup.R loads shared packages + paths (TBL_DIR points at output/tables).
source(here::here("code", "00_setup.R"))

DIAG_OUT <- here("output", "diagnostics")
dir.create(DIAG_OUT, recursive = TRUE, showWarnings = FALSE)

# The consolidated table everything is checked against.
master <- read_csv(file.path(TBL_DIR, "20_master_results.csv"),
                   show_col_types = FALSE)

# ---- Coverage per source CSV -----------------------------------------------
# The result CSVs we expect to be represented in the master. Infrastructure /
# bookkeeping CSVs are left off this list.
src_tables <- c(
  "12_anova_summary.csv",
  "12_genet_treatment_effects.csv",
  "12_r2_summary.csv",
  "13_genet_anova.csv",
  "14_cox_hazard_ratios.csv",
  "14_cox_genet_LRT.csv",
  "05_buoyant_weight_lm.csv",
  "12b_color_clmm.csv",
  "12c_morph_blme_fixed_effects.csv",
  "12c_morph_blme_anova.csv",
  "15_pca_loadings.csv",
  "15_genet_pca_displacement.csv",
  "19_genet_resilience_summary.csv",
  "19c_resilience_decomp_by_scope.csv"
)

findings <- list()

for (src in src_tables) {
  path <- file.path(TBL_DIR, src)
  # Source file doesn't exist: record it and skip — can't audit what isn't there.
  if (!file.exists(path)) {
    findings[[src]] <- tibble(source_csv = src, rows = 0,
                              note = "FILE MISSING")
    next
  }
  d <- read_csv(path, show_col_types = FALSE)
  # Count master rows whose source_artifact names this CSV. The master tags each
  # row with where it came from, so this is how many of its rows propagated.
  n_master <- sum(grepl(src, master$source_artifact, fixed = TRUE))
  findings[[src]] <- tibble(
    source_csv = src,
    src_rows   = nrow(d),
    master_rows = n_master,
    # MISSING: nothing copied. COVERED: at least half the source rows present
    # (some source rows are summarized, so 50% is the threshold, not 100%).
    # PARTIAL: present but under that threshold.
    coverage_status = if (n_master == 0) "MISSING"
                       else if (n_master >= nrow(d) * 0.5) "COVERED"
                       else "PARTIAL"
  )
}

coverage <- bind_rows(findings)

# ---- Integrity checks on the master ----------------------------------------
# Duplicate check: the same statistic should appear once. Any (domain, response,
# model_type, term, source) combination occurring >1 time is a copy/merge error.
dups <- master |>
  count(domain, response, model_type, term, source_artifact, name = "n_rows") |>
  filter(n_rows > 1)

# Domain sanity: any response named morph_* must be tagged "Morphology". Catches
# morphology results filed under Physiology.
domain_audit <- master |>
  filter(grepl("^morph_", response, ignore.case = FALSE) &
         domain != "Morphology")

# Same idea, scoped to the penalized-morphology ANOVA source: every row from that
# file should be domain == "Morphology".
miscat <- master |>
  filter(source_artifact == "output/tables/12c_morph_blme_anova.csv" &
         domain != "Morphology")

# Model-type sanity: morphology rows in the 12_anova summary come from GLMMs.
# Tally model_type — anything labeled LMM here is mislabeled.
glmm_label <- master |>
  filter(source_artifact == "output/tables/12_anova_summary.csv" &
         grepl("morph_", response, ignore.case = TRUE)) |>
  count(model_type, name = "n")

write_csv(coverage, file.path(DIAG_OUT, "H_coverage_by_source.csv"))
write_csv(dups, file.path(DIAG_OUT, "H_duplicates.csv"))

# ---- Report ----------------------------------------------------------------
# sink() captures the cat() output below into the markdown report file.
sink(file.path(DIAG_OUT, "H_coverage_report.md"))
cat("# H. Master spreadsheet coverage check\n\n")
cat("Generated:", format(Sys.time()), "\n\n")
cat("Master rows:", nrow(master), "\n\n")

cat("## Coverage per source CSV\n\n| File | Src rows | Master rows | Status |\n|---|---|---|---|\n")
for (i in seq_len(nrow(coverage))) {
  cat("| ", coverage$source_csv[i],
      " | ", coverage$src_rows[i],
      " | ", coverage$master_rows[i],
      " | ", coverage$coverage_status[i], " |\n", sep = "")
}

cat("\n## Duplicate (domain, response, model_type, term, source) rows:",
    nrow(dups), "\n")
if (nrow(dups) > 0) {
  for (i in seq_len(min(20, nrow(dups))))
    cat("- ", dups$response[i], " / ", dups$term[i],
        " (n=", dups$n_rows[i], ")\n", sep = "")
}

cat("\n## Mis-categorized morph rows (tagged Physiology):", nrow(domain_audit), "\n")
if (nrow(domain_audit) > 0)
  print(domain_audit |> select(domain, response, term, source_artifact))

cat("\n## Rows from 12c morphology that aren't Morphology:", nrow(miscat), "\n")
if (nrow(miscat) > 0)
  print(miscat |> select(domain, response, term))

cat("\n## 12_anova morph_* rows by model_type (should all be GLMM, not LMM):\n")
print(glmm_label)
sink()

cat("=== Coverage audit complete ===\n")
cat("Source CSVs:", nrow(coverage), "\n")
cat("Duplicate rows:", nrow(dups), "\n")
cat("Mis-categorized morph rows:", nrow(domain_audit), "\n")
cat("12_anova morph rows by model_type:\n")
print(glmm_label)
