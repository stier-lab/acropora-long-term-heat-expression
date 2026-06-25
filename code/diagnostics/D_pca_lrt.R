# =============================================================================
# PCA and genet-LRT diagnostic — verifies:
#   1. PCA from 15_multivariate.R (n, variance explained, Kaiser, loadings,
#      centering/scaling).
#   2. Genet LRTs from 13_genet_interaction.R (REML status, df, convergence,
#      effect direction vs. emmeans).
# Output:
#   output/diagnostics/D_pca_lrt_diagnostics.csv
#   output/diagnostics/D_pca_lrt_report.md
#   figures/diagnostics/D_pca_scree.png
#   figures/diagnostics/D_pca_loadings.png
# =============================================================================

suppressPackageStartupMessages({
  source(here::here("code", "00_setup.R"))
  library(lme4)
})

DIAG_OUT <- here::here("output", "diagnostics")
DIAG_FIG <- here::here("figures", "diagnostics")
dir.create(DIAG_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(DIAG_FIG, recursive = TRUE, showWarnings = FALSE)

results <- list()
add_check <- function(analysis, check, statistic, status, notes = "") {
  results[[length(results) + 1]] <<- tibble(
    analysis = analysis, check = check,
    statistic = as.character(statistic), status = status, notes = notes
  )
}

report_lines <- c("# PCA + Genet LRT Diagnostics",
                  sprintf("_Run: %s_\n", format(Sys.time())))

# ----------------------------------------------------------------------------
# PCA CHECKS
# ----------------------------------------------------------------------------
report_lines <- c(report_lines, "## PCA (15_multivariate.R)\n")

pam   <- readRDS(file.path(DATA_PROC, "pam_clean.rds"))
color <- readRDS(file.path(DATA_PROC, "color_clean.rds"))
bw    <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds"))
phys  <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
  filter(is.finite(cells_per_cm2))

pam_last <- pam |> group_by(id, treatment, wound, thicket) |>
  filter(day == max(day, na.rm = TRUE)) |>
  summarise(pam_end = mean(fv_fm), .groups = "drop")
color_last <- color |> group_by(id, treatment, wound, thicket) |>
  filter(day == max(day, na.rm = TRUE)) |>
  summarise(color_end = mean(color_num), .groups = "drop")
zoox_last <- phys |> group_by(id, treatment, wound, thicket) |>
  filter(biopsy_day == max(biopsy_day, na.rm = TRUE)) |>
  summarise(zoox_end = log10(mean(cells_per_cm2)), .groups = "drop")
bw_summary <- bw |> group_by(id, treatment, wound, thicket) |>
  summarise(growth_pct = pct_growth[1], .groups = "drop")

wide <- pam_last |>
  inner_join(color_last, by = c("id", "treatment", "wound", "thicket")) |>
  inner_join(bw_summary, by = c("id", "treatment", "wound", "thicket")) |>
  left_join(zoox_last,   by = c("id", "treatment", "wound", "thicket"))

pca_vars <- c("pam_end", "color_end", "growth_pct")
if (any(!is.na(wide$zoox_end))) pca_vars <- c(pca_vars, "zoox_end")
pca_input <- wide[, pca_vars, drop = FALSE] |> drop_na()
n_pca <- nrow(pca_input)

add_check("PCA", "Sample size (complete cases)", n_pca,
          if (n_pca >= 30) "PASS" else "WARN",
          sprintf("%d corals with all %d vars", n_pca, length(pca_vars)))
report_lines <- c(report_lines,
  sprintf("- n = %d corals on %d vars (%s)", n_pca, length(pca_vars),
          paste(pca_vars, collapse = ", ")))

pca <- prcomp(pca_input, center = TRUE, scale. = TRUE)
var_exp <- summary(pca)$importance[2, ] * 100
cum_var <- summary(pca)$importance[3, ] * 100
eigen <- pca$sdev^2

for (i in seq_along(var_exp)) {
  add_check("PCA", sprintf("Variance explained PC%d", i),
            sprintf("%.1f%%", var_exp[i]), "INFO",
            sprintf("cumulative = %.1f%%", cum_var[i]))
}
pc12 <- sum(var_exp[1:2])
add_check("PCA", "PC1+PC2 cumulative variance",
          sprintf("%.1f%%", pc12),
          if (pc12 >= 60) "PASS" else "WARN",
          if (pc12 >= 60) "Biplot interpretation supported"
          else "PC1+PC2 < 60% — biplot may miss meaningful structure")
report_lines <- c(report_lines,
  sprintf("- Variance explained: %s",
          paste(sprintf("PC%d=%.1f%%", seq_along(var_exp), var_exp),
                collapse = ", ")),
  sprintf("- PC1+PC2 = **%.1f%%** (%s)", pc12,
          if (pc12 >= 60) "PASS" else "WARN — below 60% threshold"))

n_kaiser <- sum(eigen > 1)
add_check("PCA", "Kaiser criterion (eigenvalues > 1)",
          n_kaiser, "INFO",
          sprintf("eigenvalues: %s",
                  paste(sprintf("%.2f", eigen), collapse = ", ")))
report_lines <- c(report_lines,
  sprintf("- Kaiser: %d PCs with eigenvalue > 1 (eigenvals %s)", n_kaiser,
          paste(sprintf("%.2f", eigen), collapse = ", ")))

# Centering & scaling
centered <- !is.null(pca$center) && length(pca$center) > 0
scaled   <- !is.null(pca$scale) && !isFALSE(pca$scale)
add_check("PCA", "Centered AND scaled",
          sprintf("center=%s, scale=%s", centered, scaled),
          if (centered && scaled) "PASS" else "FAIL",
          "Different units across PAM/color/growth/log-zoox require scaling")
report_lines <- c(report_lines,
  sprintf("- Center = %s, Scale = %s (%s)", centered, scaled,
          if (centered && scaled) "PASS" else "FAIL"))

# Loadings
load_mat <- pca$rotation
report_lines <- c(report_lines, "", "### Loadings", "")
for (pc in colnames(load_mat)) {
  ord <- load_mat[order(abs(load_mat[, pc]), decreasing = TRUE), pc]
  top_var <- names(ord)[1]
  report_lines <- c(report_lines,
    sprintf("- **%s** (%.1f%%): %s", pc, var_exp[pc],
            paste(sprintf("%s=%+.2f", names(ord), ord), collapse = ", ")))
  add_check("PCA", sprintf("Top loader on %s", pc),
            sprintf("%s (%+.2f)", top_var, ord[1]), "INFO", "")
}

# Biological interpretation of PC1
pc1 <- load_mat[, "PC1"]
hp_load <- pc1[c("pam_end", "color_end")]
hp_sym <- if (all(hp_load > 0) | all(hp_load < 0)) "PASS" else "WARN"
add_check("PCA", "PC1 biological interpretation",
          sprintf("PAM=%+.2f color=%+.2f growth=%+.2f%s",
                  pc1["pam_end"], pc1["color_end"], pc1["growth_pct"],
                  if ("zoox_end" %in% names(pc1))
                    sprintf(" zoox=%+.2f", pc1["zoox_end"]) else ""),
          hp_sym,
          paste0("Expect health vars (PAM, color, zoox) to load same sign — ",
                 if (hp_sym == "PASS") "consistent."
                 else "inconsistent signs flagged."))

# Scree + loadings figures
scree_df <- tibble(PC = factor(paste0("PC", seq_along(var_exp)),
                               levels = paste0("PC", seq_along(var_exp))),
                   variance = var_exp,
                   eigenvalue = eigen)
p_scree <- ggplot(scree_df, aes(PC, variance)) +
  geom_col(fill = "#56B4E9") +
  geom_text(aes(label = sprintf("%.0f%%", variance)), vjust = -0.4, size = 3) +
  geom_hline(yintercept = 100 / length(var_exp), linetype = "dashed",
             colour = "grey50") +
  labs(x = "Component", y = "Variance explained (%)",
       title = "PCA scree diagnostic",
       subtitle = sprintf("PC1+PC2 = %.1f%% | Kaiser PCs = %d",
                          pc12, n_kaiser)) +
  theme_pub(10)
ggsave(file.path(DIAG_FIG, "D_pca_scree.png"), p_scree,
       width = 120, height = 90, units = "mm", dpi = 200)

load_long <- as.data.frame(load_mat) |>
  tibble::rownames_to_column("variable") |>
  pivot_longer(-variable, names_to = "PC", values_to = "loading")
p_load <- ggplot(load_long, aes(variable, loading, fill = loading > 0)) +
  geom_col() +
  facet_wrap(~ PC, ncol = length(var_exp)) +
  scale_fill_manual(values = c("TRUE" = "#56B4E9", "FALSE" = "#D55E00"),
                    guide = "none") +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
  labs(x = NULL, y = "Loading",
       title = "PCA loadings diagnostic") +
  theme_pub(9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(DIAG_FIG, "D_pca_loadings.png"), p_load,
       width = 160, height = 90, units = "mm", dpi = 200)

# ----------------------------------------------------------------------------
# GENET LRT CHECKS
# ----------------------------------------------------------------------------
report_lines <- c(report_lines, "", "## Genet LRTs (13_genet_interaction.R)\n")

anova_tab <- read_csv(file.path(TBL_DIR, "13_genet_anova.csv"),
                      show_col_types = FALSE)
emm_tab   <- read_csv(file.path(TBL_DIR, "13_genet_emmeans.csv"),
                      show_col_types = FALSE)

pam_d   <- pam   |> mutate(thicket = factor(thicket))
color_d <- color |> mutate(thicket = factor(thicket))
bw_d    <- bw    |> mutate(thicket = factor(thicket))
phys_d  <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
  filter(is.finite(cells_per_cm2), cells_per_cm2 > 0) |>
  mutate(thicket = factor(thicket), biopsy_day_c = biopsy_day - 1)

check_lmer_pair <- function(name, response, data, fixed_term = "day",
                            include_id = TRUE) {
  id_term <- if (include_id) " + (1 | id)" else ""
  rhs0 <- sprintf("treatment * wound * %s + thicket + (1 | tank)%s",
                  fixed_term, id_term)
  rhs1 <- sprintf("treatment * wound * %s * thicket + (1 | tank)%s",
                  fixed_term, id_term)
  f0 <- as.formula(paste(response, "~", rhs0))
  f1 <- as.formula(paste(response, "~", rhs1))

  warns0 <- character(); warns1 <- character()
  m0 <- withCallingHandlers(
    lme4::lmer(f0, data = data, REML = FALSE),
    warning = function(w) { warns0 <<- c(warns0, conditionMessage(w))
                             invokeRestart("muffleWarning") })
  m1_path <- file.path(MOD_DIR, sprintf("13_%s_genet_lmm.rds", name))
  m1_saved <- readRDS(m1_path)
  m1 <- withCallingHandlers(
    lme4::lmer(f1, data = data, REML = FALSE),
    warning = function(w) { warns1 <<- c(warns1, conditionMessage(w))
                             invokeRestart("muffleWarning") })

  reml0 <- isREML(m0); reml1 <- isREML(m1); reml_saved <- isREML(m1_saved)
  add_check(name, "Null model REML status", sprintf("REML=%s", reml0),
            if (!reml0) "PASS" else "FAIL",
            "LRT on fixed effects requires REML=FALSE")
  add_check(name, "Full model REML status", sprintf("REML=%s", reml1),
            if (!reml1) "PASS" else "FAIL",
            "LRT on fixed effects requires REML=FALSE")
  add_check(name, "Saved model REML status (on disk)",
            sprintf("REML=%s", reml_saved),
            if (!reml_saved) "PASS" else "FAIL",
            paste("saved object at", basename(m1_path)))

  p0_fx <- length(lme4::fixef(m0)); p1_fx <- length(lme4::fixef(m1))
  # npar from anova() includes σ² + random-effect variances; the null has an
  # extra (1|thicket) variance because the full moves thicket to fixed effects.
  npar0 <- attr(logLik(m0), "df"); npar1 <- attr(logLik(m1), "df")
  ddf_total <- npar1 - npar0
  ddf_fx    <- p1_fx - p0_fx
  reported_df <- anova_tab$lrt_df[anova_tab$response == name]
  add_check(name, "LRT df matches anova(m0,m1) Df",
            sprintf("anova_df=%d, reported=%d", ddf_total, reported_df),
            if (ddf_total == reported_df) "PASS" else "FAIL",
            sprintf("fixed Δ=%d (null p_fx=%d, full p_fx=%d); total npar Δ=%d (incl. random-effect rebalance)",
                    ddf_fx, p0_fx, p1_fx, ddf_total))

  sing0 <- isSingular(m0); sing1 <- isSingular(m1)
  conv_status <- if (!sing0 && !sing1 &&
                     length(warns0) == 0 && length(warns1) == 0) "PASS"
                 else if (sing0 || sing1) "WARN" else "WARN"
  add_check(name, "Convergence / singular fits",
            sprintf("null sing=%s, full sing=%s, warns=%d",
                    sing0, sing1, length(warns0) + length(warns1)),
            conv_status,
            if (conv_status != "PASS")
              paste(c(warns0, warns1), collapse = " | ")
            else "no warnings")

  lrt <- anova(m0, m1)
  recomputed_p <- lrt$`Pr(>Chisq)`[2]
  reported_p <- anova_tab$lrt_p[anova_tab$response == name]
  p_match <- abs(recomputed_p - reported_p) < 1e-6
  add_check(name, "LRT p-value reproducibility",
            sprintf("recomputed=%.3g, reported=%.3g",
                    recomputed_p, reported_p),
            if (p_match) "PASS" else "WARN",
            "")

  report_lines <<- c(report_lines,
    sprintf("### %s", name),
    sprintf("- n_obs = %d; fixed params: null=%d, full=%d (Δfx=%d). Total npar Δ=%d (reported df=%d).",
            nrow(data), p0_fx, p1_fx, ddf_fx, ddf_total, reported_df),
    sprintf("- REML status: null=%s, full=%s, saved=%s (all should be FALSE)",
            reml0, reml1, reml_saved),
    sprintf("- Convergence: null singular=%s, full singular=%s, warnings=%d",
            sing0, sing1, length(warns0) + length(warns1)),
    sprintf("- LRT χ²(%d) = %.2f, p = %.3g (recomputed %.3g)",
            ddf_total, lrt$Chisq[2], reported_p, recomputed_p))
  list(m0 = m0, m1 = m1, ddf = ddf_total, p = recomputed_p, p_match = p_match,
       sing = c(sing0, sing1))
}

lmer_checks <- list(
  pam_fvfm     = check_lmer_pair("pam_fvfm",     "fv_fm",            pam_d),
  color_dscale = check_lmer_pair("color_dscale", "color_num",        color_d),
  log_zoox     = check_lmer_pair("log_zoox",     "log(cells_per_cm2)", phys_d,
                                  fixed_term = "biopsy_day_c",
                                  include_id = FALSE)
)

# Growth LMM
m_bw_null  <- lme4::lmer(
  areal_calc ~ treatment * wound + thicket + (1 | tank),
  data = bw_d, REML = FALSE,
  control = lme4::lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))
)
m_bw_genet <- lme4::lmer(
  areal_calc ~ treatment * wound * thicket + (1 | tank),
  data = bw_d, REML = FALSE,
  control = lme4::lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4))
)
bw_anova <- anova(m_bw_null, m_bw_genet)
p_bw0 <- length(lme4::fixef(m_bw_null)); p_bw1 <- length(lme4::fixef(m_bw_genet))
ddf_bw <- if ("Chi Df" %in% names(bw_anova)) bw_anova$`Chi Df`[2] else bw_anova$Df[2]
reported_df_bw <- anova_tab$lrt_df[anova_tab$response == "growth_areal"]
add_check("growth_areal", "Model type", "ML LMM + LRT", "INFO",
          "Growth has no time dim; tank retained as random block")
add_check("growth_areal", "df = diff in fixed-effect parameters",
          sprintf("computed=%d, reported=%d", ddf_bw, reported_df_bw),
          if (ddf_bw == reported_df_bw) "PASS" else "FAIL", "")
add_check("growth_areal", "LRT p-value",
          sprintf("p=%.3g", bw_anova$`Pr(>Chisq)`[2]),
          if (bw_anova$`Pr(>Chisq)`[2] < 0.05) "SIG" else "NS", "")
report_lines <- c(report_lines, "",
  "### growth_areal",
  sprintf("- ML LMM LRT (tank random intercept); χ²(%d) = %.2f, p = %.3g",
          ddf_bw, bw_anova$Chisq[2],
          bw_anova$`Pr(>Chisq)`[2]),
  sprintf("- Fixed params: null=%d, full=%d (Δ=%d, reported df=%d)",
          p_bw0, p_bw1, ddf_bw, reported_df_bw))

# Effect direction — does interaction make sense given emmeans?
report_lines <- c(report_lines, "", "## Effect direction vs. emmeans\n")
for (resp in unique(emm_tab$response)) {
  sub <- emm_tab |> filter(response == resp)
  delt <- sub |>
    pivot_wider(id_cols = thicket, names_from = treatment, values_from = mean) |>
    mutate(delta_31_minus_28 = `31C` - `28C`)
  signs <- sign(delt$delta_31_minus_28)
  cross <- length(unique(signs)) > 1
  range_delta <- diff(range(delt$delta_31_minus_28))
  add_check("emmeans", sprintf("%s — direction of heat effect across genets", resp),
            sprintf("Δ range = %.3f, sign mix = %s",
                    range_delta, if (cross) "crossing" else "all same"),
            "INFO",
            paste(sprintf("%s: %+.3f", delt$thicket, delt$delta_31_minus_28),
                  collapse = "; "))
  report_lines <- c(report_lines,
    sprintf("- **%s**: Δ(31C-28C) by genet = %s%s", resp,
            paste(sprintf("%s=%+.3f", delt$thicket, delt$delta_31_minus_28),
                  collapse = ", "),
            if (cross) " — **G×E (sign crossing)**"
            else " — additive (parallel)"))
}

# ----------------------------------------------------------------------------
# WRITE OUTPUTS
# ----------------------------------------------------------------------------
diag_df <- bind_rows(results)
write_csv(diag_df, file.path(DIAG_OUT, "D_pca_lrt_diagnostics.csv"))

# Overall verdicts for report header
pca_status_overall <- if (any(diag_df$status[diag_df$analysis == "PCA"] == "FAIL"))
  "FAIL" else if (any(diag_df$status[diag_df$analysis == "PCA"] == "WARN"))
  "WARN" else "PASS"
lrt_verdict <- function(name) {
  sub <- diag_df |> filter(analysis == name)
  if (any(sub$status == "FAIL")) "FAIL"
  else if (any(sub$status == "WARN")) "WARN" else "PASS"
}

report_lines <- c(report_lines, "", "## Summary verdicts\n",
  sprintf("- PCA: **%s**", pca_status_overall),
  sprintf("- pam_fvfm LRT: **%s**", lrt_verdict("pam_fvfm")),
  sprintf("- color_dscale LRT: **%s**", lrt_verdict("color_dscale")),
  sprintf("- log_zoox LRT: **%s**", lrt_verdict("log_zoox")),
  sprintf("- growth_areal LMM LRT: **%s**", lrt_verdict("growth_areal")))

writeLines(report_lines, file.path(DIAG_OUT, "D_pca_lrt_report.md"))

cat("\n=== PCA and genet-LRT diagnostics complete ===\n")
cat("CSV: ", file.path(DIAG_OUT, "D_pca_lrt_diagnostics.csv"), "\n")
cat("MD:  ", file.path(DIAG_OUT, "D_pca_lrt_report.md"), "\n")
cat("Variance explained:",
    paste(sprintf("PC%d=%.1f%%", seq_along(var_exp), var_exp),
          collapse = ", "), "\n")
cat("PCA verdict:", pca_status_overall, "\n")
cat("pam_fvfm:", lrt_verdict("pam_fvfm"),
    "| color_dscale:", lrt_verdict("color_dscale"),
    "| log_zoox:", lrt_verdict("log_zoox"),
    "| growth_areal:", lrt_verdict("growth_areal"), "\n")
