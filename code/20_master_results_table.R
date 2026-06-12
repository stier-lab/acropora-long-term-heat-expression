# =============================================================================
# Purpose: Aggregate every statistical result in the repo into one tidy
#          spreadsheet that the manuscript can pull from directly. One row
#          per (response, term, contrast, model) with a unified schema so
#          when an analysis updates, the paper text can be updated by
#          re-reading this single CSV.
#
#          Schema:
#            domain            (Physiology | Morphology | Survival | Multivariate)
#            response          PAM Fv/Fm, Color D-scale, Growth %, etc.
#            model_type        LMM | GLMM (binomial) | LM | Cox PH | PCA | LRT
#            term              Treatment, Wound, Day, Treatment:Wound, ...
#                                or contrast label (e.g. "28C - 31C | a | wounded")
#            test              ANOVA F | LRT chi-sq | Wald z | t | HR | LRT
#            statistic         numeric test stat
#            df1, df2          degrees of freedom (df1 = numerator/test df,
#                                df2 = denominator if Satterthwaite/F)
#            n                 sample size used for the test
#            estimate          model-scale estimate (effect, mean diff, log HR)
#            units             units / scale ("Fv/Fm units", "D-scale", "%",
#                                "log10 cells cm-2", "log HR", etc.)
#            pct_change        signed % change of 31C relative to 28C
#                                (where computable on natural scale)
#            ci_low, ci_high   95% CI for the estimate
#            p_value           p-value
#            qualitative       short biological-language summary
#            source_script     code file that produced the row
#            source_artifact   CSV / model RDS it came from
#
# Input:   output/tables/*.csv, output/models/*.rds
# Output:  output/tables/20_master_results.csv             — tidy spreadsheet
#          output/tables/20_master_results_paper_ready.csv — manuscript-formatted
# =============================================================================

source(here::here("code", "00_setup.R"))

# ---- Helpers --------------------------------------------------------------
fmt_p   <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}
fmt_num <- function(x, d = 2) {
  ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = d))
}

qual_dir <- function(estimate, p, response_label,
                     positive = "higher", negative = "lower",
                     ref = "31 °C vs 28 °C") {
  if (is.na(estimate) || is.na(p)) return(NA_character_)
  if (p >= 0.05) return(paste0("no detectable difference (", ref, ")"))
  direction <- if (estimate > 0) positive else negative
  paste0(response_label, " ", direction, " under ", ref)
}

response_label_map <- c(
  pam_fvfm         = "PAM Fv/Fm",
  color_dscale     = "Color (Siebeck D)",
  growth_areal     = "Areal calcification (mg cm-2 d-1)",
  growth_pct       = "Buoyant weight growth (%)",
  log_zoox_density = "log10 symbionts cm-2",
  pct_growth       = "Buoyant weight growth (%)"
)

natural_units <- c(
  pam_fvfm         = "Fv/Fm",
  color_dscale     = "D-scale units",
  growth_areal     = "mg cm-2 d-1",
  growth_pct       = "%",
  log_zoox_density = "log10 cells cm-2",
  pct_growth       = "%"
)

baseline_means <- function() {
  # Compute mean of the response at 28C (end of experiment) so the per-genet
  # treatment effects can be expressed as % change. Used for pct_change column.
  pam  <- readRDS(file.path(DATA_PROC, "pam_clean.rds"))
  col  <- readRDS(file.path(DATA_PROC, "color_clean.rds"))
  bw   <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds"))
  zx   <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds"))

  list(
    pam_fvfm     = pam |>
      filter(day == max(day, na.rm = TRUE), treatment == "28C") |>
      group_by(thicket, wound) |>
      summarise(baseline = mean(fv_fm, na.rm = TRUE), .groups = "drop"),
    color_dscale = col |>
      filter(day == max(day, na.rm = TRUE), treatment == "28C") |>
      group_by(thicket, wound) |>
      summarise(baseline = mean(color_num, na.rm = TRUE), .groups = "drop"),
    growth_areal = bw |>
      filter(treatment == "28C", is.finite(areal_calc)) |>
      group_by(thicket, wound) |>
      summarise(baseline = mean(areal_calc, na.rm = TRUE), .groups = "drop"),
    growth_pct   = bw |>
      filter(treatment == "28C") |>
      group_by(thicket, wound) |>
      summarise(baseline = mean(pct_growth, na.rm = TRUE), .groups = "drop"),
    log_zoox_density = zx |>
      filter(is.finite(cells_per_cm2), cells_per_cm2 > 0) |>
      filter(biopsy_day == max(biopsy_day, na.rm = TRUE),
             treatment == "28C") |>
      group_by(thicket, wound) |>
      summarise(baseline = log10(mean(cells_per_cm2, na.rm = TRUE)),
                .groups = "drop")
  )
}

bases <- baseline_means()

# ===========================================================================
# Block 1 — Continuous LMM ANOVA (script 12)
# ===========================================================================
# Block 1 covers continuous responses only. Morph rows in 12_anova are
# duplicates of those in 04_morphology_trait_anova_genet.csv (same model,
# different script that writes the table) — drop them here and pick up
# the morphology ANOVAs in Block 4 below.
anova12 <- read_csv(file.path(TBL_DIR, "12_anova_summary.csv"),
                    show_col_types = FALSE) |>
  filter(!grepl("^morph_", response_id),
         !grepl("^morph_", response)) |>
  transmute(
    domain          = "Physiology",
    response        = response_label_map[response] |> coalesce(response),
    model_type      = if_else(response_id %in% c("growth_pct", "growth_areal"),
                              "LM", "LMM"),
    term            = term,
    test            = if_else(is.na(`F value`), "Wald chi-sq", "ANOVA F"),
    statistic       = coalesce(`F value`, Chisq),
    df1             = coalesce(NumDF, Df),   # numerator/test df (lmerTest NumDF; car Df)
    df2             = DenDF,                 # denominator df (Satterthwaite) for F-tests
    n               = n_obs,
    estimate        = NA_real_,
    units           = NA_character_,
    pct_change      = NA_real_,
    ci_low          = NA_real_,
    ci_high         = NA_real_,
    p_value         = coalesce(`Pr(>F)`, `Pr(>Chisq)`),
    qualitative     = paste0(term, " effect on ", response,
                             if_else(coalesce(`Pr(>F)`, `Pr(>Chisq)`) < 0.05,
                                     " (significant)", " (n.s.)")),
    source_script   = "code/12_extended_stats.R",
    source_artifact = "output/tables/12_anova_summary.csv"
  )

# ===========================================================================
# Block 2 — Per-genet treatment contrasts at end of experiment (script 12)
# ===========================================================================
genet_eff <- read_csv(file.path(TBL_DIR, "12_genet_treatment_effects.csv"),
                      show_col_types = FALSE)

# Attach baseline and compute pct change
attach_pct <- function(df, response_key) {
  base_df <- bases[[response_key]]
  if (is.null(base_df)) {
    df$pct_change <- NA_real_
    return(df)
  }
  df |>
    left_join(base_df, by = c("thicket", "wound")) |>
    mutate(pct_change = if (response_key == "growth_pct") {
                          -estimate
                        } else {
                          -estimate / baseline * 100
                        })
}

genet_eff_pct <- bind_rows(
  attach_pct(filter(genet_eff, response == "pam_fvfm"),       "pam_fvfm"),
  attach_pct(filter(genet_eff, response == "color_dscale"),   "color_dscale"),
  attach_pct(filter(genet_eff, response == "growth_areal"),   "growth_areal"),
  attach_pct(filter(genet_eff, response == "log_zoox_density"),"log_zoox_density")
)

genet_rows <- genet_eff_pct |>
  transmute(
    domain          = "Physiology",
    response        = response_label_map[response] |> coalesce(response),
    model_type      = if_else(response %in% c("growth_pct", "growth_areal"),
                              "LM", "LMM"),
    term            = sprintf("contrast: 28C - 31C | genet=%s | wound=%s",
                              thicket, wound),
    test            = if_else(!is.na(z.ratio), "Wald z", "Satterthwaite t"),
    statistic       = coalesce(t.ratio, z.ratio),
    df1             = df,
    df2             = NA_real_,
    n               = NA_real_,
    estimate        = estimate,
    units           = natural_units[response] |> coalesce(""),
    pct_change      = pct_change,
    ci_low          = estimate - 1.96 * SE,
    ci_high         = estimate + 1.96 * SE,
    p_value         = p.value,
    qualitative     = sprintf(
      "Genet %s, %s: %s under 31C %s",
      thicket, if_else(wound == "yes", "wounded", "unwounded"),
      if_else(p.value < 0.05,
              if_else(estimate > 0, "lower", "higher"),
              "no change"),
      if_else(!is.na(pct_change), sprintf("(%+.1f%%)", pct_change), "")
    ),
    source_script   = "code/12_extended_stats.R",
    source_artifact = "output/tables/12_genet_treatment_effects.csv"
  )

# ===========================================================================
# Block 3 — R² per response (script 12)
# ===========================================================================
r2_rows <- read_csv(file.path(TBL_DIR, "12_r2_summary.csv"),
                    show_col_types = FALSE) |>
  pivot_longer(c(R2_marginal, R2_conditional),
               names_to = "metric", values_to = "value") |>
  transmute(
    domain          = "Physiology",
    response        = response_label_map[response] |> coalesce(response),
    model_type      = "LMM",
    term            = metric,
    test            = "R-squared",
    statistic       = value,
    df1             = NA_real_, df2 = NA_real_, n = NA_real_,
    estimate        = value, units = "proportion",
    pct_change      = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
    p_value         = NA_real_,
    qualitative     = sprintf("%s = %.3f", metric, value),
    source_script   = "code/12_extended_stats.R",
    source_artifact = "output/tables/12_r2_summary.csv"
  )

# ===========================================================================
# Block 4 — Morphology GLMM ANOVA (script 04, integrated with genet)
# ===========================================================================
morph_anova <- read_csv(file.path(TBL_DIR, "04_morphology_trait_anova_genet.csv"),
                        show_col_types = FALSE) |>
  transmute(
    domain          = "Morphology",
    response        = str_to_sentence(gsub("_", " ", trait)),
    model_type      = "GLMM (binomial)",
    term            = term,
    test            = "Wald chi-sq",
    statistic       = Chisq,
    df1             = Df,
    df2             = NA_real_, n = NA_real_,
    estimate        = NA_real_, units = NA_character_,
    pct_change      = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
    p_value         = `Pr(>Chisq)`,
    qualitative     = paste0(term, " effect on ", trait,
                             if_else(`Pr(>Chisq)` < 0.05,
                                     " (significant)", " (n.s.)")),
    source_script   = "code/04_physio_morphology.R",
    source_artifact = "output/tables/04_morphology_trait_anova_genet.csv"
  )

# Primary morphology fixed-effects: pull from the blme (Cauchy-penalized)
# refit so that Wald tests are interpretable under separation. The original
# unpenalized glmer rows are preserved as a separate block flagged
# `model_type = "GLMM (binomial, unpenalized)"`.
morph_fixed <- read_csv(
  file.path(TBL_DIR, "12c_morph_blme_fixed_effects.csv"),
  show_col_types = FALSE
) |>
  filter(term != "(Intercept)") |>
  mutate(separated = std.error > 50 | abs(estimate) > 10) |>
  transmute(
    domain          = "Morphology",
    response        = str_to_sentence(gsub("_", " ", trait)),
    model_type      = "GLMM (binomial, blme Cauchy(0,2.5))",
    term            = paste0("fixed: ", term),
    test            = "Wald z",
    statistic       = statistic,
    df1             = NA_real_, df2 = NA_real_, n = NA_real_,
    estimate        = if_else(separated, NA_real_, estimate),
    units           = if_else(separated, "log-odds [residual separation]",
                              "log-odds"),
    pct_change      = if_else(separated, NA_real_, (exp(estimate) - 1) * 100),
    ci_low          = if_else(separated, NA_real_, estimate - 1.96 * std.error),
    ci_high         = if_else(separated, NA_real_, estimate + 1.96 * std.error),
    p_value         = if_else(separated, NA_real_,
                              2 * pnorm(-abs(statistic))),
    qualitative     = if_else(
      separated,
      sprintf("%s: residual separation after Cauchy(0,2.5) prior — report omnibus chi-sq only", term),
      sprintf("%s: OR=%.2f (%s, blme penalized)", term, exp(estimate),
              if_else(2 * pnorm(-abs(statistic)) < 0.05, "sig", "n.s."))
    ),
    source_script   = "code/12c_morph_blme.R",
    source_artifact = "output/tables/12c_morph_blme_fixed_effects.csv"
  )

# ===========================================================================
# Block 5 — Cox PH hazard ratios (script 14)
# ===========================================================================
cox_rows <- read_csv(file.path(TBL_DIR, "14_cox_hazard_ratios.csv"),
                     show_col_types = FALSE) |>
  filter(is.finite(HR_31_vs28), HR_31_vs28 > 0, HR_31_vs28 < Inf) |>
  transmute(
    domain          = "Survival",
    response        = str_to_sentence(gsub("_", " ", trait)),
    model_type      = "Cox PH",
    term            = sprintf("HR 31C/28C [%s]", scope),
    test            = "Wald z",
    statistic       = z,
    df1             = 1, df2 = NA_real_, n = n,
    estimate        = HR_31_vs28,
    units           = "hazard ratio",
    pct_change      = (HR_31_vs28 - 1) * 100,
    ci_low          = HR_lo,
    ci_high         = HR_hi,
    p_value         = p,
    qualitative     = sprintf(
      "%s: %s under 31C (HR=%.2f, %s)",
      trait,
      case_when(HR_31_vs28 > 1.1 & p < 0.05 ~ "faster onset",
                HR_31_vs28 < 0.9 & p < 0.05 ~ "delayed onset",
                TRUE                         ~ "no detectable change"),
      HR_31_vs28,
      if_else(p < 0.05, "sig", "n.s.")
    ),
    source_script   = "code/14_morphology_kaplan.R",
    source_artifact = "output/tables/14_cox_hazard_ratios.csv"
  )

cox_tt <- if (file.exists(file.path(TBL_DIR, "14c_cox_tt_pigment_genetC.csv"))) {
  read_csv(file.path(TBL_DIR, "14c_cox_tt_pigment_genetC.csv"),
           show_col_types = FALSE) |>
    transmute(
      domain="Survival",
      response=str_to_sentence(gsub("_", " ", trait)),
      model_type="Cox PH (time-varying)",
      term=sprintf("tt(treatment) [%s]", scope),
      test="Wald z",
      statistic=z,
      df1=1, df2=NA_real_, n=n,
      estimate=coef,
      units="log-HR slope vs log(t+1)",
      pct_change=NA_real_,
      ci_low=coef - 1.96 * se, ci_high=coef + 1.96 * se,
      p_value=p,
      qualitative=sprintf("PH-corrected per-genet test (%s): %s. %s",
                          scope,
                          if_else(p < 0.05, "sig", "n.s. — original per-genet HR was inflated by PH violation"),
                          note),
      source_script="code/14_morphology_kaplan.R",
      source_artifact="output/tables/14c_cox_tt_pigment_genetC.csv"
    )
} else tibble()

cox_genet_lrt <- read_csv(file.path(TBL_DIR, "14_cox_genet_LRT.csv"),
                          show_col_types = FALSE)
cox_genet_rows <- bind_rows(
  cox_genet_lrt |>
    transmute(
      domain="Survival",
      response = str_to_sentence(gsub("_", " ", trait)),
      model_type="Cox PH",
      term="LRT: genet main effect",
      test="LRT chi-sq",
      statistic=chisq_genet_main,
      df1=df_genet_main, df2=NA_real_, n=NA_real_,
      estimate=NA_real_, units=NA_character_,
      pct_change=NA_real_, ci_low=NA_real_, ci_high=NA_real_,
      p_value=p_genet_main,
      qualitative=sprintf("Genet main effect on %s healing: %s", trait,
                          if_else(p_genet_main < 0.05, "sig", "n.s.")),
      source_script="code/14_morphology_kaplan.R",
      source_artifact="output/tables/14_cox_genet_LRT.csv"
    ),
  cox_genet_lrt |>
    transmute(
      domain="Survival",
      response = str_to_sentence(gsub("_", " ", trait)),
      model_type="Cox PH",
      term="LRT: genet x treatment",
      test="LRT chi-sq",
      statistic=chisq_genet_trt_int,
      df1=df_genet_trt_int, df2=NA_real_, n=NA_real_,
      estimate=NA_real_, units=NA_character_,
      pct_change=NA_real_, ci_low=NA_real_, ci_high=NA_real_,
      p_value=p_genet_trt_int,
      qualitative=sprintf("Genet x treatment on %s healing: %s", trait,
                          if_else(p_genet_trt_int < 0.05, "sig", "n.s.")),
      source_script="code/14_morphology_kaplan.R",
      source_artifact="output/tables/14_cox_genet_LRT.csv"
    )
)

# ===========================================================================
# Block 6 — Genet LRT comparisons (script 13)
# ===========================================================================
lrt13 <- read_csv(file.path(TBL_DIR, "13_genet_anova.csv"),
                  show_col_types = FALSE) |>
  transmute(
    domain          = "Physiology",
    response        = response_label_map[response] |> coalesce(response),
    model_type      = "LRT (null vs genet x trt)",
    term            = "LRT: adding genet x treatment x wound x day",
    test            = "LRT chi-sq",
    statistic       = lrt_chisq,
    df1             = lrt_df, df2 = NA_real_, n = n_obs,
    estimate        = delta_aic,
    units           = "delta AIC",
    pct_change      = NA_real_,
    ci_low          = NA_real_, ci_high = NA_real_,
    p_value         = lrt_p,
    qualitative     = sprintf(
      "Adding genet x trt improves fit by dAIC=%.1f (%s)", delta_aic,
      if_else(lrt_p < 0.05, "sig", "n.s.")
    ),
    source_script   = "code/13_genet_interaction.R",
    source_artifact = "output/tables/13_genet_anova.csv"
  )

# ===========================================================================
# Block 7 — PCA + genet displacement (script 15)
# ===========================================================================
pca_load <- read_csv(file.path(TBL_DIR, "15_pca_loadings.csv"),
                     show_col_types = FALSE) |>
  pivot_longer(starts_with("PC"), names_to = "PC", values_to = "loading") |>
  transmute(
    domain="Multivariate", response="End-of-experiment physiology PCA",
    model_type="PCA",
    term=sprintf("loading[%s] -> %s", variable, PC),
    test="PCA loading",
    statistic=loading,
    df1=NA_real_, df2=NA_real_, n=NA_real_,
    estimate=loading, units="loading",
    pct_change=NA_real_, ci_low=NA_real_, ci_high=NA_real_,
    p_value=NA_real_,
    qualitative=sprintf("%s loads %.2f on %s", variable, loading, PC),
    source_script="code/15_multivariate.R",
    source_artifact="output/tables/15_pca_loadings.csv"
  )

pca_disp <- read_csv(file.path(TBL_DIR, "15_genet_pca_displacement.csv"),
                     show_col_types = FALSE) |>
  transmute(
    domain="Multivariate", response="Per-genet thermal displacement (PCA)",
    model_type="PCA",
    term=sprintf("genet=%s", thicket),
    test="Euclidean displacement",
    statistic=displacement,
    df1=NA_real_, df2=NA_real_, n=NA_real_,
    estimate=displacement, units="PC units",
    pct_change=NA_real_, ci_low=NA_real_, ci_high=NA_real_,
    p_value=NA_real_,
    qualitative=sprintf("Genet %s shifts %.2f PC-units under heat", thicket,
                        displacement),
    source_script="code/15_multivariate.R",
    source_artifact="output/tables/15_genet_pca_displacement.csv"
  )

# ===========================================================================
# Block 8 — Buoyant weight LM (script 05 — single fixed-effect table)
# ===========================================================================
bw_lm <- read_csv(file.path(TBL_DIR, "05_buoyant_weight_lm.csv"),
                  show_col_types = FALSE) |>
  filter(term != "(Intercept)") |>
  transmute(
    domain="Physiology", response="Buoyant weight growth (%)",
    model_type="LM",
    term=paste0("fixed: ", term),
    test="t",
    statistic=statistic,
    df1=NA_real_, df2=NA_real_, n=NA_real_,
    estimate=estimate, units="%",
    pct_change=NA_real_, ci_low=conf.low, ci_high=conf.high,
    p_value=p.value,
    qualitative=sprintf("%s: est=%.2f%% (%s)", term, estimate,
                        if_else(p.value < 0.05, "sig", "n.s.")),
    source_script="code/05_buoyant_weight.R",
    source_artifact="output/tables/05_buoyant_weight_lm.csv"
  )

# ===========================================================================
# Block 9 — Color CLMM ordinal robustness (script 12b)
# ===========================================================================
clmm_rows <- if (file.exists(file.path(TBL_DIR, "12b_color_clmm.csv"))) {
  read_csv(file.path(TBL_DIR, "12b_color_clmm.csv"),
           show_col_types = FALSE) |>
    filter(!is.na(LRT)) |>
    transmute(
      domain          = "Physiology",
      response        = "Color (Siebeck D) — ordinal robustness",
      model_type      = "CLMM (cumulative-link mixed)",
      term            = term,
      test            = "LRT chi-sq",
      statistic       = LRT,
      df1             = Df, df2 = NA_real_, n = NA_real_,
      estimate        = NA_real_, units = NA_character_,
      pct_change      = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_value         = `Pr(>Chi)`,
      qualitative     = sprintf("%s (ordinal CLMM): %s", term,
                                if_else(`Pr(>Chi)` < 0.05, "sig", "n.s.")),
      source_script   = "code/12b_color_clmm_robustness.R",
      source_artifact = "output/tables/12b_color_clmm.csv"
    )
} else {
  tibble()
}

# ===========================================================================
# Block 10 — Raw summary stats for the manuscript narrative
#            (cells / growth means, so the prose cites the table, not ad-hoc numbers)
# ===========================================================================
bw_raw <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds")) |>
  group_by(treatment) |>
  summarise(mean_pct = mean(pct_growth, na.rm = TRUE),
            sd_pct   = sd(pct_growth, na.rm = TRUE),
            n        = sum(!is.na(pct_growth)),
            .groups  = "drop")

bw_means_rows <- bw_raw |>
  transmute(
    domain="Physiology", response="Buoyant weight growth (%)",
    model_type="raw summary",
    term=sprintf("mean(%s)", treatment),
    test="mean +/- SD",
    statistic=mean_pct,
    df1=NA_real_, df2=NA_real_, n=n,
    estimate=mean_pct, units="%",
    pct_change=NA_real_, ci_low=mean_pct - sd_pct, ci_high=mean_pct + sd_pct,
    p_value=NA_real_,
    qualitative=sprintf("%s mean growth = %.2f%% (SD %.2f, n=%d)",
                        treatment, mean_pct, sd_pct, n),
    source_script="code/05_buoyant_weight.R (raw)",
    source_artifact="data/processed/buoyant_weight_clean.rds"
  )

bw_pct_drop <- bw_raw |>
  summarise(
    pct_drop_31C_vs_28C = (mean_pct[treatment == "31C"] -
                           mean_pct[treatment == "28C"]) /
                           mean_pct[treatment == "28C"] * 100
  ) |>
  transmute(
    domain="Physiology", response="Buoyant weight growth (%)",
    model_type="derived",
    term="percent reduction 31C vs 28C", test="ratio",
    statistic=pct_drop_31C_vs_28C,
    df1=NA_real_, df2=NA_real_, n=NA_real_,
    estimate=pct_drop_31C_vs_28C, units="%",
    pct_change=pct_drop_31C_vs_28C, ci_low=NA_real_, ci_high=NA_real_,
    p_value=NA_real_,
    qualitative=sprintf("31C corals grew %.0f%% less mass than 28C corals",
                        abs(pct_drop_31C_vs_28C)),
    source_script="code/05_buoyant_weight.R (raw)",
    source_artifact="data/processed/buoyant_weight_clean.rds"
  )

zoox_raw <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
  filter(is.finite(cells_per_cm2), cells_per_cm2 > 0)

zoox_means_rows <- zoox_raw |>
  group_by(treatment, biopsy_day) |>
  summarise(mean_cells = mean(cells_per_cm2, na.rm = TRUE),
            n = n(), .groups = "drop") |>
  filter(biopsy_day %in% c(min(biopsy_day), max(biopsy_day))) |>
  transmute(
    domain="Physiology", response="Symbiont density (cells cm-2)",
    model_type="raw summary",
    term=sprintf("mean(%s, day %s)", treatment, biopsy_day),
    test="mean", statistic=mean_cells,
    df1=NA_real_, df2=NA_real_, n=n,
    estimate=mean_cells, units="cells cm-2",
    pct_change=NA_real_, ci_low=NA_real_, ci_high=NA_real_,
    p_value=NA_real_,
    qualitative=sprintf("%s, biopsy day %s: %.2e cells cm-2 (n=%d)",
                        treatment, biopsy_day, mean_cells, n),
    source_script="code/06_symbiont_chl.R (raw)",
    source_artifact="data/processed/symbiont_chl_clean.rds"
  )

# ===========================================================================
# Block 11 — Time-series diagnostics (script 23): autocorrelation,
#            random slope, nonlinearity-of-time for repeated-measures responses
# ===========================================================================
ts_rows <- if (file.exists(file.path(TBL_DIR, "23_timeseries_diagnostics.csv"))) {
  read_csv(file.path(TBL_DIR, "23_timeseries_diagnostics.csv"),
           show_col_types = FALSE) |>
    transmute(
      domain          = "Time-series diagnostic",
      response        = response,
      model_type      = "diagnostic",
      term            = check,
      test            = "LRT / Wald",
      statistic       = statistic,
      df1             = suppressWarnings(as.numeric(df)), df2 = NA_real_,
      n               = NA_real_,
      estimate        = NA_real_, units = NA_character_,
      pct_change      = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_value         = p_value,
      qualitative     = conclusion,
      source_script   = "code/23_timeseries_diagnostics.R",
      source_artifact = "output/tables/23_timeseries_diagnostics.csv"
    )
} else tibble()

# ===========================================================================
# Block 12 — Cox proportional-hazards tests (script 14, all overall models)
# ===========================================================================
coxph_rows <- if (file.exists(file.path(TBL_DIR, "14_cox_ph_tests.csv"))) {
  read_csv(file.path(TBL_DIR, "14_cox_ph_tests.csv"), show_col_types = FALSE) |>
    transmute(
      domain          = "Survival diagnostic",
      response        = str_to_sentence(gsub("_", " ", trait)),
      model_type      = "cox.zph (Schoenfeld)",
      term            = "proportional-hazards test (treatment)",
      test            = "cox.zph chi-sq",
      statistic       = zph_chisq,
      df1             = 1, df2 = NA_real_, n = n_event,
      estimate        = NA_real_, units = NA_character_,
      pct_change      = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_value         = zph_p,
      qualitative     = if_else(ph_ok, "PH assumption met (p>=0.05)",
                                "PH violated — see time-varying refit"),
      source_script   = "code/14_morphology_kaplan.R",
      source_artifact = "output/tables/14_cox_ph_tests.csv"
    )
} else tibble()

# ===========================================================================
# Block 13 — Thermal context vs Cunning et al. 2024 acute ED50 (script 26)
# ===========================================================================
thermal_rows <- if (file.exists(file.path(TBL_DIR, "26_thermal_context.csv"))) {
  read_csv(file.path(TBL_DIR, "26_thermal_context.csv"), show_col_types = FALSE) |>
    transmute(
      domain          = "Thermal context (external benchmark)",
      response        = "A. pulchra acute Fv/Fm ED50 (Cunning et al. 2024, Mahana)",
      model_type      = "CBASS reference",
      term            = metric,
      test            = "reference value",
      statistic       = value,
      df1 = NA_real_, df2 = NA_real_, n = NA_real_,
      estimate = value, units = "°C or count",
      pct_change = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_value = NA_real_,
      qualitative = "acute thermal-tolerance benchmark; 31C chronic is sublethal (below acute ED50)",
      source_script = "code/26_thermal_context.R",
      source_artifact = "data/external/cunning2024_apulchra_ed50.csv"
    )
} else tibble()

# ===========================================================================
# Block 14 — Inter-milestone lag (script 14): closure -> regeneration
# ===========================================================================
lag_rows <- if (file.exists(file.path(TBL_DIR, "14_milestone_lag_summary.csv"))) {
  read_csv(file.path(TBL_DIR, "14_milestone_lag_summary.csv"),
           show_col_types = FALSE) |>
    transmute(
      domain = "Healing milestone lag",
      response = pair, model_type = "event timing",
      term = paste0(treatment, " (n closed=", n_closed,
                    ", reached both=", n_reached_both, ")"),
      test = "median lag (days)", statistic = median_lag,
      df1 = NA_real_, df2 = NA_real_, n = n_reached_both,
      estimate = median_lag, units = "days",
      pct_change = pct_closed_no_regen, ci_low = iqr_low, ci_high = iqr_high,
      p_value = NA_real_,
      qualitative = sprintf("%s: %.0f%% closed but never regenerated; median lag %s d",
                            treatment, pct_closed_no_regen,
                            ifelse(is.na(median_lag), "NA", as.character(round(median_lag,1)))),
      source_script = "code/14_morphology_kaplan.R",
      source_artifact = "output/tables/14_milestone_lag_summary.csv"
    )
} else tibble()

# ===========================================================================
# Block 15 — Variance partitioning / ICC (script 27)
# ===========================================================================
icc_rows <- if (file.exists(file.path(TBL_DIR, "27_variance_partitioning.csv"))) {
  read_csv(file.path(TBL_DIR, "27_variance_partitioning.csv"),
           show_col_types = FALSE) |>
    transmute(
      domain = "Variance partitioning", response = model,
      model_type = "ICC", term = paste0("ICC[", component, "]"),
      test = "variance fraction", statistic = icc,
      df1 = NA_real_, df2 = NA_real_, n = NA_real_,
      estimate = icc, units = "fraction of variance",
      pct_change = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_value = NA_real_,
      qualitative = sprintf("%s: %.1f%% of variance", component, 100 * icc),
      source_script = "code/27_variance_partitioning.R",
      source_artifact = "output/tables/27_variance_partitioning.csv"
    )
} else tibble()

# ===========================================================================
# Block 16 — Multiple-testing sensitivity (script 28)
# ===========================================================================
mt_rows <- if (file.exists(file.path(TBL_DIR, "28_multiple_testing.csv"))) {
  read_csv(file.path(TBL_DIR, "28_multiple_testing.csv"),
           show_col_types = FALSE) |>
    transmute(
      domain = "Confirmatory / exploratory testing", response = hypothesis,
      model_type = "a priori (raw) / exploratory (BH)", term = test,
      test = "reported p", statistic = p_reported,
      df1 = NA_real_, df2 = NA_real_, n = NA_real_,
      estimate = p_reported, units = "p (raw if a priori, BH if exploratory)",
      pct_change = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_value = p_value,
      qualitative = sprintf("%s; %s (rationale: %s)",
                            hypothesis,
                            ifelse(significant, "significant", "n.s."),
                            rationale),
      source_script = "code/28_multiple_testing.R",
      source_artifact = "output/tables/28_multiple_testing.csv"
    )
} else tibble()

# ===========================================================================
# Block 17 — Morphology probability-scale contrasts (script 29)
# ===========================================================================
probc_rows <- if (file.exists(file.path(TBL_DIR, "29_morphology_prob_contrasts.csv"))) {
  read_csv(file.path(TBL_DIR, "29_morphology_prob_contrasts.csv"),
           show_col_types = FALSE) |>
    filter(is.finite(delta_prob)) |>
    transmute(
      domain = "Morphology", response = str_to_sentence(gsub("_", " ", trait)),
      model_type = "GLMM contrast (prob scale)",
      term = sprintf("28C - 31C at day %d", day),
      test = "delta probability", statistic = delta_prob,
      df1 = NA_real_, df2 = NA_real_, n = NA_real_,
      estimate = delta_prob, units = "probability",
      pct_change = 100 * delta_prob, ci_low = NA_real_, ci_high = NA_real_,
      p_value = p_value,
      qualitative = sprintf("P(28C)=%.2f vs P(31C)=%.2f; OR(31vs28)=%.2g",
                            prob_28, prob_31, or_31_vs_28),
      source_script = "code/29_morphology_prob_contrasts.R",
      source_artifact = "output/tables/29_morphology_prob_contrasts.csv"
    )
} else tibble()

# ===========================================================================
# Combine and write
# ===========================================================================
master <- bind_rows(anova12, genet_rows, r2_rows,
                    morph_anova, morph_fixed,
                    cox_rows, cox_genet_rows, cox_tt,
                    lrt13, pca_load, pca_disp, bw_lm,
                    clmm_rows, bw_means_rows, bw_pct_drop, zoox_means_rows,
                    ts_rows, coxph_rows, thermal_rows,
                    lag_rows, icc_rows, mt_rows, probc_rows) |>
  mutate(across(c(statistic, estimate, pct_change, ci_low, ci_high, p_value),
                \(x) round(x, 4))) |>
  arrange(domain, response, model_type, term) |>
  # Single human-readable description of WHAT each analysis is — the "what was
  # tested" sentence — so the table reads as: [description] | stat | df | est | p | ...
  mutate(description = sprintf("%s: %s — %s [%s, %s]",
                              domain, response, term, model_type, test)) |>
  relocate(description)

write_csv(master, file.path(TBL_DIR, "20_master_results.csv"))

# ---- Paper-ready formatted version ---------------------------------------
paper_ready <- master |>
  transmute(
    Description   = description,
    Domain        = domain,
    Response      = response,
    Model         = model_type,
    Term          = term,
    Test          = test,
    Statistic     = fmt_num(statistic, 2),
    df            = if_else(!is.na(df2),
                            sprintf("%s, %s", fmt_num(df1,0), fmt_num(df2,1)),
                            fmt_num(df1, 0)),
    N             = fmt_num(n, 0),
    Estimate      = if_else(!is.na(units),
                            paste(fmt_num(estimate, 3), units),
                            fmt_num(estimate, 3)),
    `% change`    = if_else(!is.na(pct_change),
                            sprintf("%+.1f%%", pct_change),
                            NA_character_),
    `95% CI`      = if_else(!is.na(ci_low),
                            sprintf("[%.2f, %.2f]", ci_low, ci_high),
                            NA_character_),
    p             = fmt_p(p_value),
    Interpretation = qualitative,
    Source         = source_script
  )
write_csv(paper_ready, file.path(TBL_DIR, "20_master_results_paper_ready.csv"))

cat("\n=== Master results table ===\n")
cat("Rows:", nrow(master), "\n")
cat("By domain:\n")
print(table(master$domain))
cat("\nBy model type:\n")
print(table(master$model_type))
cat("\nWrote:\n",
    " - output/tables/20_master_results.csv             (tidy, programmatic)\n",
    " - output/tables/20_master_results_paper_ready.csv (formatted)\n")
