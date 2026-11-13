# =============================================================================
# 06_genetics.R: final genetic clusters and the RNA-seq handoff table
#
# Builds one row per RNA-seq library (144 libraries; a "library" is one
# sequenced tissue sample) with its temperature, wound status, tank, day and
# symbiont density. It then joins Rachael Bay's final host clusters
# (2026-10-12; groups of fragments that are genetically similar, from DNA
# markers called SNPs) and her rough symbiont read-mapping shares (the share of
# symbiont reads that map to reference groups A-D; only group A is confirmed as
# Symbiodinium). Every join is checked, and a symbiont density model is fitted
# within clusters. That model is exploratory: clusters are unevenly spread
# across treatments.
#
# Until 2026-11-12 the library table also carried source-level resilience
# scores from a frozen exploratory table in extras/; those composites were
# retired and are no longer joined.
#
# Sections
#   1. RNA-seq library table
#   2. Final host clusters: read and check
#   3. Join clusters to the library table and audit the joins
#   4. Design balance and cluster summaries
#   5. Exploratory symbiont model within clusters
#   6. Figures
#   7. Numbers for the results text
#
# Reads:  data/raw/plate_layout/Selected_Samples.csv,
#         data/raw/rnaseq/FINAL_LTH_genoclusters_syms.csv,
#         data/raw/rnaseq/PRELIM_LTH_genoclusters.csv,
#         data/processed/*_clean.rds
# Writes: output/tables/31_rnaseq_*.csv, output/tables/43_final_snp_*.csv,
#         output/models/43_final_snp_*.rds,
#         figures/43_final_snp_{structure,design,symbiont_density,
#           symbiont_typing,diagnostics},
#         output/paper/genetics_numbers.csv
# Run:    sourced by code/run_all.R (or: Rscript code/run_all.R)
# =============================================================================

source(here::here("code", "00_setup.R"))


# ---- 1. RNA-seq library table ----------------------------------------------
# One row per gene-expression library, with the treatment information that the
# expression analysis needs as covariates. The raw lookup keeps the plate
# layout's own column names so it can be checked against the lab sheet.
sel_raw <- read_csv(file.path(DATA_RAW, "plate_layout", "Selected_Samples.csv"),
                    show_col_types = FALSE) |>
  filter(SampleType == "gene")

raw_lookup <- sel_raw |>
  transmute(library_id  = Sample_ID,
            fragment_id = Fragment_ID,
            Temp_C, Day, Wound, Tank, Genotype, Plate, WoundOrder)
write_csv(raw_lookup, file.path(TBL_DIR, "31_rnaseq_library_lookup_raw.csv"))

# Recode to the names used everywhere else in the pipeline.
libs <- sel_raw |>
  transmute(
    library_id = Sample_ID,
    id         = as.integer(Fragment_ID),
    treatment  = paste0(as.integer(Temp_C), "C"),     # 28C / 31C
    wound      = if_else(Wound == "W", "yes", "no"),  # U/W -> no/yes
    # `genet` is a legacy name: it holds the source patch, A/C/D -> a/c/d
    genet      = str_to_lower(Genotype),
    tank       = as.integer(Tank),
    day        = as.integer(Day),
    plate      = as.integer(Plate)
  )

# Add each fragment's symbiont density. distinct() keeps one row per fragment
# so the join cannot duplicate libraries.
sym <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
  distinct(id, .keep_all = TRUE) |>
  transmute(id, symbiont_cells_per_cm2 = cells_per_cm2)

covariates <- libs |>
  left_join(sym, by = "id") |>
  arrange(plate, library_id)
write_csv(covariates, file.path(TBL_DIR, "31_rnaseq_phenotype_covariates.csv"))


# ---- 2. Final host clusters: read and check --------------------------------
# Rachael Bay's final file (2026-10-12): one row per library with the host
# cluster, the first four principal-component axes of the SNP data (PCs:
# summary axes of genetic variation), sequencing coverage, and rough symbiont
# read-mapping shares for groups A-D.
# name_repair = "minimal" keeps the unnamed first column (a row index), which
# select(-1) then drops.
raw <- read_csv(file.path(DATA_RAW, "rnaseq", "FINAL_LTH_genoclusters_syms.csv"),
                show_col_types = FALSE, name_repair = "minimal") |>
  select(-1) |>
  janitor::clean_names()

# Guarantee every column used below is present in the file.
stopifnot(all(c("sample", "id", "pc1", "pc2", "pc3", "pc4", "thicket", "cov",
                "cluster", "a", "b", "c", "d", "total_ab", "dom_sym") %in% names(raw)))

# Pull the fragment id, day and wound token out of the sample labels with
# regular expressions (e.g. "^[0-9]+" = the leading digits). raw$id keeps the
# original library label, because `id` is overwritten on the line above it.
snp <- raw |>
  transmute(
    id = as.integer(str_extract(id, "^[0-9]+")),
    sample_id = as.integer(str_extract(sample, "^[0-9]+")),
    snp_library_id = raw$id, snp_sample = sample,
    snp_day = as.integer(str_match(sample, "_D([0-9]+)$")[, 2]),
    snp_label_wound_token = str_match(sample, "^[0-9]+_(.*)_D[0-9]+$")[, 2],
    snp_source_patch = thicket, host_cluster = as.character(cluster),
    snp_pc1 = pc1, snp_pc2 = pc2, snp_pc3 = pc3, snp_pc4 = pc4, snp_cov = cov,
    sym_read_A = a, sym_read_B = b, sym_read_C = c, sym_read_D = d,
    sym_total_ab_reported = total_ab, sym_dominant_reported = dom_sym,
    host_assignment_status = "final_2026-10-12", sym_typing_status = "rough_read_mapping"
  )

# Guarantee 144 libraries, one per fragment, ids that agree between the two
# label columns, no missing values, and positive coverage.
stopifnot(nrow(snp) == 144, !anyDuplicated(snp$id), all(snp$id == snp$sample_id),
          !anyNA(snp), all(snp$snp_cov > 0))

# Guarantee the read shares are proportions that sum to 1 per library, and that
# the group with the largest share (max.col = column of the row maximum)
# matches Rachael's reported dominant group.
fractions <- as.matrix(select(snp, starts_with("sym_read_")))
stopifnot(all(fractions >= 0 & fractions <= 1), all(abs(rowSums(fractions) - 1) < 1e-8),
          all(c("A", "B", "C", "D")[max.col(fractions)] == snp$sym_dominant_reported))


# ---- 3. Join clusters to the library table and audit the joins -------------
# The library table from section 1 is re-read from disk so the join uses
# exactly the file that is handed off. `genet` is renamed to source_patch here.
lib_tbl <- read_csv(file.path(TBL_DIR, "31_rnaseq_phenotype_covariates.csv"),
                show_col_types = FALSE) |>
  rename(source_patch = genet)

# Guarantee one row per fragment and the same set of fragments in both tables.
stopifnot(!anyDuplicated(lib_tbl$id), setequal(snp$id, lib_tbl$id))

# relationship = "one-to-one": the join stops with an error if any row would
# match more than one row in the other table.
j <- left_join(lib_tbl, snp, by = "id", relationship = "one-to-one")

# Guarantee the two sources agree on source patch, day and wound status
# (unwounded library labels contain "_U_id").
stopifnot(all(j$source_patch == j$snp_source_patch), all(j$day == j$snp_day),
          all(ifelse(str_detect(j$library_id, "_U_id"), "no", "yes") == j$wound))

# Cross-tab of the wound token in Rachael's labels against our wound column.
write_csv(j |> count(snp_label_wound_token, wound),
          file.path(TBL_DIR, "43_final_snp_label_token_audit.csv"))
write_csv(j, file.path(TBL_DIR, "43_final_snp_rnaseq_covariates.csv"))

# Compare final cluster labels with the preliminary ones. Cluster numbers are
# arbitrary labels, so a "changed" label need not mean a changed grouping.
old <- read_csv(file.path(DATA_RAW, "rnaseq", "PRELIM_LTH_genoclusters.csv"),
                show_col_types = FALSE) |>
  janitor::clean_names() |>
  transmute(id = as.integer(str_extract(id, "^[0-9]+")),
            previous_cluster = as.character(cluster))

# Guarantee the preliminary file covers exactly the same fragments.
stopifnot(!anyDuplicated(old$id), setequal(old$id, j$id))

comparison <- j |>
  select(id, host_cluster) |>
  left_join(old, by = "id", relationship = "one-to-one") |>
  mutate(label_changed = host_cluster != previous_cluster)
write_csv(comparison, file.path(TBL_DIR, "43_final_snp_assignment_comparison.csv"))

# How many of these fragments can be matched to each physiology data file?
files <- c("symbiont_chl_clean.rds", "pam_clean.rds", "color_clean.rds",
           "buoyant_weight_clean.rds", "physio_clean.rds")
joinability <- map_dfr(files, function(f) {
  x <- readRDS(file.path(DATA_PROC, f))
  tibble(response_file = f, n_same_fragment_ids = sum(unique(j$id) %in% unique(x$id)))
})
write_csv(joinability, file.path(TBL_DIR, "43_final_snp_response_joinability.csv"))

# One-table summary of the checks above, each value counted from the data.
agree <- j$source_patch == j$snp_source_patch & j$day == j$snp_day &
  ifelse(str_detect(j$library_id, "_U_id"), "no", "yes") == j$wound
write_csv(tibble(check = c("unique libraries", "metadata joins",
                           "patch/day/wound agreement", "fractions sum to one",
                           "dominant label agrees", "changed cluster labels"),
                 value = c(n_distinct(j$library_id), nrow(j), sum(agree),
                           sum(abs(rowSums(fractions) - 1) < 1e-8),
                           sum(c("A", "B", "C", "D")[max.col(fractions)] ==
                                 snp$sym_dominant_reported),
                           sum(comparison$label_changed))),
          file.path(TBL_DIR, "43_final_snp_join_audit.csv"))


# ---- 4. Design balance and cluster summaries -------------------------------
# Clusters were not part of the design, so check how they fall across
# temperature, wound and day. complete() adds the empty cells as n = 0.
balance <- j |>
  count(host_cluster, treatment, wound, day, name = "n") |>
  complete(host_cluster, treatment, wound, day, fill = list(n = 0))
write_csv(balance, file.path(TBL_DIR, "43_final_snp_design_balance.csv"))

cluster_summary <- j |>
  group_by(host_cluster) |>
  summarise(n = n(),
            source_patches = paste(sort(unique(source_patch)), collapse = "/"),
            n_28 = sum(treatment == "28C"), n_31 = sum(treatment == "31C"),
            n_tanks = n_distinct(tank), .groups = "drop")
write_csv(cluster_summary, file.path(TBL_DIR, "43_final_snp_cluster_summary.csv"))

# Dominant symbiont group by treatment cell, and the libraries not dominated by
# the confirmed group A.
write_csv(j |> count(treatment, wound, day, sym_dominant_reported, name = "n"),
          file.path(TBL_DIR, "43_final_snp_symbiont_typing_counts.csv"))
write_csv(j |> filter(sym_dominant_reported != "A"),
          file.path(TBL_DIR, "43_final_snp_non_A_samples.csv"))


# ---- 5. Exploratory symbiont model within clusters -------------------------
# Does heat lower symbiont density within each host cluster? Only clusters with
# at least 3 libraries at each temperature can support a heat contrast.
eligible <- cluster_summary |>
  filter(n_28 >= 3, n_31 >= 3) |>
  pull(host_cluster)
dat <- j |>
  filter(host_cluster %in% eligible) |>
  mutate(cluster = factor(host_cluster), treatment = factor(treatment, c("28C", "31C")),
         wound = factor(wound), day_f = factor(day), tank = factor(tank),
         log_density = log10(symbiont_cells_per_cm2))

# Guarantee no zero or missing densities (log10 of those is not finite).
stopifnot(all(is.finite(dat$log_density)))

## 5a. Fit the model ---------------------------------------------------------
# Linear mixed model: day is a factor (day_f), and (1 | tank) is a random
# intercept for each tank, since fragments in a tank share conditions.
# REML = TRUE: restricted maximum likelihood, the usual fit for estimating
# variance components.
form <- log_density ~ treatment * cluster + treatment * day_f + wound * day_f + (1 | tank)
fit <- function(d) lmerTest::lmer(form, data = d, REML = TRUE)
m <- fit(dat)

# Guarantee no fixed-effect columns were dropped (that would mean some
# treatment x cluster combinations cannot be estimated).
stopifnot(is.null(attr(lme4::getME(m, "X"), "col.dropped")))
saveRDS(m, file.path(MOD_DIR, "43_final_snp_symbiont_lmm.rds"))

## 5b. Heat contrasts within each cluster ------------------------------------
# emmeans gives the model's mean for each temperature within each cluster;
# the pairwise contrast is 28C minus 31C. Holm adjustment controls the chance
# of any false positive across the clusters.
heat_contrasts <- function(m, label)
  as_tibble(summary(contrast(emmeans(m, ~ treatment | cluster), "pairwise"),
                    infer = TRUE)) |>
  mutate(model = label, p_holm = p.adjust(p.value, "holm"))
ct <- heat_contrasts(m, "all eligible samples")

# Sensitivity check: refit using only libraries dominated by symbiont group A.
ma <- fit(filter(dat, sym_dominant_reported == "A"))
saveRDS(ma, file.path(MOD_DIR, "43_final_snp_symbiont_A_only_lmm.rds"))
write_csv(bind_rows(ct, heat_contrasts(ma, "dominant A only")),
          file.path(TBL_DIR, "43_final_snp_symbiont_contrasts.csv"))

## 5c. Model diagnostics -----------------------------------------------------
# Singular = a variance component estimated at zero; also record convergence
# messages, dropped columns and the largest standardised residual.
model_diag <- function(m, label)
  tibble(model = label, n = nobs(m), singular = isSingular(m),
         convergence = paste(m@optinfo$conv$lme4$messages, collapse = "; "),
         fixed_columns_dropped = length(attr(lme4::getME(m, "X"), "col.dropped")),
         max_abs_scaled_residual = max(abs(residuals(m) / sigma(m))))
write_csv(bind_rows(model_diag(m, "all eligible samples"), model_diag(ma, "dominant A only")),
          file.path(TBL_DIR, "43_final_snp_model_diagnostics.csv"))

# Leave-one-tank-out: refit without each tank in turn to see whether one tank
# drives the contrasts. tryCatch() records the error message instead of
# stopping the script if a refit fails.
loo <- map_dfr(levels(dat$tank), function(t) {
  tryCatch({
    mm <- fit(filter(dat, tank != t));
    heat_contrasts(mm, paste("omit tank", t)) |> mutate(singular = isSingular(mm))
  }, error = function(e) tibble(model = paste("omit tank", t),
                                error = conditionMessage(e)))
})
write_csv(loo, file.path(TBL_DIR, "43_final_snp_leave_one_tank_out.csv"))

# Residual spread by treatment x cluster x day: similar SDs support the
# equal-variance assumption.
res <- dat |> mutate(fitted = fitted(m), residual = residuals(m))
write_csv(res |>
            group_by(treatment, cluster, day_f) |>
            summarise(n = n(), sd_residual = sd(residual), .groups = "drop"),
          file.path(TBL_DIR, "43_final_snp_residual_spread.csv"))

# Diagnostic plots: residuals vs fitted (A) and a normal Q-Q plot (B).
q <- ggplot(res, aes(fitted, residual, colour = treatment)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  scale_colour_manual(values = PAL_TEMP) +
  theme_pub(10)
qq <- ggplot(res, aes(sample = residual)) +
  stat_qq() +
  stat_qq_line() +
  theme_pub(10)
save_fig(patch_margins(q + qq + plot_annotation(tag_levels = "A")), "43_final_snp_diagnostics",
         width = 180, height = 90)


# ---- 6. Figures ------------------------------------------------------------
# No titles or legends in the artwork; keys go in the figure captions.

## 6a. Genetic structure -----------------------------------------------------
# A: fragments on the first two SNP axes, coloured by source patch.
# B: number of fragments from each source patch in each final cluster.
counts <- j |>
  count(host_cluster, source_patch) |>
  mutate(host_cluster = factor(host_cluster, levels = as.character(1:8)))
p1 <- ggplot(j, aes(snp_pc1, snp_pc2, colour = source_patch)) +
  geom_point(size = 1.8, alpha = .7) +
  scale_colour_manual(values = PAL_GENO, name = "Source patch", labels = toupper) +
  labs(x = "DNA-marker axis 1", y = "DNA-marker axis 2") +
  theme_pub(10)
p2 <- ggplot(counts, aes(host_cluster, n, fill = source_patch)) +
  geom_col() +
  geom_text(aes(label = ifelse(n >= 6, n, "")), position = position_stack(vjust = .5),
            size = 3) +   # label only bars tall enough to hold a number
  scale_fill_manual(values = PAL_GENO, name = "Source patch", labels = toupper) +
  labs(x = "Final host cluster", y = "Fragments") +
  theme_pub(10)
save_fig(patch_margins(p1 + p2 + plot_annotation(tag_levels = "A")),
         "43_final_snp_structure", width = 190, height = 100)

## 6b. Design balance --------------------------------------------------------
# Libraries per cluster in each day x wound cell, one panel per temperature.
bp <- balance |>
  mutate(cell = factor(paste(day, wound, sep = " / "),
                       levels = c("1 / no", "1 / yes", "3 / no", "3 / yes",
                                  "10 / no", "10 / yes")))
p3 <- ggplot(bp, aes(cell, host_cluster, fill = n)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = n), size = 3) +
  facet_wrap(~ treatment) +
  scale_fill_gradient(low = "white", high = "#66B2A3") +
  labs(x = "Sampling day / clipped (no or yes)", y = "Final host cluster") +
  theme_pub(10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p3, "43_final_snp_design", width = 190, height = 110)

## 6c. Symbiont typing -------------------------------------------------------
# One stacked bar per fragment showing its read shares in groups A-D, sorted
# within each cluster by group A share. pivot_longer() turns the four share
# columns into rows so ggplot can stack them.
typ <- j |>
  arrange(host_cluster, sym_read_A, id) |>
  mutate(sample_order = row_number()) |>
  pivot_longer(starts_with("sym_read_"), names_to = "type", values_to = "fraction") |>
  mutate(type = str_remove(type, "sym_read_"))
cluster_labels <- setNames(paste0("Cluster ", cluster_summary$host_cluster), cluster_summary$host_cluster)
p4 <- ggplot(typ, aes(sample_order, fraction, fill = type)) +
  geom_col(width = 1) +
  facet_wrap(~ host_cluster, scales = "free_x", ncol = 4,
             labeller = as_labeller(cluster_labels)) +
  scale_fill_manual(values = c(A = "#009E73", B = "#CC79A7",
                               C = "#E69F00", D = "#0072B2"),
                    name = "Mapping group") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Fragment", y = "Relative mapped-read share") +
  theme_pub(10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_fig(p4, "43_final_snp_symbiont_typing", width = 190, height = 135)

## 6d. Symbiont density by cluster -------------------------------------------
# Raw densities by day, temperature and wound status, one panel per cluster.
# The jitter has a fixed seed so points land in the same place every run.
p5 <- ggplot(j, aes(factor(day), symbiont_cells_per_cm2 / 1e6, colour = treatment,
                    shape = wound)) +
  geom_point(position = position_jitterdodge(jitter.width = .12, dodge.width = .5,
                                             seed = 42),
             alpha = .75, size = 1.6) +
  facet_wrap(~ host_cluster, ncol = 4) +
  scale_colour_manual(values = PAL_TEMP) +
  scale_y_log10() +
  labs(x = "Sampling day", y = "Symbionts (million cells per sq. cm; log scale)",
       shape = "Clipped", colour = "Temperature") +
  theme_pub(10)
save_fig(p5, "43_final_snp_symbiont_density", width = 190, height = 130)

# Print the key tables to the console log.
print(cluster_summary); print(table(j$sym_dominant_reported)); print(ct)


# ---- 7. Numbers for the results text ---------------------------------------
# Symbionts by source patch. Open question on patch C: rough read shares
# (Rachael's mapping) for reference group D, in the 28 C samples, where heat has
# not shifted them. Only group A is confirmed (Symbiodinium); confirm with
# Rachael what group D is.
rows <- list()
for (p in sort(unique(j$source_patch))) {
  z <- j |> filter(treatment == "28C", source_patch == p, !is.na(sym_read_D))
  rows[[length(rows) + 1]] <- paper_row(
    paste0("sym_D_over5pct_28C_patch_", toupper(p)),
    paste0("28 C samples from patch ", toupper(p),
           " with more than 5% of symbiont reads in group D (statistic = median D share)"),
    estimate = sum(z$sym_read_D > 0.05), statistic = median(z$sym_read_D), n = nrow(z),
    source = "output/tables/43_final_snp_rnaseq_covariates.csv")
}
write_paper(rows, "genetics")
