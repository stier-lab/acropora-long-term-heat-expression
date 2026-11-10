# Final SNP integration: October 13, 2026

## Source and Status

Rachael Bay's October 12 email (22:58:48 UTC), `Re: updated LTH summary + Day 15 genotype question`, Gmail message `1a1325721b8b4c29`, supplied `FINAL_LTH_genoclusters_syms.csv`. Attachment size: 34,403 bytes. SHA-256: `d70ebfb58e4fcd02c67bdbc7f9b66468b3de888630f8e8acaaea72d7a0d0025d`.

Host genotype clusters are final according to Rachael. Symbiont typing is explicitly rough, based on mapping to four reference genomes. Expression reruns are pending at the Bay lab; no counts, raw variant calls, kinship matrix, or Day 15 genotypes were attached. The preliminary attachment remains unchanged for historical comparison. No shared Drive file or email was changed in this integration.

## Data Contract

Script 43 joins by the leading numeric fragment ID in both `sample` and `ID`, checks their agreement, requires 144 unique IDs and a one-to-one exact match to script 31, and verifies source patch and sampling day. The unnamed CSV row index is ignored. All fields used in the join and composition checks are complete.

Every supplied sequencing label includes `wound` / `W`, even the 72 experimentally unwounded controls. This was caught by the initial validation. The label cannot be used as an experimental treatment code. Wound state comes from the established script-31 metadata, not a reassignment based on this token; the 72/72 cross-tab is saved. This nomenclature should be confirmed before collaborators derive treatment variables from filenames.

All 144 host-cluster labels match the preliminary file. Final PCs differ, so retain the new numeric coordinates and do not mix old and new axes. Source C has 48 fragments in cluster 1. A spans 4/5/6/7; D spans 2/3/4/5/8. Cluster sizes: 48,17,16,52,7,2,1,1. No Day 15 phenotype IDs overlap. Same-fragment symbiont density overlaps all 144; PAM, color, growth and main morphology overlap none.

Symbiont columns A/B/C/D are nonnegative shares summing to one; `dom.sym` matches the largest share for every row. Group A corresponds to Symbiodinium in the email. Other letters are retained as supplied reference-group labels. They are not coral source labels or host-genotype numbers. Read shares are not cell fractions, and mixed mapping alone does not establish biological mixed infection. Mapping specificity, reference choice, RNA activity and read depth may affect those shares. No arbitrary purity threshold is imposed.

The meaning/denominator of `total.ab` is not documented in the attachment or email; it is retained as `sym_total_ab_reported` and excluded from biological density inference. The supplied `cov` field is preserved, not converted to a new final quality threshold.

## New Observation

138/144 libraries are A-dominated, five D-dominated, and one C-dominated. All six non-A calls occur in heated, wounded fragments. IDs 251,260,269 (host cluster 1), 186 (cluster 2) and 215 (cluster 8) are Day 10; ID 276 (cluster 1) is Day 1. There are no pre-exposure measurements of these fragments' symbiont mixtures. Do not describe this as observed symbiont switching, treatment causation, or an independent symbiont-type effect.

## Density Analysis and Diagnostics

Only host clusters with at least three fragments at each temperature enter the exploratory model: clusters 1-5, n=140. Clusters 6-8 remain separate in descriptive tables and figures; unlike the old script, they are not pooled into an artificial genotype. The eligibility threshold is a minimum replication screen, not proof of adequate support in every day-by-wound cell.

The model uses log10 symbiont density, temperature-by-cluster, temperature-by-day and wound-by-day terms, categorical day, and a tank random intercept. Source patch is not simultaneously fitted as an independent genotype effect. Heat contrasts average sampling days and wound states equally under the model; empty cells require model assumptions, so observations remain the main figure. These are descriptive/exploratory genotype comparisons, not a fully balanced genetic experiment. Contrast p-values receive Holm correction within each fitted model; intervals are pointwise.

Full-data heat contrasts (28 minus 31 C, log10 units) for clusters 1-5: 0.141,0.285,0.273,0.435,0.194. Cluster 1 remains smaller in estimated loss than clusters 2-4. This is not a formal test of every pairwise difference between heat effects. Restricting to dominant-A samples gives n=135 and contrasts 0.172,0.274,0.265,0.437,0.199. The sixth non-A sample belongs to excluded cluster 8. This sensitivity does not remove all mixed reads or demonstrate a host effect independent of symbiont composition; restricting on a possible treatment response can itself bias inference.

Both models converge, have full fixed-effect rank, and are not singular. Maximum absolute scaled residual is about 3.48. Visual residual/Q-Q inspection shows a heavier lower tail and unequal spread: goodness-of-fit remains partial, not all-pass. Leave-one-tank-out estimates and flags are saved. Main temperature replication remains eight tanks, and cluster 5 occupies only four tanks. Do not elevate exploratory p-values to the paper's main genetic claim.

## Pipeline and Handoff

Run `Rscript code/43_final_snp_integration.R`, `Rscript code/38_analysis_diagnostic_inventory.R`, then render the team-summary Rmd. On this host use `/usr/local/bin/Rscript` (R 4.5.2 with project packages). Script 43 is in the main runner; old script 32 outputs remain historical. Current combined covariates are `output/tables/43_final_snp_rnaseq_covariates.csv`; source-derived phenotype summaries retain `source_` prefixes and must not be treated as individual measurements.

For expression work: use these final cluster/PC fields, verify rank and treatment overlap for each proposed design, do not use source patch and cluster as interchangeable independent predictors, and assess rough symbiont mixture through explicitly labeled sensitivities. Do not automatically adjust the host heat effect for a symbiont measure potentially changed by treatment. The decision to proceed without new Day 15 sequencing remains unchanged.

The team summary replaces its three preliminary SNP panels with final structure, design, observed density, and rough symbiont composition plots. Other phenotype findings are unchanged. It distinguishes finalized host labels from unresolved symbiont typing, endpoint identity, and external Hauru matches.
