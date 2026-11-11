# RNA-seq / SNP Inputs

Last updated: 2026-10-13

This folder holds local RNA-seq-adjacent inputs received from the Bay lab.

## Current: `FINAL_LTH_genoclusters_syms.csv`

Rachael Bay sent this attachment on October 12, 2026, 22:58 UTC, in
`Re: updated LTH summary + Day 15 genotype question` (Gmail message
`1a1325721b8b4c29`). Host SNP clusters are final; symbiont typing is explicitly
rough. Original attachment bytes are retained. SHA-256:
`d70ebfb58e4fcd02c67bdbc7f9b66468b3de888630f8e8acaaea72d7a0d0025d`.

Use `code/analysis/43_genetic_clusters.R` and its `43_final_snp_*` outputs for current
work. All 144 memberships match the September file, but genetic PCs have changed.
The file covers Days 1/3/10, not Day 15. Symbiont A/B/C/D columns sum to one and
are read-mapping shares, not cell proportions. Preserve `total.ab` as reported;
its denominator is undocumented and it is not used as symbiont density. Every
sample label contains `wound`, including 72 controls: join by fragment ID and
use experimental metadata for wound state. See the October 13 provenance note.

## `PRELIM_LTH_genoclusters.csv`

Source: email attachment from Rachael A. Bay (`rbay@ucdavis.edu`) received
2026-09-02 14:06 PDT, subject `Re: LTH data`.

Status: **preliminary**. Rachael described these as preliminary clustering
results from an incomplete SNP set. She noted that a few individuals are assigned
to their own cluster and that some of those individuals also had lower values in
the `cov` column, so not all cluster assignments should be treated as final.
When the full SNP set arrives, replace or supersede this file rather than
treating it as confirmed genet identity.

Historical use (superseded by script 43):

- Read by `code/rnaseq/32_prelim_snp_integration.R`.
- Joined to `output/tables/31_rnaseq_phenotype_covariates.csv` by parsed coral
  fragment ID.
- Used for exploratory summaries of preliminary SNP cluster structure, design
  balance, and same-fragment symbiont density.

Important limits:

- A/C/D remain field source-patch labels, not verified genetic individuals.
- The current SNP file covers the 144 RNA-seq libraries. Those fragment IDs
  match the symbiont-density biopsy data but do not match the PAM, color-card,
  buoyant-weight, gross morphology, or microscope morphology fragment sets.
- Coverage bands created by `code/32` are descriptive only. They are not a Bay
  lab QC threshold.
