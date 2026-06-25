# External reference data

> 🗂️ **Provenance for the Cunning 2024 CBASS ED50 data** · Updated 2026-06-07 · Index: [`README.md`](../../README.md) · used by `code/sensitivity/26_thermal_context.R`.

## `cunning2024_apulchra_ed50.csv`

Per-genet acute thermal-tolerance thresholds (Fv/Fm **ED50**, °C) for *Acropora
pulchra* from **Mahana, Mo'orea**, measured by CBASS rapid acute heat-stress
assay (~18 h: ramp → 3 h hold → ramp down).

**Source:** Cunning R, Matsuda SB, Bartels E, D'Alessandro M, Detmer AR,
Harnay P, Levy J, Lirman D, Moeller HV, Muller EM, Nedimyer K, Pfab F,
Putnam HM (2024). *On the use of rapid acute heat tolerance assays to resolve
ecologically relevant differences among corals.* **Coral Reefs**.
doi:10.1007/s00338-024-02577-7. Data: github.com/jrcunning/CBASS_methods
(`data/classic_cbass/processed/ed50.csv`; collection site "mahana" per
`collection_metadata.csv`, collected 2022-12-03, CBASS 2022-12-04).

**Columns:** `geno` genet ID (1–20); `ed50`/`std.error` ED50 (°C) and SE from
the primary dose-response fit; `ed50.f`/`std.error.f` the refined fit.

**Why it's here (LTH project #17):** Same species and island as the LTH
heat × wound experiment. Used by `code/sensitivity/26_thermal_context.R` to (a) place the
LTH chronic treatments (28 °C, 31 °C) on a calibrated *A. pulchra* thermal-
tolerance axis and (b) compare the magnitude of among-genotype thermal-
tolerance variation across the acute (CBASS) and chronic (LTH) methods.

**Caveats.** ED50 is an **acute** (18 h) threshold; the LTH experiment is
**chronic** (weeks at +3 °C). The two are different exposure regimes and ED50
must not be used as the temperature axis of the chronic experiment — only as a
reference for where the chronic treatments sit relative to acute limits.
The LTH thicket labels (A, C, D) are **arbitrary** and not genotype-matched to
Cunning's numbered genets; an individual-genet correlation requires matching
LTH RNA-seq SNPs to Cunning's `genet_map` (future work).
