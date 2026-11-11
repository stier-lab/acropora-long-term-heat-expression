# Independent rerun check, 2026-11-10

`code/_run_all.R` was rerun from a clean copy of commit bd31783 on a separate Linux machine:
- R 4.6.2 with current CRAN binaries (lme4 2.1.0, ggplot2 4.0.3)
- the Mac outputs were produced with R 4.5.2 and lme4 2.0-1

The run finished in 12.7 minutes with no errors. The new outputs were compared with the committed ones; nothing in the repo was changed.

## Reproduced (relative differences of 1e-6 or smaller, usually below 1e-10)

- Physiology: PAM, buoyant weight and growth, symbiont density, Apex temperatures, YSI
- Survival and timing: Cox, Kaplan–Meier and interval-censored models; budding counts (12/12 at 28 °C vs 4/12 at 31 °C)
- Genetics (script 43), trade-off screens, PAM location audit, Figure 2 summary
- 124 tables and figures were byte-identical

## Not reproducible across machines

**Daily binary morphology GLMMs** (blme, scripts 04 and 12; tables `04_morphology_*`, `12_morph_*`, `12c_*`, and the rows of `20_master_results.csv` and `28_multiple_testing.csv` that derive from them):
- A different set of traits fit successfully on each machine. On the Mac, tip_exist and new_corallites_on_tip fit; here those two failed and were archived, while axial_polyp_formation and wound_smoothed fit instead.
- Estimates for the sparse traits swing widely. Pigment over wound changes sign on several terms. The Wald chi-square for tip extension × day is about 6×10^10 on the Mac and 35.6 here.
- The a priori tip-extension treatment test is p < 0.001 ("significant") on the Mac and p = 0.055 here.
- These models were already flagged as unstable (see the README). This run confirms that their p-values depend on the machine and package versions, so they should not be reported. The team summary does not cite them. The manuscript draft should be checked.

**Smaller shifts:**
- The colour ordinal model (`12b_color_clmm*`) changes in its third or fourth digit.
- The Wilcoxon milestone-lag test gives p = 0.061 on the Mac and 0.093 here, probably because of how the R versions handle ties.
- DHARMa checks on the PAM and colour LMMs flip from WARN to FAIL.

## Machine-specific artefacts to fix

- Several tables store absolute file paths (`37_tradeoff_*`, `42_pam_csv_parse_audit.csv`), so they change on any other computer. Write repo-relative paths instead.
- In a non-UTF-8 locale, R writes characters such as °, — and × as byte escapes. Run with a UTF-8 locale (the Mac default).

## Second run: Mac package versions (same day)

The pipeline was rerun with the modelling packages pinned to the Mac versions: lme4 2.0-1, blme 1.0-7, lmerTest 3.1-3, emmeans 2.0.1 and DHARMa 0.4.7. This run used a UTF-8 locale. R stayed at 4.6.2 and Matrix at 1.7-6 (pinning Matrix 1.7-4 crashed R). The run finished in 4.7 minutes, against 12.7 with lme4 2.1.0.

- **Closer than the first run:** the same set of morphology models fit as on the Mac. 130 outputs were byte-identical, and the encoding differences disappeared.
- **Tip extension:** the a priori treatment test is p = 0 again. Its Wald chi-square is still astronomical (21,409,117 here, 20,975 on the Mac; about 9×10^10 vs 6×10^10 for the day term). That is a separation artefact, not evidence of a heat effect.
- **The blme robustness fits and the sparse primary GLMMs still differ across platforms:**
  - Pigment over wound flips sign on several terms.
  - Some primary GLMM estimates for axial polyp formation, wound smoothed and tip exist differ by 30–170% on individual terms.
  - Polyps out matches to about 1e-4.
  - The standard tip-extension GLMM fits here but not on the Mac.
- **Unchanged from the first run:** the ordinal colour model and the milestone-lag Wilcoxon test still differ by the same amounts.

**Conclusion:** the lme4 version changes which daily morphology models converge, and how fast. Even with matched versions, the separated or sparse fits are not numerically stable across platforms. Report the event counts and the Cox, Kaplan–Meier and interval-censored timing results, which reproduce exactly. Do not report p-values from the daily binary GLMMs, including the a priori tip-extension test.
