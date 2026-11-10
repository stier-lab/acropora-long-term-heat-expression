# October 16 scoring refresh

## Sources

Molly's October 16 email in the "updated LTH summary HTML" thread confirmed field-sheet corrections. Retrieved the `data` tab from both native Google Sheets using authenticated Sheets API reads:

- Main experiment: `1_6CCmYHGpocI7q9XDzaJHya5NaR-SqArgd24aOSpuoM`, modified October 16 at 00:14:04.569 UTC.
- Separate microscope experiment: `1rB5VFJqQov0ZXsPjRKnjwVDOtiPxgbfMMDm9l_9Y-Q4`, modified October 16 at 16:27:21.528 UTC.

Current CSV exports remain at each cohort's `data/raw/.../data.csv`. Earlier exports are preserved beside them as `data_before_2026-10-16.csv`. Retrieval metadata and SHA256 hashes are in `data/metadata/scoring_refresh_2026-10-16.json`. No remote source sheet was edited.

## Changes and interpretation

The main experiment retains all 768 ID/day rows and the microscope experiment all 256. Treatment, tank, wound status, source patch, and dates are unchanged. `code/44_corrected_scoring_audit.R` reproduces the cell-level comparison in `output/tables/44_scoring_changes.csv`.

- Tip presence: ID 116 Day 9 changes from no to yes; IDs 121 and 143 Day 8 change from yes to no. The latter two now have no recorded tip formation.
- Tip extension: ID 175 Day 8 changes from yes to no. This reduces heated attainment from 11/12 to 10/12.
- Pigment: 287 cells change. Main-experiment wounded fragments have observations only on Days 3, 4, 9, 12, and 15. Other days remain missing. There are two changed observation notes and 30 added QA notes.
- Some QA notes describe a biological assumption of pigment persistence for ID 121. The actual score cells remain NA on unmeasured days. We use the cells, not these suggested assumptions, and do not carry observations forward.
- The note on ID 115 explains that budding can occur in the wound area before axial-tip extension. Stages are not forced into a strict sequence.
- Microscope IDs 5 and 7 have Day 1 tissue cover filled as yes. All 16 fragments now have observed tissue cover on Day 1.

## Analysis changes

Rerun the analysis pipeline from the revised exports. Code 14 now excludes unscored visits before deriving onset bounds: a missing observation cannot be a known absent score or extend follow-up. Code 44 verifies missing pigment remains missing and reproduces corrected attainment counts. Both cohorts remain separate.

The main budding comparison is unchanged: 12/12 at 28 degrees C and 4/12 at 31 degrees C. Tip formation is 11/12 in each treatment. Tip extension is 11/12 at 28 degrees C and 10/12 at 31 degrees C. First observed pigment occurs in 5/12 and 6/12, respectively. Sparse pigment observations do not establish equivalent biological onset times.

Code 45 creates a candidate publication Figure 2 with consistent temperature panels and wound symbols for photosynthetic efficiency, symbiont density, and endpoint growth. Means and descriptive standard errors use four tanks per temperature. It changes presentation, not the physiology models. The existing physiology metadata gives 13 unwounded and 11 wounded heated endpoint fragments; it is not silently recoded to the morphology cohort's 12/12 design.

The discrepancy is fragment 119: morphology says wounded, master physiology metadata says unwounded. Confirm with Molly before finalizing treatment-specific manuscript panels or a unified wound-status analysis. Figure 2 remains a candidate with source labels preserved.

## Limitations retained

The daily binary models include sparse or separated outcomes, and some penalized robustness fits fail or remain unstable. A saved robustness model alone is no longer enough for the diagnostic script to mark a problem handled: it must match current response data, converge, and have finite coefficient standard errors below the separation screen. Failed refits archive their earlier model artifact rather than leaving it available as a current fit.

Timing-model uncertainty does not fully account for shared tanks. The observed budding counts do not depend on those models. The source-patch heat-penalty composite remains exploratory and its pigment contribution is imprecise.

The trade-off diagnostic table contains 57 passed and seven flagged checks after rerunning. Flags concern residual uniformity/normality or influential observations; the companion robust and leave-one-out analyses remain available. This is not an all-models-passed result.

The reproducibility diagnostic now uses the same standardized day predictor as the penalized morphology fitter. Refit-based residual simulation is explicitly skipped and flagged when the base model already fails convergence or separation checks; repeated fits of an invalid base model are not treated as a remedy.

## Verification

All pipeline analysis stages were rerun. A final diagnostic refresh completed successfully; its log is `output/logs/scoring_refresh_2026-10-16_checks.log`. Code 44's regression checks confirm that unscored visits neither narrow event intervals nor extend follow-up. Code 45 checks fragment counts, duplicate keys, finite values, and four contributing tanks per group.

The HTML was rebuilt and its embedded Figure 2 was verified byte-for-byte against the current PNG. All 24 images and both videos are embedded. The figure PNGs were visually inspected; the candidate PDF is one page with embedded fonts. Browser preview of the local HTML was blocked by browser URL policy, so HTML verification used its generated contents rather than a browser screenshot.

The older `manuscript/Manuscript_LTH.docx` was not rewritten in this task. The manuscript audit correctly flags its old composite heat-penalty numbers; the team summary contains the refreshed values.
