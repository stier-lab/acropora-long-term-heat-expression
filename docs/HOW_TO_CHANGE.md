# How to change something

*Created: 2026-11-14 · Last updated: 2026-11-14*

A walkthrough for a new lab member making a change to this analysis. Read `code/README.md` first for what each script does, and `docs/GLOSSARY.md` for terms.

## Set up once

1. Install R 4.6 and RStudio, and clone the repo.
2. Open `Stier-LTH-expression-by-temperature-2025.Rproj`. renv starts automatically.
3. In the R console run `renv::restore()` and say yes. This installs the exact package versions in `renv.lock`, which matters because some model checks change between versions.
4. Run `Rscript code/run_all.R` from the terminal (or `source("code/run_all.R")` in R). It takes about a minute. Then run `git status`: if nothing but PDFs and `output/session_info.txt` changed, your setup matches.

## The rule

Every number, table and figure comes from the code. Never edit a file in `output/`, `figures/` or `data/processed/` by hand, and never edit `data/raw/`. Change a script, rerun, and commit the script and its outputs together.

## Common changes

**A raw data correction** (for example, Molly re-scores a fragment):
1. Export the corrected sheet from Google Drive into the same `data/raw/<folder>/data.csv`. Keep the old file as `data_before_YYYY-MM-DD.csv` if the change is large.
2. If the correction is a single relabel (like fragment 119), add it to `correct_wound_labels()` in `code/00_setup.R` instead, with a comment saying who decided it and when.
3. Write a line in `docs/DATA_NOTES.md`.

**A new variable or a changed calculation:** edit `code/01_clean_data.R`. That script only cleans; it fits no models. The cleaned tables go to `data/processed/` and every later script reads them from there.

**A change to a model or test:** edit the script for that results section (`03`–`06`). If the change affects a reported number, make sure the "paper numbers" section at the end of the script still writes it to `output/paper/`.

**A figure change:** figures use `theme_pub()` and the colours in `code/00_setup.R`. Keep the lab's figure rules: no titles, legends or captions inside the image (the key goes in the caption text, see `figures/FIGURE_INDEX.md`), letter tags for panels, 10 mm outer margin. Save with `save_fig()` so the PDF and PNG are written together.

## Check your change

1. Rerun everything: `Rscript code/run_all.R`. Some checks stop the run with an error, for example if growth no longer matches the metadata column `percent_growth_bw`.
2. The end of the run prints a count of PASS, HANDLED, WARN and FAIL rows; the full list is in `output/tables/18_validation_summary.csv`. A FAIL does not stop the run, but it must be fixed before committing.
3. Sort real changes from noise:
   ```bash
   python3 tools/check_output_changes.py            # lists real changes vs noise
   python3 tools/check_output_changes.py --restore  # undoes the noise-only files
   ```
4. Look at what's left with `git diff`. Every changed number should be one you expected.

## Update the write-ups

- If a reported number moved, re-render the team summary: open `docs/team_summary/LTH_results_summary.Rmd` in RStudio and click Knit. It reads its numbers from `output/paper/`.
- Add a dated line under "Reporting decisions" in `docs/OUTPUTS_STATUS.md` saying what changed and why.
- Tell Molly which numbers changed, since her Google Doc draft is typed by hand.

## Commit

Commit the code, data and outputs together, with a message that says what changed and why, for example: `Re-score fragment 121 budding per Molly (2026-11-20); budding 12/12 vs 4/12 unchanged`. Ask Adrian before pushing to GitHub. The repo is public.
