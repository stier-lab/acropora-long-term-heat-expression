# =============================================================================
# Purpose: Acropora-eating flatworm presence sanity check (3 dates: 06/07,
#          06/08, 06/12 2025). Worm pressure was a known risk for A. pulchra;
#          this script tallies worm-positive corals per tank/treatment to flag
#          any tank-level contamination that could confound treatment effects.
#
# What & why: Acropora-eating flatworms (AEFW) are a known pest of staghorn
#   corals — an infested tank can lose tissue fast, which would mimic
#   a heat/wound effect and bias the experiment. During the trial we
#   inspected every coral on three dates and recorded the worm count. This is
#   surveillance, not a hypothesis test: we confirm worms were
#   absent (or at most isolated), so any treatment differences reported later
#   reflect temperature and wounding, not a pest outbreak in one tank.
#   The output is a per-tank tally of worm-positive corals over time;
#   all zeros indicates no contamination.
# Input:   data/raw/worm_presence/Sheet1.csv
# Output:  data/processed/worm_clean.rds
#          figures/10_worm_presence.{pdf,png}
#          output/tables/10_worm_summary.csv
# =============================================================================

# 00_setup.R loads packages and shared paths (DATA_RAW, DATA_PROC, TBL_DIR, ...),
# theme_pub(), and save_fig().
source(here::here("code", "00_setup.R"))

# ---- Load ------------------------------------------------------------------
# clean_names() standardises the spreadsheet headers to snake_case so columns
# are referred to predictably below (e.g. worm-count columns become worms_*).
raw <- read_csv(file.path(DATA_RAW, "worm_presence", "Sheet1.csv"),
                show_col_types = FALSE) |>
  janitor::clean_names()

# ---- Clean + reshape -------------------------------------------------------
# Tidy the columns, tag each coral with its treatment, then pivot the three
# per-date worm-count columns into long format (one row per coral × date), which
# is what both the summary table and the faceted figure need.
worms <- raw |>
  # The thicket-ID column header varies slightly; match it by prefix and rename.
  rename(thicket = matches("^thicket")) |>
  mutate(
    thicket = str_to_lower(str_squish(thicket)),   # normalise label whitespace/case
    tank    = as.integer(tank),
    id      = as.integer(id),
    # Map tank -> treatment via the shared plumbing layout (tank_treatment() in
    # 00_setup.R; same assignment used by 08 (Apex), 09 (YSI), and 16 (Fig 1)).
    treatment = tank_treatment(tank)
  ) |>
  # Stack the per-date count columns (worms_06_07, worms_06_08, ...) into one
  # date/count pair so each row is a single coral on a single inspection date.
  pivot_longer(starts_with("worms_"),
               names_to = "date_raw", values_to = "n_worms") |>
  mutate(
    # Recover the real date from the column name (strip the "worms_" prefix,
    # parse the remaining m_d_y). mdy() is lubridate's month-day-year parser.
    date    = mdy(str_remove(date_raw, "^worms_")),
    # Collapse count -> presence/absence (1 if any worms seen, else 0). The raw
    # counts are noisy; presence is the robust indicator of tank infestation.
    present = as.integer(n_worms > 0)
  ) |>
  filter(!is.na(treatment))                          # keep experimental tanks only

saveRDS(worms, file.path(DATA_PROC, "worm_clean.rds"))

# ---- Summarise per tank × date --------------------------------------------
# Count corals inspected and corals worm-positive in each tank on each date
# (plus a percentage), giving the tank-level infestation tally we report.
summary_tbl <- worms |>
  group_by(date, treatment, tank) |>
  summarise(
    n_corals = n(),                                      # corals inspected
    n_with_worms = sum(present, na.rm = TRUE),           # of those, worm-positive
    pct      = 100 * n_with_worms / n_corals,            # infestation rate (%)
    .groups = "drop"
  )
write_csv(summary_tbl, file.path(TBL_DIR, "10_worm_summary.csv"))

# ---- Figure ----------------------------------------------------------------
# One small panel per tank; within each, bars show worm-positive coral counts by
# date, dodged and coloured by treatment. Empty panels (all zeros) indicate no
# infestation; any persistent bar flags a tank as a possible confound.
p_worm <- ggplot(summary_tbl,
                 aes(as.factor(date), n_with_worms,           # date as discrete axis (only 3 dates)
                     fill = treatment, group = treatment)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, alpha = 0.85, colour = "black", linewidth = 0.25) +
  facet_wrap(~ tank, ncol = 4, labeller = label_both) +       # label_both prints "tank: N" on each strip
  scale_fill_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),  # blue=ambient, orange=heated
                    name = "Treatment") +
  labs(x = "Inspection date", y = "Corals with flatworms",
       title = "AEFW presence by tank",
       subtitle = "All zeros = good; persistent positives flag tank contamination") +
  theme_pub(9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))     # angle dates so they don't overlap

save_fig(p_worm, "10_worm_presence", width = 170, height = 110)

cat("Wrote worm_clean.rds, 10_worm_summary.csv, 10_worm_presence.{pdf,png}\n")
