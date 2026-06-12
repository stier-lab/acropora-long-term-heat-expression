# =============================================================================
# Purpose: Build the publication Figure 1 — a four-panel narrative figure
#          for the LTH manuscript.
#          Panel A: Apex tank temperature trace (treatment validation)
#          Panel B: PAM Fv/Fm trajectory (whole-colony stress)
#          Panel C: Kaplan-Meier curves for the two diagnostic morphology
#                   traits (wound smoothed vs new corallites — closure
#                   preserved, regeneration impaired)
#          Panel D: PCA biplot (multivariate summary)
# Input:   data/processed/{apex_temperature_daily,pam_clean,physio_clean,
#                          coral_physio_wide}.rds
# Output:  figures/16_manuscript_fig1.{pdf,png}
# =============================================================================

source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages(library(survival))

apex   <- readRDS(file.path(DATA_PROC, "apex_temperature_daily.rds"))
pam    <- readRDS(file.path(DATA_PROC, "pam_clean.rds"))
physio <- readRDS(file.path(DATA_PROC, "physio_clean.rds"))
wide   <- readRDS(file.path(DATA_PROC, "coral_physio_wide.rds"))

# ---- Panel A: tank temperature --------------------------------------------
temp_pa <- apex |>
  mutate(probe = str_squish(probe)) |>
  filter(str_detect(probe, "^Temp\\d+$"),
         is.finite(value_mean), value_mean > 0) |>
  mutate(
    tank = as.integer(str_extract(probe, "\\d+")),
    value_c = if_else(value_mean > 60, (value_mean - 32) * 5 / 9, value_mean),
    treatment = case_when(
      tank %in% c(3, 6, 9, 12)  ~ "28C",
      tank %in% c(4, 5, 10, 11) ~ "31C",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(treatment), value_c > 15, value_c < 40,
         between(date, as_date("2025-05-25"), as_date("2025-06-22")))

# Direct labels for the two-group temperature treatment
last_28 <- temp_pa |> filter(treatment == "28C") |>
  group_by(tank) |> slice_max(date, n = 1) |> ungroup() |>
  slice(1)
last_31 <- temp_pa |> filter(treatment == "31C") |>
  group_by(tank) |> slice_max(date, n = 1) |> ungroup() |>
  slice(1)

pA <- ggplot(temp_pa, aes(date, value_c, group = tank, colour = treatment)) +
  geom_hline(yintercept = c(28, 31), linetype = "dashed",
             colour = "grey55", linewidth = 0.3) +
  geom_line(linewidth = 0.4, alpha = 0.85) +
  geom_point(size = 0.7, alpha = 0.7) +
  annotate("text", x = last_28$date, y = 27.5, label = "28 °C",
           hjust = 1.1, vjust = 1, size = 3, colour = "#0072B2",
           fontface = "bold") +
  annotate("text", x = last_31$date, y = 32, label = "31 °C",
           hjust = 1.1, vjust = -0.1, size = 3, colour = "#D55E00",
           fontface = "bold") +
  scale_colour_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                      name = NULL, guide = "none") +  # direct-labeled
  scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
  labs(x = NULL, y = "Tank T (°C)",
       tag = "A") +   # panel descriptions belong in the figure caption, not subtitles
  theme_pub(9)

# ---- Panel B: PAM trajectory ----------------------------------------------
pam_df <- pam |>
  group_by(day, treatment, wound) |>
  summarise(m = mean(fv_fm, na.rm = TRUE),
            se = sd(fv_fm, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

pB <- ggplot(pam_df, aes(day, m, colour = treatment,
                          fill = treatment, shape = wound)) +
  geom_ribbon(aes(ymin = m - se, ymax = m + se,
                  group = interaction(treatment, wound)),
              alpha = 0.15, colour = NA) +
  geom_line(aes(group = interaction(treatment, wound)),
            linewidth = 0.6) +
  geom_point(size = 1.6) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  # Embed small treatment labels at right edge instead of using a legend
  annotate("text", x = max(pam_df$day) + 0.5, y = 0.575,
           label = "31 °C", colour = "#D55E00", size = 2.7,
           fontface = "bold", hjust = 0) +
  annotate("text", x = max(pam_df$day) + 0.5, y = 0.675,
           label = "28 °C", colour = "#0072B2", size = 2.7,
           fontface = "bold", hjust = 0) +
  coord_cartesian(xlim = c(min(pam_df$day), max(pam_df$day) + 2),
                  clip = "off") +
  scale_colour_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                      guide = "none") +
  scale_fill_manual(values   = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                    guide = "none") +
  scale_shape_manual(values = c(no = 16, yes = 17), name = "Wound",
                     guide = "none") +  # collected from panel D
  labs(x = "Days post-wounding", y = expression(F[v]/F[m]),
       tag = "B") +
  theme_pub(9)

# ---- Panel C: KM curves for the two diagnostic traits --------------------
focus_traits <- c(wound_smoothed = "Wound smoothed (closure)",
                  new_corallites_on_tip = "New corallites on tip (regeneration)")

ph_events <- map_dfr(names(focus_traits), function(tr) {
  physio |>
    filter(wound == "yes", !is.na(day), day >= 0) |>
    mutate(y = .data[[tr]]) |>
    group_by(id, treatment) |>
    arrange(day, .by_group = TRUE) |>
    summarise(
      event_day = {
        first1 <- which(y == 1)[1]
        if (is.na(first1)) max(day, na.rm = TRUE) else day[first1]
      },
      event = as.integer(any(y == 1, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    mutate(trait = unname(focus_traits[tr]))
})

km_curves <- ph_events |>
  group_by(trait, treatment) |>
  group_modify(\(d, k) {
    d <- arrange(d, event_day)
    times <- sort(unique(d$event_day))
    n_at_risk <- sapply(times, \(t) sum(d$event_day >= t))
    n_events  <- sapply(times, \(t) sum(d$event_day == t & d$event == 1))
    surv      <- cumprod(1 - n_events / pmax(n_at_risk, 1))
    tibble(day = c(0, times), cum_event = 1 - c(1, surv))
  }) |> ungroup() |>
  mutate(trait = factor(trait, levels = unname(focus_traits)))

pC <- ggplot(km_curves, aes(day, cum_event,
                            colour = treatment, group = treatment)) +
  geom_step(linewidth = 0.8) +
  facet_wrap(~ trait, ncol = 2) +
  scale_colour_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                      name = NULL, guide = "none") +  # collected below
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(x = "Days post-wounding",
       y = "Corals expressing trait",
       tag = "C") +
  theme_pub(9)

# ---- Panel D: PCA biplot -------------------------------------------------
pca_vars <- c("pam_end", "color_end", "growth_areal", "zoox_end")
pretty   <- c(pam_end      = "Fv/Fm",
              color_end    = "Color",
              growth_areal = "Calcification",
              zoox_end     = "Symbionts")
pca_in   <- wide[, pca_vars, drop = FALSE] |> drop_na()
keep     <- complete.cases(wide[, pca_vars])
groups   <- wide[keep, c("treatment", "wound")]

pca   <- prcomp(pca_in, center = TRUE, scale. = TRUE)
var_e <- summary(pca)$importance[2, ] * 100
scores <- as.data.frame(pca$x) |> bind_cols(groups)
load_a <- as.data.frame(pca$rotation) |>
  tibble::rownames_to_column("variable") |>
  mutate(label = unname(pretty[variable]),
         across(c(PC1, PC2), \(x) x * 2.7))

pD <- ggplot(scores, aes(PC1, PC2, colour = treatment, shape = wound)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.25) +
  stat_ellipse(aes(group = treatment, fill = treatment),
               geom = "polygon", level = 0.68,
               alpha = 0.10, colour = NA, show.legend = FALSE) +
  geom_point(size = 1.8, alpha = 0.9, stroke = 0.4) +
  geom_segment(data = load_a,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(1.5, "mm")),
               colour = "grey25", linewidth = 0.3) +
  ggrepel::geom_text_repel(data = load_a,
            aes(PC1 * 1.13, PC2 * 1.13, label = label),
            inherit.aes = FALSE,
            size = 3, colour = "grey15", fontface = "bold",
            box.padding = 0.2, point.padding = 0.5,
            min.segment.length = Inf, seed = 42) +
  # Direct-label wound status inside panel
  annotate("text", x = min(scores$PC1) + 0.2, y = max(scores$PC2) - 0.3,
           label = "● unwounded   ▲ wounded", colour = "grey25",
           size = 2.6, hjust = 0) +
  scale_colour_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                      guide = "none") +
  scale_fill_manual(values = c(`28C` = "#56B4E9", `31C` = "#D55E00"),
                    guide = "none") +
  scale_shape_manual(values = c(no = 16, yes = 17), guide = "none") +
  labs(x = sprintf("PC1 (%.0f%%)", var_e[1]),
       y = sprintf("PC2 (%.0f%%)", var_e[2]),
       tag = "D") +
  theme_pub(9)

# ---- Compose ---------------------------------------------------------------
# Consolidate legends; add top margin so the title doesn't kiss the panel
# letters; use patchwork to align axes.
top_row <- pA + pB + patchwork::plot_layout(widths = c(1, 1))
bot_row <- pC + pD + patchwork::plot_layout(widths = c(1.4, 1))
fig1    <- (top_row / bot_row) +
  patchwork::plot_layout(heights = c(0.9, 1.1), guides = "collect") +
  patchwork::plot_annotation(
    title = expression("Heat compromises photochemistry, growth, and regeneration in"~italic(A.~pulchra)),
    subtitle = "n = 192 fragments, 3 genets, 8 tanks, 14 days post-wounding",
    theme = theme(
      plot.title    = element_text(size = 11, face = "bold",
                                   margin = margin(b = 4, t = 4)),
      plot.subtitle = element_text(size = 9, colour = "grey30",
                                   margin = margin(b = 10))
    )
  ) &
  theme(legend.position = "bottom",
        legend.box      = "horizontal",
        legend.margin   = margin(t = 4))

save_fig(fig1, "16_manuscript_fig1", width = 220, height = 175)

cat("Wrote 16_manuscript_fig1.{pdf,png}\n")
