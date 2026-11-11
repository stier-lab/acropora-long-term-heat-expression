# Figure 2: the same temperature columns and wound symbols for three responses.
# Descriptive uncertainty uses four tank means, not independent coral replicates.
source(here::here("code", "00_setup.R"))

pam <- readRDS(file.path(DATA_PROC, "pam_clean.rds")) |>
  transmute(id, tank, treatment, wound, day, value = fv_fm)
sym <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |>
  transmute(id, tank, treatment, wound, day = biopsy_day, value = cells_per_cm2 / 1e6)
growth <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds")) |>
  filter(biopsy_day == 15, is.finite(pct_growth)) |>
  transmute(id, tank, treatment, wound, day = 15L, value = pct_growth)
stopifnot(nrow(pam) == 336, nrow(sym) == 192, nrow(growth) == 48,
          !anyDuplicated(pam[c("id", "day")]), !anyDuplicated(sym$id),
          !anyDuplicated(growth$id))

raw <- bind_rows(PAM = pam, Symbionts = sym, Growth = growth, .id = "response")
stopifnot(all(is.finite(raw$value)), all(!is.na(raw$wound)))
tanks <- raw |>
  group_by(response, treatment, wound, day, tank) |>
  summarise(value = mean(value), n_fragments = n(), .groups = "drop")
means <- tanks |>
  group_by(response, treatment, wound, day) |>
  summarise(mean = mean(value), se = sd(value) / sqrt(n()), n_tanks = n(),
            n_fragments = sum(n_fragments), .groups = "drop")
stopifnot(all(means$n_tanks == 4), all(tanks$n_fragments >= 2),
          all(tanks$n_fragments <= 4))
write_csv(tanks, file.path(TBL_DIR, "45_figure2_tank_means.csv"))
write_csv(means, file.path(TBL_DIR, "45_figure2_summary.csv"))

common <- function(lines = TRUE) list(
  facet_wrap(~treatment, ncol = 2,
             labeller = as_labeller(c(`28C` = "28 \u00b0C", `31C` = "31 \u00b0C"))),
  scale_colour_manual(values = PAL_TEMP, guide = "none"),
  scale_shape_manual(values = c(no = 1, yes = 16), name = NULL,
                     labels = c(no = "Unwounded", yes = "Wounded")),
  if (lines) scale_linetype_manual(values = c(no = "dashed", yes = "solid"), name = NULL,
                                  labels = c(no = "Unwounded", yes = "Wounded")),
  theme_pub(11),
  theme(legend.position = "bottom", plot.title = element_text(size = 13, face = "bold"),
        panel.grid.minor = element_blank(), strip.text = element_text(face = "bold")))

trajectory <- function(response_name, title, ylab, breaks) {
  d <- filter(means, response == response_name)
  ggplot(d, aes(day, mean, colour = treatment, shape = wound,
                linetype = wound, group = wound)) +
    geom_vline(xintercept = 0, colour = "grey70", linetype = "dotted") +
    geom_line(linewidth = 0.7) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  width = 0.25, position = position_dodge(width = 0.3),
                  linewidth = 0.45, show.legend = FALSE) +
    geom_point(size = 2.6, stroke = 0.9, position = position_dodge(width = 0.3)) +
    scale_x_continuous(breaks = breaks, limits = c(-1.5, 15.5)) +
    labs(title = title, x = "Day relative to clipping (Day 0)", y = ylab) + common()
}
a <- trajectory("PAM", "A  Photosynthesis efficiency", expression(F[v]/F[m]),
                c(-1, 0, 3, 6, 9, 12, 14))
box_style <- list(
  scale_fill_manual(values = c(no = "white", yes = "grey65"), name = NULL,
                    labels = c(no = "Unwounded", yes = "Wounded")),
  common(lines = FALSE))
b <- ggplot(sym, aes(factor(day), value, colour = treatment,
                     fill = wound, group = interaction(day, wound))) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.65,
               outlier.shape = NA, linewidth = 0.5) +
  geom_point(aes(shape = wound), colour = "grey35", size = 1.2,
             position = position_jitterdodge(jitter.width = 0.12,
                                           dodge.width = 0.8, seed = 42)) +
  labs(title = "B  Algal symbiont density", x = "Sampling day",
       y = expression(10^6~cells~cm^{-2})) + box_style
c <- ggplot(growth, aes(wound, value, colour = treatment, fill = wound)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.5) +
  geom_point(aes(shape = wound), colour = "grey35", size = 1.7,
             position = position_jitter(width = 0.12, height = 0, seed = 42)) +
  scale_x_discrete(labels = c(no = "Unwounded", yes = "Wounded")) +
  labs(title = "C  Skeletal growth over 15 days", x = NULL,
       y = "Skeletal mass gain (%)") + box_style + theme(legend.position = "none")

figure2 <- (a / (b + guides(fill = "none", shape = "none")) /
              (c + guides(fill = "none", shape = "none"))) +
  plot_layout(guides = "collect", heights = c(1, 1, 0.95)) &
  theme(legend.position = "bottom")
save_fig(figure2, "45_physiology_figure2", width = 185, height = 235)
save_fig(a, "45a_photosynthesis", width = 185, height = 85)
save_fig(b, "45b_symbiont_boxplots", width = 185, height = 85)
save_fig(c, "45c_growth_boxplots", width = 185, height = 85)
stopifnot(all(raw$wound[raw$id == 119L] == "yes"),
          all(count(growth, treatment, wound)$n == 12L))
cat("Figure 2 verified: 48 repeated PAM fragments, 192 destructive symbiont samples,",
    "48 endpoint growth fragments; four tanks per temperature.\n")
