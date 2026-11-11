suppressMessages(source(here::here("code", "00_setup.R"))); suppressMessages(library(lmerTest))
ctrl <- lme4::lmerControl(check.conv.singular = .makeCC("ignore", tol = 1e-4), optimizer = "bobyqa")
pam   <- readRDS(file.path(DATA_PROC, "pam_clean.rds")) |> mutate(thicket = factor(thicket)) |> filter(day >= 0)
color <- readRDS(file.path(DATA_PROC, "color_clean.rds")) |> mutate(thicket = factor(thicket)) |> filter(day >= 0)
sym   <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |> filter(is.finite(cells_per_cm2), cells_per_cm2 > 0) |>
  mutate(thicket = factor(thicket), biopsy_day_c = biopsy_day - 1, y = log(cells_per_cm2))
show <- function(lab, m, keys) { a <- as.data.frame(anova(m, type = 3)); cat(sprintf("%-60s singular=%s ", lab, lme4::isSingular(m)));
  for (k in keys) cat(sprintf("| %s F=%.1f df=%.1f p=%.2g ", k, a[k,"F value"], a[k,"DenDF"], a[k,"Pr(>F)"])); cat("\n") }
k <- c("treatment:day", "treatment:wound:day")
show("pam (1+day||tank)+(1|id)", lmer(fv_fm ~ treatment*wound*day*thicket + (1+day||tank) + (1|id), pam, control = ctrl), k)
show("color (1+day||tank)+(1|id)", lmer(color_num ~ treatment*wound*day*thicket + (1+day||tank) + (1|id), color, control = ctrl), k)
show("color (1+day||tank)+(1+day||id)", lmer(color_num ~ treatment*wound*day*thicket + (1+day||tank) + (1+day||id), color, control = ctrl), k)
ks <- c("treatment:biopsy_day_c", "treatment:wound:biopsy_day_c")
show("sym (1+day||tank)", lmer(y ~ treatment*wound*biopsy_day_c*thicket + (1+biopsy_day_c||tank), sym, control = ctrl), ks)
# Symbiont Day-15 means by group, with and without the 5 outliers
o <- c(161,163,170,174,175)
s15 <- sym |> filter(biopsy_day == 15) |> group_by(treatment, wound) |> summarise(n = n(), geo_mean_all = exp(mean(y)), n_out = sum(id %in% o),
  geo_mean_wo = exp(mean(y[!id %in% o])), .groups = "drop"); print(as.data.frame(s15))
print(as.data.frame(sym |> filter(biopsy_day == 15, treatment == "31C") |> group_by(thicket, wound) |> summarise(n = n(), min = min(cells_per_cm2), median = median(cells_per_cm2), .groups = "drop")))
