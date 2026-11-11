suppressMessages(source(here::here("code", "00_setup.R")))
ctrl <- lme4::lmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
pam   <- readRDS(file.path(DATA_PROC, "pam_clean.rds")) |> mutate(thicket = factor(thicket)) |> filter(day >= 0)
color <- readRDS(file.path(DATA_PROC, "color_clean.rds")) |> mutate(thicket = factor(thicket)) |> filter(day >= 0)
sym   <- readRDS(file.path(DATA_PROC, "symbiont_chl_clean.rds")) |> filter(is.finite(cells_per_cm2), cells_per_cm2 > 0) |>
  mutate(thicket = factor(thicket), day = biopsy_day - 1, y = log(cells_per_cm2))
bw <- readRDS(file.path(DATA_PROC, "buoyant_weight_clean.rds")) |> mutate(thicket = factor(thicket))
lrt <- function(lab, d, y, re, f0rhs, f1rhs) {
  m0 <- suppressWarnings(suppressMessages(lme4::lmer(as.formula(paste(y, "~", f0rhs, re)), d, REML = FALSE, control = ctrl)))
  m1 <- suppressWarnings(suppressMessages(lme4::lmer(as.formula(paste(y, "~", f1rhs, re)), d, REML = FALSE, control = ctrl)))
  a <- anova(m0, m1); cat(sprintf("%-55s chisq=%6.1f df=%2d p=%.3g  sing0=%s sing1=%s\n", lab, a$Chisq[2], a$Df[2], a$`Pr(>Chisq)`[2], lme4::isSingular(m0), lme4::isSingular(m1)))
}
for (s in list(list("Fv/Fm", pam, "fv_fm"), list("color", color, "color_num"), list("symbionts", sym, "y"))) {
  re0 <- if (s[[1]] == "symbionts") "+ (1 | tank)" else "+ (1 | tank) + (1 | id)"
  re1 <- if (s[[1]] == "symbionts") "+ (1 + day || tank)" else "+ (1 + day || tank) + (1 + day || id)"
  cat("==", s[[1]], "\n")
  lrt("  REPORTED: all patch interactions, intercepts only", s[[2]], s[[3]], re0, "treatment*wound*day + thicket", "treatment*wound*day*thicket")
  lrt("  all patch interactions, day slopes", s[[2]], s[[3]], re1, "treatment*wound*day + thicket", "treatment*wound*day*thicket")
  lrt("  only patch x temperature terms, day slopes", s[[2]], s[[3]], re1,
      "treatment*wound*day + thicket*wound*day", "treatment*wound*day + thicket*wound*day + thicket:treatment + thicket:treatment:day")
}
cat("== growth\n")
lrt("  REPORTED-style: all patch interactions", bw, "pct_growth", "+ (1 | tank)", "treatment*wound + thicket", "treatment*wound*thicket")
lrt("  only patch x temperature", bw, "pct_growth", "+ (1 | tank)", "treatment*wound + thicket*wound", "treatment*wound + thicket*wound + thicket:treatment")
