PAM

# Factors: Wound (yes/no), Treatment (28/31), Thicket (A,C,D), Day (-1 to 14), 
# Tank (1-8), Location (top/bottom). We are mostly interested in FvFm by Wound and Treatment,
# however, thicket could be an issue, Location could also be an issue. We tank and id to 
# be a random factor

pam_data <- read.csv("~/Downloads/PAM log - complete - PAM data.csv")

View(pam_data)
library(dplyr)
library(lme4)
library(lmerTest)  # for p-values
library(ggplot2)

# Rename the column to FvFm
pam_data <- pam_data %>% rename(FvFm = `Fv.Fm`)

# Fixed effects: Treatment, wound
# Random effects: id (nested in tank), different wound types are included in each tank
# Each tank is a replicate of treatment (ie 3)
###############   

# Make sure coding is correct
pam_data$treatment <- factor(pam_data$treatment)
pam_data$wound     <- factor(pam_data$wound)
pam_data$tank      <- factor(pam_data$tank)
pam_data$id        <- factor(pam_data$id)
#pam_data$day       <- as.numeric(pam_data$day)   # numeric

# Run models
model_1 <- lmer(FvFm ~ treatment * wound * day + (1|id), data = pam_data)
model_2 <- lmer(FvFm ~ treatment * wound + (day | tank), data = pam_data)

# Mixed model: treatment * wound fixed, random intercepts for tanks and IDs, with repeated measures over day
model_3 <- lmer(FvFm ~ treatment * wound + day + (1 + day | tank/id), data = pam_data)
summary(model_3)
model_3a <- lmer(FvFm ~ treatment * wound + day + (1 | tank/id), data = pam_data)

## try accounting for autocorrelation by treating day as a random effect (this will lower df), maybe
## add thicket for random effect too

# random intercept only (id nested in tank)
m_int <- lmer(FvFm ~ treatment * wound + day + (1 | tank/id), data = pam_data, REML = FALSE)
# intercepts for tank + id separately (crossed instead of strict nesting)
m_int2 <- lmer(FvFm ~ treatment * wound + day + (1 | tank) + (1 | id), data = pam_data, REML = FALSE)
# random slope of day by id
m_slope <- lmer(FvFm ~ treatment * wound + day + (1 + day | id) + (1 | tank), data = pam_data, REML = FALSE)

m_slope <- lmer(
  FvFm ~ treatment * wound + day + (1 + day | id) + (1 | tank),
  data = pam_data, REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
m_slope_nc <- lmer(
  FvFm ~ treatment * wound + day + (1 + day || id) + (1 | tank),
  data = pam_data, REML = FALSE
)
anova(m_int) # Treatment by day but not wound
anova(m_int2)
anova(m_slope)
anova(m_slope_nc)
anova(m_int, m_slope_nc)  # compare intercept vs slope at id
final_model <- update(m_slope_nc, REML = TRUE)
summary(final_model)

optional_r2_pkg <- "performance"
if (requireNamespace(optional_r2_pkg, quietly = TRUE)) {
  getExportedValue(optional_r2_pkg, "r2")(final_model)
}

anova(m_int, m_slope)
anova(m_int, m_int2)
AIC(m_int, m_int2, m_slope)



summary(model_3a)
anova(model, type = 3)


# Step 1: fit linear model
model <- lmer(FvFm ~ treatment * wound * day + (1|id), data = pam_data)

# 1a: Test assumption of normality

res <- resid(model)
hist(res, main = "Residuals Histogram")
qqnorm(res); qqline(res)

# 1b:homoscedasticity
plot(fitted(model), res,
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, lty = 2)

# 1c:outliers

infl <- influence(model, obs = TRUE)
cd <- cooks.distance(infl)

# Plot Cook's distance **Haven't used this method before, is there a better way?
plot(cd, type = "h",
     main = "Cook's Distance for observations",
     ylab = "Cook's Distance")
abline(h = 4/length(cd), col = "red", lty = 2)  # common cutoff rule

# Step 2
summary(model)

# Type III ANOVA table
anova(model, type = 3)

# Rough interpretation: The rate of change in Fv/Fm is determined by treatment (28 and 31)
# rather than wound (yes and no). 31 changes more rapidly than 28.

# Summarize means and SE
pam_summary <- pam_data %>%
  group_by(day, treatment, wound) %>%
  summarise(mean_FvFm = mean(FvFm, na.rm = TRUE),
            se_FvFm = sd(FvFm, na.rm = TRUE)/sqrt(n()))

# Plot
ggplot(pam_summary, aes(x = day, y = mean_FvFm, 
                       color = factor(treatment), 
                       shape = wound, group = interaction(treatment, wound))) +
  geom_line() +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_FvFm - se_FvFm, ymax = mean_FvFm + se_FvFm), width = 0.2) +
  scale_color_manual(values = c("28" = "cyan3", "31" = "orange")) +
  scale_shape_manual(values = c("yes" = 16, "no" = 1)) +  
  labs(title = "PAM Fluorometry",
       x = "Day", y = "Fv/Fm", 
       color = "Treatment (°C)", 
       shape = "Wound") +
  theme_minimal(base_size = 14) +
  theme(
      plot.title = element_text(hjust = 0.5 ),
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 8),   
    legend.key.size = unit(0.5, "lines")   
  )

# Plot 2, look at tanks

pam_summary <- pam_data %>%
  group_by(day, treatment, wound, tank) %>%
  summarise(mean_FvFm = mean(FvFm, na.rm = TRUE),
            se_FvFm = sd(FvFm, na.rm = TRUE)/sqrt(n()))

# Plot 2
ggplot(pam_summary, aes(x = day, y = mean_FvFm, 
                        color = factor(treatment), 
                        shape = wound, group = interaction(treatment, wound))) +
  geom_line() +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_FvFm - se_FvFm, ymax = mean_FvFm + se_FvFm), width = 0.2) +
  scale_color_manual(values = c("28" = "cyan3", "31" = "orange")) +
  scale_shape_manual(values = c("yes" = 16, "no" = 1)) +  
  facet_wrap(~ tank) +
  labs(title = "PAM Fluorometry",
       x = "Day", y = "Fv/Fm", 
       color = "Treatment (°C)", 
       shape = "Wound") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5 ),
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 8),   
    legend.key.size = unit(0.5, "lines")   
  )

# Plot 3, look at top/bottom

pam_summary <- pam_data %>%
  group_by(day, treatment, wound, location) %>%
  summarise(mean_FvFm = mean(FvFm, na.rm = TRUE),
            se_FvFm = sd(FvFm, na.rm = TRUE)/sqrt(n()))

# Plot 3
ggplot(pam_summary, aes(x = day, y = mean_FvFm, 
                        color = factor(treatment), 
                        shape = wound, group = interaction(treatment, wound))) +
  geom_line() +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_FvFm - se_FvFm, ymax = mean_FvFm + se_FvFm), width = 0.2) +
  scale_color_manual(values = c("28" = "cyan3", "31" = "orange")) +
  scale_shape_manual(values = c("yes" = 16, "no" = 1)) +  
  facet_wrap(~ location) +
  labs(title = "PAM Fluorometry",
       x = "Day", y = "Fv/Fm", 
       color = "Treatment (°C)", 
       shape = "Wound") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5 ),
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 8),   
    legend.key.size = unit(0.5, "lines")   
  )

# Okay so it does look like top and bottom measurements do make a difference. Let's 
# test this statistically. Okay wait but how. Is this right?

model_loc <- lmer(
  FvFm ~ treatment * wound * location * day + (1 | id),
  data = pam_data
)

summary(model_loc)
anova(model_loc)
AIC(model_loc)


# Plot 4 Or should we look at thicket effects

pam_summary <- pam_data %>%
  group_by(day, treatment, wound, thicket) %>%
  summarise(mean_FvFm = mean(FvFm, na.rm = TRUE),
            se_FvFm = sd(FvFm, na.rm = TRUE)/sqrt(n()))
# Plot 4
ggplot(pam_summary, aes(x = day, y = mean_FvFm, 
                        color = factor(treatment), 
                        shape = wound, group = interaction(treatment, wound))) +
  geom_line() +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_FvFm - se_FvFm, ymax = mean_FvFm + se_FvFm), width = 0.2) +
  scale_color_manual(values = c("28" = "cyan3", "31" = "orange")) +
  scale_shape_manual(values = c("yes" = 16, "no" = 1)) +  
  facet_wrap(~thicket) +
  labs(title = "PAM Fluorometry",
       x = "Day", y = "Fv/Fm", 
       color = "Treatment (°C)", 
       shape = "Wound") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5 ),
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 8),   
    legend.key.size = unit(0.5, "lines")   
  )


model_rm <- lmer(
  FvFm ~ treatment * thicket * day +
    (1 | location) +
    (1 | tank) +
    (1 | id),
  data = pam_data,
  REML = TRUE
)

summary(model_rm)
anova(model_rm) #Same idea, treatment differs by time, no difference in treatment
# by wound status

model_thicket <- lmer(
  FvFm ~ treatment * wound * day +
    (1 | location) +
    (1 | tank) +
    (1 | thicket),
  data = pam_data,
  REML = TRUE
)
VarCorr(model_thicket)
summary(model_thicket)
anova(model_thicket)


#Play around with testing thicket

model_no_thicket <- lmer(
  FvFm ~ treatment * wound * day +
    (1 | id) +
    (1 | tank),
  data = pam_data,
  REML = FALSE
)

model_with_thicket <- lmer(
  FvFm ~ treatment * wound * day +
    (1 | id) +
    (1 | tank) +
    (1 | thicket),
  data = pam_data,
  REML = FALSE
)

anova(model_no_thicket, model_with_thicket) # Model with thicket improves AIC
