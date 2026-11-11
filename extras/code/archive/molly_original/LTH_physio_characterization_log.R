library(tidyverse)

phys_char <- read.csv("~/Downloads/physio characterization log - complete - data.csv")
View(phys_char)

# long format
phys_char_long <- phys_char %>%
  pivot_longer(
    cols = c(hole_in_center, polyp_in_hole, wound_smoothed,
             pigment_over_wound, tip_exist, tip_extension, new_corallites_on_tip),
    names_to = "factor",
    values_to = "status"
  )

# Convert status to logical (assuming "yes"/"no")
phys_char_long <- phys_char_long %>%
  mutate(status = tolower(status) == "yes")

# Find first day where 49% are "yes" for each factor × wound × treatment
summary_phys_char_long <- phys_char_long %>%
  group_by(factor, wounded, treatment, day) %>%
  summarise(prop_yes = mean(status), .groups = "drop") %>%
  filter(prop_yes >= 0.25) %>%
  group_by(factor, wounded, treatment) %>%
  summarise(first_day = min(day), .groups = "drop")

## My thoughts so far... we can play around with the yes:no ratio all day to make
# a nice graph, but what makes the most sense statistically? We should likely follow
# a different approach, possibly one that takes the "tank" factor into account (ie
# not just 75% of all samples yes:no, but how many tanks had samples that were >75%
# yes:no)

# put in order
summary_phys_char_long <- summary_phys_char_long %>%
  filter(factor != "pigment_over_wound") %>%
 mutate(factor = factor(factor,
                         levels = c("hole_in_center",
                                    "polyp_in_hole",
                                    "wound_smoothed",
                                    "tip_exist",
                                    "tip_extension",
                                    "new_corallites_on_tip"),
                        labels = c("Hole in Center",
                                   "Polyp in Hole",
                                   "Wound Smoothed",
                                   "Tip Exist",
                                   "Tip Extension",
                                   "New Corallites on Tip")))


ggplot(summary_phys_char_long, aes(x = first_day, y = factor,
                       color = as.factor(treatment), group = interaction(wounded, treatment))) +
  geom_line(aes(linetype = wounded), size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c("28" = "cyan3", "31" = "orange")) +
  labs(x = "Day (first ≥25% yes)", y = "Factor",
       color = "Treatment", linetype = "Wound") +
  theme_minimal()

# Let's try viewing each coral separately even if it's noisy


# Reshape to long format
phys_char_long_individual <- phys_char %>%
  pivot_longer(
    cols = c(hole_in_center, polyp_in_hole, wound_smoothed,
             tip_exist, tip_extension, new_corallites_on_tip),
    names_to = "factor",
    values_to = "status"
  ) %>%
  mutate(status = tolower(status) == "yes")

phys_char_first_yes <- phys_char_long_individual %>%
  filter(status) %>%
  group_by(id, factor, wounded, treatment) %>%
  summarise(first_day = min(day), .groups = "drop")


# order bottom-to-top
phys_char_first_yes <- phys_char_first_yes %>%
  filter(factor != "pigment_over_wound") %>%
  mutate(factor = factor(factor,
                         levels = c("hole_in_center",
                                    "polyp_in_hole",
                                    "wound_smoothed",
                                    "tip_exist",
                                    "tip_extension",
                                    "new_corallites_on_tip"),
                        labels = c("Hole in Center",
                                    "Polyp in Hole",
                                     "Wound Smoothed",
                                     "Tip Exist",
                                    "Tip Extension",
                                     "New Corallites on Tip")))


ggplot(phys_char_first_yes, aes(x = first_day, y = factor,
                                color = as.factor(treatment), shape = wounded,
                                group = id)) +
  geom_line(size = 0.8, alpha = 0.6) +
  geom_point(size = 3) +
  scale_color_manual(values = c("28" = "cyan3", "31" = "orange")) +
  labs(x = "Day (first yes)", y = "Factor",
       color = "Treatment", shape = "Wound") +
  theme_minimal() +
  facet_wrap(~ id)
#### Coral 121 and 116 look weird... do another QA/QC on this?

## Can we facet_wrap by tank to see if there are any obvious odd-balls?
## Should we also look at differences in thicket


# Calculate proportion "yes" per day × treatment × factor
df_prop <- phys_char_long_individual %>%
  group_by(factor, day, treatment, thicket) %>%
  summarise(prop_yes = mean(status, na.rm = TRUE), .groups = "drop")

# Relabel factors for plotting
df_prop <- df_prop %>%
  mutate(factor = factor(factor,
                         levels = c("hole_in_center",
                                    "polyp_in_hole",
                                    "wound_smoothed",
                                    "tip_exist",
                                    "tip_extension",
                                    "new_corallites_on_tip"),
                         labels = c("Hole in Center",
                                    "Polyp in Hole",
                                    "Wound Smoothed",
                                    "Tip Exist",
                                    "Tip Extension",
                                    "New Corallites on Tip")))

# Plot
ggplot(df_prop, aes(x = day, y = prop_yes,
                    color = as.factor(treatment), group = treatment)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("28" = "cyan3", "31" = "orange")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Characteristics by Treatment",
       x = "Day", y = "Positive for Characteristic",
       color = "Treatment (°C)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 8),        
        legend.title = element_text(size = 9)) +
  facet_wrap(~ factor + thicket, nrow = 6, ncol = 3)

