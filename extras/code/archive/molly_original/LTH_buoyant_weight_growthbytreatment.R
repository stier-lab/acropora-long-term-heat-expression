


library(tidyverse)
library(dplyr)
library(lme4)
library(lmerTest)  # for p-values
library(ggplot2)

buoy <- read.csv("~/Downloads/Buoyant_weight_calculations - standard_curve_data (use me).csv")
View(`Buoyant_weight_calculations...standard_curve_data.(use.me)`)

buoy$wound <- as.factor(buoy$wound)
buoy$treatment <- as.factor(buoy$treatment)
buoy$tank <- as.factor(buoy$tank)
buoy
str(buoy)


####################
ggplot(pam_summary, aes(x = day, y = mean_FvFm, 
                        color = factor(treatment), 
                        shape = wound, group = interaction(treatment, wound))) +
  geom_line() +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_FvFm - se_FvFm, ymax = mean_FvFm + se_FvFm), width = 0.2) +
  scale_color_manual(values = c("28" = "cyan3", "31" = "orange")) +
  scale_shape_manual(values = c("yes" = 16, "no" = 1)) +  
  facet_wrap(~thicket)
  labs(title = "PAM Fluorometry",
       x = "Day", y = "Fv/Fm", 
       color = "Treatment (°C)", 
       shape = "Wound") +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 8),  
    legend.key.size = unit(0.5, "lines")   
  )
###################

#Boxplot just to visualize a little

ggplot(buoy, aes(x = wound, y = final.initial.initial, fill = treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("28" = "cyan3", "31" = "orange")) +
    facet_wrap(~thicket)+
  labs(
    title = "Change in Growth by Wound x Treatment",
    x = "Wound",
    y = "Percent Growth"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_text(size = 10), # Makes the legend title smaller
    legend.text = element_text(size = 8),   # Makes the legend text (labels) smaller
    legend.key.size = unit(0.5, "lines")   # Makes the legend keys (symbols) smaller
  )

# Whoa, looks fairly clear that temperature has a bigger effect on growth than wound status
# Or another way to state it is that wounded corals are resilient to changes in growth under
# heat stress...

#Try plot again reordered
ggplot(buoy, aes(x = treatment, y = final.initial.initial, fill = wound)) +
  geom_boxplot() +
  scale_fill_manual(values = c("yes" = "black", "no" = "white")) +
  labs(
    title = "Change in Growth by Wound x Treatment",
    x = "Treatment",
    y = "Percent Growth"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_text(size = 10), # Makes the legend title smaller
    legend.text = element_text(size = 8),   # Makes the legend text (labels) smaller
    legend.key.size = unit(0.5, "lines")   # Makes the legend keys (symbols) smaller
  )

# Any visual difference in tanks before we do stats?
ggplot(buoy, aes(
  x = wound,
  y = final.initial.initial,
  fill = treatment
)) +
  geom_boxplot() +
  facet_wrap(~ tank) +
  labs(
    title = "Boxplot of final.initial.initial by Wound and Treatment Faceted by Tank",
    x = "Wound",
    y = "final.initial.initial"
  ) +
  theme_minimal()
#Tank 3 (ambient) is a slow grower compared to other ambient tanks. Wounding has  
# variable effects on growth, tank 11 shows negative effect and all other


