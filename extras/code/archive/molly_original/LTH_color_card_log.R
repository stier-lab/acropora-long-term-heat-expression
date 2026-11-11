library(tidyverse)


# Read the data
data <-read.csv("~/Downloads/Color card log - complete - data.csv")

# Convert "D3/D4" -> 3.5
convert_color <- function(color_str) {
  nums <- str_remove_all(color_str, "D") %>%
    str_split("/", simplify = TRUE) %>%
    as.numeric()
  mean(nums, na.rm = TRUE)
}

# 28 & 31, add numeric color
filtered_data <- data %>%
  filter(treatment %in% c(28, 31)) %>%
  mutate(
    color_numeric = map_dbl(color, convert_color),
    treatment = as.factor(treatment)
  )

# Summarize
avg_data <- filtered_data %>%
  group_by(day, treatment, thicket) %>%
  summarise(
    mean_color = mean(color_numeric, na.rm = TRUE),
    sd_color   = sd(color_numeric, na.rm = TRUE),
    n          = n(),
    se_color   = sd_color / sqrt(n),   # standard error
    .groups = "drop"
  )

# Plot average 
ggplot(avg_data, aes(x = day, y = mean_color, color = treatment, group = treatment)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_color - se_color, ymax = mean_color + se_color), width = 0.2) +
  facet_wrap(~thicket)+
  labs(
    title = "Average Color Card Score Over Time",
    x = "Day",
    y = "Average Color Score",
    color = "Treatment (°C)"
  ) +
  theme_minimal()
scale_color_manual(values = c("28" = "cyan3", "31" = "orange")
)


data <- read.csv("~/Downloads/Color card log - complete - data.csv")

# Define conversion function
convert_color <- function(color_str) {
  nums <- str_remove_all(color_str, "D") %>%
    str_split("/", simplify = TRUE) %>%
    as.numeric()
  mean(nums, na.rm = TRUE)
}

# Filter, convert, and add wound status
filtered_data <- data %>%
  filter(treatment %in% c(28, 31)) %>%
  mutate(
    color_numeric = map_dbl(color, convert_color),
    treatment = as.factor(treatment),
    wound_status = case_when(
      str_to_lower(wounded) == "yes" ~ "Wounded",
      str_to_lower(wounded) == "no"  ~ "Unwounded",
      TRUE ~ NA_character_
    )
  )

# Summarize by day × treatment × wound_status
avg_data <- filtered_data %>%
  group_by(day, treatment, wound_status, thicket) %>%
  summarise(
    mean_color = mean(color_numeric, na.rm = TRUE),
    sd_color   = sd(color_numeric, na.rm = TRUE),
    n          = n(),
    se_color   = sd_color / sqrt(n),
    .groups = "drop"
  )

# Plot: averages with se-error bar
  ggplot(avg_data, aes(x = day, y = mean_color, 
                       color = factor(treatment), 
                       shape = wound_status, 
                       group = interaction(treatment, wound_status))) +
    geom_line() +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_color - se_color, ymax = mean_color + se_color), width = 0.2) +
    scale_color_manual(values = c("28" = "cyan3", "31" = "orange")) +
    scale_shape_manual(values = c("Wounded" = 16, "Unwounded" = 1)) +  
    facet_wrap(~thicket)+
    labs(
      title = "Average Coral Color Score",
      x = "Day", 
      y = "Average Color Score", 
      color = "Treatment (°C)", 
      shape = "Wound"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.5, "lines")
    )
  
 