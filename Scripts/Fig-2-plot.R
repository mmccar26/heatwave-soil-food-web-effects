###########################################
# Predator-Prey Analysis Visualization for Ambient & Heatwave Temperature Treatments
# Code to reproduce Fig 2
###########################################

# 1. Packages ----------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)
library(emmeans)
library(multcomp)
library(multcompView)

## IMPORT DATASETS FOR AMBIENT AND HEATWAVE TEMP TREATMENTS
# set pathname
setwd("C:/Users/mattm/Dropbox/Matt's folder/Rice University/Research/
      Soil and climate Predator-Prey/Analysis/predator_prey/Data/Data")

# 2. Import data --------------------------------------------------------------

ambient_data <- read_csv(
  "ambient_temp_final.csv",
  show_col_types = FALSE
) |>
  mutate(Temperature = "Ambient")


heatwave_data <- read_csv(
  "heatwave_temp_final.csv",
  show_col_types = FALSE
) |>
  mutate(Temperature = "Heatwave")


# 3. Combine and clean data ---------------------------------------------------

all_data <- bind_rows(
  ambient_data,
  heatwave_data
) |>
  mutate(
    # Ensure variables are numeric
    harvest_week = parse_number(
      as.character(harvest_week)
    ),
    
    collembola = parse_number(
      as.character(collembola)
    ),
    
    # Standardize temperature names
    Temperature = case_when(
      str_detect(
        str_to_lower(str_trim(Temperature)),
        "ambient"
      ) ~ "Ambient",
      
      str_detect(
        str_to_lower(str_trim(Temperature)),
        "heat"
      ) ~ "Heatwave",
      
      TRUE ~ NA_character_
    ),
    
    # Standardize habitat names
    site = case_when(
      str_detect(
        str_to_lower(str_trim(site)),
        "grass"
      ) ~ "Grassland",
      
      str_detect(
        str_to_lower(str_trim(site)),
        "agric"
      ) ~ "Intensive Agriculture",
      
      str_detect(
        str_to_lower(str_trim(site)),
        "urban"
      ) ~ "Urbanized",
      
      TRUE ~ NA_character_
    ),
    
    # Standardize predation-treatment names
    predation_trt = case_when(
      str_to_upper(str_trim(predation_trt)) %in%
        c(
          "C",
          "ABSENT",
          "PREY ONLY"
        ) ~ "Absent",
      
      str_to_upper(str_trim(predation_trt)) %in%
        c(
          "CB",
          "PRESENT",
          "PREY + PREDATOR"
        ) ~ "Present",
      
      TRUE ~ NA_character_
    )
  ) |>
  
  # Remove rows with missing required values
  filter(
    !is.na(Temperature),
    !is.na(site),
    !is.na(predation_trt),
    !is.na(harvest_week),
    !is.na(collembola)
  ) |>
  
  # Set factor levels
  mutate(
    Temperature = factor(
      Temperature,
      levels = c(
        "Ambient",
        "Heatwave"
      )
    ),
    
    site = factor(
      site,
      levels = c(
        "Grassland",
        "Intensive Agriculture",
        "Urbanized"
      )
    ),
    
    predation_trt = factor(
      predation_trt,
      levels = c(
        "Absent",
        "Present"
      )
    )
  )


# 4. Retain week 3 observations ----------------------------------------------

week3 <- all_data |>
  filter(
    harvest_week == 3
  ) |>
  droplevels()


# Check number of observations in each group
dplyr::count(
  week3,
  site,
  predation_trt,
  Temperature
)


# 5. Factorial ANOVA ----------------------------------------------------------

mod5 <- aov(
  log(collembola + 1) ~
    predation_trt *
    site *
    Temperature,
  data = week3
)


# Display ANOVA results
summary(mod5)


# 6. Estimated marginal means and Tukey letters -------------------------------

prey_emmeans <- emmeans(
  mod5,
  ~ predation_trt *
    site *
    Temperature
)

prey_letters <- multcomp::cld(
  prey_emmeans,
  adjust = "tukey",
  Letters = letters,
  sort = FALSE
) |>
  as.data.frame() |>
  transmute(
    site = str_trim(
      as.character(site)
    ),
    
    predation_trt = str_trim(
      as.character(predation_trt)
    ),
    
    Temperature = str_trim(
      as.character(Temperature)
    ),
    
    tukey_letter = str_trim(
      as.character(.group)
    )
  )


# Inspect Tukey letters
print(
  prey_letters,
  row.names = FALSE
)


# 7. Calculate means and standard errors --------------------------------------

sum_week3 <- week3 |>
  mutate(
    site = as.character(site),
    predation_trt = as.character(predation_trt),
    Temperature = as.character(Temperature)
  ) |>
  
  group_by(
    site,
    predation_trt,
    Temperature
  ) |>
  
  summarise(
    n = sum(
      !is.na(collembola)
    ),
    
    mean_collembola = mean(
      collembola,
      na.rm = TRUE
    ),
    
    sd_collembola = if_else(
      n > 1,
      sd(
        collembola,
        na.rm = TRUE
      ),
      0
    ),
    
    se_collembola = if_else(
      n > 1,
      sd_collembola / sqrt(n),
      0
    ),
    
    .groups = "drop"
  ) |>
  
  mutate(
    site = str_trim(site),
    predation_trt = str_trim(predation_trt),
    Temperature = str_trim(Temperature)
  ) |>
  
  # Add Tukey letters
  left_join(
    prey_letters,
    by = c(
      "site",
      "predation_trt",
      "Temperature"
    )
  ) |>
  
  # Restore factor levels
  mutate(
    site = factor(
      site,
      levels = c(
        "Grassland",
        "Intensive Agriculture",
        "Urbanized"
      )
    ),
    
    predation_trt = factor(
      predation_trt,
      levels = c(
        "Absent",
        "Present"
      )
    ),
    
    Temperature = factor(
      Temperature,
      levels = c(
        "Ambient",
        "Heatwave"
      )
    ),
    
    # Position Tukey letters above error bars
    letter_y =
      mean_collembola +
      se_collembola +
      5
  )


# Inspect summary statistics
sum_week3 |>
  select(
    site,
    predation_trt,
    Temperature,
    n,
    mean_collembola,
    sd_collembola,
    se_collembola,
    tukey_letter
  ) |>
  print(
    n = Inf
  )


# Stop if any Tukey letters failed to join
if (any(is.na(sum_week3$tukey_letter))) {
  stop(
    paste(
      "Some Tukey letters are missing.",
      "Check treatment names in prey_letters and sum_week3."
    )
  )
}


# 8. Calculate heatwave percentage changes -----------------------------------

percent_labels <- sum_week3 |>
  dplyr::select(
    site,
    predation_trt,
    Temperature,
    mean_collembola,
    se_collembola
  ) |>
  
  tidyr::pivot_wider(
    names_from = Temperature,
    values_from = c(
      mean_collembola,
      se_collembola
    )
  ) |>
  
  dplyr::mutate(
    # Percentage change from Ambient to Heatwave
    percent_change =
      100 *
      (
        mean_collembola_Heatwave -
          mean_collembola_Ambient
      ) /
      mean_collembola_Ambient,
    
    # Format percentage labels
    percent_label = dplyr::case_when(
      is.na(percent_change) ~ NA_character_,
      
      is.infinite(percent_change) ~ NA_character_,
      
      round(percent_change, 1) == -100 ~ "-100%",
      
      TRUE ~ paste0(
        sprintf(
          "%.2f",
          percent_change
        ),
        "%"
      )
    ),
    
    # Initial vertical position
    calculated_y = pmax(
      mean_collembola_Ambient +
        se_collembola_Ambient,
      
      mean_collembola_Heatwave +
        se_collembola_Heatwave,
      
      na.rm = TRUE
    ) + 18,
    
    # Final vertical positions
    percent_y = dplyr::case_when(
      site == "Intensive Agriculture" &
        predation_trt == "Absent" ~ 70,
      
      TRUE ~ pmin(
        calculated_y,
        76
      )
    )
  ) |>
  
  dplyr::mutate(
    site = factor(
      as.character(site),
      levels = c(
        "Grassland",
        "Intensive Agriculture",
        "Urbanized"
      )
    ),
    
    predation_trt = factor(
      as.character(predation_trt),
      levels = c(
        "Absent",
        "Present"
      )
    )
  ) |>
  
  dplyr::filter(
    !is.na(percent_label),
    !is.na(percent_y)
  ) |>
  
  # Ensure one percentage per site × treatment
  dplyr::distinct(
    site,
    predation_trt,
    .keep_all = TRUE
  )

# Inspect percentage changes
percent_labels |>
  dplyr::select(
    site,
    predation_trt,
    percent_change,
    percent_label,
    percent_y
  ) |>
  print(n = Inf)

# 9. Separate the shifted percentage label -----------------------------------

# Intensive-agriculture absent label
percent_label_shifted <- percent_labels |>
  filter(
    site == "Intensive Agriculture",
    predation_trt == "Absent"
  )


# All remaining percentage labels
percent_labels_regular <- percent_labels |>
  filter(
    !(
      site == "Intensive Agriculture" &
        predation_trt == "Absent"
    )
  )


# 10. Prepare raw observations ------------------------------------------------

raw_week3 <- week3 |>
  mutate(
    site = factor(
      as.character(site),
      levels = c(
        "Grassland",
        "Intensive Agriculture",
        "Urbanized"
      )
    ),
    
    predation_trt = factor(
      as.character(predation_trt),
      levels = c(
        "Absent",
        "Present"
      )
    ),
    
    Temperature = factor(
      as.character(Temperature),
      levels = c(
        "Ambient",
        "Heatwave"
      )
    )
  ) |>
  
  filter(
    !is.na(site),
    !is.na(predation_trt),
    !is.na(Temperature),
    !is.na(collembola)
  )

# Set y-axis high enough to display all raw observations
y_upper <- max(
  raw_week3$collembola,
  sum_week3$letter_y,
  percent_labels$percent_y,
  na.rm = TRUE
) + 8

y_upper <- ceiling(y_upper / 5) * 5

# 11. Panel letters ------------------------------------------------------------

panel_letters <- tibble(
  site = factor(
    c(
      "Grassland",
      "Intensive Agriculture",
      "Urbanized"
    ),
    levels = c(
      "Grassland",
      "Intensive Agriculture",
      "Urbanized"
    )
  ),
  
  panel_letter = c(
    "A",
    "B",
    "C"
  )
)


# 12. Colors ------------------------------------------------------------------

temperature_colors <- c(
  "Ambient" = "#CC7810",
  "Heatwave" = "#EFC519"
)


# 13. Position settings -------------------------------------------------------

bar_dodge <- position_dodge(
  width = 0.82
)


raw_point_position <- position_jitterdodge(
  jitter.width = 0.12,
  jitter.height = 0,
  dodge.width = 0.82,
  seed = 123
)


# 14. Create figure -----------------------------------------------------------

week3_figure2 <- ggplot(
  sum_week3,
  aes(
    x = predation_trt,
    y = mean_collembola,
    fill = Temperature,
    group = Temperature
  )
) +
  
  # Mean abundance bars
  geom_col(
    position = bar_dodge,
    width = 0.72,
    color = "black",
    linewidth = 0.6
  ) +
  
  # Raw week-3 observations
  geom_point(
    data = raw_week3,
    aes(
      x = predation_trt,
      y = collembola,
      group = Temperature
    ),
    inherit.aes = FALSE,
    position = raw_point_position,
    shape = 21,
    size = 4,
    stroke = 0.75,
    fill = "lightblue",
    color = "black",
    alpha = 0.9,
    show.legend = FALSE
  ) +
  
  # Error bars: mean ± 1 standard error
  geom_errorbar(
    aes(
      ymin = pmax(
        mean_collembola -
          se_collembola,
        0
      ),
      
      ymax =
        mean_collembola +
        se_collembola
    ),
    position = bar_dodge,
    width = 0.28,
    linewidth = 1.2,
    color = "black"
  ) +
  
  # Tukey grouping letters
  geom_text(
    aes(
      y = letter_y,
      label = tukey_letter
    ),
    position = bar_dodge,
    size = 5,
    fontface = "bold",
    color = "black",
    show.legend = FALSE
  ) +
  
  # Regular percentage labels
  geom_text(
    data = percent_labels_regular,
    aes(
      x = predation_trt,
      y = percent_y,
      label = percent_label
    ),
    inherit.aes = FALSE,
    size = 4.2,
    fontface = "bold",
    color = "black"
  ) +
  
  # Shifted intensive-agriculture absent label
  geom_text(
    data = percent_label_shifted,
    aes(
      x = predation_trt,
      y = percent_y,
      label = percent_label
    ),
    inherit.aes = FALSE,
    nudge_x = 0.22,
    size = 4.2,
    fontface = "bold",
    color = "black"
  ) +
  
  # Panel letters
  geom_text(
    data = panel_letters,
    aes(
      x = -Inf,
      y = Inf,
      label = panel_letter
    ),
    inherit.aes = FALSE,
    hjust = -0.45,
    vjust = 1.35,
    size = 5,
    fontface = "bold",
    color = "black"
  ) +
  
  # Habitat panels
  facet_wrap(
    ~site,
    nrow = 1,
    drop = FALSE
  ) +
  
  # Temperature colors
  scale_fill_manual(
    name = NULL,
    values = temperature_colors,
    breaks = c(
      "Ambient",
      "Heatwave"
    ),
    drop = FALSE
  ) +
  
  # Legend formatting
  guides(
    fill = guide_legend(
      override.aes = list(
        color = "black",
        linewidth = 0.8
      )
    )
  ) +
  
  scale_x_discrete(
    limits = c(
      "Absent",
      "Present"
    ),
    labels = c(
      "Absent" = "Absent",
      "Present" = "Present"
    )
  ) +
  
  scale_y_continuous(
    breaks = seq(
      from = 0,
      to = ceiling(y_upper / 25) * 25,
      by = 25
    ),
    expand = expansion(
      mult = c(0, 0.03)
    )
  ) +
  
  coord_cartesian(
    ylim = c(-4, y_upper),
    clip = "on"
  ) +
  
  labs(
    x = "Predation pressure",
    y = paste0(
      "Average abundance\n",
      "(Collembola/microcosm)"
    )
  ) +
  
  theme_bw(
    base_size = 13
  ) +
  
  theme(
    # Remove grid lines
    panel.grid = element_blank(),
    
    # Panel borders
    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = 0.8
    ),
    
    # Space between panels
    panel.spacing.x = unit(
      0.15,
      "cm"
    ),
    
    # Facet-strip background
    strip.background = element_rect(
      fill = "grey75",
      color = NA
    ),
    
    # Facet-strip text
    strip.text = element_text(
      size = 13,
      face = "bold",
      color = "black",
      margin = margin(
        t = 5,
        b = 5
      )
    ),
    
    # Axis titles
    axis.title = element_text(
      size = 13,
      face = "bold",
      color = "black"
    ),
    
    # Axis tick labels
    axis.text = element_text(
      size = 11,
      face = "bold",
      color = "black"
    ),
    
    # Axis ticks
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.8
    ),
    
    axis.ticks.length = unit(
      0.15,
      "cm"
    ),
    
    # Legend
    legend.position = "right",
    
    legend.text = element_text(
      size = 13,
      face = "bold",
      color = "black"
    ),
    
    legend.key.height = unit(
      0.8,
      "cm"
    ),
    
    legend.key.width = unit(
      0.8,
      "cm"
    ),
    
    # Figure margins
    plot.margin = margin(
      t = 8,
      r = 8,
      b = 8,
      l = 8
    )
  )


# 15. Display figure ----------------------------------------------------------

week3_figure2


# 16. Save high-resolution PNG ------------------------------------------------

ggsave(
  filename = "week3_collembola_abundance_raw_data.png",
  plot = week3_figure2,
  width = 10,
  height = 5.5,
  units = "in",
  dpi = 600,
  bg = "white"
)