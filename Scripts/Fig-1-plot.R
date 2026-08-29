###########################################
# Predator-Prey Analysis Visualization for Ambient & Heatwave Temperature Treatments
# Code to reproduce Fig 1
###########################################


#load libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)

# 1. Read and standardize the datasets ------------------------------------

ambient <- read_csv(
  "ambient_temp_final.csv",
  show_col_types = FALSE
) |>
  mutate(temperature = "Ambient")

heatwave <- read_csv(
  "heatwave_temp_final.csv",
  show_col_types = FALSE
) |>
  mutate(temperature = "Heatwave")

# 3. Combine and standardize data ---------------------------------------------

raw_data <- dplyr::bind_rows(
  ambient,
  heatwave
) |>
  dplyr::mutate(
    # Ensure numeric variables are numeric
    harvest_week = readr::parse_number(
      as.character(harvest_week)
    ),
    
    collembola = readr::parse_number(
      as.character(collembola)
    ),
    
    staphylinidae = readr::parse_number(
      as.character(staphylinidae)
    ),
    
    staphylinidae_larval = readr::parse_number(
      as.character(staphylinidae_larval)
    ),
    
    # Standardize site names
    site = dplyr::case_when(
      stringr::str_detect(
        stringr::str_to_lower(site),
        "grass"
      ) ~ "Grassland",
      
      stringr::str_detect(
        stringr::str_to_lower(site),
        "agric"
      ) ~ "Intensive Agriculture",
      
      stringr::str_detect(
        stringr::str_to_lower(site),
        "urban"
      ) ~ "Urbanized",
      
      TRUE ~ NA_character_
    ),
    
    # Standardize predation-treatment names
    predation_trt = stringr::str_to_upper(
      stringr::str_trim(predation_trt)
    )
  ) |>
  dplyr::filter(
    !is.na(temperature),
    !is.na(site),
    !is.na(harvest_week)
  )


# Check observations by temperature, habitat, and week
dplyr::count(
  raw_data,
  temperature,
  site,
  harvest_week
)


# 4. Create the three arthropod groups ----------------------------------------

plot_long <- dplyr::bind_rows(
  
  # Prey-only treatment
  raw_data |>
    dplyr::filter(
      predation_trt == "C"
    ) |>
    dplyr::transmute(
      temperature,
      site,
      harvest_week,
      rep,
      group = "Collembola without beetles",
      abundance = collembola
    ),
  
  # Prey in the predator-addition treatment
  raw_data |>
    dplyr::filter(
      predation_trt == "CB"
    ) |>
    dplyr::transmute(
      temperature,
      site,
      harvest_week,
      rep,
      group = "Collembola with beetles",
      abundance = collembola
    ),
  
  # Predators in the predator-addition treatment
  raw_data |>
    dplyr::filter(
      predation_trt == "CB"
    ) |>
    dplyr::transmute(
      temperature,
      site,
      harvest_week,
      rep,
      group = "Beetles",
      
      abundance =
        staphylinidae +
        staphylinidae_larval
    )
) |>
  dplyr::mutate(
    group = factor(
      group,
      levels = c(
        "Collembola without beetles",
        "Collembola with beetles",
        "Beetles"
      )
    )
  )


# 5. Assign panels to the observed data ---------------------------------------

plot_long <- plot_long |>
  dplyr::mutate(
    panel_id = dplyr::case_when(
      temperature == "Ambient" &
        site == "Grassland" ~ "A",
      
      temperature == "Ambient" &
        site == "Intensive Agriculture" ~ "B",
      
      temperature == "Ambient" &
        site == "Urbanized" ~ "C",
      
      temperature == "Heatwave" &
        site == "Grassland" ~ "D",
      
      temperature == "Heatwave" &
        site == "Intensive Agriculture" ~ "E",
      
      temperature == "Heatwave" &
        site == "Urbanized" ~ "F",
      
      TRUE ~ NA_character_
    ),
    
    panel_id = factor(
      panel_id,
      levels = LETTERS[1:6]
    )
  ) |>
  dplyr::filter(
    !is.na(panel_id)
  )


# 6. Calculate observed means and standard errors -----------------------------

observed_summary <- plot_long |>
  dplyr::filter(
    !is.na(temperature),
    !is.na(site),
    !is.na(harvest_week),
    !is.na(abundance)
  ) |>
  dplyr::group_by(
    temperature,
    site,
    panel_id,
    harvest_week,
    group
  ) |>
  dplyr::summarise(
    n = sum(
      !is.na(abundance)
    ),
    
    mean_abundance = mean(
      abundance,
      na.rm = TRUE
    ),
    
    se_abundance = dplyr::if_else(
      n > 1,
      stats::sd(
        abundance,
        na.rm = TRUE
      ) / sqrt(n),
      0
    ),
    
    .groups = "drop"
  )


# 7. Add the initial mean abundances ------------------------------------------

initial_data <- tidyr::expand_grid(
  temperature = c(
    "Ambient",
    "Heatwave"
  ),
  
  site = c(
    "Grassland",
    "Intensive Agriculture",
    "Urbanized"
  ),
  
  group = c(
    "Collembola without beetles",
    "Collembola with beetles",
    "Beetles"
  )
) |>
  dplyr::mutate(
    harvest_week = 0,
    
    mean_abundance = dplyr::if_else(
      group == "Beetles",
      5,
      20
    ),
    
    se_abundance = 0,
    n = NA_integer_,
    
    panel_id = dplyr::case_when(
      temperature == "Ambient" &
        site == "Grassland" ~ "A",
      
      temperature == "Ambient" &
        site == "Intensive Agriculture" ~ "B",
      
      temperature == "Ambient" &
        site == "Urbanized" ~ "C",
      
      temperature == "Heatwave" &
        site == "Grassland" ~ "D",
      
      temperature == "Heatwave" &
        site == "Intensive Agriculture" ~ "E",
      
      temperature == "Heatwave" &
        site == "Urbanized" ~ "F"
    ),
    
    panel_id = factor(
      panel_id,
      levels = LETTERS[1:6]
    ),
    
    group = factor(
      group,
      levels = c(
        "Collembola without beetles",
        "Collembola with beetles",
        "Beetles"
      )
    )
  )


# 8. Combine initial and observed means ---------------------------------------

plot_summary <- dplyr::bind_rows(
  initial_data,
  observed_summary
) |>
  dplyr::mutate(
    harvest_week = as.numeric(
      harvest_week
    ),
    
    panel_id = factor(
      panel_id,
      levels = LETTERS[1:6]
    ),
    
    group = factor(
      group,
      levels = c(
        "Collembola without beetles",
        "Collembola with beetles",
        "Beetles"
      )
    )
  ) |>
  dplyr::arrange(
    panel_id,
    group,
    harvest_week
  )


# Confirm that weeks 0–3 occur in every panel
dplyr::count(
  plot_summary,
  panel_id,
  harvest_week
)


# 9. Prepare observed raw data from weeks 1–3 ---------------------------------

plot_raw_observed <- plot_long |>
  dplyr::filter(
    !is.na(temperature),
    !is.na(site),
    !is.na(panel_id),
    !is.na(harvest_week),
    !is.na(abundance)
  ) |>
  dplyr::mutate(
    panel_id = factor(
      panel_id,
      levels = LETTERS[1:6]
    ),
    
    group = factor(
      group,
      levels = c(
        "Collembola without beetles",
        "Collembola with beetles",
        "Beetles"
      )
    )
  )


# 10. Create initial observations for every microcosm -------------------------

# Each experimental unit began with:
# 20 Collembola or 5 beetles.
#
# scheduled_harvest_week preserves separate experimental units that were
# designated for different harvest weeks.

initial_raw <- plot_raw_observed |>
  dplyr::distinct(
    temperature,
    site,
    panel_id,
    rep,
    group,
    harvest_week
  ) |>
  dplyr::rename(
    scheduled_harvest_week = harvest_week
  ) |>
  dplyr::mutate(
    harvest_week = 0,
    
    abundance = dplyr::if_else(
      group == "Beetles",
      5,
      20
    )
  )


# 11. Combine initial and observed raw data -----------------------------------

plot_raw <- dplyr::bind_rows(
  initial_raw,
  plot_raw_observed
) |>
  dplyr::mutate(
    panel_id = factor(
      panel_id,
      levels = LETTERS[1:6]
    ),
    
    group = factor(
      group,
      levels = c(
        "Collembola without beetles",
        "Collembola with beetles",
        "Beetles"
      )
    )
  ) |>
  dplyr::arrange(
    panel_id,
    group,
    harvest_week,
    rep
  )


# Check number of raw points
dplyr::count(
  plot_raw,
  panel_id,
  group,
  harvest_week
)


# Confirm the largest raw values
plot_raw |>
  dplyr::arrange(
    dplyr::desc(abundance)
  ) |>
  dplyr::select(
    temperature,
    site,
    harvest_week,
    group,
    abundance
  ) |>
  head(10)


# 12. Facet labels and panel letters ------------------------------------------

panel_names <- c(
  "A" = "Grassland",
  "B" = "Intensive Agriculture",
  "C" = "Urbanized",
  "D" = "Grassland",
  "E" = "Intensive Agriculture",
  "F" = "Urbanized"
)


panel_letters <- tibble::tibble(
  panel_id = factor(
    LETTERS[1:6],
    levels = LETTERS[1:6]
  ),
  
  letter = LETTERS[1:6]
)


# 13. Treatment colors --------------------------------------------------------

group_colors <- c(
  "Collembola without beetles" = "#79AFA5",
  "Collembola with beetles" = "#2F6475",
  "Beetles" = "#DEB35B"
)


# 14. Calculate y-axis limit --------------------------------------------------

# This ensures that all raw observations, including the value of 122,
# are visible.

y_upper <- max(
  plot_raw$abundance,
  plot_summary$mean_abundance +
    plot_summary$se_abundance,
  na.rm = TRUE
) + 8


# Round upward to the nearest multiple of 5
y_upper <- ceiling(
  y_upper / 5
) * 5


# 15. Create figure -----------------------------------------------------------

arthropod_plot <- ggplot(
  plot_summary,
  aes(
    x = harvest_week,
    y = mean_abundance,
    color = group,
    group = group
  )
) +
  
  # Standard-error bars
  geom_errorbar(
    aes(
      ymin = pmax(
        mean_abundance -
          se_abundance,
        0
      ),
      
      ymax =
        mean_abundance +
        se_abundance
    ),
    width = 0,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  # Lines connecting mean abundances
  geom_line(
    linewidth = 0.9
  ) +
  
  # Solid-colored mean points
  geom_point(
    size = 3.2
  ) +
  
  # Treatment-colored raw observations
  #
  # This layer appears after the means so the raw points are drawn on top.
  geom_point(
    data = plot_raw,
    aes(
      x = harvest_week,
      y = abundance,
      color = group
    ),
    inherit.aes = FALSE,
    
    position = position_jitter(
      width = 0.18,
      height = 0,
      seed = 123
    ),
    
    shape = 16,
    size = 1.9,
    alpha = 0.75,
    show.legend = FALSE
  ) +
  
  # Six panels
  facet_wrap(
    ~panel_id,
    ncol = 3,
    labeller = as_labeller(
      panel_names
    ),
    drop = FALSE
  ) +
  
  # Panel letters A–F
  geom_text(
    data = panel_letters,
    aes(
      x = -Inf,
      y = Inf,
      label = letter
    ),
    inherit.aes = FALSE,
    hjust = -0.65,
    vjust = 1.35,
    size = 5.5,
    fontface = "bold",
    color = "black"
  ) +
  
  # Treatment colors and legend labels
  scale_color_manual(
    name = NULL,
    
    values = group_colors,
    
    breaks = c(
      "Collembola without beetles",
      "Collembola with beetles",
      "Beetles"
    ),
    
    labels = c(
      "Prey only",
      "Prey + predator",
      "Predator"
    ),
    
    drop = FALSE
  ) +
  
  # Display week 0 but label only weeks 1–3
  scale_x_continuous(
    breaks = 1:3,
    
    expand = expansion(
      mult = c(
        0.03,
        0.03
      )
    )
  ) +
  
  # Dynamic y-axis includes all raw observations
  scale_y_continuous(
    breaks = seq(
      from = 0,
      to = ceiling(y_upper / 25) * 25,
      by = 25
    ),
    
    expand = expansion(
      mult = c(
        0.02,
        0.03
      )
    )
  ) +
  
  coord_cartesian(
    # Extra horizontal space prevents jittered week-0 points from clipping
    xlim = c(
      -0.25,
      3.25
    ),
    
    ylim = c(
      -5,
      y_upper
    ),
    
    clip = "on"
  ) +
  
  labs(
    x = "Microcosm harvest week",
    
    y = paste0(
      "Average abundance\n",
      "(arthropods/microcosm)"
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
      color = "grey25",
      fill = NA,
      linewidth = 0.8
    ),
    
    # Panel spacing
    panel.spacing.x = unit(
      0.15,
      "cm"
    ),
    
    panel.spacing.y = unit(
      0.30,
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
        t = 4,
        b = 4
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
      linewidth = 0.7
    ),
    
    axis.ticks.length = unit(
      0.15,
      "cm"
    ),
    
    # Legend inside upper-right panel
    legend.position = "inside",
    
    legend.position.inside = c(
      0.98,
      0.98
    ),
    
    legend.justification = c(
      1,
      1
    ),
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 11,
      color = "black"
    ),
    
    legend.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    legend.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 5
    ),
    
    # Figure margins
    plot.margin = margin(
      t = 8,
      r = 10,
      b = 8,
      l = 8
    )
  )


# 16. Display figure ----------------------------------------------------------

arthropod_plot


# 17. Save publication-quality PNG --------------------------------------------

ggsave(
  filename = "arthropod_abundance_with_raw_data.png",
  plot = arthropod_plot,
  width = 10,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "white"
)