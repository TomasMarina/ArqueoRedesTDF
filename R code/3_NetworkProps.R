# R script for "Ancient consumer-resource networks in Tierra del Fuego"
# Authors: Fernando C. Santiago, Ulises Balza & Tomás I. Marina
# Script #3: Calculate network properties

# Load packages -----------------------------------------------------------

library(tidyverse)
library(purrr)
library(vegan)


# Load data ---------------------------------------------------------------
# This file should contain the 'interactions', 'ages', and 'traits' objects
load("results/int_model_outputs_180725.Rdata")


# Calculate network props -------------------------------------------------

network_properties <- purrr::map_dfr(
  .x = results,  # The list of interaction matrices
  .f = ~{
    # Ensure the element is a matrix and handle empty cases
    if (!is.matrix(.x) || nrow(.x) == 0 || ncol(.x) == 0) {
      return(tibble()) # Return an empty tibble to skip this locality
    }
    
    # Calculate all properties for the current matrix
    ## Predator/Prey Ratio
    ppr                <- nrow(.x) / ncol(.x)
    ## Quantitative connectance
    connectance        <- sum(.x) / (nrow(.x) * ncol(.x))
    shannon_diversity  <- vegan::diversity(as.numeric(.x), index = "shannon")
    mean_degree_carn   <- mean(rowSums(.x))
    
    # Return all properties as a single-row tibble
    tibble(
      PPR = ppr,
      Connectance = connectance,
      Diversity = shannon_diversity,
      Mean_Degree_Carn = mean_degree_carn
    )
  },
  .id = "Locality" # Creates the 'Locality' column from the list names
)


# Merge with age data ----------------------------------------------------
# This joins the calculated properties with the locality age information
# all_properties <- network_properties %>%
#   left_join(ages, by = "Locality") %>%
#   # Calculate mean age and place it after the 'Locality' column
#   mutate(
#     Age = (MIN_age + MAX_age) / 2,
#     .after = Locality
#   ) %>%
#   # Optional: remove redundant age columns
#   select(-MIN_age, -MAX_age)


# Save results ------------------------------------------------------------

save(network_properties,
     file = "results/net_props_180725.Rdata")

write.csv(network_properties, file = "results/NetworkProperties_180725.csv")
