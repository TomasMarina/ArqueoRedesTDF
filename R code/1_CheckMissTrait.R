# R script for "Ancient consumer-resource networks in Tierra del Fuego"
# Authors: Fernando C. Santiago, Ulises Balza & Tomás I. Marina
# Script #1.1: Check missing traits


# Load packages -----------------------------------------------------------

library(tidyverse)
library(writexl)


# Loada data --------------------------------------------------------------

load("data - Santiago2025/tidy_data_180725.Rdata")


# Identify missing species ------------------------------------------------
# 1. Identify species present in communities but missing from the traits file
species_in_communities <- communities %>%
  rename(TrophicSpecies = 1) %>%
  pivot_longer(cols = -TrophicSpecies, values_to = "Abundance") %>%
  filter(Abundance > 0) %>%
  distinct(TrophicSpecies) %>%
  pull(TrophicSpecies)

species_with_traits <- traits %>%
  distinct(TrophicSpecies) %>%
  pull(TrophicSpecies)

missing_species_names <- setdiff(species_in_communities, species_with_traits)

# 2. Append the missing species to the original traits data frame
if (length(missing_species_names) > 0) {
  # Create a new data frame with just the names of the missing species
  missing_species_df <- tibble(TrophicSpecies = missing_species_names)
  
  # Bind the missing species to the end of your original traits data frame
  # All other trait columns will be automatically filled with NA
  traits_updated <- bind_rows(traits, missing_species_df)
  
# 3. Save the updated data frame to a new Excel file
  write_xlsx(
    traits_updated,
    path = "data - Santiago2025/Species_traits_miss.xlsx"
  )
  
  print(paste(
    "Success! An updated file named 'Species_traits_miss.xlsx' has been saved.",
    "It contains", length(missing_species_names), "new species rows for you to complete. ✍️"
  ))
  
} else {
  print("No missing species were found. Your traits file is already complete! ✅")
}
