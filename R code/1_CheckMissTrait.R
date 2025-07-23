# Load necessary packages
library(tidyverse)

# Load your initial data
# Make sure the file path is correct
load("data - Santiago2025/tidy_data_180725.Rdata")

# --- Data Validation Step ---

# 1. Get a unique list of all species with an abundance > 0 across all localities
species_in_communities <- communities %>%
  # Ensure the first column is named 'TrophicSpecies' for consistency
  rename(TrophicSpecies = 1) %>%
  # Pivot to a long format to easily filter
  pivot_longer(
    cols = -TrophicSpecies,
    names_to = "Locality",
    values_to = "Abundance"
  ) %>%
  filter(Abundance > 0) %>%
  distinct(TrophicSpecies) %>% # Get a unique list of species names
  pull(TrophicSpecies)         # Convert to a simple vector

# 2. Get a unique list of all species that have trait data
species_with_traits <- traits %>%
  distinct(TrophicSpecies) %>%
  pull(TrophicSpecies)

# 3. Compare the lists to find which species are missing from the traits file
missing_species <- data.frame(setdiff(species_in_communities, species_with_traits))

# 4. Display the results
if (length(missing_species) > 0) {
  print("Warning: The following species are present in communities but are missing trait data:")
  print(missing_species)
} else {
  print("Excellent! All species present in the communities have corresponding trait data. ✅")
}
