# R script for "Ancient consumer-resource networks in Tierra del Fuego"
# Authors: Fernando C. Santiago, Ulises Balza & Tomás I. Marina
# Script #2: Estimate consumer-resource interactions

# Load packages -----------------------------------------------------------

library(tidyverse)
library(purrr)


# Load data ---------------------------------------------------------------

load("data - Santiago2025/tidy_data_290825.Rdata")


# Community data frames ---------------------------------------------------
## Discard localities with NA values (more robust method)
# This keeps only the columns that are complete cases.
comm_full <- communities %>%
  select(where(~!any(is.na(.))))

## Classify species by feeding strategy
carn_species <- traits %>%
  filter(FeedingStrategy == "Carnivore") %>%
  pull(TrophicSpecies) # pull() is a clean way to get a single column as a vector

herb_species <- traits %>%
  filter(FeedingStrategy != "Carnivore") %>%
  pull(TrophicSpecies)

## Create community data frames based on the filtered species lists
# The first column is 'TrophicSpecies', so we keep it.
occ_carn <- comm_full %>% filter(TrophicSpecies %in% carn_species)
occ_herb <- comm_full %>% filter(TrophicSpecies %in% herb_species)

# Get the list of locality names to iterate over
locality_names <- names(comm_full)[-1] # Exclude the 'TrophicSpecies' column

## Create list of carnivore species present at each locality
communities_carn <- map(locality_names, ~{
  occ_carn %>%
    filter(.data[[.x]] > 0) %>% # Use .data[[]] to access column by name
    pull(TrophicSpecies)
}) %>%
  set_names(locality_names) # Set the names of the list elements

## Create list of herbivore species present at each locality
communities_herb <- map(locality_names, ~{
  occ_herb %>%
    filter(.data[[.x]] > 0) %>%
    pull(TrophicSpecies)
}) %>%
  set_names(locality_names)

## Include trait data for each locality's herbivore list
herb_traits <- map(communities_herb, ~{
  tibble(TrophicSpecies = .x) %>% # Create a tibble from the species vector
    left_join(traits, by = "TrophicSpecies") %>%
    select(TrophicSpecies, BodySize_max, FeedingStrategy)
})

## Include trait data for each locality's carnivore list
carn_traits <- map(communities_carn, ~{
  tibble(TrophicSpecies = .x) %>%
    left_join(traits, by = "TrophicSpecies") %>%
    select(TrophicSpecies, BodySize_max, FeedingStrategy)
})


# Estimate interaction prob -----------------------------------------------
## Log-Ratio Model ----
## Prey/Predator body size (in kg) ratio function
LRM <- function(mass_C, mass_H) {
  m <- nrow(mass_C)
  n <- nrow(mass_H)
  P <- matrix(NA, m, n)
  
  # Sample model parameters from parameter ranges
  # Carnivores
  alphaP <- runif(1, min = 1, max = 2)
  betaP  <- runif(1, min = -2, max = -1)
  gammaP <- runif(1, min = -1, max = -0.5)
  # Omnivores
  alphaO <- runif(1, min = -6, max = -4)
  betaO  <- runif(1, min = -3, max = -2)
  gammaO <- 0
  
  for (i in 1:m) {
    # Select parameters based on feeding strategy
    if (mass_C$FeedingStrategy[i] != "Omnivore") {
      alpha <- alphaP
      beta  <- betaP
      gamma <- gammaP
    } else {
      alpha <- alphaO
      beta  <- betaO
      gamma <- gammaO
    }
    
    for (j in 1:n) {
      log_ratio <- log(mass_H$BodySize_max[j] / mass_C$BodySize_max[i])
      term <- exp(alpha + (beta * log_ratio) + (gamma * (log_ratio^2)))
      P[i, j] <- term / (1 + term) # Probability of interaction
    }
    
    if (mass_C$FeedingStrategy[i] == "Omnivore") {
      P[i, ] <- P[i, ] * 0.5 # Rescaling predation probability
    }
  }
  
  row.names(P) <- mass_C$TrophicSpecies
  colnames(P) <- mass_H$TrophicSpecies
  return(P)
}


## Calculate probability ----
# Determine which localities have both carnivores and herbivores with complete data
localities_to_run <- names(which(
  sapply(carn_traits, nrow) > 0 & sapply(herb_traits, nrow) > 0
))

# Set number of replicates
N_REPLICATES <- 100

# Run simulation for each valid locality
# `map` is a cleaner alternative to a for-loop
results <- purrr::map(localities_to_run, ~{
  
  # Get trait data for the current locality
  carn_data <- carn_traits[[.x]]
  herb_data <- herb_traits[[.x]]
  
  # Replicate the LRM simulation N times
  replicated_matrices <- replicate(N_REPLICATES, LRM(carn_data, herb_data), simplify = FALSE)
  
  # Apply probability threshold correctly to each matrix
  # Set any probability < 0.05 to 0
  threshold_matrices <- purrr::map(replicated_matrices, ~ .x * (.x >= 0.05))
  
  # Average the thresholded matrices to get the final probability matrix
  # `simplify2array` converts the list to a 3D array for easy averaging
  final_matrix <- apply(simplify2array(threshold_matrices), 1:2, mean)
  
  return(round(final_matrix, 5))
  
}) %>%
  # Name the final list elements with their locality names
  set_names(localities_to_run)

# You can also store the raw replicas if needed
# replicas <- ... # (This part of the logic can be added here if you need the raw data)

## Bind all interactions ----
all_interactions <- purrr::map_dfr(
  .x = results,          # The list of matrices to iterate over
  .f = function(matrix) {
    # For each matrix, convert it to a tidy (long) data frame
    matrix %>%
      as.data.frame() %>%
      rownames_to_column(var = "Predator") %>%
      pivot_longer(
        cols = -Predator,
        names_to = "Prey",
        values_to = "Int_prob"
      )
  },
  .id = "Locality"  # Creates a 'Locality' column from the names of the 'results' list
)

# Optional: Filter out the non-interactions to keep the data frame smaller
all_interactions_filtered <- all_interactions %>%
  filter(Int_prob > 0)


# Save results ------------------------------------------------------------

save(results, all_interactions_filtered, traits, ages,
     file = "results/int_model_outputs_180725.Rdata")

write.csv(all_interactions, file = "results/InteractionProbability_180725.csv")
