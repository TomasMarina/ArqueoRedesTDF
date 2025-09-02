# -----------------------------------------------------------------------------
#
# R SCRIPT FOR RECONSTRUCTING ANCIENT CONSUMER-RESOURCE NETWORKS
# A complete and optimized workflow from raw data to final interaction list,
# now including spatio-temporal context with a focus on regional grouping.
#
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# PART 1: SETUP
# Load all necessary packages for the analysis.
# -----------------------------------------------------------------------------

# The 'tidyverse' is a collection of R packages designed for data science.
# It includes 'dplyr' for data manipulation, 'tidyr' for tidying data,
# and 'purrr' for functional programming.
library(tidyverse)

# The 'readxl' package is used to read data directly from Excel files.
library(readxl)


# -----------------------------------------------------------------------------
# PART 2: DATA LOADING AND INITIAL CLEANING
# Read the raw data files and perform essential cleaning steps.
# ASSUMPTION: The script assumes your .xlsx files are in a 'data' sub-folder.
# Please adjust the file paths if your files are located elsewhere.
# -----------------------------------------------------------------------------

# --- Load Species Localities Data ---
# This file contains the abundance of each species at each archaeological site.
localities_raw <- read_excel("data/Species localities - TDF_290825.xlsx")

# --- Load Species Trait Data ---
# This file contains biological traits for each species, like body size.
traits_raw <- read_excel("data/Species traits - TDF_290825.xlsx", sheet = "traits")

# --- Load Age and Spatial Data ---
# This file contains the age and regional context for each locality.
ages_raw <- read_excel("data/Age localities - TDF.xlsx")


# -----------------------------------------------------------------------------
# PART 3: DATA TIDYING AND PREPARATION
# This is a critical step to transform the raw data into a structured,
# analysis-ready format. We follow the "tidy data" principles.
# -----------------------------------------------------------------------------

# --- Tidy the Localities Data ---
# The goal is to convert the "wide" localities data into a "long" format.
localities_df <- localities_raw %>%
  # Before pivoting, we must ensure all locality columns have the same data type.
  mutate(across(-TrophicSpecies, as.numeric)) %>%
  # Use pivot_longer() to transform the data from wide to long.
  pivot_longer(
    cols = -TrophicSpecies,
    names_to = "Locality",
    values_to = "Abundance"
  ) %>%
  # Replace all NAs with 0, assuming NA means the species was absent.
  mutate(Abundance = replace_na(Abundance, 0)) %>%
  # Filter out rows where the species was not present.
  filter(Abundance > 0) %>%
  select(Locality, TrophicSpecies)

# --- Tidy the Traits Data ---
# The goal is to create a clean table of species traits.
traits_df <- traits_raw %>%
  select(TrophicSpecies, BodySize_min, BodySize_max, FeedingStrategy) %>%
  rowwise() %>%
  mutate(
    BodyMass_g = exp(mean(log(c(BodySize_min, BodySize_max)), na.rm = TRUE)) * 1000
  ) %>%
  ungroup() %>%
  mutate(
    Guild = if_else(
      FeedingStrategy %in% c("Carnivore", "Omnivore", "Piscivore"),
      "Consumer",
      "Resource"
    )
  ) %>%
  select(TrophicSpecies, BodyMass_g, Guild) %>%
  filter(!is.na(BodyMass_g))

# --- Tidy the Age and Spatial Data ---
# Create a clean lookup table for locality age and regional context.
ages_df <- ages_raw %>%
  # Select the key columns for joining and analysis, using 'Region' as requested.
  select(Locality, Age_Average, Region)


# -----------------------------------------------------------------------------
# PART 4: THE LOG-RATIO MODEL (LRM) FUNCTION
# This section defines the core function to calculate interaction probabilities.
# This version is "vectorized" for speed and remains unchanged.
# -----------------------------------------------------------------------------

#' Calculate Interaction Probabilities Using a Vectorized Log-Ratio Model
#'
#' @param consumers A data frame of consumer species with a 'BodyMass_g' column.
#' @param resources A data frame of resource species with a 'BodyMass_g' column.
#' @return A matrix of interaction probabilities.
calculate_lrm_matrix <- function(consumers, resources) {
  params <- list(
    alphaP = 1.5, betaP = -1.5, gammaP = -0.75,
    alphaO = -5.0, betaO = -2.5, gammaO = 0
  )
  consumer_masses <- consumers$BodyMass_g
  resource_masses <- resources$BodyMass_g
  log_ratio_matrix <- outer(resource_masses, consumer_masses, function(R, C) log(R / C))
  logit_p_matrix <- matrix(0, nrow = nrow(resources), ncol = nrow(consumers))
  is_omnivore <- consumers$Guild == "Omnivore"
  
  if (any(!is_omnivore)) {
    logit_p_matrix[, !is_omnivore] <- params$alphaP +
      (params$betaP * log_ratio_matrix[, !is_omnivore]) +
      (params$gammaP * (log_ratio_matrix[, !is_omnivore]^2))
  }
  if (any(is_omnivore)) {
    logit_p_matrix[, is_omnivore] <- params$alphaO +
      (params$betaO * log_ratio_matrix[, is_omnivore]) +
      (params$gammaO * (log_ratio_matrix[, is_omnivore]^2))
  }
  
  prob_matrix <- plogis(logit_p_matrix)
  if (any(is_omnivore)) {
    prob_matrix[, is_omnivore] <- prob_matrix[, is_omnivore] * 0.5
  }
  
  final_matrix <- t(prob_matrix)
  rownames(final_matrix) <- consumers$TrophicSpecies
  colnames(final_matrix) <- resources$TrophicSpecies
  return(final_matrix)
}


# -----------------------------------------------------------------------------
# PART 5: EXECUTION OF THE SIMULATION WORKFLOW
# This is the main part of the script where we apply our function to each site.
# -----------------------------------------------------------------------------

sites_list <- localities_df %>%
  inner_join(traits_df, by = "TrophicSpecies") %>%
  group_by(Locality) %>%
  group_split()

N_REPLICATES <- 100

interaction_list <- sites_list %>%
  set_names(map_chr(., ~ .x$Locality[1])) %>%
  map(~ {
    site_data <- .x
    consumers <- site_data %>% filter(Guild == "Consumer")
    resources <- site_data
    if (nrow(consumers) == 0 || nrow(resources) == 0) {
      return(NULL)
    }
    replicated_matrices <- replicate(
      N_REPLICATES,
      calculate_lrm_matrix(consumers, resources),
      simplify = FALSE
    )
    threshold_matrices <- map(replicated_matrices, ~ .x * (.x >= 0.05))
    final_prob_matrix <- apply(simplify2array(threshold_matrices), 1:2, mean)
    as.data.frame(final_prob_matrix) %>%
      rownames_to_column(var = "Predator") %>%
      pivot_longer(
        cols = -Predator,
        names_to = "Prey",
        values_to = "Int_prob"
      )
  }) %>%
  list_rbind(names_to = "Locality")


# -----------------------------------------------------------------------------
# PART 6: FINAL OUTPUT AND DISPLAY
# Filter, join spatio-temporal data, and arrange the final results.
# -----------------------------------------------------------------------------

final_interaction_df <- interaction_list %>%
  # Join the interaction data with the age/regional context data.
  # This adds the new columns based on the matching 'Locality' name.
  left_join(ages_df, by = "Locality") %>%
  # Filter out the zero-probability interactions.
  filter(Int_prob > 0) %>%
  # Arrange the columns in a logical order, using 'Region' as the primary spatial variable.
  select(Locality, Region, Age_Average, Predator, Prey, Int_prob) %>%
  # Sort the results for easier reading, now prioritizing the spatial 'Region'.
  arrange(Region, Age_Average, Locality, Predator, desc(Int_prob))

# --- Display Results ---
# Print the first few rows of the final interaction list to the console.
cat("Successfully generated and integrated spatio-temporal data.\n")
cat("Total interactions with probability > 0:", nrow(final_interaction_df), "\n\n")
print(head(final_interaction_df))

# --- Optional: Save Results ---
# You can uncomment the line below to save the final data frame to a CSV file.
# write_csv(final_interaction_df, "ancient_foodweb_interactions_region.csv")


