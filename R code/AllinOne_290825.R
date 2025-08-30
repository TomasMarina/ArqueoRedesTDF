# -----------------------------------------------------------------------------
#
# R SCRIPT FOR RECONSTRUCTING ANCIENT CONSUMER-RESOURCE NETWORKS
# A complete and optimized workflow from raw data to final interaction list.
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
# We read it directly from the Excel file.
localities_raw <- read_excel("data - Santiago2025/Species localities - TDF_290825.xlsx")

# --- Load Species Trait Data ---
# This file contains biological traits for each species, like body size.
# We specify the sheet name 'traits' to ensure we read the correct data.
traits_raw <- read_excel("data - Santiago2025/Species traits - TDF_290825.xlsx", sheet = "traits")


# -----------------------------------------------------------------------------
# PART 3: DATA TIDYING AND PREPARATION
# This is a critical step to transform the raw data into a structured,
# analysis-ready format. We follow the "tidy data" principles.
# -----------------------------------------------------------------------------

# --- Tidy the Localities Data ---
# The goal is to convert the "wide" localities data into a "long" format.
# A long format has one row per observation (e.g., one row for each species
# found at a specific site), which is much easier to work with.
localities_df <- localities_raw %>%
  # Before pivoting, we must ensure all locality columns have the same data type.
  # We convert them all to numeric. This will force non-numeric values (like 'x')
  # into NA, which we will handle after pivoting.
  mutate(across(-TrophicSpecies, as.numeric)) %>%
  
  # Use pivot_longer() to transform the data from wide to long.
  # We are "un-pivoting" all columns EXCEPT for 'TrophicSpecies'.
  pivot_longer(
    cols = -TrophicSpecies,
    names_to = "Locality",      # The new column for the site names.
    values_to = "Abundance"     # The new column for the abundance counts.
  ) %>%
  # The 'Abundance' column now contains numbers and NAs. We replace all NAs
  # with 0, assuming NA means the species was absent at that site.
  mutate(Abundance = replace_na(Abundance, 0)) %>%
  # Filter out rows where the species was not present (Abundance == 0).
  # This makes our dataset smaller and more efficient to work with.
  filter(Abundance > 0) %>%
  # Select only the columns we need for the analysis.
  select(Locality, TrophicSpecies)


# --- Tidy the Traits Data ---
# The goal is to create a clean table of species traits, focusing on body mass
# and feeding guild, which are essential for the model.
traits_df <- traits_raw %>%
  # Select only the columns relevant to our model.
  select(TrophicSpecies, BodySize_min, BodySize_max, FeedingStrategy) %>%
  # The model requires a single body mass value per species.
  # We calculate the geometric mean of the min and max body size, as body mass
  # often scales logarithmically. We convert from kg to grams.
  # `rowwise()` is used to perform the calculation on each row individually.
  rowwise() %>%
  mutate(
    BodyMass_g = exp(mean(log(c(BodySize_min, BodySize_max)), na.rm = TRUE)) * 1000
  ) %>%
  ungroup() %>% # Always ungroup after a rowwise operation.
  # Define a 'Guild' for each species. We classify any species that eats animals
  # as a "Consumer" for the purpose of our model.
  mutate(
    Guild = if_else(
      FeedingStrategy %in% c("Carnivore", "Omnivore", "Piscivore"),
      "Consumer",
      "Resource"
    )
  ) %>%
  # Keep only the essential columns and filter out any species that still
  # don't have a body mass after our calculations (e.g., if both min and max were NA).
  select(TrophicSpecies, BodyMass_g, Guild) %>%
  filter(!is.na(BodyMass_g))


# -----------------------------------------------------------------------------
# PART 4: THE LOG-RATIO MODEL (LRM) FUNCTION
# This section defines the core function to calculate interaction probabilities.
# This version is "vectorized" for speed, avoiding slow R loops.
# -----------------------------------------------------------------------------

#' Calculate Interaction Probabilities Using a Vectorized Log-Ratio Model
#'
#' This function takes two data frames (consumers and resources) and calculates
#' the probability of interaction for every possible pair based on their body mass.
#'
#' @param consumers A data frame of consumer species with a 'BodyMass_g' column.
#' @param resources A data frame of resource species with a 'BodyMass_g' column.
#' @return A matrix of interaction probabilities, with consumers as rows and
#'   resources as columns.
calculate_lrm_matrix <- function(consumers, resources) {
  
  # --- Model Parameters ---
  # These parameters are derived from previous studies (e.g., Nascimento et al. 2024)
  # and define the shape of the body-size niche.
  params <- list(
    # Parameters for Carnivores
    alphaP = 1.5, betaP = -1.5, gammaP = -0.75,
    # Parameters for Omnivores (assuming omnivores are less predatory)
    alphaO = -5.0, betaO = -2.5, gammaO = 0
  )
  
  # --- Vectorized Calculation ---
  # Create vectors of consumer and resource masses.
  consumer_masses <- consumers$BodyMass_g
  resource_masses <- resources$BodyMass_g
  
  # Create a matrix of all log-ratios in a single operation.
  # `outer()` applies a function to every combination of two vectors.
  log_ratio_matrix <- outer(resource_masses, consumer_masses, function(R, C) log(R / C))
  
  # Initialize the matrix for the logit(p) values.
  logit_p_matrix <- matrix(0, nrow = nrow(resources), ncol = nrow(consumers))
  
  # Identify which consumers are omnivores (we'll apply different parameters to them).
  # This creates a logical vector (TRUE/FALSE).
  is_omnivore <- consumers$Guild == "Omnivore"
  
  # Calculate logit(p) for carnivores (where is_omnivore is FALSE).
  if (any(!is_omnivore)) {
    logit_p_matrix[, !is_omnivore] <- params$alphaP +
      (params$betaP * log_ratio_matrix[, !is_omnivore]) +
      (params$gammaP * (log_ratio_matrix[, !is_omnivore]^2))
  }
  
  # Calculate logit(p) for omnivores (where is_omnivore is TRUE).
  if (any(is_omnivore)) {
    logit_p_matrix[, is_omnivore] <- params$alphaO +
      (params$betaO * log_ratio_matrix[, is_omnivore]) +
      (params$gammaO * (log_ratio_matrix[, is_omnivore]^2))
  }
  
  # Convert from logit scale to probability scale using the logistic function (`plogis`).
  prob_matrix <- plogis(logit_p_matrix)
  
  # Rescale predation probability for omnivores, as in the original script.
  if (any(is_omnivore)) {
    prob_matrix[, is_omnivore] <- prob_matrix[, is_omnivore] * 0.5
  }
  
  # Transpose the matrix to have consumers as rows and resources as columns.
  # Assign species names to the dimensions for clarity.
  final_matrix <- t(prob_matrix)
  rownames(final_matrix) <- consumers$TrophicSpecies
  colnames(final_matrix) <- resources$TrophicSpecies
  
  return(final_matrix)
}


# -----------------------------------------------------------------------------
# PART 5: EXECUTION OF THE SIMULATION WORKFLOW
# This is the main part of the script where we apply our function to each site.
# -----------------------------------------------------------------------------

# --- Prepare Data for Iteration ---
# We join the locality and trait data, then split the result into a list,
# with one data frame for each locality. This is a standard "split-apply-combine"
# strategy, implemented functionally.
sites_list <- localities_df %>%
  # Use an inner_join to keep only species that are present in both datasets.
  inner_join(traits_df, by = "TrophicSpecies") %>%
  # Group by locality and then split into a list of data frames.
  group_by(Locality) %>%
  group_split()

# --- Run Simulation Across All Sites ---
# We use purrr::map to apply our simulation logic to each site's data frame.

# Set the number of simulation replicates.
N_REPLICATES <- 100

interaction_list <- sites_list %>%
  # Set the names of the list elements to be the locality name.
  # This name will be preserved in the final output.
  set_names(map_chr(., ~ .x$Locality[1])) %>%
  # map() iterates over each site's data frame in the 'sites_list'.
  map(~ {
    # For each site, define the set of consumers and resources present.
    site_data <- .x
    consumers <- site_data %>% filter(Guild == "Consumer")
    resources <- site_data # All species can potentially be a resource.
    
    # A crucial check: if a site has no consumers or no resources, we can't
    # build a network. We return NULL, and `list_rbind` will skip it.
    if (nrow(consumers) == 0 || nrow(resources) == 0) {
      return(NULL)
    }
    
    # Replicate the LRM simulation N_REPLICATES times.
    # 'replicate' returns a list of probability matrices.
    replicated_matrices <- replicate(
      N_REPLICATES,
      calculate_lrm_matrix(consumers, resources),
      simplify = FALSE
    )
    
    # Apply a probability threshold (>= 0.05) to each matrix, setting lower
    # values to 0. This filters out very unlikely interactions.
    threshold_matrices <- map(replicated_matrices, ~ .x * (.x >= 0.05))
    
    # Average the thresholded matrices to get the final probability matrix.
    # 'simplify2array' converts the list to a 3D array for easy averaging.
    final_prob_matrix <- apply(simplify2array(threshold_matrices), 1:2, mean)
    
    # Convert the final matrix to a tidy (long) data frame.
    as.data.frame(final_prob_matrix) %>%
      rownames_to_column(var = "Predator") %>%
      pivot_longer(
        cols = -Predator,
        names_to = "Prey",
        values_to = "Int_prob"
      )
  }) %>%
  # Combine the list of site-specific data frames into one large data frame.
  # `list_rbind` is the modern, robust way to do this.
  # The 'names_to' argument creates a column from the list names (our localities).
  list_rbind(names_to = "Locality")


# -----------------------------------------------------------------------------
# PART 6: FINAL OUTPUT AND DISPLAY
# Filter and arrange the final results for interpretation.
# -----------------------------------------------------------------------------

final_interaction_df <- interaction_list %>%
  # Filter out the zero-probability interactions to keep the data frame clean.
  filter(Int_prob > 0) %>%
  # Arrange the columns in a logical order.
  select(Locality, Predator, Prey, Int_prob) %>%
  # Sort the results for easier reading.
  arrange(Locality, Predator, desc(Int_prob))

# --- Display Results ---
# Print the first few rows of the final interaction list to the console.
cat("Successfully generated interaction list.\n")
cat("Total interactions with probability > 0:", nrow(final_interaction_df), "\n\n")
print(head(final_interaction_df))

# --- Optional: Save Results ---
# You can uncomment the line below to save the final data frame to a CSV file.
# write_csv(final_interaction_df, "ancient_foodweb_interactions.csv")
