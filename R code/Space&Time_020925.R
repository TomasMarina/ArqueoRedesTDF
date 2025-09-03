# -----------------------------------------------------------------------------
#
# R SCRIPT FOR RECONSTRUCTING ANCIENT CONSUMER-RESOURCE NETWORKS
# Final version: Includes spatio-temporal context (using categorical time),
# habitat filtering, removal of cannibalism, and saves final plots to a
# multi-page PDF file.
#
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# PART 1: SETUP
# Load all necessary packages for the analysis.
# -----------------------------------------------------------------------------

# The 'tidyverse' is a collection of R packages for data science.
library(tidyverse)
# The 'readxl' package is used to read data directly from Excel files.
library(readxl)
# The 'igraph' package is the primary tool for network analysis in R.
library(igraph)
# The 'ggraph' package is an extension of ggplot2 for plotting networks.
library(ggraph)
# The 'ggrepel' package is used for non-overlapping text labels.
library(ggrepel)


# -----------------------------------------------------------------------------
# PART 2: DATA LOADING AND INITIAL CLEANING
# Read the updated raw data files.
# -----------------------------------------------------------------------------

# ASSUMPTION: The script assumes your .xlsx files are in a 'data' sub-folder.
# To make this script run, please create a folder named 'data' in the same
# directory as the script, and place your Excel files inside it.
localities_raw <- read_excel("data - Santiago2025/Species localities - TDF - 030925.xlsx")
traits_raw <- read_excel("data - Santiago2025/Species traits - TDF - 030925.xlsx", sheet = "traits")
ages_raw <- read_excel("data - Santiago2025/Age localities - TDF - 030925.xlsx")


# -----------------------------------------------------------------------------
# PART 3: DATA TIDYING AND PREPARATION
# Transform the raw data into a structured, analysis-ready format.
# -----------------------------------------------------------------------------

# --- Tidy the Localities Data ---
localities_df <- localities_raw %>%
  mutate(across(-TrophicSpecies, as.numeric)) %>%
  pivot_longer(
    cols = -TrophicSpecies,
    names_to = "Locality",
    values_to = "Abundance"
  ) %>%
  mutate(Abundance = replace_na(Abundance, 0)) %>%
  filter(Abundance > 0) %>%
  select(Locality, TrophicSpecies)

# --- Tidy the Traits Data ---
traits_df <- traits_raw %>%
  select(TrophicSpecies, BodySize_min, BodySize_max, FeedingStrategy, FeedingHabitat) %>%
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
  select(TrophicSpecies, BodyMass_g, Guild, FeedingHabitat) %>%
  filter(!is.na(BodyMass_g) & !is.na(FeedingHabitat))

# --- Tidy the Age and Spatial Data ---
ages_df <- ages_raw %>%
  select(Locality, Geological_ages, Region)


# -----------------------------------------------------------------------------
# PART 4: THE LOG-RATIO MODEL (LRM) FUNCTION WITH HABITAT FILTER
# The core function to calculate interaction probabilities.
# -----------------------------------------------------------------------------

calculate_lrm_matrix <- function(consumers, resources) {
  habitat_compatibility_matrix <- outer(
    consumers$FeedingHabitat,
    resources$FeedingHabitat,
    FUN = function(consumer_hab, resource_hab) {
      is_incompatible <- (consumer_hab == 'Marino' & resource_hab == 'Terrestre') |
        (consumer_hab == 'Terrestre' & resource_hab == 'Marino')
      ifelse(is_incompatible, 0, 1)
    }
  )
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
  prob_matrix_transposed <- t(prob_matrix)
  final_matrix <- prob_matrix_transposed * habitat_compatibility_matrix
  rownames(final_matrix) <- consumers$TrophicSpecies
  colnames(final_matrix) <- resources$TrophicSpecies
  return(final_matrix)
}


# -----------------------------------------------------------------------------
# PART 5: EXECUTION OF THE SIMULATION WORKFLOW
# Apply the LRM function to each site.
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
    if (nrow(consumers) == 0 || nrow(resources) == 0) return(NULL)
    replicated_matrices <- replicate(N_REPLICATES, calculate_lrm_matrix(consumers, resources), simplify = FALSE)
    threshold_matrices <- map(replicated_matrices, ~ .x * (.x >= 0.05))
    final_prob_matrix <- apply(simplify2array(threshold_matrices), 1:2, mean)
    as.data.frame(final_prob_matrix) %>%
      rownames_to_column(var = "Predator") %>%
      pivot_longer(cols = -Predator, names_to = "Prey", values_to = "Int_prob")
  }) %>%
  list_rbind(names_to = "Locality")


# -----------------------------------------------------------------------------
# PART 6: CREATE FINAL INTERACTION DATA FRAME
# Filter, join data, remove self-loops, and arrange the results.
# -----------------------------------------------------------------------------

final_interaction_df <- interaction_list %>%
  left_join(ages_df, by = "Locality") %>%
  filter(Predator != Prey) %>%
  filter(Int_prob > 0) %>%
  select(Locality, Region, Geological_ages, Predator, Prey, Int_prob) %>%
  arrange(Region, Geological_ages, Locality, Predator, desc(Int_prob))


# -----------------------------------------------------------------------------
# PART 7: VISUALIZATION OF NETWORK EVOLUTION AND EXPORT TO PDF
# Create a separate, chronologically-ordered plot for each region and
# save them all to a single PDF file.
# -----------------------------------------------------------------------------

# --- Step 1: Prepare data for plotting ---

# Define the correct chronological order for the geological ages.
age_levels <- c("Pleistoceno final", "Holoceno medio", "Holoceno final", "Histórico")

# Prepare the final data frame for plotting.
plotting_df <- final_interaction_df %>%
  filter(!is.na(Region) & !is.na(Geological_ages)) %>%
  mutate(Geological_ages = factor(Geological_ages, levels = age_levels))

# --- Step 2: Open PDF device and loop through regions to create plots ---

# Open the PDF file device. All subsequent plots will be saved here.
# We use a landscape orientation (width > height) which is good for faceted plots.
pdf("Network_Evolution_by_Region.pdf", width = 11, height = 8.5)

# Get the unique regions to loop over.
regions_to_plot <- unique(plotting_df$Region)

for (current_region in regions_to_plot) {
  
  region_data <- plotting_df %>%
    filter(Region == current_region)
  
  if (nrow(region_data) > 0) {
    
    edges_for_plotting <- region_data %>%
      select(from = Predator, to = Prey, weight = Int_prob, Geological_ages)
    
    vertices_for_plotting <- tibble(name = unique(c(region_data$Predator, region_data$Prey))) %>%
      left_join(traits_df %>% select(TrophicSpecies, Guild), by = c("name" = "TrophicSpecies")) %>%
      mutate(type = Guild == "Consumer")
    
    region_graph <- graph_from_data_frame(
      d = edges_for_plotting,
      vertices = vertices_for_plotting,
      directed = TRUE
    )
    
    evolution_plot <- ggraph(region_graph, layout = 'bipartite') +
      geom_edge_fan(aes(alpha = weight), show.legend = FALSE) +
      geom_node_point(aes(color = Guild), size = 4) +
      geom_node_text(aes(label = name), repel = TRUE, size = 2.5, max.overlaps = 15, bg.colour = "white", segment.color = 'grey50') +
      scale_color_manual(values = c("Consumer" = "tomato", "Resource" = "skyblue")) +
      theme_graph(base_family = 'sans', background = 'white') +
      labs(
        title = paste("Network Evolution in Region:", current_region),
        subtitle = "Networks are shown chronologically by geological age.",
        caption = "Consumers (predators) are on top; Resources (prey) are on the bottom."
      ) +
      facet_edges(~ Geological_ages)
    
    # This print command now sends the plot to the open PDF file, creating a new page.
    print(evolution_plot)
  }
}

# --- Step 3: Close the PDF device ---
# This is a crucial step to finalize and save the PDF file correctly.
dev.off()

# A message to let you know the script has finished and where the file is.
cat("\nPDF file 'Network_Evolution_by_Region.pdf' has been created in your working directory.\n")

