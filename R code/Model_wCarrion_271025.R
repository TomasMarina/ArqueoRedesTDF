# ----------------------------------------------------------------------------
#
# R SCRIPT FOR RECONSTRUCTING ANCIENT CONSUMER-RESOURCE NETWORKS
# (Tierra del Fuego Archaeological Study)
#
# Version: 241025
# This version is updated to handle specific scavenging interactions
# based on 'carrion_only' and 'facultative_scavenger' traits.
#
# ----------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# PART 1: SETUP
# Load all necessary packages for the analysis.
# ----------------------------------------------------------------------------

# Ensure all packages are installed:
# install.packages(c("tidyverse", "readxl", "janitor", "bipartite", "igraph", 
#                    "tidygraph", "ggraph", "ggrepel", "patchwork", "viridis"))

library(tidyverse)
library(readxl)
library(janitor)       # For cleaning column names
library(bipartite)
library(igraph)
library(tidygraph)
library(ggraph)
library(ggrepel)
library(patchwork)
library(viridis)

set.seed(123) # for reproducibility


# ----------------------------------------------------------------------------
# PART 2: DATA LOADING AND INITIAL CLEANING
# Read the raw data files.
# *** UPDATED ***: Using user-specified paths with read_excel.
# ----------------------------------------------------------------------------

print("Loading data...")

# Use file paths from the original script
localities_raw <- read_excel("data - Santiago2025/Species localities - TDF - 030925.xlsx")
traits_raw <- read_excel("data - Santiago2025/Species traits - TDF - 241025.xlsx", sheet = "traits")
ages_raw <- read_excel("data - Santiago2025/Age localities - TDF - 030925.xlsx")


# ----------------------------------------------------------------------------
# PART 3: DATA PREPARATION
# Clean and prepare data frames for analysis.
# ----------------------------------------------------------------------------

print("Preparing data...")

# --- Clean Traits Data ---
traits <- traits_raw %>%
  clean_names() %>%
  # Ensure TrophicSpecies is a character
  mutate(trophic_species = as.character(trophic_species)) %>%
  # Calculate average log10 body mass (in kg)
  mutate(log_mass = log10((body_size_min + body_size_max) / 2)) %>%
  # Handle NAs in new scavenging columns
  # *** NEW ***: Convert NAs to 0 for logical comparisons
  mutate(
    carrion_only = ifelse(is.na(carrion_only), 0, carrion_only),
    facultative_scavenger = ifelse(is.na(facultative_scavenger), 0, facultative_scavenger)
  ) %>%
  # Select relevant columns
  select(
    trophic_species, class, order, family, 
    log_mass, feeding_strategy, feeding_habitat,
    carrion_only, facultative_scavenger, human_resource
  ) %>%
  # Ensure no duplicate species
  distinct(trophic_species, .keep_all = TRUE)

# --- Clean Localities Data ---
# *** UPDATED ***: Added 'values_transform' to handle mixed data types.
localities <- localities_raw %>%
  clean_names() %>%
  pivot_longer(
    cols = -trophic_species,
    names_to = "locality",
    values_to = "nisp_raw", # Store the raw value first
    # This forces all values (numeric or text) into a character string
    values_transform = list(nisp_raw = as.character)
  ) %>%
  # Now, safely convert the character strings to numbers.
  # Any text (like 'x' or 'present') will become NA.
  mutate(nisp = as.numeric(nisp_raw)) %>%
  # Filter for presence (NISP > 0) and ignore NAs
  filter(!is.na(nisp) & nisp > 0) %>%
  # Select only the columns we need
  select(locality, trophic_species) %>%
  # Merge with traits to ensure we only keep species with trait data
  inner_join(traits %>% select(trophic_species), by = "trophic_species")

# --- Clean Ages Data ---
# *** UPDATED ***: Renamed 'geological_ages' to 'time_bin' to match script
ages <- ages_raw %>%
  clean_names() %>%
  # The column in your file is "Geological_ages"
  select(locality, time_bin = geological_ages) %>%
  distinct()

# --- Create final species-by-locality dataset ---
species_by_locality <- localities %>%
  left_join(ages, by = "locality") %>%
  filter(!is.na(time_bin)) %>%
  select(time_bin, locality, trophic_species) %>%
  distinct()

# --- Define plotting colors for taxa ---
taxa_colors <- c(
  "Mammalia" = "#a6611a",  # Brown
  "Aves" = "#018571",      # Teal
  "Actinopterygii" = "#4575b4", # Blue
  "Chondrichthyes" = "#2c7fb8",
  "Mollusca" = "#7fbc41",     # Green
  "Malacostraca" = "#d01c8b",  # Magenta
  "Echinoidea" = "#f1b6da"   # Pink
)


# ----------------------------------------------------------------------------
# PART 4: HELPER FUNCTIONS
# Define functions for network reconstruction and analysis.
# ----------------------------------------------------------------------------

#' Check Habitat Overlap
#'
#' Checks if two species share at least one feeding habitat.
#' Habitats are provided as strings, potentially multi-valued (e.g., "Ter-Mar").
#'
#' @param hab1 Habitat(s) for species 1 (e.g., "Terrestrial", "Ter-Mar")
#' @param hab2 Habitat(s) for species 2 (e.g., "Marine", "Ter-Mar")
#' @return Logical (TRUE if overlap, FALSE otherwise)
check_habitat_overlap <- function(hab1, hab2) {
  # Return FALSE if either habitat is unknown
  if (is.na(hab1) | is.na(hab2)) {
    return(FALSE)
  }
  
  # Split multi-habitat strings into vectors
  hab1_list <- str_split(hab1, "-")[[1]]
  hab2_list <- str_split(hab2, "-")[[1]]
  
  # Check for any common elements
  return(any(hab1_list %in% hab2_list))
}


#' Calculate Interaction Probabilities
#'
#' Reconstructs a probabilistic food web for a given list of species.
#'
#' @param species_list A character vector of 'trophic_species' present.
#' @param traits_df The main 'traits' data frame.
#' @param a_param Intercept parameter from Nascimento et al. (2024).
#' @param b_param Slope parameter from Nascimento et al. (2024).
#' @return A data frame (tibble) of edges (from, to, Int_prob), or NULL.
calculate_interaction_probabilities <- function(species_list, traits_df, a_param, b_param) {
  
  traits_filtered <- traits_df %>%
    filter(trophic_species %in% species_list)
  
  # 1. Define potential consumers and resources
  consumers <- traits_filtered %>%
    filter(feeding_strategy %in% c("Carnivore", "Omnivore"))
  
  # *** UPDATED ***: Resources are herbivores, omnivores, OR 'carrion_only'
  resources <- traits_filtered %>%
    filter(feeding_strategy %in% c("Herbivore", "Omnivore") | carrion_only == 1)
  
  # Exit if no consumers or resources are present
  if (nrow(consumers) == 0 | nrow(resources) == 0) {
    return(NULL)
  }
  
  # 2. Create all possible consumer-resource pairs
  all_pairs <- expand_grid(
    Consumer = consumers$trophic_species,
    Resource = resources$trophic_species
  ) %>%
    # A species cannot consume itself
    filter(Consumer != Resource)
  
  if (nrow(all_pairs) == 0) {
    return(NULL)
  }
  
  # 3. Join trait data for pairs
  all_pairs <- all_pairs %>%
    left_join(
      traits_df %>% select(
        trophic_species, Consumer_BodySize = log_mass, 
        Consumer_Habitat = feeding_habitat, 
        Consumer_Scavenger = facultative_scavenger # *** NEW ***
      ),
      by = c("Consumer" = "trophic_species")
    ) %>%
    left_join(
      traits_df %>% select(
        trophic_species, Resource_BodySize = log_mass, 
        Resource_Habitat = feeding_habitat, 
        Resource_Carrion = carrion_only # *** NEW ***
      ),
      by = c("Resource" = "trophic_species")
    )
  
  # 4. Filter by Habitat Overlap
  # (This was the step you asked about - it's included here)
  all_pairs$Habitat_Overlap <- mapply(
    check_habitat_overlap, 
    all_pairs$Consumer_Habitat, 
    all_pairs$Resource_Habitat
  )
  
  interaction_edges <- all_pairs %>%
    filter(Habitat_Overlap == TRUE)
  
  if (nrow(interaction_edges) == 0) {
    return(NULL)
  }
  
  # 5. Calculate Interaction Probabilities
  interaction_edges <- interaction_edges %>%
    mutate(
      # Standard log-mass ratio
      log_ratio = Consumer_BodySize - Resource_BodySize,
      
      # Standard model probability (Logistic regression from Nascimento et al.)
      prob_model = exp(a_param + b_param * log_ratio) / (1 + exp(a_param + b_param * log_ratio)),
      
      # *** NEW SCAVENGING RULE ***
      # If consumer is a scavenger AND resource is carrion, force prob to 1.0
      # Otherwise, use the standard body-mass model probability.
      Int_prob = case_when(
        Consumer_Scavenger == 1 & Resource_Carrion == 1 ~ 1.0,
        TRUE ~ prob_model
      ),
      
      # Remove biologically impossible interactions (prob < 0.01)
      Int_prob = ifelse(Int_prob < 0.01, 0, Int_prob)
    ) %>%
    filter(Int_prob > 0) %>%
    select(from = Consumer, to = Resource, Int_prob)
  
  if (nrow(interaction_edges) == 0) {
    return(NULL)
  }
  
  return(interaction_edges)
}


#' Calculate Network Metrics
#'
#' Calculates a set of network metrics for a single binary adjacency matrix.
#'
#' @param web_matrix A binary adjacency matrix (rows=resources, cols=consumers).
#' @return A tibble with one row of network-level metrics.
calculate_network_metrics <- function(web_matrix) {
  
  # Ensure matrix has at least 2 rows and 2 columns
  if (nrow(web_matrix) < 2 | ncol(web_matrix) < 2) {
    return(NULL) 
  }
  
  # Remove empty rows/columns (species that are isolated)
  web_matrix <- web_matrix[rowSums(web_matrix) > 0, colSums(web_matrix) > 0, drop = FALSE]
  
  # Check again after pruning
  if (nrow(web_matrix) < 2 | ncol(web_matrix) < 2) {
    return(NULL)
  }
  
  tryCatch({
    # Calculate network-level metrics
    metrics <- networklevel(web_matrix, index = c(
      "connectance", "links per species", "cluster coefficient", 
      "nestedness", "robustness"
    ))
    
    # Calculate node-level metrics to get generality/vulnerability
    node_metrics <- specieslevel(web_matrix, index = c("generality", "vulnerability"))
    
    # Create the result tibble
    result <- tibble(
      n_consumers = ncol(web_matrix),
      n_resources = nrow(web_matrix),
      n_links = sum(web_matrix),
      connectance = metrics["connectance"],
      links_per_species = metrics["links per species"],
      cluster_coefficient = metrics["cluster coefficient"],
      nestedness = metrics["nestedness"],
      robustness_consumers = metrics["robustness.high"], # Robustness to consumer loss
      robustness_resources = metrics["robustness.low"],  # Robustness to resource loss
      mean_generality = mean(node_metrics$generality, na.rm = TRUE),
      mean_vulnerability = mean(node_metrics$vulnerability, na.rm = TRUE)
    )
    return(result)
    
  }, error = function(e) {
    print(paste("Error in metric calculation:", e$message))
    return(NULL)
  })
}


# ----------------------------------------------------------------------------
# PART 5: RECONSTRUCTION PARAMETERS
# Define model parameters and number of replicates.
# ----------------------------------------------------------------------------

# Parameters from Nascimento et al. (2024, Fig. 2b)
# logit(P(i,j)) = a + b * (log10(Mass_i) - log10(Mass_j))
PARAM_A <- -0.34
PARAM_B <- 1.63

# Number of replicated networks to generate per site/time bin
N_REPS <- 100


# ----------------------------------------------------------------------------
# PART 6: MAIN ANALYSIS LOOP
# Reconstruct networks for each locality and time bin.
# ----------------------------------------------------------------------------

print("Starting network reconstruction loop...")

# Get all unique localities to loop through
all_localities <- species_by_locality %>%
  select(time_bin, locality) %>%
  distinct()

# List to store all results
all_metrics_results <- list()
all_network_plots <- list()

for (i in 1:nrow(all_localities)) {
  
  current_locality <- all_localities$locality[i]
  current_time_bin <- all_localities$time_bin[i]
  
  cat(paste0("\nProcessing: ", current_locality, " (", current_time_bin, ") ...\n"))
  
  # 1. Get species for this locality
  current_species_list <- species_by_locality %>%
    filter(locality == current_locality, time_bin == current_time_bin) %>%
    pull(trophic_species)
  
  # Get all nodes present (for plotting)
  all_nodes <- traits %>%
    filter(trophic_species %in% current_species_list) %>%
    mutate(
      short_name = str_replace(trophic_species, "_", ". "),
      Category = class, # Use 'class' for color mapping
      is_human_node = (trophic_species == "Homo_sapiens")
    )
  
  # 2. Calculate interaction probabilities
  prob_edges <- calculate_interaction_probabilities(
    species_list = current_species_list,
    traits_df = traits,
    a_param = PARAM_A,
    b_param = PARAM_B
  )
  
  if (is.null(prob_edges) || nrow(prob_edges) == 0) {
    cat("  -> Skipping (no interactions predicted).\n")
    next
  }
  
  # 3. Run Replicates
  replicate_metrics <- list()
  for (j in 1:N_REPS) {
    
    # Generate a binary web based on probabilities
    binary_edges <- prob_edges %>%
      mutate(interacts = rbinom(n(), 1, Int_prob)) %>%
      filter(interacts == 1)
    
    if (nrow(binary_edges) == 0) {
      next
    }
    
    # Convert to adjacency matrix
    web_graph <- igraph::graph_from_data_frame(binary_edges, directed = TRUE, vertices = all_nodes)
    web_matrix <- igraph::as_adjacency_matrix(web_graph, sparse = FALSE, attr = NULL)
    
    # Bipartite functions expect (rows=resources, cols=consumers)
    consumers_names <- all_nodes %>% 
      filter(feeding_strategy %in% c("Carnivore", "Omnivore")) %>% pull(trophic_species)
    resources_names <- all_nodes %>% 
      filter(feeding_strategy %in% c("Herbivore", "Omnivore") | carrion_only == 1) %>% pull(trophic_species)
    
    # Get intersection of names that are in the matrix
    consumers_in_matrix <- intersect(consumers_names, colnames(web_matrix))
    resources_in_matrix <- intersect(resources_names, rownames(web_matrix))
    
    # Subset matrix
    bipartite_matrix <- web_matrix[resources_in_matrix, consumers_in_matrix, drop = FALSE]
    
    # Calculate metrics for this replicate
    metrics_rep <- calculate_network_metrics(bipartite_matrix)
    
    if (!is.null(metrics_rep)) {
      replicate_metrics[[j]] <- metrics_rep
    }
  }
  
  if (length(replicate_metrics) == 0) {
    cat("  -> Skipping (no links in replicates).\n")
    next
  }
  
  # 4. Aggregate metrics for this locality
  locality_metrics <- bind_rows(replicate_metrics) %>%
    mutate(
      locality = current_locality,
      time_bin = current_time_bin,
      .before = 1
    )
  
  all_metrics_results[[i]] <- locality_metrics
  
  # 5. Generate and save average network plot
  cat("  -> Generating network plot.\n")
  
  # Use tidygraph for plotting the *probabilistic* web
  plot_nodes <- all_nodes
  plot_edges <- prob_edges %>%
    mutate(
      # Highlight human interactions
      is_human_interaction = (from == "Homo_sapiens" | to == "Homo_sapiens")
    )
  
  # Check if there are any nodes left after filtering
  nodes_in_edges <- unique(c(plot_edges$from, plot_edges$to))
  plot_nodes_filtered <- plot_nodes %>% filter(trophic_species %in% nodes_in_edges)
  
  if (nrow(plot_nodes_filtered) == 0 || nrow(plot_edges) == 0) {
    cat("  -> Skipping plot (no nodes or edges after filtering).\n")
    next
  }
  
  tidy_web <- tbl_graph(nodes = plot_nodes_filtered, edges = plot_edges, directed = TRUE)
  
  # Create the ggraph plot
  net_plot <- ggraph(tidy_web, layout = 'kk') +
    geom_edge_fan(aes(alpha = Int_prob, color = Int_prob, filter = !is_human_interaction), width = 0.5) +
    # Highlight human interactions
    geom_edge_fan(aes(alpha = Int_prob, filter = is_human_interaction), color = "red", width = 1.2) +
    geom_node_point(aes(fill = Category, filter = !is_human_node), shape = 21, color = "black", size = 5, stroke = 0.5) +
    # Highlight human node
    geom_node_point(aes(fill = Category, filter = is_human_node), shape = 21, color = "black", size = 7, stroke = 1.5) +
    geom_node_text(aes(label = short_name), repel = TRUE, size = 2.5, max.overlaps = 35, bg.colour = "white", segment.color = 'grey50') +
    scale_fill_manual(values = taxa_colors, name = "Taxonomic Category", na.value = "grey50", limits = names(taxa_colors), drop = FALSE) +
    scale_edge_color_gradientn(colors = viridis::magma(256, direction = -1), name = "Interaction Prob.", limits = c(0,1)) +
    scale_edge_alpha(guide = 'none') +
    scale_y_continuous(expand = expansion(mult = 0.4)) + 
    theme_graph(base_family = 'sans') +
    theme(legend.position = "bottom", plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
    labs(
      title = paste("Probabilistic Food Web:", current_locality),
      subtitle = paste("Time Bin:", current_time_bin)
    )
  
  all_network_plots[[current_locality]] <- net_plot
  
  # Save the plot
  ggsave(
    filename = paste0("output/network_plot_", make_clean_names(current_locality), ".png"),
    plot = net_plot,
    width = 10,
    height = 10,
    dpi = 300
  )
}

print("...Reconstruction loop finished.")


# ----------------------------------------------------------------------------
# PART 7: AGGREGATE AND SAVE RESULTS
# Combine all replicate metrics into final summary tables.
# ----------------------------------------------------------------------------

print("Aggregating results...")

# Ensure 'output' directory exists
if (!dir.exists("output")) {
  dir.create("output")
}

# Bind all locality results together
final_metrics_all_reps <- bind_rows(all_metrics_results)

# Save the raw replicate data
write_csv(final_metrics_all_reps, "output/TDF_all_replicates_metrics.csv")

# Create a summary table (mean, sd, 95% CI) for each metric
metrics_summary <- final_metrics_all_reps %>%
  pivot_longer(
    cols = -(c(locality, time_bin)),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(time_bin, locality, metric) %>%
  summarise(
    n = n(),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    q025 = quantile(value, 0.025, na.rm = TRUE),
    q975 = quantile(value, 0.975, na.rm = TRUE)
  ) %>%
  ungroup()

# Save the summary table
write_csv(metrics_summary, "output/TDF_summary_metrics.csv")

print("...Results saved to 'output' folder.")


# ----------------------------------------------------------------------------
# PART 8: VISUALIZATION (SUMMARY PLOTS)
# Create summary boxplots comparing metrics across time bins.
# ----------------------------------------------------------------------------

print("Generating summary plots...")

# --- Plot 1: Connectance ---
p_conn <- final_metrics_all_reps %>%
  ggplot(aes(x = time_bin, y = connectance, fill = time_bin)) +
  geom_boxplot() +
  labs(title = "Network Connectance by Time Bin", x = "Time Bin", y = "Connectance") +
  theme_minimal() +
  theme(legend.position = "none")

# --- Plot 2: Links per Species ---
p_links <- final_metrics_all_reps %>%
  ggplot(aes(x = time_bin, y = links_per_species, fill = time_bin)) +
  geom_boxplot() +
  labs(title = "Links per Species by Time Bin", x = "Time Bin", y = "Links / Species") +
  theme_minimal() +
  theme(legend.position = "none")

# --- Plot 3: Robustness ---
p_robust <- final_metrics_all_reps %>%
  pivot_longer(
    cols = c(robustness_consumers, robustness_resources),
    names_to = "robustness_type",
    values_to = "value"
  ) %>%
  mutate(
    robustness_type = ifelse(
      robustness_type == "robustness_consumers", 
      "To Consumer Loss", "To Resource Loss"
    )
  ) %>%
  ggplot(aes(x = time_bin, y = value, fill = robustness_type)) +
  geom_boxplot() +
  labs(title = "Network Robustness by Time Bin", x = "Time Bin", y = "Robustness") +
  theme_minimal() +
  theme(legend.position = "bottom")

# --- Combine plots ---
summary_plot <- (p_conn | p_links) / p_robust +
  plot_annotation(
    title = "TDF Food Web Metrics Comparison",
    tag_levels = 'A'
  )

# Save the summary plot
ggsave(
  filename = "output/TDF_summary_plots.png",
  plot = summary_plot,
  width = 12,
  height = 8,
  dpi = 300
)

print("--- Analysis Complete ---")



