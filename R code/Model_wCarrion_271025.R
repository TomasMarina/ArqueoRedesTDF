# ----------------------------------------------------------------------------
#
# R SCRIPT FOR RECONSTRUCTING ANCIENT CONSUMER-RESOURCE NETWORKS
# (Tierra del Fuego Archaeological Study)
#
# Version: 301025-CorrectedPlots
#
# *** SCRIPT CORRECTIONS ***
# 1. PART 9 (Boxplots): Now plots the raw 'final_metrics_all_reps' data (all 100
#    replicates) and uses 'facet_grid' to match the desired PDF layout.
# 2. PART 10 (Network Plots): Now correctly creates a true bipartite graph
#    by DUPLICATING nodes (e.g., omnivores) that exist in both the
#    consumer and resource layers, matching the desired PDF layout.
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
Sys.setlocale("LC_TIME", "en_US.UTF-8") # Ensure consistent date formatting

# Get today's date for filenames (DDMMYY format)
today_date <- format(Sys.Date(), "%d%m%y")


# ----------------------------------------------------------------------------
# PART 2: DATA LOADING AND INITIAL CLEANING
# Read the raw data files.
# ----------------------------------------------------------------------------

print("Loading data...")

# Use file paths from the original script
localities_raw <- read_excel("data - Santiago2025/Species localities - TDF - 030925.xlsx")
traits_raw <- read_excel("data - Santiago2025/Species traits - TDF - 241025.xlsx", sheet = "traits")
ages_raw <- read_excel("data - Santiago2025/Age localities - TDF - 030925.xlsx")


# ----------------------------------------------------------------------------
# PART 3: DATA PREPARATION
# Clean and prepare the three data sources.
# ----------------------------------------------------------------------------

print("Preparing data...")

# --- Clean Traits Data ---
traits <- traits_raw %>%
  clean_names() %>%
  # Handle potential NAs in new columns
  mutate(
    facultative_scavenger = ifelse(is.na(facultative_scavenger), 0, facultative_scavenger),
    carrion_only = ifelse(is.na(carrion_only), 0, carrion_only)
  ) %>%
  # Calculate log10 body size
  mutate(
    # Use geometric mean if both min and max are available
    avg_bodysize_kg = exp((log(body_size_min) + log(body_size_max)) / 2),
    log10_bodysize = log10(avg_bodysize_kg)
  ) %>%
  # Select only the columns we need
  select(
    trophic_species, class, order, family,
    log10_bodysize, feeding_strategy,
    feeding_habitat, facultative_scavenger, carrion_only
  ) %>%
  # Ensure no duplicate species
  distinct(trophic_species, .keep_all = TRUE)

# --- Clean Localities Data ---
# This block pivots the "wide" format of the localities file.
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
ages <- ages_raw %>%
  clean_names() %>%
  # Fix case-sensitivity to allow joining with pivoted locality names
  mutate(locality = tolower(locality)) %>%
  # The column in your file is "Geological_ages"
  select(locality, time_bin = geological_ages, region_col) %>%
  distinct()

# *** KEY CHANGE HERE ***
# --- Create final species-by-AGGREGATED-unit dataset ---
# This merges all localities that share the same time_bin and region_col.
splocs_agg <- localities %>%
  inner_join(ages, by = "locality") %>%
  # Group by the desired spatio-temporal bins, NOT by locality
  group_by(time_bin, region_col) %>%
  # Summarise the species list, collecting all unique species from all localities in that bin
  summarise(
    species_list = list(unique(trophic_species)),
    # Optional: keep a list of the localities that were merged
    localities_merged_count = n_distinct(locality)
  ) %>%
  ungroup() %>%
  # Handle potential NAs in time_bin or region_col if they exist
  filter(!is.na(time_bin), !is.na(region_col))

print("Data prepared.")
print(paste(nrow(splocs_agg), "unique aggregated time-bin/region combinations found."))
print(splocs_agg)


# ----------------------------------------------------------------------------
# PART 4: DEFINE MODEL PARAMETERS
# This function generates the STOCHASTIC (randomized) parameters
# from your 'Space&Time_11025.R' script.
# It will be called ONCE PER REPLICATE.
# ----------------------------------------------------------------------------

get_stochastic_params <- function() {
  params <- list(
    # Parameters for Carnivores (P)
    alphaP = runif(1, 1, 2),
    betaP = runif(1, -2, -1),
    gammaP = runif(1, -1, -0.5),
    
    # Parameters for Omnivores (O)
    alphaO = runif(1, -6, -4),
    betaO = runif(1, -3, -2),
    gammaO = 0 # As specified in your original script
  )
  return(params)
}

print("Stochastic parameter function (get_stochastic_params) loaded.")


# ----------------------------------------------------------------------------
# PART 5: RECONSTRUCTION FUNCTION
# This function reconstructs a single probabilistic network.
# ----------------------------------------------------------------------------

reconstruct_network <- function(sp_traits, model_params, focal_habitat = "Ter-Mar") {
  
  # --- 1. Define Consumer and Resource Pools ---
  
  # Consumers are Carnivores and Omnivores
  consumers <- sp_traits %>%
    filter(feeding_strategy %in% c("Carnivore", "Omnivore"))
  
  # Resources are Herbivores, Omnivores, AND carrion-only species
  resources <- sp_traits %>%
    filter(feeding_strategy %in% c("Herbivore", "Omnivore") | carrion_only == 1)
  
  # Filter by habitat if specified
  if (focal_habitat != "Ter-Mar") {
    consumers <- consumers %>% filter(feeding_habitat == focal_habitat | feeding_habitat == "Ter-Mar")
    resources <- resources %>% filter(feeding_habitat == focal_habitat | feeding_habitat == "Ter-Mar")
  }
  
  # Create all possible consumer-resource pairs
  all_pairs <- expand_grid(
    consumer_species = consumers$trophic_species,
    resource_species = resources$trophic_species
  ) %>%
    # Remove self-loops (e.g., Omnivore-Omnivore)
    filter(consumer_species != resource_species) %>%
    # Join traits
    left_join(consumers %>% select(consumer_species = trophic_species, consumer = feeding_strategy, log10_c = log10_bodysize, habitat_c = feeding_habitat, scav_c = facultative_scavenger), by = "consumer_species") %>%
    left_join(resources %>% select(resource_species = trophic_species, resource = feeding_strategy, log10_r = log10_bodysize, habitat_r = feeding_habitat, carrion_r = carrion_only), by = "resource_species")
  
  if (nrow(all_pairs) == 0) {
    return(data.frame(Consumer = character(), Resource = character(), Int_prob = numeric()))
  }
  
  # --- 2. Calculate Interaction Probabilities ---
  
  # This part iterates row-by-row to apply the complex rules
  results_list <- list()
  
  for (i in 1:nrow(all_pairs)) {
    pair <- all_pairs[i, ]
    
    # --- A. Check Habitat Overlap ---
    habitat_overlap <- if_else(
      pair$habitat_c == pair$habitat_r |
        pair$habitat_c == "Ter-Mar" |
        pair$habitat_r == "Ter-Mar",
      1.0,
      0.0
    )
    
    # If no habitat overlap, skip to next pair
    if (habitat_overlap == 0) {
      results_list[[i]] <- data.frame(
        Consumer = pair$consumer_species,
        Resource = pair$resource_species,
        Int_prob = 0.0
      )
      next
    }
    
    # --- B. Apply Interaction Model (Habitat Overlaps) ---
    
    prob <- 0.0 # Initialize probability
    
    # *** MERGED MODEL LOGIC ***
    
    # Rule 1: Special scavenging interaction (overrides body mass)
    if (pair$scav_c == 1 && pair$carrion_r == 1) {
      
      prob <- 1.0 # Force this interaction
      
      # Rule 2: Check for missing body size data *before* comparison
      # If either body size is NA, we can't use the model, so prob = 0
    } else if (is.na(pair$log10_c) || is.na(pair$log10_r)) {
      
      prob <- 0.0 # Cannot model interaction if body size is missing
      
      # Rule 3: Standard interaction (falls back to STOCHASTIC body mass model)
      # This block is now safe because we've already checked for NAs
    } else if (pair$log10_c > pair$log10_r) {
      
      log_ratio <- pair$log10_c - pair$log10_r
      
      # Use the correct parameters based on consumer strategy
      if (pair$consumer == "Carnivore") {
        # Use Carnivore (P) parameters
        z <- model_params$alphaP + model_params$betaP * log_ratio + model_params$gammaP * (log_ratio^2)
        prob <- exp(z) / (1 + exp(z))
        
      } else { # This handles "Omnivore"
        # Use Omnivore (O) parameters
        z <- model_params$alphaO + model_params$betaO * log_ratio + model_params$gammaO * (log_ratio^2)
        prob <- exp(z) / (1 + exp(z))
      }
      
      # Rule 4: No interaction (consumer < resource or not meeting other rules)
    } else {
      prob <- 0.0
    }
    
    # --- C. Store Result ---
    results_list[[i]] <- data.frame(
      Consumer = pair$consumer_species,
      Resource = pair$resource_species,
      Int_prob = prob # Final probability (habitat was already checked)
    )
  }
  
  # Combine all results
  final_links <- bind_rows(results_list) %>%
    filter(Int_prob > 0)
  
  return(final_links)
}


# ----------------------------------------------------------------------------
# PART 6: METRIC CALCULATION
# Helper functions to calculate metrics for a single *weighted* network.
# ----------------------------------------------------------------------------

# Function to calculate network metrics using 'bipartite'
calculate_metrics <- function(adj_matrix) {
  
  # *** UPDATED ***: NO binarization. Use the raw probability matrix.
  
  # Remove empty rows/columns
  adj_matrix_weighted <- adj_matrix[rowSums(adj_matrix) > 0, , drop = FALSE]
  adj_matrix_weighted <- adj_matrix_weighted[, colSums(adj_matrix_weighted) > 0, drop = FALSE]
  
  if (nrow(adj_matrix_weighted) < 2 | ncol(adj_matrix_weighted) < 2) {
    # Network is too small
    return(
      tibble(
        n_consumers = nrow(adj_matrix_weighted),
        n_resources = ncol(adj_matrix_weighted),
        n_nodes = nrow(adj_matrix_weighted) + ncol(adj_matrix_weighted),
        connectance = NA, links_per_species = NA, 
        weighted_nestedness = NA, modularity = NA, 
        generality = NA, vulnerability = NA
      )
    )
  }
  
  # Calculate metrics using bipartite
  # Note: Modularity calculation can sometimes fail
  mod <- try(computeModules(adj_matrix_weighted), silent = TRUE)
  mod_value <- if (inherits(mod, "try-error")) NA else mod@likelihood
  
  # *** UPDATED ***: Corrected metric calls
  metrics <- networklevel(
    adj_matrix_weighted,
    index = c("connectance", "links per species", "weighted nestedness")
  )
  
  # *** UPDATED ***: Generality and Vulnerability are separate functions
  # We take the mean to get a network-level value
  gen <- try(mean(generality.HL(adj_matrix_weighted), na.rm = TRUE), silent = TRUE)
  vul <- try(mean(vulnerability.LL(adj_matrix_weighted), na.rm = TRUE), silent = TRUE)
  
  gen_value <- if (inherits(gen, "try-error")) NA else gen
  vul_value <- if (inherits(vul, "try-error")) NA else vul
  
  # Return metrics as a tibble
  return(tibble(
    n_consumers = nrow(adj_matrix_weighted),
    n_resources = ncol(adj_matrix_weighted),
    n_nodes = nrow(adj_matrix_weighted) + ncol(adj_matrix_weighted),
    connectance = metrics["connectance"],
    links_per_species = metrics["links per species"],
    weighted_nestedness = metrics["weighted nestedness"],
    modularity = mod_value,
    generality = gen_value,
    vulnerability = vul_value
  ))
}


# ----------------------------------------------------------------------------
# PART 7: MAIN ANALYSIS LOOP
# Run the analysis 'n_reps' times.
# Loops over aggregated bins, not localities.
# ----------------------------------------------------------------------------

n_reps <- 100 # Number of replicates (set to 100)

# This function runs the full analysis for ONE replicate
run_analysis <- function(splocs_table, traits) {
  
  # Get ONE set of random parameters for this ENTIRE replicate
  model_params <- get_stochastic_params()
  
  # This iterates over each row of the AGGREGATED splocs_agg tibble
  metrics_per_bin <- map_dfr(1:nrow(splocs_table), function(i) {
    
    # Get the time_bin, region, and species list for this row
    focal_time_bin <- splocs_table$time_bin[i]
    focal_region <- splocs_table$region_col[i]
    species_names <- splocs_table$species_list[[i]]
    
    # Get the trait data for *only* the species in this aggregated bin
    sp_traits <- traits %>%
      filter(trophic_species %in% species_names)
    
    # Reconstruct the network for this bin
    # We pass the STOCHASTIC model_params for this replicate
    net_links <- reconstruct_network(sp_traits, model_params, focal_habitat = "Ter-Mar")
    
    # Convert to an adjacency matrix for bipartite
    consumers <- unique(net_links$Consumer)
    resources <- unique(net_links$Resource)
    
    adj_matrix <- matrix(0, 
                         nrow = length(consumers), 
                         ncol = length(resources),
                         dimnames = list(consumers, resources))
    
    # Fill the matrix with probabilities
    if(nrow(net_links) > 0) {
      for (k in 1:nrow(net_links)) {
        adj_matrix[net_links$Consumer[k], net_links$Resource[k]] <- net_links$Int_prob[k]
      }
    }
    
    # Calculate metrics on the weighted matrix
    metrics_tibble <- calculate_metrics(adj_matrix)
    
    # Add identifying info and return
    return(
      metrics_tibble %>%
        mutate(
          time_bin = focal_time_bin,
          region_col = focal_region
        )
    )
  })
  
  return(metrics_per_bin)
}

# --- Run the main loop ---
print(paste("Running analysis for", n_reps, "replicates..."))
start_time <- Sys.time()

# Use map_dfr to run the analysis 'n_reps' times
# We pass the AGGREGATED 'splocs_agg' tibble
final_metrics_all_reps <- map_dfr(
  1:n_reps,
  ~run_analysis(splocs_agg, traits), # Pass the aggregated splocs
  .id = "rep"
)

end_time <- Sys.time()
print(paste("Analysis complete. Time taken:", round(difftime(end_time, start_time, units = "secs"), 1), "seconds"))

# ----------------------------------------------------------------------------
# PART 8: SUMMARIZE METRICS
# Aggregate the results from all replicates.
# Summarizes by bin, not locality.
# ----------------------------------------------------------------------------

print("Summarizing metrics...")

# This block calculates the mean, sd, and 95% CIs for all metrics,
# grouped by time_bin and region_col.
metrics_summary <- final_metrics_all_reps %>%
  pivot_longer(
    cols = -(c(rep, time_bin, region_col)), 
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(time_bin, region_col, metric) %>% 
  summarise(
    n = n(),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    q025 = quantile(value, 0.025, na.rm = TRUE),
    q975 = quantile(value, 0.975, na.rm = TRUE)
  ) %>%
  ungroup()

# --- Save Results ---
output_dir <- "results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

write_csv(
  final_metrics_all_reps, 
  file.path(output_dir, paste0(today_date, "_TDF_all_replicates_metrics_AGGREGATED.csv"))
)
write_csv(
  metrics_summary, 
  file.path(output_dir, paste0(today_date, "_TDF_summary_metrics_AGGREGATED.csv"))
)

print(paste("Results saved to", output_dir, "directory."))


# ----------------------------------------------------------------------------
# PART 9: PLOT METRICS
# *** COMPLETELY REVISED ***
# This part now plots all 100 replicates from 'final_metrics_all_reps'
# and uses 'facet_grid' to create a plot matrix similar to the desired PDF.
# ----------------------------------------------------------------------------

print("Generating plots...")

# --- 1. Define the correct order for time bins ---
time_order <- c("Pleistoceno final", "Holoceno medio", "Holoceno final", "Histórico")

# --- 2. Prepare the data for plotting ---
plot_data_metrics <- final_metrics_all_reps %>%
  # Pivot all metrics into a single column
  pivot_longer(
    cols = c(connectance, links_per_species, weighted_nestedness, modularity, generality, vulnerability),
    names_to = "metric",
    values_to = "value"
  ) %>%
  # Ensure factors are in the correct order for plotting
  mutate(
    time_bin = factor(time_bin, levels = time_order),
    metric = factor(metric, levels = c("connectance", "links_per_species", "weighted_nestedness", "modularity", "generality", "vulnerability"))
  ) %>%
  # Remove any NA values that would break plotting
  filter(!is.na(value))

# --- 3. Create the multi-faceted plot ---
combined_plot <- ggplot(plot_data_metrics, aes(x = time_bin, y = value, fill = region_col)) +
  # Use geom_boxplot, which will show the distribution of the 100 replicates
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  
  # *** This is the key change ***
  # Use facet_grid to create rows of metrics and columns of regions
  # scales = "free_y" allows each metric to have its own y-axis
  facet_grid(metric ~ region_col, scales = "free_y", switch = "y") +
  
  # Add labels and titles
  labs(
    title = "Tierra del Fuego Aggregated Network Metrics",
    subtitle = "Boxplots show distribution of 100 stochastic replicates",
    x = "Geological Age",
    y = "Metric Value",
    fill = "Region"
  ) +
  
  # Apply themes
  theme_minimal() +
  theme(
    legend.position = "bottom",
    # Rotate x-axis labels
    axis.text.x = element_text(angle = 45, hjust = 1),
    # Move the metric labels (on the left) to be outside the plot
    strip.placement = "outside",
    strip.text.y = element_text(angle = 0, face = "bold"),
    strip.text.x = element_text(face = "bold"),
    # Add a border
    panel.border = element_rect(color = "grey80", fill = NA)
  )

# --- 4. Save the plot ---
ggsave(
  file.path(output_dir, paste0(today_date, "_TDF_metrics_boxplots_AGGREGATED.png")),
  combined_plot,
  width = 8, # Narrower plot is better for this grid
  height = 12, # Taller plot to accommodate all metric rows
  dpi = 300,
  bg = "white"
)

print(paste("Metric boxplots saved to", output_dir, "directory."))


# ----------------------------------------------------------------------------
# PART 10: GENERATE INDIVIDUAL NETWORK PLOTS
# *** COMPLETELY REVISED ***
# This part now generates a true bipartite plot by duplicating nodes
# that are both consumers and resources (e.g., omnivores).
# ----------------------------------------------------------------------------

print("Generating individual network plots...")

# --- 1. Helper function to shorten names ---
shorten_name <- function(species_name) {
  parts <- str_split(species_name, "_")[[1]]
  if (length(parts) == 2) {
    # e.g., "Homo_sapiens" -> "H. sapiens"
    return(paste0(str_sub(parts[1], 1, 1), ". ", parts[2]))
  } else {
    # e.g., "Cricetidae" -> "Cricetidae"
    return(species_name)
  }
}

# --- 2. Define color palette ---
taxa_colors <- c(
  "Mammalia" = "#7f3b08",
  "Aves" = "#b35806",
  "Actinopterygii" = "#e08214",
  "Chondrichthyes" = "#fdb863",
  "Malacostraca" = "#fee0b6",
  "Bivalvia" = "#d8daeb",
  "Gastropoda" = "#b2abd2",
  "Cephalopoda" = "#8073ac",
  "Echinoidea" = "#542788",
  "Other" = "grey50"
)

# --- 3. Define mean parameters for representative plot ---
mean_params <- list(
  alphaP = mean(c(1, 2)),
  betaP = mean(c(-2, -1)),
  gammaP = mean(c(-1, -0.5)),
  alphaO = mean(c(-6, -4)),
  betaO = mean(c(-3, -2)),
  gammaO = 0
)

# --- 4. Define the main plotting function ---
# This function is now designed for a true bipartite layout
plot_network_bipartite <- function(graph, title) {
  
  # Create the bipartite layout
  layout <- create_layout(graph, layout = 'bipartite')
  
  # Manually set y-coordinates: resources (type=FALSE) at y=0, consumers (type=TRUE) at y=1
  layout$y <- ifelse(layout$type, 1, 0)
  
  ggraph(layout) +
    # Draw edges
    # Non-human interactions
    geom_edge_fan(aes(alpha = Int_prob, color = Int_prob, filter = !is_human_interaction), width = 0.5) +
    # Human interactions (plotted on top)
    geom_edge_fan(aes(alpha = Int_prob, filter = is_human_interaction), color = "red", width = 1.0) +
    
    # Draw nodes
    # Non-human nodes
    geom_node_point(aes(fill = Category, filter = !is_human_node), shape = 21, color = "black", size = 5, stroke = 0.5) +
    # Human nodes (both consumer and resource)
    geom_node_point(aes(fill = Category, filter = is_human_node), shape = 21, color = "red", size = 7, stroke = 1.5) +
    
    # Draw labels
    geom_node_text(aes(label = short_name), repel = TRUE, size = 2.5, max.overlaps = 35, 
                   bg.colour = "white", bg.r = 0.1, segment.color = 'grey50') +
    
    # Scales
    scale_fill_manual(values = taxa_colors, name = "Taxonomic Category", na.value = "grey50", limits = names(taxa_colors), drop = FALSE) +
    scale_edge_color_gradientn(colors = viridis::magma(256, direction = -1, begin = 0.1), name = "Interaction Prob.", limits = c(0,1)) +
    scale_edge_alpha(guide = 'none') +
    scale_y_continuous(expand = expansion(mult = 0.4)) + # Add space for labels
    
    # Theme
    theme_graph(base_family = 'sans') +
    theme(
      legend.position = "bottom",
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      plot.title = element_text(hjust = 0.5, size = 16)
    ) +
    labs(title = title)
}


# --- 5. Prepare data and loop through plots ---

# Create a data frame with all data needed for plotting
plot_data_prep <- splocs_agg %>%
  mutate(
    # Get the species traits for each bin
    sp_traits = map(species_list, ~traits %>% filter(trophic_species %in% .x)),
    
    # --- *** CRITICAL NEW LOGIC FOR BIPARTITE PLOTS *** ---
    
    # Generate the single representative link list
    links_df = map(sp_traits, ~reconstruct_network(.x, mean_params) %>%
                     filter(Int_prob > 0) %>%
                     mutate(
                       # Create 'from' and 'to' names for the duplicated nodes
                       from = paste0(Consumer, "_C"), # e.g., "Homo_sapiens_C"
                       to = paste0(Resource, "_R"),   # e.g., "Lama_guanicoe_R"
                       # Identify human interactions for special plotting
                       is_human_interaction = (Consumer == "Homo_sapiens" | Resource == "Homo_sapiens")
                     )
    ),
    
    # Create the nodes data frame by duplicating omnivores
    nodes_df = map2(sp_traits, links_df, function(st, links) {
      
      # Define consumer and resource names *from the links*
      consumer_names <- unique(links$Consumer)
      resource_names <- unique(links$Resource)
      
      # 1. Create consumer nodes (top layer)
      nodes_C <- st %>%
        filter(trophic_species %in% consumer_names) %>%
        mutate(
          name = paste0(trophic_species, "_C"), # Duplicated name
          type = TRUE # Bipartite layout type (top layer)
        )
      
      # 2. Create resource nodes (bottom layer)
      nodes_R <- st %>%
        filter(trophic_species %in% resource_names) %>%
        mutate(
          name = paste0(trophic_species, "_R"), # Duplicated name
          type = FALSE # Bipartite layout type (bottom layer)
        )
      
      # 3. Combine them. Omnivores will now appear in both dataframes.
      bind_rows(nodes_C, nodes_R) %>%
        mutate(
          # Create 'Category' for colors
          Category = if_else(class %in% names(taxa_colors), class, "Other"),
          # Create short names for labels (from original species name)
          short_name = map_chr(trophic_species, shorten_name),
          # Identify human nodes
          is_human_node = (trophic_species == "Homo_sapiens")
        ) %>%
        # Ensure 'name' is the first column for tbl_graph
        select(name, short_name, type, Category, is_human_node, everything())
    }),
    
    # --- End of new logic ---
    
    # Create title and filename
    plot_title = paste0("Aggregated Network: ", time_bin, " - ", region_col),
    file_name = file.path(output_dir, paste0(today_date, "_NetworkPlot_", time_bin, "_", region_col, ".png"))
  )

# Use pwalk to iterate and save each plot
pwalk(plot_data_prep, function(nodes_df, links_df, plot_title, file_name, ...) {
  
  print(paste("Plotting:", plot_title))
  
  # Check if there are any nodes or links
  if (nrow(nodes_df) == 0 || nrow(links_df) == 0) {
    print(paste("Skipping plot for", plot_title, "due to no nodes or links."))
    return()
  }
  
  # Create graph object
  g <- tbl_graph(nodes = nodes_df, edges = links_df, directed = TRUE)
  
  # Generate the plot using the new bipartite function
  final_plot <- plot_network_bipartite(g, plot_title)
  
  # Save the plot
  ggsave(
    filename = file_name,
    plot = final_plot,
    width = 14,
    height = 10,
    dpi = 300,
    bg = "white"
  )
})

print(paste("Network plots saved to", output_dir, "directory."))
print("--- SCRIPT FINISHED ---")
