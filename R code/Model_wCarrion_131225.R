# ----------------------------------------------------------------------------
#
# R SCRIPT FOR RECONSTRUCTING ANCIENT CONSUMER-RESOURCE NETWORKS
# (Tierra del Fuego Archaeological Study)
#
# Version: 131225-Final-BiplotFix-Complete
#
# *** SCRIPT CORRECTIONS ***
# 1. PART 1 (SETUP): Added 'tidytext' package for reorder_within function.
# 2. PART 3 (DATA PREP): Robust cleaning for 1/0/NA columns (Kept).
# 3. PART 5 (reconstruct_network):
#    - Filter 'Int_prob > 0'.
#    - Filter to prevent humans as prey.
#    - Filter to prevent 'carrion_only' species from being consumers.
#    - Updated human predation logic: 'human_r > 0' captures categories 1 and 2.
# 4. PART 6 (calculate_metrics): Correct calculation using pivoted matrix.
# 5. PART 10 (Network Plots): Correct colors and layout.
# 6. PART 12 (Biplot): 
#    - Correctly sorts bars High -> Low within each facet using tidytext.
#    - Homo sapiens colored #33CCFF, others steelblue.
#
# ----------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# PART 1: SETUP
# Load all necessary packages for the analysis.
# ----------------------------------------------------------------------------

# Ensure all packages are installed:
# install.packages(c("tidyverse", "readxl", "janitor", "bipartite", "igraph", 
#                    "tidygraph", "ggraph", "ggrepel", "patchwork", "viridis", 
#                    "ggnewscale", "tidytext"))

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
library(ggnewscale)
library(tidytext)      # *** ADDED for reorder_within ***

set.seed(123) # for reproducibility
Sys.setlocale("LC_TIME", "en_US.UTF-8") # Ensure consistent date formatting

# Get today's date for filenames (DDMMYY format)
today_date <- format(Sys.Date(), "%d%m%y")


# ----------------------------------------------------------------------------
# PART 2: DATA LOADING AND INITIAL CLEANING
# Read the raw data files.
# ----------------------------------------------------------------------------

print("Loading data...")

# Use file paths exactly as provided by user
localities_raw <- read_excel("data - Santiago2025/Species localities - TDF  - 131225.xlsx")
traits_raw <- read_excel("data - Santiago2025/Species traits - TDF - 131225.xlsx", sheet = "traits")
ages_raw <- read_excel("data - Santiago2025/Age localities - TDF  - 081125.xlsx")


# ----------------------------------------------------------------------------
# PART 3: DATA PREPARATION
# Clean and prepare the three data sources.
# ----------------------------------------------------------------------------

print("Preparing data...")

# --- Clean Traits Data ---
traits <- traits_raw %>%
  clean_names() %>%
  # *** ROBUST CLEANING FOR 1/0/BLANK COLUMNS ***
  mutate(
    # Force to numeric first. This turns "1" -> 1 and ""/blanks -> NA
    facultative_scavenger = as.numeric(facultative_scavenger),
    carrion_only = as.numeric(carrion_only),
    human_resource = as.numeric(human_resource),
    
    # Now, replace the resulting NAs with 0
    facultative_scavenger = ifelse(is.na(facultative_scavenger), 0, facultative_scavenger),
    carrion_only = ifelse(is.na(carrion_only), 0, carrion_only),
    human_resource = ifelse(is.na(human_resource), 0, human_resource)
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
    feeding_habitat, facultative_scavenger, carrion_only, human_resource
  ) %>%
  # Ensure no duplicate species
  distinct(trophic_species, .keep_all = TRUE)

# --- Clean Localities Data ---
# This block pivots the "wide" format of the localities file.
# *** UPDATED: Robust handling of site names and values to keep rare species ***
localities <- localities_raw %>%
  clean_names() %>%
  pivot_longer(
    cols = -trophic_species,
    names_to = "locality",
    values_to = "nisp_raw", # Store the raw value first
    values_transform = list(nisp_raw = as.character)
  ) %>%
  mutate(
    # Safely convert to numeric, handling potential whitespace
    nisp = as.numeric(str_trim(nisp_raw)),
    # Clean locality names to ensure matching (remove extra spaces if any)
    locality = str_trim(locality) 
  ) %>%
  # Filter: Keep if NISP is a number AND greater than 0
  filter(!is.na(nisp) & nisp > 0) %>%
  # Select only the columns we need
  select(locality, trophic_species) %>%
  # Merge with traits to ensure we only keep species with trait data
  inner_join(traits %>% select(trophic_species), by = "trophic_species")

# --- Clean Ages Data ---
ages <- ages_raw %>%
  clean_names() %>%
  # Fix case-sensitivity and whitespace to allow joining with pivoted locality names
  mutate(
    # The locality names in ages file might need cleaning too
    locality = tolower(str_trim(locality))
  ) %>%
  # The column in your file is "Geological_ages"
  select(locality, time_bin = geological_ages, region_col) %>%
  distinct()

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
  filter(!is.na(time_bin), !is.na(region_col)) %>%
  
  # *** THE KEY FIX ***
  # Manually add "Homo_sapiens" to every aggregated species list
  mutate(
    species_list = map(species_list, ~ unique(c(.x, "Homo_sapiens")))
  )

print("Data prepared.")
print(paste(nrow(splocs_agg), "unique aggregated time-bin/region combinations found."))


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
  # Explicitly exclude 'carrion_only' species from being consumers
  consumers <- sp_traits %>%
    filter(
      feeding_strategy %in% c("Carnivore", "Omnivore") &
        carrion_only == 0 # Ensure carrion species (e.g. whales) don't act as consumers
    )
  
  # Resources are Herbivores, Omnivores, AND carrion-only species
  # Filter to explicitly remove Homo_sapiens from being a resource.
  resources <- sp_traits %>%
    filter(
      (feeding_strategy %in% c("Herbivore", "Omnivore") | carrion_only == 1) &
        trophic_species != "Homo_sapiens"
    )
  
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
    # Join traits for resources (including human_resource flag)
    left_join(resources %>% select(resource_species = trophic_species, resource = feeding_strategy, log10_r = log10_bodysize, habitat_r = feeding_habitat, carrion_r = carrion_only, human_r = human_resource), by = "resource_species")
  
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
    
    # --- *** NEW REVISED MODEL LOGIC *** ---
    
    # Rule 1: Special human scavenging interaction (overrides body mass)
    # This now correctly checks numeric 1
    if (pair$consumer_species == "Homo_sapiens" && pair$scav_c == 1 && pair$carrion_r == 1) {
      
      prob <- 1.0 # Force this interaction
      
      # Rule 2: Special human PREDATION interaction (overrides body mass CHECK)
      # *** UPDATED: Checks for human_r > 0 (i.e., 1 or 2) ***
    } else if (pair$consumer_species == "Homo_sapiens" && pair$human_r > 0) {
      
      # This is human hunting. Ignore mass check and run the model.
      # Check for NA mass first
      if (is.na(pair$log10_c) || is.na(pair$log10_r)) {
        prob <- 0.0
      } else {
        log_ratio <- pair$log10_c - pair$log10_r
        # Use Omnivore (O) parameters for humans
        z <- model_params$alphaO + model_params$betaO * log_ratio + model_params$gammaO * (log_ratio^2)
        prob <- exp(z) / (1 + exp(z))
      }
      
      # Rule 3: Check for missing body size data *before* comparison
    } else if (is.na(pair$log10_c) || is.na(pair$log10_r)) {
      
      prob <- 0.0 # Cannot model interaction if body size is missing
      
      # Rule 4: Standard non-human interaction (falls back to body mass model)
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
      
      # Rule 5: No interaction (consumer < resource or not meeting other rules)
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
    # *** Reverted to prob > 0, as in old script
    filter(Int_prob > 0)
  
  return(final_links)
}


# ----------------------------------------------------------------------------
# PART 6: METRIC CALCULATION
# Helper functions to calculate metrics for a single *weighted* network.
# ----------------------------------------------------------------------------

# Function to calculate network metrics using 'bipartite'
calculate_metrics <- function(links_df) {
  
  # 1. Build the incidence matrix directly from the edge list
  #    Using quantitative values (Int_prob)
  incidence_matrix <- links_df %>%
    select(Consumer, Resource, Int_prob) %>%
    # Ensure uniqueness in case of dupes (shouldn't happen, but safe)
    group_by(Consumer, Resource) %>%
    summarise(Int_prob = max(Int_prob), .groups = "drop") %>%
    pivot_wider(
      names_from = Resource,
      values_from = Int_prob,
      values_fill = 0
    ) %>%
    column_to_rownames("Consumer") %>%
    as.matrix()
  
  # 2. Check dimensions (must have at least 2 rows and 2 cols)
  if (nrow(incidence_matrix) < 2 || ncol(incidence_matrix) < 2) {
    return(
      tibble(
        n_consumers = nrow(incidence_matrix),
        n_resources = ncol(incidence_matrix),
        n_nodes = nrow(incidence_matrix) + ncol(incidence_matrix),
        connectance = NA, links_per_species = NA, 
        weighted_nestedness = NA, modularity = NA, 
        generality = NA, vulnerability = NA
      )
    )
  }
  
  # 3. Transpose: bipartite expects (Lower Level / Resources) x (Higher Level / Consumers)
  #    Our incidence matrix is (Consumers x Resources), so we transpose it.
  transposed_matrix <- t(incidence_matrix)
  
  # 4. Calculate Basic Counts
  num_consumers <- ncol(transposed_matrix)
  num_resources <- nrow(transposed_matrix)
  total_species <- num_consumers + num_resources
  
  # 5. Calculate Metrics using networklevel()
  #    Using the specific indices requested
  level_metrics <- tryCatch(
    networklevel(transposed_matrix, 
                 index = c("connectance", "links per species", "weighted nestedness", "generality", "vulnerability")),
    error = function(e) return(NULL)
  )
  
  # 6. Calculate Modularity
  modules <- tryCatch(computeModules(transposed_matrix), error = function(e) NULL)
  modularity_q <- if (!is.null(modules)) modules@likelihood else NA
  
  # 7. Extract values safely
  if (is.null(level_metrics)) {
    return(tibble(
      n_consumers = num_consumers, n_resources = num_resources, n_nodes = total_species,
      connectance = NA, links_per_species = NA, weighted_nestedness = NA, modularity = modularity_q, 
      generality = NA, vulnerability = NA
    ))
  }
  
  return(tibble(
    n_consumers = num_consumers,
    n_resources = num_resources,
    n_nodes = total_species,
    connectance = level_metrics["connectance"],
    links_per_species = level_metrics["links per species"],
    weighted_nestedness = level_metrics["weighted nestedness"],
    modularity = modularity_q,
    # Note: networklevel returns "generality.HL" and "vulnerability.LL" when called this way
    generality = level_metrics["generality.HL"],     # Mean generality for Higher Level
    vulnerability = level_metrics["vulnerability.LL"] # Mean vulnerability for Lower Level
  ))
}


# ----------------------------------------------------------------------------
# PART 7: MAIN ANALYSIS LOOP
# ----------------------------------------------------------------------------

n_reps <- 10 # Number of replicates (set to 100)

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
    net_links <- reconstruct_network(sp_traits, model_params, focal_habitat = "Ter-Mar")
    
    # Pass the 'net_links' dataframe DIRECTLY to calculate_metrics
    metrics_tibble <- calculate_metrics(net_links)
    
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

final_metrics_all_reps <- map_dfr(
  1:n_reps,
  ~run_analysis(splocs_agg, traits), 
  .id = "rep"
)

end_time <- Sys.time()
print(paste("Analysis complete. Time taken:", round(difftime(end_time, start_time, units = "secs"), 1), "seconds"))

# ----------------------------------------------------------------------------
# PART 8: SUMMARIZE METRICS
# ----------------------------------------------------------------------------

print("Summarizing metrics...")

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
# PART 9: PLOT METRICS (BOXPLOTS)
# ----------------------------------------------------------------------------

print("Generating boxplots...")

# --- 1. Define the correct order for time bins ---
time_order <- c("Pleistoceno final", "Holoceno medio", "Holoceno final", "Histórico")

# --- 2. Prepare the data for plotting ---
plot_data_metrics <- final_metrics_all_reps %>%
  pivot_longer(
    cols = c(connectance, links_per_species, weighted_nestedness, modularity, generality, vulnerability),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    time_bin = factor(time_bin, levels = time_order),
    metric = factor(metric, levels = c("connectance", "links_per_species", "weighted_nestedness", "modularity", "generality", "vulnerability"))
  ) %>%
  filter(!is.na(value))

# --- 3. Create the multi-faceted plot ---
combined_plot <- ggplot(plot_data_metrics, aes(x = time_bin, y = value)) +
  geom_boxplot(aes(fill = region_col), 
               position = position_dodge(width = 0.9), 
               width = 0.8, 
               outlier.shape = NA) + 
  facet_wrap(~ metric, scales = "free_y", ncol = 2, strip.position = "left") +
  labs(
    title = "Tierra del Fuego Aggregated Network Metrics",
    subtitle = "Boxplots show distribution of 100 stochastic replicates. (Prob > 0)",
    x = "Geological Age",
    y = "Metric Value",
    fill = "Region"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Bosque" = "#1B9E77", "Estepa" = "#D95F02")) + 
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    strip.placement = "outside",
    strip.text.y = element_text(angle = 0, face = "bold", size = 10),
    strip.text.x = element_blank(), 
    panel.border = element_rect(color = "grey80", fill = NA),
    panel.spacing = unit(1, "lines")
  )

ggsave(
  file.path(output_dir, paste0(today_date, "_TDF_metrics_boxplots_AGGREGATED.png")),
  combined_plot,
  width = 10, height = 10, dpi = 300, bg = "white"
)


# ----------------------------------------------------------------------------
# PART 10: GENERATE INDIVIDUAL NETWORK PLOTS
# ----------------------------------------------------------------------------

print("Generating individual network plots...")

shorten_name <- function(species_name) {
  parts <- str_split(species_name, "_")[[1]]
  if (length(parts) == 2) {
    return(paste0(str_sub(parts[1], 1, 1), ". ", parts[2]))
  } else {
    return(species_name)
  }
}

taxa_colors <- c(
  "Homo_sapiens" = "#33CCFF", # Bright blue for humans
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

mean_params <- list(
  alphaP = mean(c(1, 2)),
  betaP = mean(c(-2, -1)),
  gammaP = mean(c(-1, -0.5)),
  alphaO = mean(c(-6, -4)),
  betaO = mean(c(-3, -2)),
  gammaO = 0
)

plot_network_bipartite <- function(graph, title) {
  
  layout <- create_layout(graph, layout = 'bipartite')
  layout$y <- ifelse(layout$type, 1, 0)
  
  # Base plot with non-human links first
  p <- ggraph(layout) +
    # 1. Non-Human interactions (Blue-Red gradient)
    geom_edge_link(aes(alpha = Int_prob, color = Int_prob, 
                       filter = human_link_type == "Non-Human"), 
                   width = 0.5) +
    scale_edge_color_gradient(low = "#3B4CC0", high = "#B40426", name = "Interaction Prob.", limits = c(0,1)) +
    
    # 2. Human links (Manual colors via a separate scale trick or direct color assignment)
    # Since we can't easily have two color scales without ggnewscale (which failed),
    # we will use a workaround: Plot human edges as separate layers with FIXED colors.
    
    # Human Predation (Red)
    geom_edge_link(aes(alpha = Int_prob,
                       filter = human_link_type == "Human-Predation"), 
                   color = "red", width = 1.0) +
    
    # Human Scavenging (Brown)
    geom_edge_link(aes(alpha = Int_prob,
                       filter = human_link_type == "Human-Scavenging"), 
                   color = "#8B4513", width = 1.0) +
    
    # Nodes
    geom_node_point(aes(fill = Category, filter = !is_human_node), shape = 21, color = "black", size = 5, stroke = 0.5) +
    geom_node_point(aes(fill = Category, filter = is_human_node), shape = 21, color = "black", size = 7, stroke = 1.5) +
    
    geom_node_text(aes(label = short_name), repel = TRUE, size = 2.5, max.overlaps = 35, 
                   bg.colour = "white", bg.r = 0.1, segment.color = 'grey50') +
    
    scale_fill_manual(values = taxa_colors, name = "Taxonomic Category", na.value = "grey50", limits = names(taxa_colors), drop = FALSE) +
    scale_edge_alpha(guide = 'none') +
    scale_y_continuous(expand = expansion(mult = 0.4)) + 
    theme_graph(base_family = 'sans') +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      plot.title = element_text(hjust = 0.5, size = 16)
    ) +
    labs(title = title)
  
  return(p)
}


# --- Prepare data and loop ---
plot_data_prep <- splocs_agg %>%
  mutate(
    sp_traits = map(species_list, ~traits %>% filter(trophic_species %in% .x)),
    
    # Generate links (Representative)
    links_df = map(sp_traits, function(st) {
      links <- reconstruct_network(st, mean_params) %>%
        filter(Int_prob > 0)
      
      links_with_traits <- links %>%
        left_join(st %>% select(trophic_species, carrion_r = carrion_only), 
                  by = c("Resource" = "trophic_species"))
      
      links_with_traits %>%
        mutate(
          from = paste0(Consumer, "_C"),
          to = paste0(Resource, "_R"),
          human_link_type = case_when(
            Consumer == "Homo_sapiens" & carrion_r == 1 ~ "Human-Scavenging",
            Consumer == "Homo_sapiens" ~ "Human-Predation",
            TRUE ~ "Non-Human"
          )
        )
    }),
    
    # Generate Nodes
    nodes_df = map2(sp_traits, links_df, function(st, links) {
      consumer_names <- unique(links$Consumer)
      resource_names <- unique(links$Resource) 
      
      nodes_C <- st %>%
        filter(trophic_species %in% consumer_names) %>%
        mutate(name = paste0(trophic_species, "_C"), type = TRUE)
      
      nodes_R <- st %>%
        filter(trophic_species %in% resource_names) %>%
        mutate(name = paste0(trophic_species, "_R"), type = FALSE)
      
      bind_rows(nodes_C, nodes_R) %>%
        mutate(
          is_human_node = (trophic_species == "Homo_sapiens"),
          Category = case_when(
            is_human_node ~ "Homo_sapiens",
            class %in% names(taxa_colors) ~ class,
            TRUE ~ "Other"
          ),
          short_name = map_chr(trophic_species, shorten_name)
        ) %>%
        select(name, short_name, type, Category, is_human_node, everything()) %>%
        distinct(name, .keep_all = TRUE)
    }),
    
    plot_title = paste0("Aggregated Network: ", time_bin, " - ", region_col),
    file_name = file.path(output_dir, paste0(today_date, "_NetworkPlot_", time_bin, "_", region_col, ".png"))
  )

pwalk(plot_data_prep, function(nodes_df, links_df, plot_title, file_name, ...) {
  if (nrow(nodes_df) == 0 || nrow(links_df) == 0) return()
  
  all_nodes_in_links <- unique(c(links_df$from, links_df$to))
  nodes_df <- nodes_df %>% filter(name %in% all_nodes_in_links)
  
  if (nrow(nodes_df) == 0) return()
  
  g <- tbl_graph(nodes = nodes_df, edges = links_df, directed = TRUE)
  final_plot <- plot_network_bipartite(g, title = plot_title)
  
  ggsave(filename = file_name, plot = final_plot, width = 14, height = 10, dpi = 300, bg = "white")
})


# ----------------------------------------------------------------------------
# PART 11: EXPORT MASTER INTERACTION LIST
# ----------------------------------------------------------------------------

all_interactions_list <- plot_data_prep %>%
  select(time_bin, region_col, links_df) %>%
  unnest(links_df) %>%
  select(
    food_web_time = time_bin,
    food_web_region = region_col,
    Consumer,
    Resource,
    Int_prob,
    human_link_type
  ) %>%
  arrange(food_web_time, food_web_region, desc(Int_prob))

write_csv(
  all_interactions_list,
  file.path(output_dir, paste0(today_date, "_TDF_all_interactions_list.csv"))
)


# ----------------------------------------------------------------------------
# PART 12: NEW BIPLOT (INTERACTION PROBABILITY BAR PLOT)
# *** NEW SECTION ***
# ----------------------------------------------------------------------------

print("Generating Interaction Probability Bar Plots...")

# Prepare data for the bar plot
# We sum the interaction probabilities for each species in each network
int_prob_summary <- all_interactions_list %>%
  # Gather all species (both consumers and resources involved in interactions)
  pivot_longer(
    cols = c(Consumer, Resource),
    names_to = "Role",
    values_to = "Species"
  ) %>%
  # Group by Network and Species
  group_by(food_web_time, food_web_region, Species) %>%
  # Calculate total interaction probability
  summarise(Total_Int_Prob = sum(Int_prob, na.rm = TRUE), .groups = "drop") %>%
  # Create a unique network identifier for faceting
  mutate(
    Network = paste0(food_web_time, "\n", food_web_region),
    # Simplify species names for plotting
    Short_Name = map_chr(Species, shorten_name),
    # Identify Homo sapiens for coloring
    Is_Human = Species == "Homo_sapiens"
  ) %>%
  # *** SORTING FIX: Order species by Total_Int_Prob for better visualization ***
  group_by(Network) %>% 
  mutate(Short_Name = reorder_within(Short_Name, Total_Int_Prob, Network)) %>%
  ungroup()

# Create the plot
p_biplot <- ggplot(int_prob_summary, aes(x = Short_Name, y = Total_Int_Prob, fill = Is_Human)) +
  geom_bar(stat = "identity") +
  # Facet by Network
  facet_wrap(~ Network, scales = "free") +
  scale_x_reordered() +
  coord_flip() + # Horizontal bars are better for species names
  # Use specific colors: Blue for humans, Steelblue for others
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "#33CCFF")) +
  labs(
    title = "Total Interaction Probability per Species",
    x = "Species",
    y = "Total Probability Sum"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    strip.text = element_text(face = "bold"),
    legend.position = "none" # Hide the legend as it's self-explanatory
  )

# Save the plot
ggsave(
  file.path(output_dir, paste0(today_date, "_TDF_interaction_prob_barplot.png")),
  p_biplot,
  width = 15,
  height = 12,
  dpi = 300,
  bg = "white"
)

print(paste("All plots and data saved to", output_dir, "directory."))
print("--- SCRIPT FINISHED (ALL PARTS) ---")
