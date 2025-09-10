# -----------------------------------------------------------------------------
#
# R SCRIPT FOR RECONSTRUCTING ANCIENT CONSUMER-RESOURCE NETWORKS
# Definitive Final Version: This version fixes the legend display by creating
# a single, unified legend per page and adding Homo sapiens as a distinct
# category in the node legend.
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
# The 'tidygraph' package provides tools to work with graph objects cleanly.
library(tidygraph)
# The 'ggraph' package is an extension of ggplot2 for plotting networks.
library(ggraph)
# The 'ggrepel' package is used for non-overlapping text labels.
library(ggrepel)
# The 'patchwork' package is used to combine separate ggplots into one figure.
library(patchwork)
# The 'viridis' package is required to generate the viridis color palette.
library(viridis)


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
localities_df_base <- localities_raw %>%
  mutate(across(-TrophicSpecies, as.numeric)) %>%
  pivot_longer(
    cols = -TrophicSpecies,
    names_to = "Locality",
    values_to = "Abundance"
  ) %>%
  mutate(Abundance = replace_na(Abundance, 0)) %>%
  filter(Abundance > 0) %>%
  select(Locality, TrophicSpecies)

# --- Inject Homo sapiens as a Ubiquitous Consumer ---
all_localities <- unique(localities_df_base$Locality)
homo_sapiens_presence <- tibble(
  Locality = all_localities,
  TrophicSpecies = "Homo_sapiens"
)
localities_df <- bind_rows(localities_df_base, homo_sapiens_presence) %>%
  distinct()


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
      FeedingStrategy %in% c("Carnivore", "Piscivore"), # "Omnivore", "Insectivore"
      "Consumer",
      "Resource"
    )
  ) %>%
  select(TrophicSpecies, BodyMass_g, Guild, FeedingHabitat) %>%
  filter(!is.na(BodyMass_g) & !is.na(FeedingHabitat))

# --- Tidy the Age and Spatial Data ---
ages_df <- ages_raw %>%
  select(Locality, Geological_ages, Biome = Region_col)


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

N_REPLICATES <- 1000

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
  filter(Prey != "Homo_sapiens") %>%
  filter(Int_prob > 0) %>%
  select(Locality, Biome, Geological_ages, Predator, Prey, Int_prob) %>%
  arrange(Biome, Geological_ages, Locality, Predator, desc(Int_prob))


# -----------------------------------------------------------------------------
# PART 7: VISUALIZATION OF NETWORK EVOLUTION BY BIOME
# Creates a final, polished plot for each biome.
# -----------------------------------------------------------------------------

# --- Step 1: Prepare data for plotting ---
age_levels <- c("Pleistoceno final", "Holoceno medio", "Holoceno final", "Histórico")
plotting_df <- final_interaction_df %>%
  filter(!is.na(Biome) & !is.na(Geological_ages)) %>%
  group_by(Biome, Geological_ages, Predator, Prey) %>%
  summarise(Int_prob = mean(Int_prob, na.rm = TRUE), .groups = 'drop') %>%
  mutate(Geological_ages = factor(Geological_ages, levels = age_levels))

# --- Step 2: Open PDF device and loop through biomes to create composite plots ---
pdf("results/Network_Evolution_by_Biome_Strict.pdf", width = 15, height = 8.5)
biomes_to_plot <- unique(plotting_df$Biome)

for (current_biome in biomes_to_plot) {
  
  biome_data <- plotting_df %>%
    filter(Biome == current_biome)
  
  plot_list <- age_levels %>%
    map(~ {
      age_data <- biome_data %>% filter(Geological_ages == .x)
      
      if (nrow(age_data) == 0) {
        return(ggplot() + theme_void() + ggtitle(.x))
      }
      
      edges_for_graph <- age_data %>%
        mutate(
          from = paste0(Predator, "_C"), 
          to = paste0(Prey, "_R"),
          is_human_interaction = (Predator == "Homo_sapiens" | Prey == "Homo_sapiens")
        )
      
      # MODIFIED: Create a new 'NodeType' column for the legend
      vertices <- tibble(name = unique(c(edges_for_graph$from, edges_for_graph$to))) %>%
        mutate(
          short_name = str_remove(name, "_[CR]$"),
          Guild = if_else(str_ends(name, "_C$"), "Consumer", "Resource"),
          type = Guild == "Consumer",
          is_human_node = (short_name == "Homo_sapiens"),
          NodeType = case_when(
            is_human_node ~ "Homo sapiens",
            TRUE ~ Guild
          )
        )
      
      clean_edges <- edges_for_graph %>% select(from, to)
      graph <- graph_from_data_frame(d = clean_edges, vertices = vertices, directed = TRUE)
      
      E(graph)$Int_prob <- edges_for_graph$Int_prob
      E(graph)$is_human_interaction <- edges_for_graph$is_human_interaction
      
      ggraph(graph, layout = 'bipartite') +
        geom_edge_fan(aes(alpha = Int_prob, color = Int_prob, 
                          filter = !is_human_interaction), width = 0.5) +
        geom_edge_fan(aes(alpha = Int_prob, filter = is_human_interaction), 
                      color = "red", width = 1.2) +
        # MODIFIED: Map color and fill aesthetics to the new 'NodeType' column
        geom_node_point(aes(color = NodeType, filter = !is_human_node), size = 4) +
        geom_node_point(aes(fill = NodeType, filter = is_human_node), 
                        color = "black", size = 6, shape = 21, stroke = 1.2) +
        geom_node_text(aes(label = short_name), repel = TRUE, size = 2.5, max.overlaps = 15, bg.colour = "white", segment.color = 'grey50') +
        # MODIFIED: Update the manual scale with three categories
        scale_color_manual(
          values = c("Consumer" = "tomato", "Resource" = "skyblue", "Homo sapiens" = "tomato"), 
          name = "Node Type", 
          aesthetics = c("color", "fill")
        ) +
        # MODIFIED: Add limits to the edge scale to ensure it's collectible
        scale_edge_color_gradientn(
          colors = viridis::viridis(256), 
          name = "Interaction Prob.",
          limits = c(0, 1) 
        ) +
        scale_edge_alpha(guide = 'none') +
        theme_graph(base_family = 'sans') +
        ggtitle(.x)
    })
  
  composite_plot <- wrap_plots(plot_list, ncol = 4) +
    plot_layout(guides = 'collect') +
    plot_annotation(
      title = paste("Network Evolution in Biome:", current_biome),
      subtitle = "Networks are shown chronologically by geological age.",
      caption = "Edge color indicates interaction probability (yellow=high). Red lines and bordered node highlight Homo sapiens."
    )
  
  print(composite_plot)
}

# --- Step 3: Close the PDF device ---
dev.off()

cat("\nPDF file 'Network_Evolution_by_Biome_Final.pdf' has been created in your working directory.\n")

