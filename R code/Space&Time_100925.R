# -----------------------------------------------------------------------------
#
# R SCRIPT FOR RECONSTRUCTING ANCIENT CONSUMER-RESOURCE NETWORKS
# Definitive Final Version: This version implements the final, definitive fix
# for the plotting error by building a complete graph object upfront and using
# ggraph's 'filter' aesthetic, which is the robust and correct method.
#
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# PART 1: SETUP
# Load all necessary packages for the analysis.
# -----------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(igraph)
library(tidygraph)
library(ggraph)
library(ggrepel)
library(patchwork)
library(viridis)


# -----------------------------------------------------------------------------
# PART 2: DATA LOADING AND INITIAL CLEANING
# Read the updated raw data files.
# -----------------------------------------------------------------------------

# ASSUMPTION: The script assumes your .xlsx files are in a 'data' sub-folder.
# To make this script run, please create a folder named 'data' in the same
# directory as the script, and place your Excel files inside it.
localities_raw <- read_excel("data - Santiago2025/Species localities - TDF - 030925.xlsx")
traits_raw <- read_excel("data - Santiago2025/Species traits - TDF - 100925.xlsx", sheet = "traits")
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
  select(Class, TrophicSpecies, BodySize_min, BodySize_max, FeedingStrategy, FeedingHabitat) %>%
  filter(!is.na(FeedingStrategy), FeedingStrategy != "??", !is.na(Class)) %>%
  mutate(
    Category = case_when(
      Class == "Aves" ~ "Aves",
      Class == "Mammalia" ~ "Mamiferos",
      Class %in% c("Teleostei", "Actinopterygii", "Chondrichthyes") ~ "Peces",
      TRUE ~ "Invertebrados"
    )
  ) %>%
  rowwise() %>%
  mutate(
    BodyMass_g = exp(mean(log(c(BodySize_min, BodySize_max)), na.rm = TRUE)) * 1000
  ) %>%
  ungroup() %>%
  mutate(
    Guild = if_else(
      FeedingStrategy %in% c("Carnivore", "Piscivore", "Omnivore", "Insectivore"),
      "Consumer",
      "Resource"
    )
  ) %>%
  select(Category, TrophicSpecies, BodyMass_g, Guild, FeedingHabitat) %>%
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
# PART 6: CREATE FINAL INTERACTION DATA FRAME WITH TAXONOMIC CATEGORIES
# Filter, join data, and add the new Predator_cat and Prey_cat columns.
# -----------------------------------------------------------------------------

traits_for_join <- traits_df %>% select(TrophicSpecies, Category)

final_interaction_df <- interaction_list %>%
  left_join(ages_df, by = "Locality") %>%
  filter(Predator != Prey) %>%
  filter(Prey != "Homo_sapiens") %>%
  filter(Int_prob > 0) %>%
  left_join(traits_for_join, by = c("Predator" = "TrophicSpecies")) %>%
  rename(Predator_cat = Category) %>%
  left_join(traits_for_join, by = c("Prey" = "TrophicSpecies")) %>%
  rename(Prey_cat = Category) %>%
  mutate(
    Predator_cat = if_else(Predator == "Homo_sapiens", "Humano", Predator_cat)
  ) %>%
  select(Locality, Biome, Geological_ages, Predator, Predator_cat, Prey, Prey_cat, Int_prob) %>%
  arrange(Biome, Geological_ages, Locality, Predator, desc(Int_prob))


# -----------------------------------------------------------------------------
# PART 7: VISUALIZATION OF NETWORK EVOLUTION BY BIOME
# Creates a final, polished plot using the new taxonomic categories.
# -----------------------------------------------------------------------------

# --- Step 1: Prepare data and define color palette ---
age_levels <- c("Pleistoceno final", "Holoceno medio", "Holoceno final", "Histórico")

taxa_colors <- c(
  "Aves" = rgb(221, 165, 2, maxColorValue = 255),
  "Mamiferos" = rgb(205, 58, 46, maxColorValue = 255),
  "Peces" = rgb(47, 99, 184, maxColorValue = 255),
  "Invertebrados" = rgb(192, 80, 1, maxColorValue = 255),
  "Humano" = rgb(0, 158, 115, maxColorValue = 255)
)

plotting_df <- final_interaction_df %>%
  filter(!is.na(Biome) & !is.na(Geological_ages)) %>%
  group_by(Biome, Geological_ages, Predator, Predator_cat, Prey, Prey_cat) %>%
  summarise(Int_prob = mean(Int_prob, na.rm = TRUE), .groups = 'drop') %>%
  mutate(Geological_ages = factor(Geological_ages, levels = age_levels))

# --- Step 2: Open PDF device and loop through biomes to create plots ---
pdf("Network_Evolution_by_Biome_Final_Categorized.pdf", width = 15, height = 8.5)
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
      
      # DEFINITIVE FIX: Build the graph with all attributes from the start
      
      edges_for_graph <- age_data %>%
        mutate(
          from = paste0(Predator, "_C"), 
          to = paste0(Prey, "_R"),
          is_human_interaction = (Predator == "Homo_sapiens")
        )
      
      vertices <- tibble(
        short_name = c(age_data$Predator, age_data$Prey),
        Category = c(age_data$Predator_cat, age_data$Prey_cat),
        Guild = c(rep("Consumer", nrow(age_data)), rep("Resource", nrow(age_data)))
      ) %>%
        mutate(
          name = if_else(Guild == "Consumer", paste0(short_name, "_C"), paste0(short_name, "_R")),
          type = Guild == "Consumer",
          is_human_node = (short_name == "Homo_sapiens")
        ) %>%
        distinct(name, .keep_all = TRUE) %>%
        select(name, everything())
      
      # Create the graph object with all attributes included in the edge data frame
      graph <- graph_from_data_frame(
        d = edges_for_graph %>% select(from, to, Int_prob, is_human_interaction), 
        vertices = vertices, 
        directed = TRUE
      )
      
      # Convert to tbl_graph for ggraph compatibility
      tbl_graph <- as_tbl_graph(graph)
      
      # Generate the plot for this specific age, using the 'filter' aesthetic
      ggraph(tbl_graph, layout = 'bipartite') +
        # Layer for standard edges
        geom_edge_fan(aes(alpha = Int_prob, color = Int_prob, filter = !is_human_interaction), 
                      width = 0.5) +
        # Layer for highlighted human edges
        geom_edge_fan(aes(alpha = Int_prob, filter = is_human_interaction), 
                      color = "red", width = 1.2) +
        # Layer for standard nodes
        geom_node_point(aes(fill = Category, filter = !is_human_node), 
                        shape = 21, color = "black", size = 5, stroke = 0.5) +
        # Layer for highlighted human node
        geom_node_point(aes(fill = Category, filter = is_human_node), 
                        shape = 21, color = "black", size = 7, stroke = 1.5) +
        geom_node_text(aes(label = short_name), repel = TRUE, size = 2.5, max.overlaps = 15, bg.colour = "white", segment.color = 'grey50') +
        scale_fill_manual(values = taxa_colors, name = "Taxonomic Category", na.value = "grey50") +
        scale_edge_color_gradient(low = "grey85", high = "black", name = "Interaction Prob.", limits = c(0,1)) +
        scale_edge_alpha(guide = 'none') +
        theme_graph(base_family = 'sans') +
        ggtitle(.x)
    })
  
  composite_plot <- wrap_plots(plot_list, ncol = 4) +
    plot_layout(guides = 'collect') +
    plot_annotation(
      title = paste("Network Evolution in Biome:", current_biome),
      subtitle = "Networks are shown chronologically by geological age.",
      caption = "Node color indicates taxonomic category. Edge color indicates probability (black=high). Red lines & bordered node highlight Homo sapiens."
    )
  
  print(composite_plot)
}

# --- Step 3: Close the PDF device ---
dev.off()

cat("\nPDF file 'Network_Evolution_by_Biome_Final_Categorized.pdf' has been created in your working directory.\n")
