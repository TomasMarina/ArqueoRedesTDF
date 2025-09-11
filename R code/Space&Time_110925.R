# -----------------------------------------------------------------------------
#
# R SCRIPT FOR RECONSTRUCTING ANCIENT CONSUMER-RESOURCE NETWORKS
# Definitive Final Version: This script performs the full replicated analysis,
# calculates network metrics, produces a summary table and boxplots, AND
# generates the individual, polished network plots for each food web.
#
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# PART 1: SETUP
# Load all necessary packages for the analysis.
# -----------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(bipartite)
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
localities_raw <- read_excel("data/Species localities - TDF - 030925.xlsx")
traits_raw <- read_excel("data/Species traits - TDF - 100925.xlsx", sheet = "traits")
ages_raw <- read_excel("data/Age localities - TDF - 030925.xlsx")


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


# --- Tidy the Traits Data with Taxonomic Categories ---
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
# PART 4: THE LOG-RATIO MODEL (LRM) FUNCTION WITH PARAMETER SAMPLING
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
    alphaP = runif(1, 1, 2),
    betaP = runif(1, -2, -1),
    gammaP = runif(1, -1, -0.5),
    alphaO = runif(1, -6, -4),
    betaO = runif(1, -3, -2),
    gammaO = 0
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
# PART 5: REVISED SIMULATION WORKFLOW
# Generate and store all unique replicates for each site.
# -----------------------------------------------------------------------------

N_REPLICATES <- 1000

replicated_interaction_list <- localities_df %>%
  inner_join(traits_df, by = "TrophicSpecies") %>%
  group_by(Locality) %>%
  summarise(
    replicates = list(
      replicate(N_REPLICATES, {
        consumers <- cur_data() %>% filter(Guild == "Consumer")
        resources <- cur_data()
        if (nrow(consumers) == 0 || nrow(resources) == 0) return(NULL)
        calculate_lrm_matrix(consumers, resources)
      }, simplify = FALSE)
    ),
    .groups = "drop"
  ) %>%
  unnest_longer(replicates, indices_to = "Replicate") %>%
  filter(!is.null(replicates)) %>%
  mutate(
    tidy_replicate = map(replicates, ~ .x %>%
                           as.data.frame() %>%
                           rownames_to_column(var = "Predator") %>%
                           pivot_longer(cols = -Predator, names_to = "Prey", values_to = "Int_prob")
    )
  ) %>%
  select(Replicate, Locality, tidy_replicate) %>%
  unnest(tidy_replicate)


# -----------------------------------------------------------------------------
# PART 6: CREATE FINAL REPLICATED INTERACTION DATA FRAME
# Apply the standard filters to the full replicated dataset.
# -----------------------------------------------------------------------------

full_replicated_df <- replicated_interaction_list %>%
  left_join(ages_df, by = "Locality") %>%
  filter(Predator != Prey) %>%
  filter(Prey != "Homo_sapiens") %>%
  filter(Int_prob >= 0.1) %>%
  select(Replicate, Biome, Geological_ages, Predator, Prey)


# -----------------------------------------------------------------------------
# PART 7: REPLICATED NETWORK METRICS CALCULATION
# Calculate metrics for each replicate of each unique network.
# -----------------------------------------------------------------------------

age_levels <- c("Pleistoceno final", "Holoceno medio", "Holoceno final", "Histórico")

replicated_metrics_df <- full_replicated_df %>%
  filter(!is.na(Biome) & !is.na(Geological_ages)) %>%
  group_by(Replicate, Biome, Geological_ages) %>%
  group_map(~ {
    
    network_data <- .x
    group_key <- .y
    
    incidence_matrix <- network_data %>%
      distinct(Predator, Prey) %>%
      mutate(value = 1) %>%
      pivot_wider(
        id_cols = Predator,
        names_from = Prey,
        values_from = value,
        values_fill = 0
      ) %>%
      column_to_rownames("Predator") %>%
      as.matrix()
    
    if (nrow(incidence_matrix) < 2 || ncol(incidence_matrix) < 2) {
      return(NULL)
    }
    
    transposed_matrix <- t(incidence_matrix)
    
    num_consumers <- ncol(transposed_matrix)
    num_resources <- nrow(transposed_matrix)
    total_species <- num_consumers + num_resources
    num_links <- sum(transposed_matrix)
    
    level_metrics <- networklevel(transposed_matrix, index = c("connectance", "generality", "vulnerability"))
    
    modules <- tryCatch(computeModules(transposed_matrix), error = function(e) NULL)
    modularity_q <- if (!is.null(modules)) modules@likelihood else NA
    
    tibble(
      Replicate = group_key$Replicate,
      Biome = group_key$Biome,
      Geological_ages = group_key$Geological_ages,
      Num_Consumers = num_consumers,
      Num_Resources = num_resources,
      Total_Species = total_species,
      Num_Links = num_links,
      Connectance = level_metrics["connectance"],
      Generality = level_metrics["generality.HL"],
      Vulnerability = level_metrics["vulnerability.LL"],
      Modularity = modularity_q
    )
  }) %>%
  list_rbind() %>%
  mutate(Geological_ages = factor(Geological_ages, levels = age_levels))

# --- Display and save a summary of the metrics with final formatting ---
cat("\n--- Summary of Replicated Network Metrics ---\n")

structural_summary <- replicated_metrics_df %>%
  group_by(Biome, Geological_ages) %>%
  summarise(
    across(
      c(Num_Consumers, Num_Resources, Total_Species, Num_Links),
      list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Consumers = sprintf("%.3f (%.3f)", Num_Consumers_mean, Num_Consumers_sd),
    Resources = sprintf("%.3f (%.3f)", Num_Resources_mean, Num_Resources_sd),
    `Total Species` = sprintf("%.3f (%.3f)", Total_Species_mean, Total_Species_sd),
    Links = sprintf("%.3f (%.3f)", Num_Links_mean, Num_Links_sd)
  ) %>%
  select(Biome, Geological_ages, Consumers, Resources, `Total Species`, Links)

metrics_summary <- replicated_metrics_df %>%
  group_by(Biome, Geological_ages) %>%
  summarise(
    across(
      Connectance:Modularity, 
      list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE))
    ),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

summary_metrics <- left_join(structural_summary, metrics_summary, by = c("Biome", "Geological_ages"))

print(summary_metrics)
write_csv(summary_metrics, "network_metrics_summary.csv")
cat("\nSummary of network metrics saved to 'network_metrics_summary.csv'\n\n")


# -----------------------------------------------------------------------------
# PART 8: VISUALIZATION SETUP
# Define common elements for all plots.
# -----------------------------------------------------------------------------

taxa_colors <- c(
  "Aves" = rgb(221, 165, 2, maxColorValue = 255),
  "Mamiferos" = rgb(205, 58, 46, maxColorValue = 255),
  "Peces" = rgb(47, 99, 184, maxColorValue = 255),
  "Invertebrados" = rgb(192, 80, 1, maxColorValue = 255),
  "Humano" = rgb(0, 158, 115, maxColorValue = 255)
)


# -----------------------------------------------------------------------------
# PART 9: VISUALIZE METRIC DISTRIBUTIONS AS BOXPLOTS
# -----------------------------------------------------------------------------

cat("Creating PDF of metric boxplots...\n")
pdf("Network_Metrics_Boxplots.pdf", width = 12, height = 10)

plotting_metrics <- replicated_metrics_df %>%
  pivot_longer(
    cols = Connectance:Modularity,
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value)) %>%
  mutate(Metric = factor(Metric, levels = c("Connectance", "Generality", "Vulnerability", "Modularity")))

biome_colors <- c(
  "Estepa" = "#D55E00",
  "Bosque" = "#0072B2"
)

metrics_boxplot <- ggplot(plotting_metrics, aes(x = Geological_ages, y = Value, fill = Biome)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, outlier.shape = NA) +
  facet_wrap(
    ~ Metric,
    scales = "free_y",
    ncol = 2
  ) +
  scale_fill_manual(values = biome_colors) +
  labs(
    title = "Comparison of Network Metrics Across Time and Biome",
    subtitle = paste("Boxplots show the distribution of each metric from", N_REPLICATES, "replicates."),
    x = "Geological Age",
    y = "Metric Value",
    fill = "Biome"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top"
  )

print(metrics_boxplot)
dev.off()
cat("PDF file 'Network_Metrics_Boxplots.pdf' has been created.\n")


# -----------------------------------------------------------------------------
# PART 10: VISUALIZATION: SAVE EACH NETWORK TO A SEPARATE PDF
# -----------------------------------------------------------------------------

cat("\nCreating individual network plots...\n")

traits_for_join_viz <- traits_df %>% select(TrophicSpecies, Category)

plotting_df_viz <- replicated_interaction_list %>%
  filter(Int_prob >= 0.1) %>%
  left_join(ages_df, by = "Locality") %>%
  filter(Predator != Prey) %>%
  filter(Prey != "Homo_sapiens") %>%
  left_join(traits_for_join_viz, by = c("Predator" = "TrophicSpecies")) %>%
  rename(Predator_cat = Category) %>%
  left_join(traits_for_join_viz, by = c("Prey" = "TrophicSpecies")) %>%
  rename(Prey_cat = Category) %>%
  mutate(Predator_cat = if_else(Predator == "Homo_sapiens", "Humano", Predator_cat)) %>%
  filter(!is.na(Biome) & !is.na(Geological_ages) & !is.na(Predator_cat) & !is.na(Prey_cat)) %>%
  group_by(Biome, Geological_ages, Predator, Predator_cat, Prey, Prey_cat) %>%
  summarise(Int_prob = mean(Int_prob, na.rm = TRUE), .groups = 'drop') %>%
  mutate(Geological_ages = factor(Geological_ages, levels = age_levels))

if (!dir.exists("individual_networks")) {
  dir.create("individual_networks")
}

plotting_df_viz %>%
  group_by(Biome, Geological_ages) %>%
  group_walk(~ {
    age_data <- .x
    group_key <- .y
    
    filename <- paste0("individual_networks/Network_", group_key$Biome, "_", str_replace_all(group_key$Geological_ages, " ", "_"), ".pdf")
    cat("Creating PDF:", filename, "\n")
    
    pdf(filename, width = 11, height = 8.5)
    
    edges_for_graph <- age_data %>%
      mutate(from = paste0(Predator, "_C"), to = paste0(Prey, "_R"), is_human_interaction = (Predator == "Homo_sapiens"))
    
    vertices_for_graph <- tibble(
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
    
    graph <- graph_from_data_frame(d = edges_for_graph %>% select(from, to, Int_prob, is_human_interaction), vertices = vertices_for_graph, directed = TRUE)
    
    tbl_graph <- as_tbl_graph(graph)
    
    layout <- create_layout(tbl_graph, layout = 'fr')
    layout$y <- ifelse(layout$Guild == 'Consumer', 2, 0)
    layout$y <- layout$y + runif(nrow(layout), min = -0.3, max = 0.3)
    
    final_plot <- ggraph(layout) +
      geom_edge_fan(aes(alpha = Int_prob, color = Int_prob, filter = !is_human_interaction), width = 0.5) +
      geom_edge_fan(aes(alpha = Int_prob, filter = is_human_interaction), color = "red", width = 1.2) +
      geom_node_point(aes(fill = Category, filter = !is_human_node), shape = 21, color = "black", size = 5, stroke = 0.5) +
      geom_node_point(aes(fill = Category, filter = is_human_node), shape = 21, color = "black", size = 7, stroke = 1.5) +
      geom_node_text(aes(label = short_name), repel = TRUE, size = 2.5, max.overlaps = 35, bg.colour = "white", segment.color = 'grey50') +
      scale_fill_manual(values = taxa_colors, name = "Taxonomic Category", na.value = "grey50", limits = names(taxa_colors), drop = FALSE) +
      scale_edge_color_gradientn(colors = viridis::magma(256, direction = -1), name = "Interaction Prob.", limits = c(0,1)) +
      scale_edge_alpha(guide = 'none') +
      scale_y_continuous(expand = expansion(mult = 0.4)) + 
      theme_graph(base_family = 'sans') +
      theme(legend.position = "bottom", plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
      labs(
        title = paste("Food Web for Biome:", group_key$Biome),
        subtitle = paste("Geological Age:", group_key$Geological_ages),
        caption = "Node color indicates taxonomic category. Edge color indicates probability (yellow=high). Red lines & bordered node highlight Homo sapiens."
      )
    
    print(final_plot)
    dev.off()
  })

cat("\nAll individual network plots have been created.\n")

