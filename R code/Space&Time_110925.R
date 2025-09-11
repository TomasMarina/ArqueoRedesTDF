# -----------------------------------------------------------------------------
#
# R SCRIPT FOR RECONSTRUCTING ANCIENT CONSUMER-RESOURCE NETWORKS
# Definitive Final Version: This version implements the final, definitive fix
# for the network metrics calculation by correctly transposing the incidence
# matrix and updates the interaction threshold to 0.1 as requested.
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
library(bipartite)


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
# PART 5: REVISED SIMULATION WORKFLOW
# Generate and store all replicates for each site.
# -----------------------------------------------------------------------------

# For quick testing, you can reduce this number. For final results, use 1000.
N_REPLICATES <- 100 

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
  # MODIFIED: Use the new, stricter threshold of 0.1
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
    
    if (nrow(incidence_matrix) < 1 || ncol(incidence_matrix) < 1) {
      return(NULL)
    }
    
    # DEFINITIVE FIX: Transpose the matrix to have HIGHER level (consumers) as COLUMNS
    transposed_matrix <- t(incidence_matrix)
    
    # Calculate metrics on the correctly oriented binary matrix
    level_metrics <- networklevel(transposed_matrix, index = c("connectance", "generality", "vulnerability"))
    
    modules <- tryCatch(computeModules(transposed_matrix), error = function(e) NULL)
    modularity_q <- if (!is.null(modules)) modules@likelihood else NA
    
    # Combine all metrics into a single row, extracting the correct level-specific values
    tibble(
      Replicate = group_key$Replicate,
      Biome = group_key$Biome,
      Geological_ages = group_key$Geological_ages,
      Connectance = level_metrics["connectance"],
      # Generality is a property of the HIGHER level (now the columns)
      Generality = level_metrics["generality.HL"],
      # Vulnerability is a property of the LOWER level (now the rows)
      Vulnerability = level_metrics["vulnerability.LL"],
      Modularity = modularity_q
    )
  }) %>%
  list_rbind() %>%
  mutate(Geological_ages = factor(Geological_ages, levels = age_levels))

# --- Display and save a summary of the metrics ---
cat("\n--- Summary of Replicated Network Metrics ---\n")

summary_metrics <- replicated_metrics_df %>%
  group_by(Biome, Geological_ages) %>%
  summarise(
    across(
      Connectance:Modularity, 
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd = ~sd(.x, na.rm = TRUE)
      )
    ),
    .groups = "drop"
  )

print(summary_metrics)
write_csv(summary_metrics, "network_metrics_summary.csv")
cat("\nSummary of network metrics saved to 'network_metrics_summary.csv'\n\n")


# -----------------------------------------------------------------------------
# PART 8: VISUALIZE METRIC DISTRIBUTIONS
# Create histograms for each metric, faceted by network.
# -----------------------------------------------------------------------------

cat("Creating PDF of metric distributions...\n")
pdf("Network_Metrics_Distributions.pdf", width = 12, height = 8)

plotting_metrics <- replicated_metrics_df %>%
  pivot_longer(
    cols = Connectance:Modularity,
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))

metrics_plot <- ggplot(plotting_metrics, aes(x = Value)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  facet_grid(
    Metric ~ Biome + Geological_ages,
    scales = "free"
  ) +
  labs(
    title = "Distributions of Network Metrics Across Replicates",
    subtitle = "Each panel represents a unique food web, showing the distribution of the metric.",
    x = "Metric Value",
    y = "Frequency"
  ) +
  theme_bw() +
  theme(
    strip.text.x = element_text(size = 8),
    strip.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(metrics_plot)
dev.off()
cat("PDF file 'Network_Metrics_Distributions.pdf' has been created.\n")

