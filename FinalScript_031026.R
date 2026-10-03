# ============================================================================
# TITLE: Reconstructing Ancient Ecosystems: Guanaco and Marine Carrion as 
#        Drivers of Predator-Prey Networks in Southern South America
# AUTHOR: Tomás I. Marina
# EMAIL: tomasimarina@gmail.com
# ============================================================================

# ----------------------------------------------------------------------------
# R SCRIPT: COMPLETE SPATIO-TEMPORAL ASSEMBLAGE 
# ----------------------------------------------------------------------------

# Load standard libraries required for data manipulation and visualization
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("readxl")) install.packages("readxl")

library(tidyverse)
library(readxl)

print("Starting structuring of the global matrix...")

# 1. DATA LOADING
# Read the functional traits, localities (NISP), and chronological ages datasets
traits_raw <- read_excel("../data/Species_traits_TDF_270426.xlsx", sheet = "traits")
localities_raw <- read_excel("../data/Species_localities_TDF_270526.xlsx")
ages_raw <- read_excel("../data/Age_localities_TDF_080926.xlsx")

# 2. FILTERING AND CLEANING OF FUNCTIONAL TRAITS
# Standardize taxonomic resolution, format categorical variables, and keep unique species
traits_clean <- traits_raw %>%
  mutate(
    # Assign taxonomic resolution based on the presence of genus or specific epithet markers
    Resolution = case_when(
      is.na(Genus) | Genus == "" ~ "Indeterminate/Broad",
      str_detect(TrophicSpecies, "_| ") ~ "Species",
      TRUE ~ "Genus"
    ),
    # Fill missing values with "no" for specific binary ecological traits
    Carrion = ifelse(is.na(Carrion), "no", Carrion),
    HumanResource = ifelse(is.na(HumanResource), "no", HumanResource)
  ) %>%
  select(
    TrophicSpecies, Class, Resolution, Bodymass_min = BodySize_min,
    Bodymass_max = BodySize_max, FeedingStrategy, FeedingHabitat,
    Carrion, HumanResource
  ) %>%
  distinct(TrophicSpecies, .keep_all = TRUE)

# 3. PROCESSING OF LOCALITIES AND AGES
# Transform the localities matrix into a long format and filter out zero-counts
localities_long <- localities_raw %>%
  pivot_longer(
    cols = -TrophicSpecies, 
    names_to = "Locality", 
    values_to = "NISP"
  ) %>%
  mutate(NISP = as.numeric(NISP)) %>%
  filter(!is.na(NISP) & NISP > 0)

# Extract unique chronological and spatial assignments for each locality
ages_clean <- ages_raw %>%
  select(Locality, Time = Geological_ages, Region = Region_col) %>%
  distinct()

# 4. CREATION OF THE FINAL SPATIO-TEMPORAL ASSEMBLAGE MATRIX
# Merge localities, ages, and traits to build the foundational faunal matrix
assemblage_table <- localities_long %>%
  inner_join(ages_clean, by = "Locality") %>%
  inner_join(traits_clean, by = "TrophicSpecies") %>%
  # Group by temporal, spatial, and biological attributes to aggregate site statistics
  group_by(
    Region, Time, TrophicSpecies, Class, Resolution, 
    Bodymass_min, Bodymass_max, FeedingStrategy, 
    FeedingHabitat, Carrion, HumanResource
  ) %>%
  # Compute site presence count (N_Sites) and dispersion statistics (NISP)
  summarise(
    N_Sites = n_distinct(Locality),
    Sum_NISP = sum(NISP, na.rm = TRUE),
    Min_NISP = min(NISP, na.rm = TRUE),
    Max_NISP = max(NISP, na.rm = TRUE),
    Median_NISP = median(NISP, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Compute the raw relative importance of the taxon across its specific sites
  mutate(
    Relative_Importance = Sum_NISP / N_Sites
  ) %>%
  # 5. CALCULATION OF RELATIVE IMPORTANCE BY TIME AND SPACE
  # Group strictly by spatial-temporal block to derive community-level proportions
  group_by(Region, Time) %>%
  mutate(
    # Compute the total NISP and Total Relative Importance for the current spatial-temporal block
    Total_TimeSpace_NISP = sum(Sum_NISP, na.rm = TRUE),
    Total_TimeSpace_RelImp = sum(Relative_Importance, na.rm = TRUE),
    
    # Derive intra-assemblage percentage proportions
    Spatiotemporal_NISP_Pct = (Sum_NISP / Total_TimeSpace_NISP) * 100,
    Spatiotemporal_RelImp_Pct = (Relative_Importance / Total_TimeSpace_RelImp) * 100
  ) %>%
  ungroup() %>%
  # Remove auxiliary total columns to keep the matrix clean
  select(-Total_TimeSpace_NISP, -Total_TimeSpace_RelImp) %>%
  # Sort hierarchically prioritizing the relative weight of the taxon in its specific block
  arrange(Region, Time, desc(Spatiotemporal_RelImp_Pct))

print("Matrix 'assemblage_table' with spatio-temporal standardization calculated successfully.")


# ----------------------------------------------------------------------------
# 6. DATASET EXPORTS
# ----------------------------------------------------------------------------
print("Generating requested CSV files...")

# Ensure results directory exists
if (!dir.exists("results")) dir.create("results")

# PRE-FILTERING: Discard "Historic" period to strictly evaluate the prehistoric trajectory
ages_filtered <- ages_clean %>%
  filter(Time != "Historic")

# 6.1 Locality Space-Time (Export)
locality_space_time_export <- ages_filtered %>%
  select(Locality = Locality, Period = Time, Region = Region)
write_csv(locality_space_time_export, "results/Locality space-time.csv")

# 6.2 Locality Species Matrix (Export)
# Apply the filtered ages to eliminate historic localities and their occurrences
localities_long_filtered <- localities_long %>%
  semi_join(ages_filtered, by = "Locality") 

locality_species_export <- localities_long_filtered %>%
  pivot_wider(
    names_from = Locality, 
    values_from = NISP, 
    values_fill = 0
  ) %>%
  rename(`Trophic species` = TrophicSpecies)
write_csv(locality_species_export, "results/Locality species.csv")

# 6.3 Species Traits (Export)
# Ensure species existing *only* in Historic localities are dropped from the traits table
species_traits_export <- traits_clean %>%
  semi_join(localities_long_filtered, by = "TrophicSpecies") %>%
  select(
    `Trophic species` = TrophicSpecies,
    Resolution = Resolution,
    `Bodymass (kg)` = Bodymass_max, 
    `Feeding habitat` = FeedingHabitat,
    `Feeding strategy` = FeedingStrategy,
    `Human resource` = HumanResource
  )
write_csv(species_traits_export, "results/Species traits.csv")

print("Files successfully saved to the 'results' folder: Species traits.csv, Locality species.csv, Locality space-time.csv")


# ----------------------------------------------------------------------------
# R SCRIPT: STOCHASTIC TROPHIC NETWORKS (Universal Carnivore Parameters)
# ----------------------------------------------------------------------------
print("Starting network reconstruction...")

# Define stochastic parameters for the allometric probabilistic model (1000 simulations)
n_sims <- 1000
set.seed(123) 

params_sim <- tibble(
  Sim_ID = 1:n_sims,
  alpha = runif(n_sims, min = 1, max = 2),
  beta  = runif(n_sims, min = -2, max = -1),
  gamma = runif(n_sims, min = -1, max = -0.5)
)

prob_cutoff <- 0.0 

# Define network estimation function combining allometric scaling and ecological constraints
estimate_interactions_variance <- function(assemblage_data, human_scenario, scenario_label) {
  
  final_list <- list()
  variability_list <- list()
  strata <- assemblage_data %>% distinct(Region, Time)
  
  # Iterate through every spatial-temporal block
  for (i in 1:nrow(strata)) {
    curr_region <- strata$Region[i]
    curr_time <- strata$Time[i]
    
    local_sp <- assemblage_data %>% filter(Region == curr_region, Time == curr_time)
    
    # Establish effective mass based on human cooperative hunting scenarios
    if (human_scenario == "cooperative") {
      curr_human_mass <- 325 
    } else {
      curr_human_mass <- 65  
    }
    
    # Convert individual body mass into aggregate effective biomass for colonial mollusks
    local_sp <- local_sp %>%
      mutate(
        Bodymass_max = case_when(
          Class == "Bivalvia"   ~ (Sum_NISP / 2) * Bodymass_max,
          Class == "Gastropoda" ~ Sum_NISP * Bodymass_max,
          TRUE                  ~ Bodymass_max
        )
      )
    
    # Isolate functional consumers (Carnivores/Omnivores) while excluding humans initially
    consumers <- local_sp %>%
      filter(FeedingStrategy %in% c("Carnivore", "Omnivore") & 
               FeedingStrategy != "Piscivore" & 
               !(FeedingStrategy == "Carnivore" & FeedingHabitat == "Marine") &
               Carrion == "no" & 
               TrophicSpecies != "Homo_sapiens") %>%
      bind_rows(tibble(
        TrophicSpecies = "Homo_sapiens", 
        Bodymass_max = curr_human_mass, 
        FeedingHabitat = "Ter-Mar", 
        FeedingStrategy = "Omnivore" 
      ))
    
    # Isolate functional resources (Herbivores/Omnivores/Carrion) available for consumption
    resources <- local_sp %>%
      filter(
        (FeedingStrategy %in% c("Herbivore", "Omnivore") | Carrion == "yes" |
           (FeedingStrategy == "Piscivore" & HumanResource %in% c("yes", "potential")) |
           (FeedingStrategy == "Carnivore" & FeedingHabitat == "Marine" & HumanResource %in% c("yes", "potential"))) & 
          TrophicSpecies != "Homo_sapiens"
      ) %>% distinct()
    
    # Generate all combinatoric pairs and compute their stochastic interaction probabilities
    all_pairs <- expand_grid(Consumer = consumers$TrophicSpecies, Resource = resources$TrophicSpecies) %>%
      left_join(consumers %>% select(Consumer = TrophicSpecies, mass_c = Bodymass_max, Habitat_c = FeedingHabitat), by = "Consumer") %>%
      left_join(resources %>% select(Resource = TrophicSpecies, mass_r = Bodymass_max, Carrion_r = Carrion, Human_r = HumanResource, Habitat_r = FeedingHabitat, FeedingStrategy_r = FeedingStrategy), by = "Resource")
    
    all_pairs_sims <- all_pairs %>%
      mutate(Sim_ID = list(1:n_sims)) %>%
      unnest(Sim_ID) %>%
      left_join(params_sim, by = "Sim_ID") %>%
      mutate(
        log_ratio = log10(mass_c) - log10(mass_r),
        z_sim = alpha + beta * log_ratio + gamma * (log_ratio^2),
        prob_sim = exp(z_sim) / (1 + exp(z_sim)),
        
        # Apply strict ecological constraints overriding standard probabilities
        Interaction_prob_sim = case_when(
          Habitat_c != Habitat_r & Habitat_c != "Ter-Mar" & Habitat_r != "Ter-Mar" ~ 0.0,
          Habitat_c == "Marine" & Habitat_r == "Marine" ~ 0.0,
          Consumer != "Homo_sapiens" & (FeedingStrategy_r == "Piscivore" | (FeedingStrategy_r == "Carnivore" & Habitat_r == "Marine")) ~ 0.0,
          Consumer == "Homo_sapiens" & Carrion_r == "yes" ~ 1.0, # Human scavenging on marine necromass = 1.0
          Consumer == "Homo_sapiens" & Human_r %in% c("yes", "potential") & mass_c > mass_r ~ prob_sim,
          Consumer != "Homo_sapiens" & mass_c > mass_r ~ prob_sim,
          TRUE ~ 0.0
        ),
        
        # Classify the type of interaction
        Interaction_type = case_when(
          Consumer == "Homo_sapiens" & Carrion_r == "yes" ~ "Human-scavenging",
          Consumer == "Homo_sapiens" & Human_r %in% c("yes", "potential") ~ "Human-predation",
          Consumer != "Homo_sapiens" ~ "Natural",
          TRUE ~ "Excluded"
        )
      )
    
    # Calculate the central tendency (mean/median) for each interaction across all 1000 simulations
    aggregated_pairs <- all_pairs_sims %>%
      group_by(Consumer, Resource, Interaction_type) %>%
      summarise(sim_vals = list(Interaction_prob_sim), .groups = "drop") %>%
      rowwise() %>%
      mutate(
        Interaction_probability = {
          vals <- unlist(sim_vals)
          if (var(vals) < 1e-8) {
            mean(vals)
          } else {
            pval <- tryCatch(shapiro.test(vals)$p.value, error = function(e) 0)
            if (pval > 0.05) mean(vals) else median(vals)
          }
        }
      ) %>%
      ungroup() %>%
      filter(Interaction_probability > prob_cutoff) %>%
      mutate(Time = curr_time, Region = curr_region, Scenario = scenario_label)
    
    # Retain the raw variability dataset for subsequent statistical testing
    valid_variability <- all_pairs_sims %>%
      inner_join(aggregated_pairs %>% select(Consumer, Resource), by = c("Consumer", "Resource")) %>%
      mutate(Time = curr_time, Region = curr_region, Scenario = scenario_label) %>%
      select(Scenario, Time, Region, Resource, Consumer, Sim_ID, Interaction_prob_sim, Interaction_type)
    
    final_list[[i]] <- aggregated_pairs %>% select(Scenario, Time, Region, Resource, Consumer, Interaction_probability, Interaction_type)
    variability_list[[i]] <- valid_variability
  }
  
  return(list(final = bind_rows(final_list), variability = bind_rows(variability_list)))
}

# Execute simulations and extract results
print("Simulating Scenario 1 (Basal anatomy 65 kg)...")
results_65kg <- estimate_interactions_variance(assemblage_table, "anatomical", "1 (Human 65 kg)")

print("Simulating Scenario 2 (Cooperative Hunting Unit 325 kg and Malacological Biomass)...")
results_cooperative <- estimate_interactions_variance(assemblage_table, "cooperative", "2 (Human 325 kg)")

# Define the cooperative scenario as the definitive baseline for subsequent analyses
variability_trophic_interactions <- results_cooperative$variability
final_trophic_interactions <- results_cooperative$final


# ----------------------------------------------------------------------------
# 8. PREPARATION FOR EXPERT VALIDATION (MANUAL STEP)
# ----------------------------------------------------------------------------
# Note: This section isolates all unique modeled interactions and creates an empty 
# 'Expert_validation' column.
# In a real-world workflow, this file is exported, evaluated offline 
# by an expert (who fills the column with "Valid" or "Invalid" based on empirical 
# ecological/archaeological evidence), and then re-imported in the following step.

print("Generating pre-validation interaction file for expert review...")

interactions_to_validate <- variability_trophic_interactions %>%
  distinct(Time, Region, Consumer, Resource, Interaction_type) %>%
  mutate(Expert_validation = "") # Empty column ready for manual input

write_csv(interactions_to_validate, "results/TDF_Scenario2_Interactions_ToBeValidated.csv")
print("File for manual expert validation saved as: results/TDF_Scenario2_Interactions_ToBeValidated.csv")


# ----------------------------------------------------------------------------
# 9. TOPOLOGICAL FILTERING BY EXCLUSION (DISCARDING INVALID LINKS)
# ----------------------------------------------------------------------------
# NOTE: The following lines load the manually validated file. For reproducibility,
# we assume the file "TDF_Scenario2_Interactions_Validated_Apr2026.csv" has been 
# populated by an expert and placed in the 'results' directory.

print("Loading expert-validated interactions...")
expert_validation_df <- read_csv("results/TDF_Scenario2_Interactions_Validated_Apr2026.csv")

# Extract invalid interaction keys to create an exclusion list
invalid_keys <- expert_validation_df %>%
  filter(Expert_validation == "Invalid") %>%
  select(Time, Region, Consumer, Resource, Interaction_type) %>%
  distinct()

# Apply the exclusion filter using anti_join
validated_variability <- variability_trophic_interactions %>%
  anti_join(invalid_keys, by = c("Time", "Region", "Consumer", "Resource", "Interaction_type"))


# ----------------------------------------------------------------------------
# 10. TOPOLOGICAL CONSENSUS FILTERING
# ----------------------------------------------------------------------------
print("Applying topological consensus filter based on stochastic variability...")

# Evaluate structural stability across the 1000 stochastic runs
consensus_variability <- validated_variability %>%
  group_by(Scenario, Time, Region, Consumer, Resource, Interaction_type) %>%
  mutate(
    Total_Sims = n(),
    Sims_Below_Threshold = sum(Interaction_prob_sim < 0.01, na.rm = TRUE),
    Prop_Below = Sims_Below_Threshold / Total_Sims
  ) %>%
  ungroup() %>%
  # Discard the entire interaction link if > 50% of its simulations fall below 0.01
  filter(Prop_Below <= 0.50) %>%
  select(-Total_Sims, -Sims_Below_Threshold, -Prop_Below)

# Calculate final aggregated parameters for the remaining robust interactions
final_validated_interactions <- consensus_variability %>%
  group_by(Scenario, Time, Region, Consumer, Resource, Interaction_type) %>%
  summarise(
    Median_Probability = median(Interaction_prob_sim, na.rm = TRUE),
    Mean_Probability = mean(Interaction_prob_sim, na.rm = TRUE),
    n_sims_retained = n(),
    .groups = "drop"
  )

write_csv(final_validated_interactions, "results/Definitive_Consensus_Validated_Networks.csv")


# ----------------------------------------------------------------------------
# 11. BIPARTITE NETWORKS VISUALIZATION
# ----------------------------------------------------------------------------
# Load graphic network libraries
if (!require("tidygraph")) install.packages("tidygraph")
if (!require("ggraph")) install.packages("ggraph")
if (!require("patchwork")) install.packages("patchwork")
if (!require("cowplot")) install.packages("cowplot")

library(tidygraph)
library(ggraph)
library(patchwork)
library(cowplot)

print("Starting reconstruction with unified legend extracted from assemblage_table...")

# Extract median probabilities for the baseline bipartite visual structure
median_networks <- final_validated_interactions %>%
  group_by(Scenario, Time, Region, Resource, Consumer, Interaction_type) %>%
  summarise(
    Median_Prob = median(Mean_Probability, na.rm = TRUE), 
    .groups = "drop"
  )

# Extract habitat information for node styling
traits_habitat <- assemblage_table %>%
  select(TrophicSpecies, FeedingHabitat) %>%
  distinct(TrophicSpecies, .keep_all = TRUE)

# Define universal categorical palette for nodes
habitat_colors <- c(
  "Homo sapiens" = "#33CCFF",   
  "Lama guanicoe" = "#B40426",  
  "Marine" = "#3182bd",         
  "Terrestrial" = "#2ca25f",    
  "Ter-Mar" = "grey50"          
)

# Build a phantom plot strictly to extract a universally standardized legend
dummy_df <- tibble(
  Category = factor(names(habitat_colors), levels = names(habitat_colors))
)

dummy_plot <- ggplot(dummy_df, aes(x = 1, y = 1, fill = Category)) +
  geom_point(shape = 21, size = 5, color = "black") +
  scale_fill_manual(
    name = "Species / Habitat", 
    values = habitat_colors,
    drop = FALSE
  ) +
  theme_void() +
  theme(legend.position = "right", 
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 11))

master_category_legend <- get_legend(dummy_plot)

# Function to generate individual bipartite plots
build_bipartite_plot <- function(net_data, custom_title) {
  
  # Distinguish consumers and resources to form bipartite layers
  consumers <- net_data %>% 
    distinct(Consumer) %>% 
    transmute(name = paste0(Consumer, "_C"), type = TRUE, species = Consumer)
  
  resources <- net_data %>% 
    distinct(Resource) %>% 
    transmute(name = paste0(Resource, "_R"), type = FALSE, species = Resource)
  
  # Assign ecological categories for nodal coloring
  nodes_df <- bind_rows(consumers, resources) %>%
    left_join(traits_habitat, by = c("species" = "TrophicSpecies")) %>%
    mutate(
      Category = case_when(
        species == "Homo_sapiens" ~ "Homo sapiens",
        species == "Lama_guanicoe" ~ "Lama guanicoe",
        grepl("ter-mar", FeedingHabitat, ignore.case = TRUE) ~ "Ter-Mar",
        grepl("marine", FeedingHabitat, ignore.case = TRUE) ~ "Marine",
        grepl("terrestrial", FeedingHabitat, ignore.case = TRUE) ~ "Terrestrial",
        TRUE ~ "Ter-Mar" 
      ),
      Category = factor(Category, levels = names(habitat_colors))
    )
  
  # Format edges and attach interaction probabilities
  edges_df <- net_data %>% 
    transmute(
      from = paste0(Consumer, "_C"), 
      to = paste0(Resource, "_R"),
      Median_Prob = Median_Prob
    )
  
  # Construct and layout the graph
  g <- tbl_graph(nodes = nodes_df, edges = edges_df, directed = TRUE)
  layout <- create_layout(g, layout = 'bipartite')
  layout$y <- ifelse(layout$type, 1, 0)
  
  # Render the plot
  p <- ggraph(layout) +
    geom_edge_link(aes(color = Median_Prob, edge_alpha = Median_Prob, edge_width = Median_Prob)) +
    geom_node_point(aes(fill = Category), shape = 21, size = 5, color = "black", stroke = 0.5) +
    scale_edge_color_gradient(
      low = "blue", high = "red",
      name = "Interaction probability",
      limits = c(0.01, 1.0),
      guide = guide_edge_colourbar(barwidth = 1.5, barheight = 10)
    ) +
    scale_edge_width_continuous(range = c(0.4, 2.0), guide = "none") +
    scale_edge_alpha_continuous(range = c(0.3, 1), guide = "none") +
    scale_fill_manual(values = habitat_colors, guide = "none") +
    labs(title = custom_title) +
    theme_graph(base_family = "sans") +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
    )
  
  return(p)
}

# Define spatial-temporal strata constraints for visualization loop
strata_definitions <- tibble(
  Time = c("Pleistocene final", "Holocene mid", "Holocene final", "Holocene final"),
  Region = c("Steppe", "Steppe", "Steppe", "Forest"),
  Custom_Title = c("late Pleistocene - Steppe", 
                   "middle Holocene - Steppe",
                   "late Holocene - Steppe", 
                   "late Holocene - Forest")
)

plot_list <- list()

# Generate and store each block's network
for (i in 1:nrow(strata_definitions)) {
  curr_net <- median_networks %>% 
    filter(Time == strata_definitions$Time[i], Region == strata_definitions$Region[i])
  
  if (nrow(curr_net) > 0) {
    plot_list[[i]] <- build_bipartite_plot(curr_net, strata_definitions$Custom_Title[i])
  }
}

# Assemble final figure matrix with shared master legend
combined_networks <- wrap_plots(plot_list, ncol = 2) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

final_figure <- plot_grid(
  combined_networks,
  master_category_legend,
  ncol = 2,
  rel_widths = c(0.85, 0.15) 
)

output_pdf <- "results/TierraDelFuego_Bipartite_Networks.pdf"
ggsave(output_pdf, final_figure, width = 18, height = 14, device = "pdf")
ggsave("results/TierraDelFuego_Bipartite_Networks.png", final_figure, width = 15, height = 12, bg = "white")

# Export standardized interactions matrix as Supplementary Material
median_networks_supmat <- final_validated_interactions %>%
  filter(Time != "Historic") %>%
  group_by(
    Period = Time, 
    Region, 
    Resource, 
    Consumer, 
    `Interaction type` = Interaction_type
  ) %>%
  summarise(
    `Interaction median` = round(median(Mean_Probability, na.rm = TRUE), 3),
    .groups = "drop"
  )
write_csv(median_networks_supmat, "results/TDF_Validated_Interactions_SupMat.csv")


# ----------------------------------------------------------------------------
# 12. CALCULATION OF STRUCTURAL PROPERTIES AND COMPLEXITY METRICS
# ----------------------------------------------------------------------------
if (!require("bipartite")) install.packages("bipartite")
library(bipartite)

print("Starting topological calculation with scaling factor (E)...")

# Calculate comprehensive bipartite metrics for every valid stochastic iteration
network_metrics_sim <- consensus_variability %>%
  mutate(Int_prob = ifelse(Interaction_prob_sim < 0.01, 0, Interaction_prob_sim)) %>%
  filter(Int_prob > 0) %>%
  group_by(Time, Region, Sim_ID) %>%
  group_modify(~ {
    
    # Construct incidence matrix for networklevel calculations
    inc_mat <- .x %>%
      select(Resource, Consumer, Int_prob) %>%
      pivot_wider(names_from = Consumer, values_from = Int_prob, values_fill = 0) %>%
      column_to_rownames("Resource") %>%
      as.matrix()
    
    # Assert minimum viable dimensions (2x2) for stable topology metrics
    if(nrow(inc_mat) < 2 || ncol(inc_mat) < 2) {
      return(tibble(
        Connectance = NA_real_, Vulnerability = NA_real_, Generality = NA_real_, 
        Nestedness = NA_real_, Interaction_Diversity = NA_real_, 
        Functional_Complementarity_Predators = NA_real_, 
        Functional_Complementarity_Prey = NA_real_,
        Predator_Prey_Ratio = NA_real_,
        Total_Resources = NA_integer_,
        Total_Consumers = NA_integer_,
        Total_Interactions = NA_integer_
      ))
    }
    
    # Extract raw counts and ratio
    num_prey <- nrow(inc_mat)
    num_predators <- ncol(inc_mat)
    num_interactions <- sum(inc_mat > 0) 
    ratio_predator_prey <- num_predators / num_prey
    
    # Apply absolute scaling factor (E) to normalize incidence matrix dimensions
    E_factor <- num_predators * num_prey * 3
    inc_mat_scaled <- inc_mat * E_factor
    
    # Compute classic bipartite indices
    mets <- tryCatch({
      networklevel(inc_mat_scaled, index = c("connectance", "vulnerability", "generality", 
                                             "weighted nestedness", "Shannon diversity"))
    }, error = function(e) {
      c("connectance"=NA, "vulnerability.LL"=NA, "generality.HL"=NA, 
        "weighted nestedness"=NA, "Shannon diversity"=NA)
    })
    
    # Compute functional complementarity metrics
    fc_prey <- tryCatch(as.numeric(fc(inc_mat_scaled)), error = function(e) NA_real_)
    fc_predators <- tryCatch(as.numeric(fc(t(inc_mat_scaled))), error = function(e) NA_real_)
    
    # Compile metrics array for current stochastic run
    tibble(
      Connectance = as.numeric(mets["connectance"]),
      Vulnerability = as.numeric(mets["vulnerability.LL"]),
      Generality = as.numeric(mets["generality.HL"]),
      Nestedness = as.numeric(mets["weighted nestedness"]),
      Interaction_Diversity = as.numeric(mets["Shannon diversity"]),
      Functional_Complementarity_Predators = fc_predators,
      Functional_Complementarity_Prey = fc_prey,
      Predator_Prey_Ratio = ratio_predator_prey,
      Total_Resources = num_prey,
      Total_Consumers = num_predators,
      Total_Interactions = num_interactions
    )
  }) %>%
  ungroup() %>%
  filter(!is.na(Connectance))

print(paste("Calculation completed. Retained stochastic networks:", nrow(network_metrics_sim)))
print("Preparing graphical architecture for network properties (excluding absolute counts)...")

# Structure dataset into tidy long format to facilitate faceted ggplotting
metrics_long <- network_metrics_sim %>%
  filter(Time != "Historic") %>%
  pivot_longer(
    cols = c(
      Connectance, Vulnerability, Generality, Nestedness, 
      Interaction_Diversity, Functional_Complementarity_Predators, 
      Functional_Complementarity_Prey, Predator_Prey_Ratio,
      Total_Resources, Total_Consumers, Total_Interactions
    ),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    # Human-readable metric formatting
    Metric = case_when(
      Metric == "Total_Resources" ~ "Total Resources",
      Metric == "Total_Consumers" ~ "Total Consumers",
      Metric == "Total_Interactions" ~ "Total Interactions",
      Metric == "Connectance" ~ "Connectance",
      Metric == "Interaction_Diversity" ~ "Interaction diversity",
      Metric == "Nestedness" ~ "Nestedness",
      Metric == "Predator_Prey_Ratio" ~ "Consumer-Resource ratio",
      Metric == "Generality" ~ "Generality",
      Metric == "Vulnerability" ~ "Vulnerability",
      Metric == "Functional_Complementarity_Predators" ~ "Func. comp. Consumer",
      Metric == "Functional_Complementarity_Prey" ~ "Func. comp. Resource",
      TRUE ~ Metric
    ),
    # Strict factorization to control facet rendering order
    Metric = factor(Metric, levels = c(
      "Total Consumers", "Total Resources", "Total Interactions", "Connectance", 
      "Interaction diversity", "Nestedness", "Consumer-Resource ratio", "Generality", 
      "Vulnerability", "Func. comp. Consumer", "Func. comp. Resource"
    )),
    Time = case_when(
      Time == "Pleistocene final" ~ "late Pleistocene",
      Time == "Holocene mid" ~ "mid Holocene",
      Time == "Holocene final" ~ "late Holocene",
      TRUE ~ Time
    ),
    Time = factor(Time, levels = c("late Pleistocene", "mid Holocene", "late Holocene"))
  )

# Extract only relative metrics for boxplot visualization (ignoring absolute count metrics)
plot_data <- metrics_long %>%
  filter(!Metric %in% c("Total Resources", "Total Consumers", "Total Interactions")) %>%
  mutate(Metric = droplevels(Metric)) 

p_network_properties <- ggplot(plot_data, aes(x = Time, y = Value, fill = Region)) +
  geom_boxplot(alpha = 0.85, outlier.size = 0.8, outlier.alpha = 0.4, color = "black") +
  facet_wrap(~ Metric, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Steppe" = "#E3D3B5", "Forest" = "#2D5F43")) +
  labs(
    x = "Period",
    y = "Variability of network property",
    fill = "Habitat"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )

output_file <- "results/TierraDelFuego_Network_Properties_Scaled_Extended.pdf"
ggsave(output_file, p_network_properties, width = 18, height = 10, device = "pdf")
ggsave("results/TierraDelFuego_Network_Properties_Scaled_Extended.png", p_network_properties, width=15, height=12, bg="white")

print(paste("Network properties figure successfully generated as:", output_file))

# Summarize central tendency (median) of raw counts for explicit reporting
print("Generating summary table for network complexity medians...")

complexity_summary <- network_metrics_sim %>%
  filter(Time != "Historic") %>%
  group_by(Time, Region) %>%
  summarise(
    `Total Consumers` = round(median(Total_Consumers, na.rm = TRUE), 1),
    `Total Resources` = round(median(Total_Resources, na.rm = TRUE), 1),
    `Total Interactions` = round(median(Total_Interactions, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  mutate(
    Time = case_when(
      Time == "Pleistocene final" ~ "late Pleistocene",
      Time == "Holocene mid" ~ "mid Holocene",
      Time == "Holocene final" ~ "late Holocene",
      TRUE ~ Time
    ),
    Time = factor(Time, levels = c("late Pleistocene", "mid Holocene", "late Holocene"))
  ) %>%
  arrange(Time, desc(Region))

output_csv_complexity <- "results/TierraDelFuego_Network_Complexity_Medians.csv"
write_csv(complexity_summary, output_csv_complexity)

print(paste("Complexity metrics summary successfully saved as:", output_csv_complexity))


# ----------------------------------------------------------------------------
# 13. HUMAN EGO-NETWORK ANALYSIS (DIET COMPOSITION)
# ----------------------------------------------------------------------------
if (!require("patchwork")) install.packages("patchwork")
library(patchwork)

print("Starting Human ego-network analysis (Proportional Diet Composition)...")

# Define taxa representing strictly marine carrion subsidy
cetacean_taxa <- c("Physeter_macrocephalus", "Balaenoptera_borealis", "Cetacea", "Grampus_griseus", "Phocoena")

# Map defined taxa uniformly as 'Necromass' and aggregate their interaction probabilities
merged_interactions <- final_validated_interactions %>%
  mutate(
    Resource = ifelse(Resource %in% cetacean_taxa, "Necromass", Resource)
  ) %>%
  group_by(Scenario, Time, Region, Consumer, Resource, Interaction_type) %>%
  summarise(
    Median_Probability = max(Median_Probability, na.rm = TRUE),
    Mean_Probability = max(Mean_Probability, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate theoretical diet proportional contribution including Necromass
human_diet_with_necro <- merged_interactions %>%
  filter(Time != "Historic", Consumer == "Homo_sapiens") %>%
  group_by(Time, Region, Resource) %>%
  summarise(Interaction_Strength = sum(Median_Probability, na.rm = TRUE), .groups = "drop") %>%
  group_by(Time, Region) %>%
  mutate(
    Total_Human_Strength = sum(Interaction_Strength, na.rm = TRUE),
    Proportion = (Interaction_Strength / Total_Human_Strength) * 100
  ) %>%
  ungroup() %>%
  mutate(Scenario = "Including Marine Necromass")

# Calculate theoretical diet proportional contribution strictly ignoring Necromass
human_diet_no_necro <- merged_interactions %>%
  filter(Time != "Historic", Consumer == "Homo_sapiens", Resource != "Necromass") %>%
  group_by(Time, Region, Resource) %>%
  summarise(Interaction_Strength = sum(Median_Probability, na.rm = TRUE), .groups = "drop") %>%
  group_by(Time, Region) %>%
  mutate(
    Total_Human_Strength = sum(Interaction_Strength, na.rm = TRUE),
    Proportion = (Interaction_Strength / Total_Human_Strength) * 100
  ) %>%
  ungroup() %>%
  mutate(Scenario = "Excluding Marine Necromass")

# Combine both conceptual models and structure formatting
network_levels <- c(
  "late Pleistocene\nSteppe", 
  "mid Holocene\nSteppe", 
  "late Holocene\nSteppe", 
  "late Holocene\nForest"
)

combined_human_diet <- bind_rows(human_diet_with_necro, human_diet_no_necro) %>%
  mutate(
    Resource = str_replace_all(Resource, "_", " "),
    
    Period = case_when(
      Time == "Pleistocene final" ~ "late Pleistocene",
      Time == "Holocene mid" ~ "mid Holocene",
      Time == "Holocene final" ~ "late Holocene",
      TRUE ~ Time
    ),
    Period = factor(Period, levels = c("late Pleistocene", "mid Holocene", "late Holocene")),
    Network_Label = factor(paste0(Period, "\n", Region), levels = network_levels),
    
    # Reclassify marginal resources grouping those under 5% to streamline the plot legend
    Plot_Resource = case_when(
      Resource == "Lama guanicoe" ~ "Lama guanicoe",
      Resource == "Necromass" ~ "Marine Necromass",
      Proportion >= 5.0 ~ Resource, 
      TRUE ~ "Other Resources (<5%)"
    )
  )

# Export calculated diet arrays for external inspection
output_csv <- "results/TDF_Human_Diet_Proportions.csv"
write_csv(combined_human_diet %>% select(Scenario, Period, Region, Resource, Interaction_Strength, Proportion), output_csv)
print(paste("Human diet proportion dataset saved as:", output_csv))

# Generate categorical palette establishing focal interest colors
unique_resources <- unique(combined_human_diet$Plot_Resource)
base_palette <- scales::hue_pal()(length(unique_resources))
names(base_palette) <- unique_resources

base_palette["Lama guanicoe"] <- "#B40426"    
base_palette["Marine Necromass"] <- "#8B4513" 
if("Other Resources (<5%)" %in% names(base_palette)) {
  base_palette["Other Resources (<5%)"] <- "grey80"
}

combined_human_diet <- combined_human_diet %>%
  mutate(Plot_Resource = factor(Plot_Resource, levels = names(base_palette)))

# Construct final 100% Stacked Bar Visualization
print("Generating Human Diet Composition stacked bar chart...")

p_diet_bar <- ggplot(combined_human_diet, aes(x = Network_Label, y = Proportion, fill = Plot_Resource)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = base_palette, name = "Human Resource") +
  facet_wrap(~ Scenario, ncol = 2) +
  labs(
    x = "Period-Habitat",
    y = "Proportional Contribution"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "right",
    panel.grid.major.x = element_blank()
  )

ggsave("results/TierraDelFuego_Human_Diet_Proportions_Bar.pdf", p_diet_bar, width = 12, height = 8, device = "pdf")
ggsave("results/TierraDelFuego_Human_Diet_Proportions_Bar.png", p_diet_bar, width = 12, height = 8, bg = "white")

print("All tasks completed successfully.")
