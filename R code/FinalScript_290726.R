# ----------------------------------------------------------------------------
# R SCRIPT: COMPLETE SPATIO-TEMPORAL ASSEMBLAGE 
# ----------------------------------------------------------------------------

# Load standard libraries
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("readxl")) install.packages("readxl")

library(tidyverse)
library(readxl)

print("Starting structuring of the global matrix...")

# 1. DATA LOADING
traits_raw <- read_excel("data/Species_traits_TDF_270426.xlsx", sheet = "traits")
localities_raw <- read_excel("data/Species_localities_TDF_270526.xlsx")
ages_raw <- read_excel("data/Age_localities_TDF_270526.xlsx")

# 2. FILTERING AND CLEANING OF FUNCTIONAL TRAITS
traits_clean <- traits_raw %>%
  # The exclusion filter for indeterminate taxa has been removed
  mutate(
    # Taxonomic resolution assignment considering broad taxa
    Resolution = case_when(
      is.na(Genus) | Genus == "" ~ "Indeterminate/Broad",
      str_detect(TrophicSpecies, "_| ") ~ "Species",
      TRUE ~ "Genus"
    ),
    Carrion = ifelse(is.na(Carrion), "no", Carrion),
    HumanResource = ifelse(is.na(HumanResource), "no", HumanResource)
  ) %>%
  select(
    TrophicSpecies, Resolution, Bodymass_min = BodySize_min,
    Bodymass_max = BodySize_max, FeedingStrategy, FeedingHabitat,
    Carrion, HumanResource
  ) %>%
  distinct(TrophicSpecies, .keep_all = TRUE)

# 3. PROCESSING OF LOCALITIES AND AGES
localities_long <- localities_raw %>%
  pivot_longer(
    cols = -TrophicSpecies, 
    names_to = "Locality", 
    values_to = "NISP"
  ) %>%
  mutate(NISP = as.numeric(NISP)) %>%
  filter(!is.na(NISP) & NISP > 0)

ages_clean <- ages_raw %>%
  select(Locality, Time = Geological_ages, Region = Region_col) %>%
  distinct()

# 4. CREATION OF THE FINAL SPATIO-TEMPORAL ASSEMBLAGE MATRIX
assemblage_table <- localities_long %>%
  inner_join(ages_clean, by = "Locality") %>%
  inner_join(traits_clean, by = "TrophicSpecies") %>%
  # Grouping by temporal, spatial, and biological matrices
  group_by(
    Region, Time, TrophicSpecies, Resolution, 
    Bodymass_min, Bodymass_max, FeedingStrategy, 
    FeedingHabitat, Carrion, HumanResource
  ) %>%
  # Compute site count (N_Sites) and dispersion statistics (NISP)
  summarise(
    N_Sites = n_distinct(Locality),
    Sum_NISP = sum(NISP, na.rm = TRUE),
    Min_NISP = min(NISP, na.rm = TRUE),
    Max_NISP = max(NISP, na.rm = TRUE),
    Median_NISP = median(NISP, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Compute the raw relative importance of the taxon
  mutate(
    Relative_Importance = Sum_NISP / N_Sites
  ) %>%
  # 5. CALCULATION OF RELATIVE IMPORTANCE BY TIME AND SPACE
  # Group again, this time exclusively by spatio-temporal stratum
  group_by(Region, Time) %>%
  mutate(
    # Compute the total sum of NISP and Relative Importance for the current Region/Time
    Total_TimeSpace_NISP = sum(Sum_NISP, na.rm = TRUE),
    Total_TimeSpace_RelImp = sum(Relative_Importance, na.rm = TRUE),
    
    # Derive intra-assemblage percentage proportions (analogous to NISP %)
    Spatiotemporal_NISP_Pct = (Sum_NISP / Total_TimeSpace_NISP) * 100,
    Spatiotemporal_RelImp_Pct = (Relative_Importance / Total_TimeSpace_RelImp) * 100
  ) %>%
  # Ungroup to avoid subsequent matrix conflicts
  ungroup() %>%
  # Clean auxiliary total columns to keep the matrix clean (Optional)
  select(-Total_TimeSpace_NISP, -Total_TimeSpace_RelImp) %>%
  # Sort hierarchically prioritizing the relative weight of the taxon in its specific stratum
  arrange(Region, Time, desc(Spatiotemporal_RelImp_Pct))

# Save or visualize the exported table
print("Matrix 'assemblage_table' with spatio-temporal standardization calculated successfully.")


# ----------------------------------------------------------------------------
# 6. DATASET EXPORTS
# ----------------------------------------------------------------------------
print("Generating requested CSV files...")

# Ensure results directory exists
if (!dir.exists("results")) dir.create("results")

# PRE-FILTERING: Discard "Historic" period
ages_filtered <- ages_clean %>%
  filter(Time != "Historic")

# 6.3 Locality Space-Time
locality_space_time_export <- ages_filtered %>%
  select(
    Locality = Locality,
    Period = Time,
    Region = Region
  )
write_csv(locality_space_time_export, "results/Locality space-time.csv")

# 6.2 Locality Species (Matrix)
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

# 6.1 Species Traits
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
# CÓDIGO DE DIAGNÓSTICO: DETECTAR DISCREPANCIAS ENTRE FAUNA Y CRONOLOGÍA
# ----------------------------------------------------------------------------

sitios_en_fauna <- colnames(localities_raw)[-1]
sitios_en_cronologia <- unique(ages_raw$Locality) 

missing_in_fauna <- setdiff(sitios_en_cronologia, sitios_en_fauna)
missing_in_cronologia <- setdiff(sitios_en_fauna, sitios_en_cronologia)

cat("\n=== RESULTADOS DEL DIAGNÓSTICO DE SITIOS ===\n")

if(length(missing_in_fauna) > 0) {
  cat("\n[ALERTA] Sitios presentes en CRONOLOGÍA pero AUSENTES en columnas de Fauna:\n")
  print(missing_in_fauna)
} else {
  cat("\n[OK] Todos los sitios de la cronología existen en la matriz de fauna.\n")
}

if(length(missing_in_cronologia) > 0) {
  cat("\n[ALERTA] Sitios presentes en FAUNA pero AUSENTES en la lista de Cronología:\n")
  print(missing_in_cronologia)
} else {
  cat("\n[OK] Todos los sitios de la fauna están registrados en la cronología.\n")
}


# ----------------------------------------------------------------------------
# R SCRIPT: STOCHASTIC TROPHIC NETWORKS (Universal Carnivore Parameters)
# ----------------------------------------------------------------------------

print("Starting network reconstruction...")

# 1. UPDATING THE ASSEMBLAGE MATRIX
traits_clean <- traits_raw %>%
  mutate(
    Resolution = case_when(
      is.na(Genus) | Genus == "" ~ "Indeterminate/Broad",
      str_detect(TrophicSpecies, "_| ") ~ "Species",
      TRUE ~ "Genus"
    ),
    Carrion = ifelse(is.na(Carrion), "no", Carrion),
    HumanResource = ifelse(is.na(HumanResource), "no", HumanResource)
  ) %>%
  select(
    TrophicSpecies, Class, Resolution, Bodymass_min = BodySize_min,
    Bodymass_max = BodySize_max, FeedingStrategy, FeedingHabitat,
    Carrion, HumanResource
  ) %>%
  distinct(TrophicSpecies, .keep_all = TRUE)

assemblage_table <- localities_long %>%
  inner_join(ages_clean, by = "Locality") %>%
  inner_join(traits_clean, by = "TrophicSpecies") %>%
  group_by(
    Region, Time, TrophicSpecies, Class, Resolution, 
    Bodymass_min, Bodymass_max, FeedingStrategy, 
    FeedingHabitat, Carrion, HumanResource
  ) %>%
  summarise(
    N_Sites = n_distinct(Locality),
    Sum_NISP = sum(NISP, na.rm = TRUE),
    Min_NISP = min(NISP, na.rm = TRUE),
    Max_NISP = max(NISP, na.rm = TRUE),
    Median_NISP = median(NISP, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Relative_Importance = Sum_NISP / N_Sites) %>%
  group_by(Region, Time) %>%
  mutate(
    Total_TimeSpace_NISP = sum(Sum_NISP, na.rm = TRUE),
    Total_TimeSpace_RelImp = sum(Relative_Importance, na.rm = TRUE),
    Spatiotemporal_NISP_Pct = (Sum_NISP / Total_TimeSpace_NISP) * 100,
    Spatiotemporal_RelImp_Pct = (Relative_Importance / Total_TimeSpace_RelImp) * 100
  ) %>%
  ungroup() %>%
  select(-Total_TimeSpace_NISP, -Total_TimeSpace_RelImp)

# 2. DEFINITION OF STOCHASTIC PARAMETERS
n_sims <- 1000
set.seed(123) 

params_sim <- tibble(
  Sim_ID = 1:n_sims,
  alpha = runif(n_sims, min = 1, max = 2),
  beta  = runif(n_sims, min = -2, max = -1),
  gamma = runif(n_sims, min = -1, max = -0.5)
)

prob_cutoff <- 0.0 

# 3. ADJUSTED NETWORK ESTIMATION FUNCTION
estimate_interactions_variance <- function(assemblage_data, human_scenario, scenario_label) {
  
  final_list <- list()
  variability_list <- list()
  strata <- assemblage_data %>% distinct(Region, Time)
  
  for (i in 1:nrow(strata)) {
    curr_region <- strata$Region[i]
    curr_time <- strata$Time[i]
    
    local_sp <- assemblage_data %>% filter(Region == curr_region, Time == curr_time)
    
    if (human_scenario == "cooperative") {
      curr_human_mass <- 325 
    } else {
      curr_human_mass <- 65  
    }
    
    local_sp <- local_sp %>%
      mutate(
        Bodymass_max = case_when(
          Class == "Bivalvia"   ~ (Sum_NISP / 2) * Bodymass_max,
          Class == "Gastropoda" ~ Sum_NISP * Bodymass_max,
          TRUE                  ~ Bodymass_max
        )
      )
    
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
    
    resources <- local_sp %>%
      filter(
        (FeedingStrategy %in% c("Herbivore", "Omnivore") | Carrion == "yes" |
           (FeedingStrategy == "Piscivore" & HumanResource %in% c("yes", "potential")) |
           (FeedingStrategy == "Carnivore" & FeedingHabitat == "Marine" & HumanResource %in% c("yes", "potential"))) & 
          TrophicSpecies != "Homo_sapiens"
      ) %>% distinct()
    
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
        
        Interaction_prob_sim = case_when(
          Habitat_c != Habitat_r & Habitat_c != "Ter-Mar" & Habitat_r != "Ter-Mar" ~ 0.0,
          Habitat_c == "Marine" & Habitat_r == "Marine" ~ 0.0,
          Consumer != "Homo_sapiens" & (FeedingStrategy_r == "Piscivore" | (FeedingStrategy_r == "Carnivore" & Habitat_r == "Marine")) ~ 0.0,
          Consumer == "Homo_sapiens" & Carrion_r == "yes" ~ 1.0,
          Consumer == "Homo_sapiens" & Human_r %in% c("yes", "potential") & mass_c > mass_r ~ prob_sim,
          Consumer != "Homo_sapiens" & mass_c > mass_r ~ prob_sim,
          TRUE ~ 0.0
        ),
        
        Interaction_type = case_when(
          Consumer == "Homo_sapiens" & Carrion_r == "yes" ~ "Human-scavenging",
          Consumer == "Homo_sapiens" & Human_r %in% c("yes", "potential") ~ "Human-predation",
          Consumer != "Homo_sapiens" ~ "Natural",
          TRUE ~ "Excluded"
        )
      )
    
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
    
    valid_variability <- all_pairs_sims %>%
      inner_join(aggregated_pairs %>% select(Consumer, Resource), by = c("Consumer", "Resource")) %>%
      mutate(Time = curr_time, Region = curr_region, Scenario = scenario_label) %>%
      select(Scenario, Time, Region, Resource, Consumer, Sim_ID, Interaction_prob_sim, Interaction_type)
    
    final_list[[i]] <- aggregated_pairs %>% select(Scenario, Time, Region, Resource, Consumer, Interaction_probability, Interaction_type)
    variability_list[[i]] <- valid_variability
  }
  
  return(list(final = bind_rows(final_list), variability = bind_rows(variability_list)))
}

print("Simulating Scenario 1 (Basal anatomy 65 kg)...")
results_65kg <- estimate_interactions_variance(assemblage_table, "anatomical", "1 (Human 65 kg)")

print("Simulating Scenario 2 (Cooperative Hunting Unit 325 kg and Malacological Biomass)...")
results_cooperative <- estimate_interactions_variance(assemblage_table, "cooperative", "2 (Human 325 kg)")

variability_trophic_interactions <- results_cooperative$variability
final_trophic_interactions <- results_cooperative$final

# ----------------------------------------------------------------------------
# R SCRIPT: TOPOLOGICAL FILTERING BY EXCLUSION (DISCARDING INVALID LINKS)
# ----------------------------------------------------------------------------

expert_validation_df <- read_csv("results/TDF_Scenario2_Interactions_Validated_Apr2026.csv")

invalid_keys <- expert_validation_df %>%
  filter(Expert_validation == "Invalid") %>%
  select(Time, Region, Consumer, Resource, Interaction_type) %>%
  distinct()

validated_variability <- variability_trophic_interactions %>%
  anti_join(invalid_keys, by = c("Time", "Region", "Consumer", "Resource", "Interaction_type"))

# ----------------------------------------------------------------------------
# R SCRIPT: TOPOLOGICAL CONSENSUS FILTERING
# ----------------------------------------------------------------------------

print("Applying topological consensus filter based on stochastic variability...")

consensus_variability <- validated_variability %>%
  group_by(Scenario, Time, Region, Consumer, Resource, Interaction_type) %>%
  mutate(
    Total_Sims = n(),
    Sims_Below_Threshold = sum(Interaction_prob_sim < 0.01, na.rm = TRUE),
    Prop_Below = Sims_Below_Threshold / Total_Sims
  ) %>%
  ungroup() %>%
  filter(Prop_Below <= 0.50) %>%
  select(-Total_Sims, -Sims_Below_Threshold, -Prop_Below)

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
# R SCRIPT: BIPARTITE NETWORKS
# ----------------------------------------------------------------------------

# Load required libraries
if (!require("tidygraph")) install.packages("tidygraph")
if (!require("ggraph")) install.packages("ggraph")
if (!require("patchwork")) install.packages("patchwork")
if (!require("cowplot")) install.packages("cowplot")

library(tidygraph)
library(ggraph)
library(patchwork)
library(cowplot)

print("Starting reconstruction with unified legend extracted from assemblage_table...")

median_networks <- final_validated_interactions %>%
  group_by(Scenario, Time, Region, Resource, Consumer, Interaction_type) %>%
  summarise(
    Median_Prob = median(Mean_Probability, na.rm = TRUE), 
    .groups = "drop"
  )

traits_habitat <- assemblage_table %>%
  select(TrophicSpecies, FeedingHabitat) %>%
  distinct(TrophicSpecies, .keep_all = TRUE)

habitat_colors <- c(
  "Homo sapiens" = "#33CCFF",   
  "Lama guanicoe" = "#B40426",  
  "Marine" = "#3182bd",         
  "Terrestrial" = "#2ca25f",    
  "Ter-Mar" = "grey50"          
)

dummy_df <- tibble(
  Category = factor(names(habitat_colors), levels = names(habitat_colors))
)

dummy_plot <- ggplot(dummy_df, aes(x = 1, y = 1, fill = Category)) +
  geom_point(shape = 21, size = 5, color = "black") +
  # --- CHANGED: Legend title updated to "Species / Habitat"
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

build_bipartite_plot <- function(net_data, custom_title) {
  
  consumers <- net_data %>% 
    distinct(Consumer) %>% 
    transmute(name = paste0(Consumer, "_C"), type = TRUE, species = Consumer)
  
  resources <- net_data %>% 
    distinct(Resource) %>% 
    transmute(name = paste0(Resource, "_R"), type = FALSE, species = Resource)
  
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
  
  edges_df <- net_data %>% 
    transmute(
      from = paste0(Consumer, "_C"), 
      to = paste0(Resource, "_R"),
      Median_Prob = Median_Prob
    )
  
  g <- tbl_graph(nodes = nodes_df, edges = edges_df, directed = TRUE)
  layout <- create_layout(g, layout = 'bipartite')
  layout$y <- ifelse(layout$type, 1, 0)
  
  p <- ggraph(layout) +
    geom_edge_link(aes(color = Median_Prob, edge_alpha = Median_Prob, edge_width = Median_Prob)) +
    geom_node_point(aes(fill = Category), shape = 21, size = 5, color = "black", stroke = 0.5) +
    # --- CHANGED: Legend title updated to "Interaction probability"
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

strata_definitions <- tibble(
  Time = c("Pleistocene final", "Holocene mid", "Holocene final", "Holocene final"),
  Region = c("Steppe", "Steppe", "Steppe", "Forest"),
  Custom_Title = c("late Pleistocene - Steppe", 
                   "middle Holocene - Steppe",
                   "late Holocene - Steppe", 
                   "late Holocene - Forest")
)

plot_list <- list()

for (i in 1:nrow(strata_definitions)) {
  curr_net <- median_networks %>% 
    filter(Time == strata_definitions$Time[i], Region == strata_definitions$Region[i])
  
  if (nrow(curr_net) > 0) {
    plot_list[[i]] <- build_bipartite_plot(curr_net, strata_definitions$Custom_Title[i])
  }
}

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

# 1. Topological Aggregation: Supplementary Material
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
# R SCRIPT: CALCULATION OF STRUCTURAL PROPERTIES WITH ABSOLUTE SCALING (E)
# AND INCLUSION OF THE PREDATOR/PREY RATIO + INTERACTIONS
# ----------------------------------------------------------------------------

if (!require("bipartite")) install.packages("bipartite")
library(bipartite)

print("Starting topological calculation with scaling factor (E)...")

# 1. Matrix Transformation and Calculation of Bipartite Metrics
network_metrics_sim <- consensus_variability %>%
  mutate(Int_prob = ifelse(Interaction_prob_sim < 0.01, 0, Interaction_prob_sim)) %>%
  filter(Int_prob > 0) %>%
  group_by(Time, Region, Sim_ID) %>%
  group_modify(~ {
    
    inc_mat <- .x %>%
      select(Resource, Consumer, Int_prob) %>%
      pivot_wider(names_from = Consumer, values_from = Int_prob, values_fill = 0) %>%
      column_to_rownames("Resource") %>%
      as.matrix()
    
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
    
    num_prey <- nrow(inc_mat)
    num_predators <- ncol(inc_mat)
    num_interactions <- sum(inc_mat > 0) 
    
    ratio_predator_prey <- num_predators / num_prey
    
    E_factor <- num_predators * num_prey * 3
    inc_mat_scaled <- inc_mat * E_factor
    
    mets <- tryCatch({
      networklevel(inc_mat_scaled, index = c("connectance", "vulnerability", "generality", 
                                             "weighted nestedness", "Shannon diversity"))
    }, error = function(e) {
      c("connectance"=NA, "vulnerability.LL"=NA, "generality.HL"=NA, 
        "weighted nestedness"=NA, "Shannon diversity"=NA)
    })
    
    fc_prey <- tryCatch(as.numeric(fc(inc_mat_scaled)), error = function(e) NA_real_)
    fc_predators <- tryCatch(as.numeric(fc(t(inc_mat_scaled))), error = function(e) NA_real_)
    
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

# 2. Data Restructuring for Faceting (Tidy Format)
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

plot_data <- metrics_long %>%
  filter(!Metric %in% c("Total Resources", "Total Consumers", "Total Interactions")) %>%
  mutate(Metric = droplevels(Metric)) 

p_network_properties <- ggplot(plot_data, aes(x = Time, y = Value, fill = Region)) +
  geom_boxplot(alpha = 0.85, outlier.size = 0.8, outlier.alpha = 0.4, color = "black") +
  facet_wrap(~ Metric, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Steppe" = "#E3D3B5", "Forest" = "#2D5F43")) +
  # --- CHANGED: Legend title updated to "Habitat"
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

# 4. Exporting Results
output_file <- "results/TierraDelFuego_Network_Properties_Scaled_Extended.pdf"
ggsave(output_file, p_network_properties, width = 18, height = 10, device = "pdf")
ggsave("results/TierraDelFuego_Network_Properties_Scaled_Extended.png", p_network_properties, width=15, height=12, bg="white")

print(paste("Network properties figure successfully generated as:", output_file))

# ----------------------------------------------------------------------------
# 5. EXTRACTION AND EXPORT OF MEDIAN COMPLEXITY METRICS TABLE
# ----------------------------------------------------------------------------
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
# R SCRIPT: HUMAN EGO-NETWORK ANALYSIS (DIET COMPOSITION)
# ----------------------------------------------------------------------------

# Load required libraries
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("patchwork")) install.packages("patchwork")

library(tidyverse)
library(patchwork)

print("Starting Human ego-network analysis (Proportional Diet Composition)...")

# ============================================================================
# PART 0: NECROMASS AGGREGATION
# ============================================================================

# Define the taxa that constitute the marine carrion subsidy
cetacean_taxa <- c("Physeter_macrocephalus", "Balaenoptera_borealis", "Cetacea", "Grampus_griseus", "Phocoena")

# Create a new analytical dataframe where these specific resources are merged
# Note: They are only merged if they act as a Resource.
merged_interactions <- final_validated_interactions %>%
  mutate(
    Resource = ifelse(Resource %in% cetacean_taxa, "Necromass", Resource)
  ) %>%
  # Because multiple merged species might now share the same Consumer, we must aggregate
  # We take the max probability to prevent artificial inflation (e.g., 1.0 + 1.0 = 2.0)
  group_by(Scenario, Time, Region, Consumer, Resource, Interaction_type) %>%
  summarise(
    Median_Probability = max(Median_Probability, na.rm = TRUE),
    Mean_Probability = max(Mean_Probability, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================================
# PART 1: DATA PREPARATION (WITH NECROMASS)
# ============================================================================

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
  # --- CHANGED: Renamed Scenario to "Including Marine Necromass"
  mutate(Scenario = "Including Marine Necromass")

# ============================================================================
# PART 2: DATA PREPARATION (WITHOUT NECROMASS)
# ============================================================================

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
  # --- CHANGED: Renamed Scenario to "Excluding Marine Necromass"
  mutate(Scenario = "Excluding Marine Necromass")

# ============================================================================
# PART 3: FORMATTING FOR VISUALIZATION
# ============================================================================

network_levels <- c(
  "late Pleistocene\nSteppe", 
  "mid Holocene\nSteppe", 
  "late Holocene\nSteppe", 
  "late Holocene\nForest"
)

combined_human_diet <- bind_rows(human_diet_with_necro, human_diet_no_necro) %>%
  mutate(
    # --- CHANGED: Remove underscores from all scientific names globally
    Resource = str_replace_all(Resource, "_", " "),
    
    Period = case_when(
      Time == "Pleistocene final" ~ "late Pleistocene",
      Time == "Holocene mid" ~ "mid Holocene",
      Time == "Holocene final" ~ "late Holocene",
      TRUE ~ Time
    ),
    Period = factor(Period, levels = c("late Pleistocene", "mid Holocene", "late Holocene")),
    
    Network_Label = factor(paste0(Period, "\n", Region), levels = network_levels),
    
    # Adjust logic now that Guanaco is already separated by a space
    Plot_Resource = case_when(
      Resource == "Lama guanicoe" ~ "Lama guanicoe",
      Resource == "Necromass" ~ "Marine Necromass",
      Proportion >= 5.0 ~ Resource, 
      TRUE ~ "Other Resources (<5%)"
    )
  )

output_csv <- "results/TDF_Human_Diet_Proportions.csv"
write_csv(combined_human_diet %>% select(Scenario, Period, Region, Resource, Interaction_Strength, Proportion), output_csv)
print(paste("Human diet proportion dataset saved as:", output_csv))

# ============================================================================
# PART 4: VISUALIZATION (100% STACKED BAR CHART)
# ============================================================================

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

print("Generating Human Diet Composition stacked bar chart...")

p_diet_bar <- ggplot(combined_human_diet, aes(x = Network_Label, y = Proportion, fill = Plot_Resource)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.2) +
  
  scale_y_continuous(labels = scales::percent_format()) +
  # --- CHANGED: Legend title updated to "Human Resource"
  scale_fill_manual(values = base_palette, name = "Human Resource") +
  
  facet_wrap(~ Scenario, ncol = 2) +
  
  # --- CHANGED: X-axis title updated to "Period-Habitat"
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
