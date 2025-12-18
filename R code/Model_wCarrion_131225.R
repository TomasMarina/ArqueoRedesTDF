# ----------------------------------------------------------------------------
#
# R SCRIPT FOR RECONSTRUCTING ANCIENT CONSUMER-RESOURCE NETWORKS
# (Tierra del Fuego Archaeological Study)
#
# Version: 171225-Final-AbundanceWeighted
#
# *** SCRIPT CORRECTIONS & UPDATES ***
# 1. PART 3 (DATA PREP):
#    - Now calculates 'nisp' for each species within each aggregated bin.
#    - Sums NISP if a species appears in multiple localities within the same bin.
# 2. PART 5 (reconstruct_network):
#    - NOW INCORPORATES ABUNDANCE (NISP).
#    - Methodology: "Mass Action" Bias on Logit.
#      a. Calculate 'abund_score' = log10(NISP_consumer * NISP_resource).
#      b. Normalize 'abund_score' to 0-1 range within the current network.
#      c. Add this score to the logit (z) of the body-mass model.
#      d. High abundance pairs get a boost in probability; low abundance pairs stay similar.
#    - Humans (Homo_sapiens) are assigned the MAX NISP of the network to represent
#      their ubiquity/agency, preventing them from having low interaction probs due to absence in bone counts.
#
# ----------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# PART 1: SETUP
# ----------------------------------------------------------------------------

if (!require("tidyverse")) install.packages("tidyverse")
if (!require("readxl")) install.packages("readxl")
if (!require("janitor")) install.packages("janitor")
if (!require("bipartite")) install.packages("bipartite")
if (!require("igraph")) install.packages("igraph")
if (!require("tidygraph")) install.packages("tidygraph")
if (!require("ggraph")) install.packages("ggraph")
if (!require("ggrepel")) install.packages("ggrepel")
if (!require("patchwork")) install.packages("patchwork")
if (!require("viridis")) install.packages("viridis")
if (!require("ggnewscale")) install.packages("ggnewscale")
if (!require("tidytext")) install.packages("tidytext")

library(tidyverse)
library(readxl)
library(janitor)
library(bipartite)
library(igraph)
library(tidygraph)
library(ggraph)
library(ggrepel)
library(patchwork)
library(viridis)
library(ggnewscale)
library(tidytext)

set.seed(123)
Sys.setlocale("LC_TIME", "en_US.UTF-8")
today_date <- format(Sys.Date(), "%d%m%y")


# ----------------------------------------------------------------------------
# PART 2: DATA LOADING
# ----------------------------------------------------------------------------

print("Loading data...")
localities_raw <- read_excel("data - Santiago2025/Species localities - TDF  - 131225.xlsx")
traits_raw <- read_excel("data - Santiago2025/Species traits - TDF - 131225.xlsx", sheet = "traits")
ages_raw <- read_excel("data - Santiago2025/Age localities - TDF  - 081125.xlsx")


# ----------------------------------------------------------------------------
# PART 3: DATA PREPARATION
# ----------------------------------------------------------------------------

print("Preparing data...")

# --- Clean Traits ---
traits <- traits_raw %>%
  clean_names() %>%
  mutate(
    facultative_scavenger = ifelse(is.na(as.numeric(facultative_scavenger)), 0, as.numeric(facultative_scavenger)),
    carrion_only = ifelse(is.na(as.numeric(carrion_only)), 0, as.numeric(carrion_only)),
    human_resource = ifelse(is.na(as.numeric(human_resource)), 0, as.numeric(human_resource))
  ) %>%
  mutate(
    avg_bodysize_kg = exp((log(body_size_min) + log(body_size_max)) / 2),
    log10_bodysize = log10(avg_bodysize_kg)
  ) %>%
  select(
    trophic_species, class, order, family,
    log10_bodysize, feeding_strategy,
    feeding_habitat, facultative_scavenger, carrion_only, human_resource
  ) %>%
  distinct(trophic_species, .keep_all = TRUE)

# --- Clean Localities & Merge Ages ---
# We now keep NISP to calculate abundance at the network level
localities <- localities_raw %>%
  clean_names() %>%
  pivot_longer(
    cols = -trophic_species,
    names_to = "locality",
    values_to = "nisp_raw",
    values_transform = list(nisp_raw = as.character)
  ) %>%
  mutate(
    nisp = as.numeric(str_trim(nisp_raw)),
    locality = str_trim(locality)
  ) %>%
  filter(!is.na(nisp)) %>% # Keep rare species (NISP > 0 handled implicitly or can add & nisp > 0)
  inner_join(traits %>% select(trophic_species), by = "trophic_species")

ages <- ages_raw %>%
  clean_names() %>%
  mutate(locality = tolower(str_trim(locality))) %>%
  select(locality, time_bin = geological_ages, region_col) %>%
  distinct()

# --- Create Aggregated Networks with ABUNDANCE ---
splocs_agg <- localities %>%
  inner_join(ages, by = "locality") %>%
  group_by(time_bin, region_col) %>%
  summarise(
    # Create a nested dataframe of species AND their total NISP in this bin
    species_data = list(
      cur_data() %>% 
        group_by(trophic_species) %>% 
        summarise(total_nisp = sum(nisp, na.rm = TRUE), .groups = "drop")
    ),
    .groups = "drop"
  ) %>%
  filter(!is.na(time_bin), !is.na(region_col)) %>%
  # *** FORCE HUMANS INTO EVERY NETWORK ***
  # We assign Humans the MAXIMUM NISP found in that network to represent high presence
  mutate(
    species_data = map(species_data, function(df) {
      if (!"Homo_sapiens" %in% df$trophic_species) {
        max_nisp <- if(nrow(df) > 0) max(df$total_nisp, na.rm=TRUE) else 100
        bind_rows(df, tibble(trophic_species = "Homo_sapiens", total_nisp = max_nisp))
      } else {
        df
      }
    })
  )

print(paste(nrow(splocs_agg), "aggregated networks prepared."))


# ----------------------------------------------------------------------------
# PART 4: DEFINE MODEL PARAMETERS
# ----------------------------------------------------------------------------

get_stochastic_params <- function() {
  list(
    alphaP = runif(1, 1, 2), betaP = runif(1, -2, -1), gammaP = runif(1, -1, -0.5),
    alphaO = runif(1, -6, -4), betaO = runif(1, -3, -2), gammaO = 0
  )
}
mean_params <- list(
  alphaP = mean(c(1, 2)), betaP = mean(c(-2, -1)), gammaP = mean(c(-1, -0.5)),
  alphaO = mean(c(-6, -4)), betaO = mean(c(-3, -2)), gammaO = 0
)


# ----------------------------------------------------------------------------
# PART 5: RECONSTRUCTION FUNCTION (with Abundance)
# ----------------------------------------------------------------------------

reconstruct_network <- function(sp_data, traits_df, model_params, focal_habitat = "Ter-Mar") {
  
  # sp_data contains 'trophic_species' and 'total_nisp'
  # traits_df is the global traits file
  
  # Join traits to the local species list
  local_sp <- sp_data %>%
    left_join(traits_df, by = "trophic_species")
  
  # 1. Pools
  consumers <- local_sp %>% filter(feeding_strategy %in% c("Carnivore", "Omnivore") & carrion_only == 0)
  resources <- local_sp %>% filter((feeding_strategy %in% c("Herbivore", "Omnivore") | carrion_only == 1) & trophic_species != "Homo_sapiens")
  
  if (focal_habitat != "Ter-Mar") {
    consumers <- consumers %>% filter(feeding_habitat %in% c(focal_habitat, "Ter-Mar"))
    resources <- resources %>% filter(feeding_habitat %in% c(focal_habitat, "Ter-Mar"))
  }
  
  all_pairs <- expand_grid(consumer_species = consumers$trophic_species, 
                           resource_species = resources$trophic_species) %>%
    filter(consumer_species != resource_species) %>%
    # Join Consumer Traits & NISP
    left_join(consumers %>% select(consumer_species = trophic_species, consumer = feeding_strategy, log10_c = log10_bodysize, habitat_c = feeding_habitat, scav_c = facultative_scavenger, nisp_c = total_nisp), by = "consumer_species") %>%
    # Join Resource Traits & NISP
    left_join(resources %>% select(resource_species = trophic_species, resource = feeding_strategy, log10_r = log10_bodysize, habitat_r = feeding_habitat, carrion_r = carrion_only, human_r = human_resource, nisp_r = total_nisp), by = "resource_species")
  
  if (nrow(all_pairs) == 0) return(tibble(Consumer=character(), Resource=character(), Int_prob=numeric()))
  
  # --- CALCULATE ABUNDANCE SCORE ---
  # We define abundance impact as log(NISP_c * NISP_r). 
  # This creates a "mass action" score.
  all_pairs <- all_pairs %>%
    mutate(
      mass_action = log10(nisp_c * nisp_r + 1), # +1 to avoid log(0)
      # Normalize to 0-1 range within this specific network
      abund_score = (mass_action - min(mass_action)) / (max(mass_action) - min(mass_action) + 1e-6)
    )
  
  # Calculate Probabilities
  results_list <- list()
  for (i in 1:nrow(all_pairs)) {
    pair <- all_pairs[i, ]
    
    habitat_overlap <- if_else(pair$habitat_c == pair$habitat_r | pair$habitat_c == "Ter-Mar" | pair$habitat_r == "Ter-Mar", 1.0, 0.0)
    
    if (habitat_overlap == 0) {
      results_list[[i]] <- tibble(Consumer = pair$consumer_species, Resource = pair$resource_species, Int_prob = 0.0)
      next
    }
    
    prob <- 0.0
    z <- NA # Logit
    
    # Rule 1: Human Scavenging (Always 1.0)
    if (pair$consumer_species == "Homo_sapiens" && pair$scav_c == 1 && pair$carrion_r == 1) {
      prob <- 1.0 
      
      # Rule 2: Human Predation
    } else if (pair$consumer_species == "Homo_sapiens" && pair$human_r > 0) {
      if (!is.na(pair$log10_c) && !is.na(pair$log10_r)) {
        log_ratio <- pair$log10_c - pair$log10_r
        z <- model_params$alphaO + model_params$betaO * log_ratio + model_params$gammaO * (log_ratio^2)
      }
      
      # Rule 3: Standard Predation
    } else if (!is.na(pair$log10_c) && !is.na(pair$log10_r) && pair$log10_c > pair$log10_r) {
      log_ratio <- pair$log10_c - pair$log10_r
      if (pair$consumer == "Carnivore") {
        z <- model_params$alphaP + model_params$betaP * log_ratio + model_params$gammaP * (log_ratio^2)
      } else {
        z <- model_params$alphaO + model_params$betaO * log_ratio + model_params$gammaO * (log_ratio^2)
      }
    }
    
    # --- APPLY ABUNDANCE BOOST ---
    # If a valid 'z' was calculated, we boost it by the abundance score.
    # We add (abund_score * 2.0) to z. This shifts the logit.
    # E.g., if prob was 0.1 (z=-2.2), and abund_score is 1.0 (max), 
    # z becomes -0.2 (prob ~0.45). Significant but not overwhelming.
    if (!is.na(z)) {
      z_final <- z + (pair$abund_score * 2.0) 
      prob <- exp(z_final) / (1 + exp(z_final))
    }
    
    results_list[[i]] <- tibble(Consumer = pair$consumer_species, Resource = pair$resource_species, Int_prob = prob)
  }
  
  bind_rows(results_list) %>% filter(Int_prob > 0)
}


# ----------------------------------------------------------------------------
# PART 6: METRICS
# ----------------------------------------------------------------------------
calculate_metrics <- function(links_df) {
  inc_mat <- links_df %>%
    group_by(Consumer, Resource) %>% summarise(Int_prob = max(Int_prob), .groups = "drop") %>%
    pivot_wider(names_from = Resource, values_from = Int_prob, values_fill = 0) %>%
    column_to_rownames("Consumer") %>% as.matrix()
  
  if (nrow(inc_mat) < 2 || ncol(inc_mat) < 2) return(tibble(n_consumers=nrow(inc_mat), n_resources=ncol(inc_mat), connectance=NA, generality=NA, vulnerability=NA, modularity=NA))
  
  t_mat <- t(inc_mat)
  mets <- tryCatch(networklevel(t_mat, index = c("connectance", "links per species", "weighted nestedness", "generality", "vulnerability")), error = function(e) NULL)
  mods <- tryCatch(computeModules(t_mat), error = function(e) NULL)
  mod_q <- if (!is.null(mods)) mods@likelihood else NA
  
  if (is.null(mets)) return(tibble(n_consumers=nrow(inc_mat), n_resources=ncol(inc_mat), connectance=NA, generality=NA, vulnerability=NA, modularity=mod_q))
  
  tibble(n_consumers = nrow(inc_mat), n_resources = ncol(inc_mat), connectance = mets["connectance"], links_per_species = mets["links per species"], weighted_nestedness = mets["weighted nestedness"], modularity = mod_q, generality = mets["generality.HL"], vulnerability = mets["vulnerability.LL"])
}


# ----------------------------------------------------------------------------
# PART 7 & 8: RUN & SAVE
# ----------------------------------------------------------------------------
run_analysis <- function(splocs_table, traits) {
  model_params <- get_stochastic_params()
  map_dfr(1:nrow(splocs_table), function(i) {
    # Note: sp_traits now passed as the 'species_data' tibble which includes NISP
    net_links <- reconstruct_network(splocs_table$species_data[[i]], traits, model_params, "Ter-Mar")
    metrics <- calculate_metrics(net_links)
    metrics %>% mutate(time_bin = splocs_table$time_bin[i], region_col = splocs_table$region_col[i])
  })
}

print("Running 100 replicates...")
final_metrics <- map_dfr(1:10, ~run_analysis(splocs_agg, traits), .id = "rep")

if (!dir.exists("results")) dir.create("results")
write_csv(final_metrics, file.path("results", paste0(today_date, "_TDF_metrics.csv")))

# Representative Links
plot_data <- splocs_agg %>%
  mutate(
    links_df = map(species_data, function(sp_dat) {
      links <- reconstruct_network(sp_dat, traits, mean_params)
      links %>% left_join(traits %>% select(trophic_species, carrion_r = carrion_only), by = c("Resource"="trophic_species")) %>%
        mutate(human_link_type = case_when(
          Consumer == "Homo_sapiens" & carrion_r == 1 ~ "Human-Scavenging",
          Consumer == "Homo_sapiens" ~ "Human-Predation",
          TRUE ~ "Non-Human"
        ))
    })
  )
all_links <- plot_data %>% select(time_bin, region_col, links_df) %>% unnest(links_df)
write_csv(all_links, file.path("results", paste0(today_date, "_TDF_interaction_list.csv")))


# ----------------------------------------------------------------------------
# PART 10: NETWORK PLOTS
# ----------------------------------------------------------------------------
shorten_name <- function(species_name) {
  parts <- str_split(species_name, "_")[[1]]
  if (length(parts) == 2) return(paste0(str_sub(parts[1], 1, 1), ". ", parts[2])) else return(species_name)
}
taxa_colors <- c("Homo_sapiens"="#33CCFF", "Mammalia"="#7f3b08", "Aves"="#b35806", "Actinopterygii"="#e08214", "Chondrichthyes"="#fdb863", "Malacostraca"="#fee0b6", "Bivalvia"="#d8daeb", "Gastropoda"="#b2abd2", "Cephalopoda"="#8073ac", "Echinoidea"="#542788", "Other"="grey50")

pwalk(plot_data, function(species_data, links_df, time_bin, region_col, ...) {
  if (nrow(links_df) == 0) return()
  
  nodes_df <- bind_rows(
    links_df %>% distinct(Consumer) %>% transmute(name=paste0(Consumer, "_C"), short_name=map_chr(Consumer, shorten_name), type=TRUE, species=Consumer),
    links_df %>% distinct(Resource) %>% transmute(name=paste0(Resource, "_R"), short_name=map_chr(Resource, shorten_name), type=FALSE, species=Resource)
  ) %>%
    left_join(traits %>% select(trophic_species, class), by=c("species"="trophic_species")) %>%
    mutate(Category = ifelse(species=="Homo_sapiens", "Homo_sapiens", ifelse(class %in% names(taxa_colors), class, "Other")),
           is_human_node = species=="Homo_sapiens") %>%
    mutate(from = name) # Dummy for filter
  
  links_plot <- links_df %>% mutate(from=paste0(Consumer, "_C"), to=paste0(Resource, "_R"))
  
  g <- tbl_graph(nodes=nodes_df, edges=links_plot, directed=TRUE)
  layout <- create_layout(g, layout='bipartite')
  layout$y <- ifelse(layout$type, 1, 0)
  
  p <- ggraph(layout) +
    geom_edge_link(aes(alpha=Int_prob, color=Int_prob, filter=human_link_type=="Non-Human"), width=0.5) +
    scale_edge_color_gradient(low="#3B4CC0", high="#B40426", name="Prob") +
    geom_edge_link(aes(alpha=Int_prob, filter=human_link_type=="Human-Predation"), color="red", width=1.0) +
    geom_edge_link(aes(alpha=Int_prob, filter=human_link_type=="Human-Scavenging"), color="#8B4513", width=1.0) +
    geom_node_point(aes(fill=Category), shape=21, color="black", size=6) +
    geom_node_text(aes(label=short_name), repel=TRUE, bg.colour="white") +
    scale_fill_manual(values=taxa_colors) +
    theme_graph(base_family='sans') + labs(title=paste(time_bin, "-", region_col))
  
  ggsave(file.path("results", paste0(today_date, "_Network_", time_bin, "_", region_col, ".png")), p, width=12, height=10, bg="white")
})


# ----------------------------------------------------------------------------
# PART 12: BIPLOT (RELATIVE)
# ----------------------------------------------------------------------------
prob_summary <- all_links %>%
  pivot_longer(cols = c(Consumer, Resource), values_to = "Species") %>%
  group_by(time_bin, region_col) %>% mutate(Network_Total = sum(Int_prob)) %>% ungroup() %>%
  group_by(time_bin, region_col, Species, Network_Total) %>% summarise(Sp_Total = sum(Int_prob), .groups="drop") %>%
  mutate(Pct = (Sp_Total/Network_Total)*100, Network = paste0(time_bin, "\n", region_col),
         Short_Name = map_chr(Species, shorten_name), Is_Human = Species == "Homo_sapiens") %>%
  group_by(Network) %>% mutate(Short_Name = reorder_within(Short_Name, Pct, Network))

p_biplot <- ggplot(prob_summary, aes(x=Short_Name, y=Pct, fill=Is_Human)) +
  geom_bar(stat="identity") + facet_wrap(~Network, scales="free") + scale_x_reordered() + coord_flip() +
  scale_fill_manual(values=c("FALSE"="steelblue", "TRUE"="#33CCFF")) +
  theme_minimal() + labs(y="% Total Prob", x="Species", title="Relative Importance")

ggsave(file.path("results", paste0(today_date, "_TDF_biplot_relative.png")), p_biplot, width=15, height=12, bg="white")

print("Finished.")