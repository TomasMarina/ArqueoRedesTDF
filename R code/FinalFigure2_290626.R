# ----------------------------------------------------------------------------
# R SCRIPT: COMPLETE SPATIO-TEMPORAL ASSEMBLAGE 
# ----------------------------------------------------------------------------

# Load standard libraries
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("readxl")) install.packages("readxl")
if (!require("patchwork")) install.packages("patchwork")

library(tidyverse)
library(readxl)
library(patchwork)

print("Starting structuring of the global matrix...")

# 1. DATA LOADING
traits_raw <- read_excel("data - Santiago2025/Species_traits_TDF_270426.xlsx", sheet = "traits")
localities_raw <- read_excel("data - Santiago2025/Species_localities_TDF_270526.xlsx")
ages_raw <- read_excel("data - Santiago2025/Age_localities_TDF_270526.xlsx")

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
    TrophicSpecies, Class, Resolution, Bodymass_min = BodySize_min,
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
    Region, Time, TrophicSpecies, Class, Resolution, 
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
# 7. COMPOUND MANUSCRIPT FIGURE (PIE CHARTS & MAMMALIAN ASSEMBLAGE)
# ----------------------------------------------------------------------------
print("Generating compound figure A (Pie charts) and B (Mammalian bar charts)...")

# ============================================================================
# PANEL A: MAJOR TAXONOMIC GROUPS PIE CHARTS
# ============================================================================

# 7.1 Aggregate Spatiotemporal Relative Importance into the 4 Major Categories
fauna_pies <- assemblage_table %>%
  filter(Time != "Historic") %>%
  mutate(
    Period = case_when(
      Time == "Pleistocene final" ~ "late Pleistocene",
      Time == "Holocene mid" ~ "mid Holocene",
      Time == "Holocene final" ~ "late Holocene",
      TRUE ~ Time
    ),
    Period = factor(Period, levels = c("late Pleistocene", "mid Holocene", "late Holocene")),
    Block = paste(Period, "-", Region),
    Block = factor(Block, levels = c(
      "late Pleistocene - Steppe",
      "mid Holocene - Steppe",
      "late Holocene - Steppe",
      "late Holocene - Forest"
    )),
    
    # Consolidate taxonomic classes into requested major groups
    Major_Category = case_when(
      Class %in% c("Actinopterygii", "Chondrichthyes", "Teleostei") ~ "Fishes",
      Class == "Aves" ~ "Aves",
      Class %in% c("Bivalvia", "Echinoidea", "Gastropoda", "Polyplacophora", "Thecostraca") ~ "Invertebrates",
      Class == "Mammalia" ~ "Mammalia",
      TRUE ~ "Other"
    )
  ) %>%
  # Group by the spatial-temporal block and sum the relative importance
  group_by(Block, Major_Category) %>%
  summarise(Total_RelImp = sum(Spatiotemporal_RelImp_Pct, na.rm = TRUE), .groups = "drop") %>%
  # Ensure clean ordering in the legend
  mutate(Major_Category = factor(Major_Category, levels = c("Mammalia", "Fishes", "Invertebrates", "Aves", "Other")))

# Set color palette for major categories (updated to prevent clashing with Forest/Steppe colors)
category_colors <- c(
  "Mammalia" = "#D55E00",       # Red-Orange
  "Fishes" = "#56B4E9",        # Light Blue
  "Invertebrates" = "#CC79A7", # Reddish-Purple (distinct from the new Forest green)
  "Aves" = "#0072B2",          # Dark Blue (distinct from the new Steppe beige)
  "Other" = "grey50"
)

# Render Panel A
p_pies <- ggplot(fauna_pies, aes(x = "", y = Total_RelImp, fill = Major_Category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  facet_wrap(~ Block, ncol = 4) +
  scale_fill_manual(values = category_colors, name = "Taxonomic group") +
  theme_void(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 12, margin = margin(b = 10)),
    plot.margin = margin(10, 10, 20, 10)
  )


# ============================================================================
# PANEL B: MAMMALIAN ASSEMBLAGE SPATIO-TEMPORAL VARIATION
# ============================================================================

# 7.2 Construct the summary dataset for mammals
mammal_summary <- tibble(
  Period = rep(c("late Pleistocene", "mid Holocene", "late Holocene"), each = 8),
  Habitat = rep(rep(c("Steppe", "Forest"), each = 4), 3),
  Taxon = rep(c("Lama guanicoe", "Rodents", "Other terrestrial", "Marine mammals"), 6),
  Percentage = c(
    # late Pleistocene
    29.8, 31.6, 38.6, 0.0,   # Steppe
    NA,   NA,   NA,   NA,    # Forest (No Network)
    # mid Holocene
    51.1, 32.5, 6.8,  9.6,   # Steppe
    NA,   NA,   NA,   NA,    # Forest (No Network)
    # late Holocene
    57.9, 37.2, 1.3,  3.6,   # Steppe
    52.8, 0.0,  0.6,  46.6   # Forest
  )
) %>%
  mutate(
    Period = factor(Period, levels = c("late Pleistocene", "mid Holocene", "late Holocene")),
    Habitat = factor(Habitat, levels = c("Steppe", "Forest")),
    # Using plotmath syntax to properly italicize the scientific name in facet strips
    Taxon_Label = case_when(
      Taxon == "Lama guanicoe" ~ "italic('Lama guanicoe')",
      TRUE ~ paste0("'", Taxon, "'")
    ),
    Taxon_Label = factor(Taxon_Label, levels = c("italic('Lama guanicoe')", "'Rodents'", "'Other terrestrial'", "'Marine mammals'")),
    Proportion = Percentage / 100
  )

# Render Panel B
p_mammals <- ggplot(mammal_summary, aes(x = Period, y = Proportion, fill = Habitat)) +
  # NA padding naturally forces Steppe bars to the left when Forest is absent
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", alpha = 0.85, width = 0.7) +
  facet_wrap(~ Taxon_Label, scales = "fixed", ncol = 4, labeller = label_parsed) +
  # Applying the specific hex codes requested for Steppe and Forest
  scale_fill_manual(values = c("Steppe" = "#E3D3B5", "Forest" = "#2D5F43")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.7)) +
  labs(
    x = "Period",
    y = "Relative Abundance (% NISP)",
    fill = "Habitat"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(10, 10, 10, 10)
  )


# ============================================================================
# PANEL C: COMPOUND ASSEMBLY & EXPORT
# ============================================================================

# Assemble using patchwork. Heights (1, 1.5) give the bar charts slightly more vertical space.
compound_fig <- p_pies / p_mammals +
  plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(1, 1.5)) 

output_pdf_compound <- "results/TierraDelFuego_Compound_Assemblage.pdf"
ggsave(output_pdf_compound, compound_fig, width = 14, height = 11, device = "pdf")
ggsave("results/TierraDelFuego_Compound_Assemblage.png", compound_fig, width = 14, height = 11, bg = "white")

print(paste("Compound figure successfully generated as:", output_pdf_compound))
