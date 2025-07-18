# R script for "Ancient consumer-resource networks in Tierra del Fuego"
# Authors: Fernando C. Santiago, Ulises Balza & Tomás I. Marina
# Script #1: Load & tidy data

# Load packages -----------------------------------------------------------

library(tidyverse)
library(readxl)


# Load data ---------------------------------------------------------------
## Species localities ----
communities <- readxl::read_excel("data - Santiago2025/Species localities - TDF.xlsx")

# Convert all species columns to numeric to prevent errors
communities <- communities %>%
  mutate(across(2:65, as.numeric))

# Calculate abundance for each site
abundance <- as.data.frame(colSums(communities[,2:65], na.rm = TRUE))
names(abundance) <- "Abundance"
abundance$Site <- colnames(communities[,2:65])

## Species traits ----
traits <- readxl::read_excel("data - Santiago2025/Species traits - TDF.xlsx")

## Age localities ----
ages <- readxl::read_excel("data - Santiago2025/Age localities - TDF.xlsx")


# Save results ------------------------------------------------------------

save(communities, abundance, traits, ages,
     file = "data - Santiago2025/tidy_data_180725.Rdata")
