# R script for "Ancient consumer-resource networks in Tierra del Fuego"
# Authors: Fernando C. Santiago, Ulises Balza & Tomás I. Marina
# Script #1: Load data


# Load packages -----------------------------------------------------------
library(tidyverse)
library(readxl)


# Load data ---------------------------------------------------------------
## Species localities
communities <- readxl::read_excel("datos/Species localities - TDF.xlsx")
abundance <- as.data.frame(colSums(communities[,2:143]))
names(abundance) <- "Abundance"
abundance$Site <- rownames(abundance)

## Species traits
traits <- readxl::read_excel("datos/Species traits - TDF.xlsx")

## Age localities
ages <- readxl::read_excel("datos/Age localities - TDF.xlsx")


# Save results ------------------------------------------------------------
save(communities, abundance, traits, ages,
     file = "datos/tidy_data.Rdata")

