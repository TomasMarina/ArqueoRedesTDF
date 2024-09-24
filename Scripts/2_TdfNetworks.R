# R script for "Ancient consumer-resource networks in Tierra del Fuego"
# Authors: Ulises Balza, Fernando Santiago & Tomás I. Marina
# Script #2: Inferring consumer-resource networks


# Load packages -----------------------------------------------------------
library(tidyverse)


# Load data ---------------------------------------------------------------
load("datos/tidy_data.Rdata")


# Community data frames ---------------------------------------------------
## Discard uncompleted localities
comm_full <- communities %>% 
  dplyr::select(where(~!any(is.na(.))))

## Classify species by feeding strategy
carn <- traits[traits$FeedingStrategy == "Carnivore",]
herb <- traits[traits$FeedingStrategy == "Herbivore",]

carn_species <- traits[traits$FeedingStrategy == "Carnivore",]$TrophicSpecies
herb_species <- traits[traits$FeedingStrategy == "Herbivore",]$TrophicSpecies

## Create community data frames by locality
occ_carn = communities[communities$TrophicSpecies %in% carn_species,]
occ_herb = communities[communities$TrophicSpecies %in% herb_species,]

communities_carn = list()
for (i in 1:142){ #nrow(communities)-1
  k = occ_carn[occ_carn[,i+1] > 1,]
  communities_carn[[i]] = k$TrophicSpecies
}
communities_carn
names(communities_carn) = names(occ_carn[2:143])

communities_herb = list()
for (i in 1:142){
  k = occ_herb[occ_herb[,i+1] > 1,]
  communities_herb[[i]] = k$TrophicSpecies
}
communities_herb
names(communities_herb) = names(occ_herb[2:143])

## Include trait data in community by locality
herb_traits = list()  # herbivores
for (i in 1:length(communities_herb)){
  k = as.data.frame(communities_herb[[i]])
  names(k)="TrophicSpecies"
  data = merge(k, traits)
  data = data[c("TrophicSpecies","BodySize_min","BodySize_max","FeedingStrategy")]
  names(data) = c("TrophicSpecies","BodySize","FeedingStrategy")
  herb_traits[[i]] = data
}
names(herb_traits) = names(communities_herb)

carn_traits = list()  # carnivores
for (i in 1:length(communities_carn)){
  k = as.data.frame(communities_carn[[i]])
  names(k) = "TrophicSpecies"
  data = merge(k, traits)
  data = data[c("TrophicSpecies","BodySize_min","BodySize_max","FeedingStrategy")]
  #data = data[data$locomotion != "aquatic",]
  #data = data[data$body.size > 10, ]
  names(data) = c("TrophicSpecies","BodySize_min","BodySize_max","FeedingStrategy")
  carn_traits[[i]] = data
}
names(carn_traits) = names(communities_carn)


# Network model -----------------------------------------------------------
## Prey/Predator body size ratio function








