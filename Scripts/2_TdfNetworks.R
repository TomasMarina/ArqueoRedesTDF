# R script for "Ancient consumer-resource networks in Tierra del Fuego"
# Authors: Ulises Balza, Fernando Santiago & Tomás I. Marina
# Script #2: Inferring consumer-resource networks


# Load packages -----------------------------------------------------------
library(tidyverse)
library(igraph)
library(ggplot2)
library(bipartite)
require(ggpubr)

# Load data ---------------------------------------------------------------
load("datos/tidy_data.Rdata")


# Community data frames ---------------------------------------------------
## Discard uncompleted localities
comm_full <- communities %>% 
  dplyr::select(where(~!any(is.na(.))))

## Classify species by feeding strategy
carn <- traits[traits$FeedingStrategy == "Carnivore",]
herb <- traits[traits$FeedingStrategy != "Carnivore",]

carn_species <- traits[traits$FeedingStrategy == "Carnivore",]$TrophicSpecies
herb_species <- traits[traits$FeedingStrategy != "Carnivore",]$TrophicSpecies

## Create community data frames by locality
occ_carn <- communities[communities$TrophicSpecies %in% carn_species,]
occ_herb <- communities[communities$TrophicSpecies %in% herb_species,]

communities_carn = list()
for (i in 1:142){  #length(communities)-1
  k = occ_carn[occ_carn[,i+1] > 0,]
  communities_carn[[i]] = k$TrophicSpecies
}
communities_carn
names(communities_carn) = names(occ_carn[2:143])

communities_herb = list()
for (i in 1:142){
  k = occ_herb[occ_herb[,i+1] > 0,]
  communities_herb[[i]] = k$TrophicSpecies
}
communities_herb
names(communities_herb) = names(communities_carn)  #names(occ_herb[2:143])

## Include trait data in community by locality
herb_traits = list()  # herbivores
for (i in 1:length(communities_herb)){
  k = as.data.frame(communities_herb[[i]])
  names(k)="TrophicSpecies"
  data = merge(k, traits)
  data = data[c("TrophicSpecies","BodySize_max","FeedingStrategy")]
  #names(data) = c("TrophicSpecies","BodySize_max","FeedingStrategy")
  herb_traits[[i]] = data
}
names(herb_traits) = names(communities_herb)

carn_traits = list()  # carnivores
for (i in 1:length(communities_carn)){
  k = as.data.frame(communities_carn[[i]])
  names(k) = "TrophicSpecies"
  data = merge(k, traits)
  data = data[c("TrophicSpecies","BodySize_max","FeedingStrategy")]
  #data = data[data$locomotion != "aquatic",]
  #data = data[data$body.size > 10, ]
  #names(data) = c("TrophicSpecies","BodySize_max","FeedingStrategy")
  carn_traits[[i]] = data
}
names(carn_traits) = names(communities_carn)


# Network model -----------------------------------------------------------
## Prey/Predator body size (in kg) ratio function. Log-Ratio Model
LRM <- function(mass_C = NULL, mass_H = NULL){ 
  m <- nrow(mass_C)
  n <- nrow (mass_H)
  P <- matrix(NA, m, n)
  
  # sample model parameters from parameter ranges 
  # carnivores
  alphaP<-runif(1, min = 1, max = 2)
  betaP<-runif(1, min = -2, max = -1)
  gammaP<-runif(1, min =-1, max = -0.5)
  # omnivores
  alphaO <- runif(1, min = -6 , max = -4)
  betaO <- runif(1, min = -3, max = -2)
  gammaO <- 0  #runif(1, min = -0.5, max = -0.2)
  
  for (i in 1:m){
    
    if(mass_C$FeedingStrategy[i]!="Omnivore"){
      alpha <- alphaP
      beta <- betaP
      gamma <- gammaP
    }else{
      alpha <- alphaO
      beta <- betaO
      gamma <- gammaO
    }
    
    for(j in 1:n){
      
      term <- (exp(alpha+(beta*log((mass_H$BodySize_max[j]/mass_C$BodySize_max[i]))) + (gamma*(log((mass_H$BodySize_max[j]/mass_C$BodySize_max[i]))^2))))
      p <- term/(1+term)  #probability of interaction
      P[i,j] <- p  #filling probability matrix	
      #if(mass_H$locomotion[j]=="arboreal"){P[i,j]=P[i,j]*0.1}  #reducing interaction probability btw arboreal vs terrestrial species
      #if(mass_C$trophic[i]=="carnivore_invert" & mass_H$trophic[j]!="invertebrate"){P[i,j]=0}  #not used in main analyses  
      #if(mass_H$trophic[j]=="invertebrate" & mass_C$trophic[i]!="carnivore_invert"){P[i,j]=0}  #not used in main analyses 
      
    }
    if(mass_C$FeedingStrategy[i]=="Omnnivore"){P[i,]=P[i,]*0.5}  #rescaling predation probability for omnivores
    
  } 
  
  row.names(P) <- as.character(mass_C$TrophicSpecies)
  colnames(P) <- as.character(mass_H$TrophicSpecies)
  return(P)
}

## Estimate interaction probability
### Tidy needed: select completed localities
### Filter spp with abundance > 0 & feeding strategy data
comm_full_pivot <- comm_full %>% 
  tidyr::pivot_longer(!TrophicSpecies, names_to = "Locality", values_to = "Abundance") %>% 
  dplyr::filter(Abundance > 0)
comm_spp_trait <- comm_full_pivot %>% 
  dplyr::left_join(traits) %>% 
  drop_na(FeedingStrategy, BodySize_max)

### Filter localities with spp abundance > 0
loc_comp <- unique(comm_spp_trait$Locality)
#loc_comp <- "RC1_6284"
carn_comp <- carn_traits[loc_comp]
herb_comp <- herb_traits[loc_comp]

results = list()
replicas = list()
for (i in 1:length(loc_comp)){  #running for all localities
  if (nrow(carn_comp[[i]]) > 0 & nrow(herb_comp[[i]]) > 0){
    if (nrow(carn_comp[[i]]) > 1 | nrow(herb_comp[[i]]) > 1){
      N = 100
      P.list <- replicate(N, LRM(carn_comp[[i]], herb_comp[[i]]))
      P.list[P.list < 0.05] <- 0
      replicas[[i]] = P.list
      results[[i]] = apply(simplify2array(P.list), 1:2, mean)
      results[[i]] = round(results[[i]], 5)
      
    }
  }
}
names(results) = names(carn_comp)
# Shows interaction matrix, predators in row, prey in col
results[["CN26_4000"]]

## Interaction list ----
### Tidy data
results_df <- lapply(results, function(x) as.data.frame(x))
my_list <- Map(cbind, results_df, Locality = names(results_df))
int_df3 <- my_list[["CN29_4700"]] %>% 
  tibble::rownames_to_column(var = "Predator") %>% 
  tidyr::pivot_longer(cols = 2:length(my_list[["CN29_4700"]]), names_to = "Prey", values_to = "Int_prob") %>%
  dplyr::select(Prey, Predator, Int_prob, Locality)
df_list <- list(int_df, int_df2, int_df3)

all_int <- bind_rows(df_list)


# Network properties ------------------------------------------------------
## Predator-prey ratio ----
prop <- matrix(ncol=2,nrow=length(loc_comp))
for (i in 1:length(loc_comp)){
  prop[i,1]=names(results[i])
  if (!is.null(results[[i]])){
    prop[i,2]=nrow(results[[i]])/ncol(results[[i]])
  }
}

prop <- as.data.frame(prop)
names(prop) <- c("Locality", "PPR")
### Include time data
# prop = merge(prop,ages)
# prop$age = (prop$MIN.age+prop$MAX.age)/2
# prop = merge(prop,riq)
# prop = merge(prop,modules, all.x=T)

## Connectance ----
conec <- matrix(nrow=length(loc_comp), ncol=2)
for (i in 1:length(loc_comp)){
  conec[i,1] = names(results[i])
  if (!is.null(results[[i]])){
    f = as.matrix(results[[i]])
    if (nrow(f)>1&ncol(f)>1){
      conec[i,2] = conec[i,2]=sum(results[[i]]/(ncol(results[[i]])*nrow(results[[i]])))
    }
  }
}

conec <- as.data.frame(conec)
names(conec) <- c("Locality", "Connectance")
conec$Connectance <- as.numeric(conec$Connectance)
### Include time data
# conec = merge(conec,ages)
# conec$age = (conec$MIN.age+conec$MAX.age)/2
# conec = merge(conec,riq)
# conec = merge(conec,modules, all.x=T)

## Interaction diversity ----
div <- matrix(nrow=length(loc_comp), ncol=2)
for (i in 1:length(loc_comp)){
  div[i,1] = names(results[i])
  if (!is.null(results[[i]])){
    f = as.matrix(results[[i]])
    if (nrow(f)>1&ncol(f)>1){
      
      div[i,2] = diversity(as.numeric(results[[i]],index="shannon"))
    }
  }
}

div <- as.data.frame(div)
names(div) <- c("Locality", "Diversity")
div$Diversity <- as.numeric(div$Diversity)
### Include time data
# div=merge(div,ages)
# div$age=(div$MIN.age+div$MAX.age)/2
# div=merge(div,riq)
# div=merge(div,modules, all.x=T)

## Mean degree ----
mean_degree_carn <- matrix(ncol=2,nrow=length(loc_comp))
for (i in 1:length(loc_comp)){
  mean_degree_carn[i,1]=names(results)[i]
  if (!is.null(results[[i]])){
    mean_degree_carn[i,2]=mean(rowSums(results[[i]]))}
}
mean_degree_carn <- as.data.frame(mean_degree_carn)
names(mean_degree_carn) <- c("Locality", "mean_degree_carn")
mean_degree_carn$mean_degree_carn <- as.numeric(mean_degree_carn$mean_degree_carn)
### Include time data
# mean_degree_carn = merge(mean_degree_carn,ages)
# mean_degree_carn$age = (mean_degree_carn$MIN.age+mean_degree_carn$MAX.age)/2
# mean_degree_carn = merge(mean_degree_carn,riq)
# mean_degree_carn = merge(mean_degree_carn,modules, all.x=T)

## All ----
all_props <- prop %>% 
  dplyr::left_join(conec) %>% 
  dplyr::left_join(div) %>% 
  dplyr::left_join(mean_degree_carn)


# Save results ------------------------------------------------------------
save(comm_spp_trait, all_int, all_props,
     file = "results/preliminaryres_250924.Rdata")

write.csv(all_int, file = "results/InteractionProbability.csv")
write.csv(all_props, file = "results/NetworkProperties.csv")
