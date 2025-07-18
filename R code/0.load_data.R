#Nascimento et al. 2024 - "The reorganization of predator-prey networks over 20  million years explains extinction patterns of mammalian carnivores"
#Loading data
#============================

##Loading packages
library(igraph)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

##Reading files with species occurences for each community
communities=read.csv("data/communities_iberian.txt", h=T,sep="\t",encoding="latin1")
riq=as.data.frame(colSums(communities[,2:170]))
names(riq)="Richness"
riq$Site=rownames(riq)

##Reading information of the module to which each locality was assigned to
modules=read.table("data/modules.txt", h=T)
modules$module=as.character(modules$module)

##Reading information of the ages of each site
ages=read.table("data/age_sites.txt", h=T, sep="\t",encoding="latin1")
ages=ages[1:169,]
ages$Site=gsub(" ",".",ages$Site)

##Loading and modifying the functional data of each species
func=read.table("data/functional_data.txt", h=T, sep="\t")
func$body.size=gsub(">1000",1000,func$body.size)
func$body.size=gsub("360-1000",680,func$body.size)
func$body.size=gsub("180-360",270,func$body.size)
func$body.size=gsub("90-180",135,func$body.size)
func$body.size=gsub("45-90",67.5,func$body.size)
func$body.size=gsub("10-45",27.5,func$body.size)
func$body.size=gsub("1-10",5.5,func$body.size)
func$body.size=gsub("< 1",0.5,func$body.size)
func$body.size=as.numeric(func$body.size)
