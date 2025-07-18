#Nascimento et al. 2024 - "The reorganization of predator-prey networks over 20  million years explains extinction patterns of mammalian carnivores"
#Constructing predator-prey networks
#============================

##loading packages
library(igraph)
library(dplyr)
library(tidyr)
library(ggplot2)
library(bipartite)
require(ggpubr)


#loading data
source("0.load_data.r")


#renaming taxon - merging genus and sp.
carn = func[func$order == "Carnivora",]
carn$species2 = paste(carn$genus, carn$species, sep="_")
herb = func[func$order != "Carnivora",]
herb$species2 = paste(herb$genus, herb$species, sep="_")

func$species = paste(func$genus,func$species, sep="_" )
herb_species = func[func$order != "Carnivora",]$species
carn_species = func[func$order == "Carnivora",]$species
aquatic = func[func$locomotion == "aquatic",]$species

#changing pack hunting species body size to account for the mode of hunting 
social=c("Canis_etruscus","Canis_lupus","Canis_mosbachensis","Cuon_alpinus","Eucyon_debonisi","Eucyon_monticinensis","Lycaon_falconeri")
func[func$species %in% social,]$body.size=4*func[func$species %in% social,]$body.size


#creating community data frames for each locality
occ_carn = communities[communities$species %in% carn_species,]
occ_herb = communities[communities$species %in% herb_species,]

communities_carn = list()
for (i in 1:169){
  k = occ_carn[occ_carn[,i+1] == 1,]
  communities_carn[[i]] = k$species
}
communities_carn
names(communities_carn) = names(occ_carn[2:170])

communities_herb = list()
for (i in 1:169){
  k = occ_herb[occ_herb[,i+1] == 1,]
  communities_herb[[i]] = k$species
}
communities_herb
names(communities_herb) = names(communities_carn)

#Including functional data in community data frame for each locality
herb_func = list() #herbivores
for (i in 1:length(communities_herb)){
  k = as.data.frame(communities_herb[[i]])
  names(k)="species"
  data = merge(k,func)
  data = data[c("species","body.size","diet","locomotion")]
  names(data) = c("species","mass","trophic","locomotion")
  herb_func[[i]] = data
}
names(herb_func) = names(communities_herb)

carn_func = list()#carnivores
for (i in 1:length(communities_carn)){
  k = as.data.frame(communities_carn[[i]])
  names(k) = "species"
  data = merge(k,func)
  data = data[c("species","body.size","diet","locomotion")]
  data = data[data$locomotion != "aquatic",]
  data = data[data$body.size > 10, ]
  names(data) = c("species","mass","trophic","locomotion")
  carn_func[[i]] = data
}
names(carn_func) = names(communities_carn)

#Log-ratio model function, the model described in the methods sections of the paper
LRM <- function(mass_C = NULL, mass_H = NULL){ 
  m <- nrow(mass_C)
  n <- nrow (mass_H)
  P<-matrix(NA,m,n)
  
  # sampling model parameters from parameter ranges 
  # Carnivores and Omnivores have different parameterization
  alphaP<-runif(1,min= 1, max= 2)
  betaP<-runif(1,min= -2,max= -1)
  gammaP<-runif(1,min= -1,max= -0.5)

  alphaO<-runif(1,min= -6 , max= -4)
  betaO<-runif(1,min= -3,max= -2)
  gammaO<- 0#runif(1,min= -0.5,max= -0.2)
  
  
  for (i in 1:m){
    
    if(mass_C$trophic[i]!="omnivore"){
      alpha <- alphaP
      beta <- betaP
      gamma <- gammaP
    }else{
      alpha <- alphaO
      beta <- betaO
      gamma <- gammaO
    }
    
    for(j in 1:n){
      
      term<-(exp(alpha+(beta*log((mass_H$mass[j]/mass_C$mass[i])))+(gamma*(log((mass_H$mass[j]/mass_C$mass[i]))^2))))
      p<-term/(1+term) #probability of interaction
      P[i,j]<-p #filling probability matrix	
      if(mass_H$locomotion[j]=="arboreal"){P[i,j]=P[i,j]*0.1} #reducing interaction probability btw arboreal vs terrestrial species
      #if(mass_C$trophic[i]=="carnivore_invert" & mass_H$trophic[j]!="invertebrate"){P[i,j]=0}  #not used in main analyses  
      #if(mass_H$trophic[j]=="invertebrate" & mass_C$trophic[i]!="carnivore_invert"){P[i,j]=0}  #not used in main analyses 
      
    }
    if(mass_C$trophic[i]=="omnnivore"){P[i,]=P[i,]*0.5} #rescaling predation probability for omnivores
    
  } 
  
  row.names(P) <- as.character(mass_C$species)
  colnames(P) <- as.character(mass_H$species)
  return(P)
}

########Calculating the networks for each locality using the log-ratio model 
########Results are the interaction probability between a pair of species for a given locality

results = list()
replicas = list()
for (i in 1:169){ #running for all localities
  if (nrow(carn_func[[i]])>0 & nrow(herb_func[[i]])>0){
    if (nrow(carn_func[[i]])>1 | nrow(herb_func[[i]])>1){
      N = 100
      P.list <- replicate(N,LRM(carn_func[[i]], herb_func[[i]]))
      P.list[P.list < 0.05] <- 0
      replicas[[i]]=P.list
      results[[i]]=apply(simplify2array(P.list), 1:2, mean)
      results[[i]]=round(results[[i]],5)
      
    }
  }
}
names(results)=names(carn_func)

###Calculating the metrics we describe in the paper: predator-prey ratio, connectance, interaction diversity and species degree

##predator-prey ratio
prop = matrix(ncol=2,nrow=169)
for (i in 1:169){
  prop[i,1]=names(results[i])
  if (!is.null(results[[i]])){
    prop[i,2]=nrow(results[[i]])/ncol(results[[i]])
  }
}
prop = as.data.frame(prop)
names(prop) = c("Site", "Prop")
prop = merge(prop,ages)
prop$age = (prop$MIN.age+prop$MAX.age)/2
prop = merge(prop,riq)
prop = merge(prop,modules, all.x=T)

##Connectance
conec = matrix(nrow=169, ncol=2)
for (i in 1:169){
  conec[i,1] = names(results[i])
  if (!is.null(results[[i]])){
    f = as.matrix(results[[i]])
    if (nrow(f)>1&ncol(f)>1){
    conec[i,2] = conec[i,2]=sum(results[[i]]/(ncol(results[[i]])*nrow(results[[i]])))
    }
  }
}
    
conec = as.data.frame(conec)
names(conec) = c("Site", "Connectance")
conec = merge(conec,ages)
conec$age = (conec$MIN.age+conec$MAX.age)/2
conec = merge(conec,riq)
conec = merge(conec,modules, all.x=T)
conec$Connectance = as.numeric(conec$Connectance)


##Interaction diversity
div = matrix(nrow=169, ncol=2)
for (i in 1:169){
  div[i,1] = names(results[i])
  if (!is.null(results[[i]])){
    f = as.matrix(results[[i]])
    if (nrow(f)>1&ncol(f)>1){
      
      div[i,2] = diversity(as.numeric(results[[i]],index="shannon"))
    }
  }
}

div = as.data.frame(div)
names(div)=c("Site", "Diversity")
div=merge(div,ages)
div$age=(div$MIN.age+div$MAX.age)/2
div=merge(div,riq)
div=merge(div,modules, all.x=T)
div$Diversity=as.numeric(div$Diversity)

##Mean degree
mean_degree_carn=matrix(ncol=2,nrow=169)
for (i in 1:169){
  mean_degree_carn[i,1]=names(results)[i]
  if (!is.null(results[[i]])){
    mean_degree_carn[i,2]=mean(rowSums(results[[i]]))}
}
mean_degree_carn=as.data.frame(mean_degree_carn
)
names(mean_degree_carn)=c("Site","mean_degree_carn")
mean_degree_carn$mean_degree_carn=as.numeric(mean_degree_carn$mean_degree_carn)
mean_degree_carn = merge(mean_degree_carn,ages)
mean_degree_carn$age = (mean_degree_carn$MIN.age+mean_degree_carn$MAX.age)/2
mean_degree_carn = merge(mean_degree_carn,riq)
mean_degree_carn = merge(mean_degree_carn,modules, all.x=T)

#Saving results
write.table(div, "results/metrics/Diversity_mean.txt", row.names=F, quote=F)
write.table(conec, "results/metrics/Connectance_mean.txt", row.names=F, quote=F)
write.table(prop, "results/metrics/Predator_prey_ratio.txt", row.names=F, quote=F)
write.table(mean_degree_carn, "results/metrics/Mean_degree_carnivore.txt", row.names=F)

###################Plotting networks################
##restricting analysis to the 3 functional faunas identified by Blanco et al. (2021)
fes=c(3,4,8)
conec=conec[conec$module %in% fes,]
div=div[div$module %in% fes,]
prop=prop[prop$module %in% fes,]
mean_degree_carn=mean_degree_carn[mean_degree_carn$module %in% fes,]
prop$module=as.character(prop$module)
prop$Prop=as.numeric(prop$Prop)


###Changing modules names

conec$module=gsub("3","FF1",conec$module)
conec$module=gsub("4","FF2",conec$module)
conec$module=gsub("8","FF3",conec$module)

div$module=gsub("3","FF1",div$module)
div$module=gsub("4","FF2",div$module)
div$module=gsub("8","FF3",div$module)

prop$module=gsub("3","FF1",prop$module)
prop$module=gsub("4","FF2",prop$module)
prop$module=gsub("8","FF3",prop$module)

mean_degree_carn$module=gsub("3","FF1",mean_degree_carn$module)
mean_degree_carn$module=gsub("4","FF2",mean_degree_carn$module)
mean_degree_carn$module=gsub("8","FF3",mean_degree_carn$module)


#Plotting results
m1=ggplot(conec, aes(x=-age,y=Connectance, col=module))+geom_point(size=3)+theme_classic()+scale_color_manual(values = c("#723294","#E0A263","#76C299"))+ylab("Connectance")+geom_smooth(aes(x=-age,y=Connectance),method="lm",inherit.aes=F,col="gray50")+xlab("age (Ma)")+labs(color='Functional fauna')+labs(color='Functional fauna')+scale_x_continuous(labels = abs)
m2=ggplot(div, aes(x=-age,y=exp(Diversity), col=module))+geom_point(size=3)+theme_classic()+scale_color_manual(values = c("#723294","#E0A263","#76C299"))+ylab("Interaction diversity")+geom_smooth(aes(x=-age,y=exp(Diversity)),method="lm",inherit.aes=F,col="gray50")+xlab("age (Ma)")+labs(color='Functional fauna')+labs(color='Functional fauna')+scale_x_continuous(labels = abs)
m3=ggplot(prop, aes(x=-age,y=Prop, col=module))+geom_point(size=3)+theme_classic()+scale_color_manual(values = c("#723294","#E0A263","#76C299"))+ylab("Predator/Prey Ratio")+ylim(0,1.8)+geom_smooth(aes(x=-age,y=Prop),method="lm",inherit.aes=F,col="gray50")+xlab("age (Ma)")+labs(color='Functional fauna')+labs(color='Functional fauna')+scale_x_continuous(labels = abs)
m4=ggplot(mean_degree_carn, aes(x=-age,y=mean_degree_carn, col=module))+geom_point(size=3)+theme_classic()+scale_color_manual(values = c("#723294","#E0A263","#76C299"))+ylab("Mean degree")+geom_smooth(aes(x=-age,y=mean_degree_carn),method="lm",inherit.aes=F,col="gray50")+xlab("age (Ma)")+labs(color='Functional fauna')+labs(color='Functional fauna')+scale_x_continuous(labels = abs)

#Plotting panels
ggarrange(m3,m1,m2,m4, labels=c("(a)", "(b)", "(c)","(d)"))

#Saving results for posterior use
saveRDS(results, "results/networks.RData")
