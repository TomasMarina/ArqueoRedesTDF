#Nascimento et al. 2024 - "The reorganization of predator-prey networks over 20  million years explains extinction patterns of mammalian carnivores"
#Computing the relationship between species longevity and predictors
#====================================================================

#loading packages
library(igraph)
library(dplyr)
library(tidyr)
library(ggplot2)
library(bipartite)
require(ggpubr)

#loading data
source("0.load_data.r")

#Creating a column with species name
func$species=paste(func$genus,func$species, sep="_" )

#Creating files with mass and diet of each species
mass=func[c("species","body.size")]
diet=func[c("species","diet")]

#Read netowkors we calculated using script 1. Iberica_networks
results=readRDS("results/networks.RData")

###Calculating carnivorans degree (sum of the probabilities of interaction) for each species in each locality it occurs###
carnivora=func[func$order == "Carnivora",]
carnivora=carnivora[carnivora$locomotion != "aquatic",]
species_carnivora=carnivora$species

degree_carnivora=matrix(ncol=157,nrow=169)
degree_carnivora[,1]=names(results)
for (j in 1:nrow(carnivora)){
  for (i in 1:length(results)){
    p=results[[i]]
    if (!is.null(p)){
      s=sum(p[row.names(p) == carnivora$species[j],])
      if (dim(as.data.frame(p[row.names(p) == carnivora$species[j],]))[1]==0)
      {degree_carnivora[i,j+1]="NA"}else{
        degree_carnivora[i,j+1]=s
      }
    }
  }
}


#Merging with information for the sites
degree_carnivora=as.data.frame(degree_carnivora)
names(degree_carnivora)=c("Site",carnivora$species)
degree_carnivora=merge(degree_carnivora,ages)
degree_carnivora$age=(degree_carnivora$MIN.age+degree_carnivora$MAX.age)/2
degree_carnivora=merge(degree_carnivora,riq)
degree_carnivora=merge(degree_carnivora,modules, all.x=T)

degree_carnivora[,2:157]=as.numeric(unlist(degree_carnivora[,2:157]))

#Changing format for posterior use
v_carnivora=pivot_longer(degree_carnivora, cols=2:157)
names(v_carnivora)=c(names(v_carnivora)[1:7],"species","value")

v_carnivora=merge(v_carnivora,mass)
v_carnivora=v_carnivora[!is.na(v_carnivora$value),]
v_carnivora$body.size=as.factor(v_carnivora$body.size)

#saving results
write.table(v_carnivora, "results/degree_carnivora.txt",sep="\t",row.names=F,quote=F)

####Calculating carnivora relative degree (the proportion of the total number of prey species available that can be consumed)
degree_carnivora_relat=matrix(ncol=157,nrow=169)
degree_carnivora_relat[,1]=names(results)
for (j in 1:nrow(carnivora)){
  for (i in 1:length(results)){
    p=results[[i]]
    if (!is.null(p)){
      s=sum(p[row.names(p) == carnivora$species[j],])/ncol(p)
      if (dim(as.data.frame(p[row.names(p) == carnivora$species[j],]))[1]==0)
      {degree_carnivora_relat[i,j+1]="NA"}else{
        degree_carnivora_relat[i,j+1]=s
      }
    }
  }
}

#Merging with information of each site
degree_carnivora_relat=as.data.frame(degree_carnivora_relat)
names(degree_carnivora_relat)=c("Site",carnivora$species)
degree_carnivora_relat=merge(degree_carnivora_relat,ages)
degree_carnivora_relat$age=(degree_carnivora_relat$MIN.age+degree_carnivora_relat$MAX.age)/2
degree_carnivora_relat=merge(degree_carnivora_relat,riq)
degree_carnivora_relat=merge(degree_carnivora_relat,modules, all.x=T)

degree_carnivora_relat[,2:157]=as.numeric(unlist(degree_carnivora_relat[,2:157]))


#Changing format to longer for posterior use
v_carnivora_relat=pivot_longer(degree_carnivora_relat, cols=2:157)
names(v_carnivora_relat)=c(names(v_carnivora_relat)[1:7],"species","value")

v_carnivora_relat=merge(v_carnivora_relat,mass)
v_carnivora_relat=v_carnivora_relat[!is.na(v_carnivora_relat$value),]
v_carnivora_relat$body.size=as.factor(v_carnivora_relat$body.size)

#saving results
write.table(v_carnivora_relat, "results/degree_carnivora_relat.txt",sep="\t",row.names=F,quote=F)


#reloading the results we calculated in script 1. Iberian networks
carnivora_degree=read.table("results/degree_carnivora.txt",h=T)
carnivora=func[func$order == "Carnivora",]
carnivora=carnivora[carnivora$body.size > 5.5,]
carnivora=carnivora[carnivora$diet != "carnivore_invert",]

#Common species
#Obtaining number of occurrences per species
ocor_species=cbind(communities$species, apply(communities[,-1],1,sum))
ocor_species=as.data.frame(ocor_species)
ocor_species$V2=as.numeric(ocor_species$V2)
ocor_species=ocor_species[order(ocor_species$V2, decreasing = T),]
ocor_carn=ocor_species[ocor_species$V1 %in% carnivora$species,]
common_species=ocor_carn[ocor_carn$V2 >4,]$V1

###Calculating carnivores longevity and mean degree of each species###
longevity=matrix(ncol=5,nrow=length(common_species))
for (i in 1:length(common_species)){
  v=common_species[i]
  k=carnivora_degree %>% filter(species == v)
  longevity[i,1]=v
  longevity[i,2]=max(k$age)
  longevity[i,3]=min(k$age)
  longevity[i,4]=max(k$age)-min(k$age)
  longevity[i,5]=mean(k$value)
}

#Renaming columns
longevity=as.data.frame(longevity)
names(longevity)=c("species","max_age","min_age","longevity","mean_degree")

#Changing format and excluding living species
longevity=merge(longevity,func,by="species")
longevity$mean_degree=as.numeric(longevity$mean_degree)
longevity$longevity=as.numeric(longevity$longevity)
longevity=longevity[longevity$min_age != 0,]

#Linear model between longevity and mean degree (Figure 3b)
lm_long_degree=lm(longevity$longevity ~longevity$mean_degree)
lm_long_degree

rsqd_mean=round(summary(lm_long_degree)$r.squared,3)

#plotting results
p4=ggplot(longevity,aes(x=mean_degree,y=longevity))+geom_point(size=2)+geom_smooth(method="lm")+theme_classic()+
  xlab("Mean degree")+ylab("Longevity")+annotate(geom="text",label=paste0("R?: ",rsqd_mean),x=0.5,y=6)
p4

#Saving results
write.table(longevity, "results/longevity.txt", row.names=F)

#Calculating relative degree
longevity_relat=matrix(ncol=5,nrow=length(common_species))
for (i in 1:length(common_species)){
  v=common_species[i]
  k=carnivora_degree_relat%>% filter(species == v)
  longevity_relat[i,1]=v
  longevity_relat[i,2]=max(k$MIN.age)
  longevity_relat[i,3]=min(k$MAX.age)
  longevity_relat[i,4]=max(k$MIN.age)-min(k$MAX.age)
  longevity_relat[i,5]=mean(k$value)
}

#changing column names
longevity_relat=as.data.frame(longevity_relat)
names(longevity_relat)=c("species","max_age","min_age","longevity","mean_degree")

#changing format and excluding living species
longevity_relat=merge(longevity_relat,func,by="species")
longevity_relat$mean_degree=as.numeric(longevity_relat$mean_degree)
longevity_relat$longevity=as.numeric(longevity_relat$longevity)
longevity_relat=longevity_relat[longevity_relat$min_age != 0,]

#Saving for posterior use
write.table(longevity_relat, "results/longevity_relat.txt", row.names=F)

####Calculating the propotion of disconnected predators (Figure 3a)#####
fes=c(3,4,8)
carnivora_degree=carnivora_degree[carnivora_degree$module %in% fes,]

desconec=matrix(ncol=3, nrow=169)
for (i in 1:169){
  v=names(results)[i]
  desconec[i,1]=v
  f=carnivora_degree[carnivora_degree$Site == v ,]
  desconec[i,2]=nrow(f[f$value < 1,])/nrow(f)
  desconec[i,3]=nrow(f[f$value < 2,])/nrow(f)
  
  
}

#renaming columns
desconec=as.data.frame(desconec)
names(desconec)=c("Site","Ratio1","Ratio2")

#Merging with information for each locality
desconec=merge(desconec, modules, all.x = T)
desconec=merge(desconec, ages)
desconec$age=(desconec$MIN.age+desconec$MAX.age)/2
desconec=merge(desconec,riq)
desconec$Ratio1=as.numeric(desconec$Ratio1)
desconec$Ratio2=as.numeric(desconec$Ratio2)
desconec=desconec[desconec$module %in% fes,]


#Plotting results
d1=ggplot(desconec, aes(x=-age,y=Ratio1))+geom_point(size=3)+theme_classic()+scale_color_brewer(palette = "Set1")+geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), col="gray50")+ylim(c(0,1))+xlab("age(Ma)")+ylab("Low connected predators (%)")+scale_x_continuous(labels=abs)
d2=ggplot(desconec, aes(x=-age,y=Ratio2))+geom_point(size=3)+theme_classic()+scale_color_brewer(palette = "Set1")+geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), col="gray50")+ylim(c(0,1))+xlab("age(Ma)")+ylab("Low connected predators (%)")+scale_x_continuous(labels=abs)

grid.arrange(d1,d2,ncol=2)
ggarrange(d1,d2, labels=c("(a)","(b)"), label.x=-0.02)

#Summary of the linear models between age and disconnected predators
summary(lm(desconec$Ratio1 ~desconec$age))
summary(lm(desconec$Ratio2 ~desconec$age))

#Plotting results
ggplot(desconec, aes(x=-age,y=Ratio2, col=module))+scale_color_manual(values = c("#723294","#E0A263","#76C299"))+geom_point(size=3)+theme_classic()+geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), col="gray50")+ylim(c(0,1))+xlab("age(Ma)")+ylab("Low connected predators (%)")+scale_x_continuous(labels=abs)
