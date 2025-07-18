#Nascimento et al. 2024 - "The reorganization of predator-prey networks over 20  million years explains extinction patterns of mammalian carnivores"
#Constructing predator-prey networks
#============================

#loading packages
require(msm)
require(survival)
require(ggfortify)
require(ggpubr)


#loading data
source("0.load_data.r")

#Creating a column with species names
func$species=paste(func$genus,func$species,sep="_")

#Reloading the results we obtained in script 2. longevity_and_extinction
carnivora_degree=read.table("results/degree_carnivora.txt",h=T)
carnivora_degree_relat=read.table("results/degree_carnivora_relat.txt",h=T)

longevity=read.table("results/longevity.txt",h=T)
longevity_relat=read.table("results/longevity_relat.txt",h=T)

#Identifying common species (>5 occurrences)
ocor_species=data.frame(species = as.character(communities[,1]),
                        occ = apply(communities[,-1],1,sum))
#ocor_species=as.data.frame(ocor_species)
ocor_species$occ=as.numeric(ocor_species$occ)
ocor_species=ocor_species[order(ocor_species$occ, decreasing = T),]

carn_species=func[func$order == "Carnivora",]$species
ocor_carn=ocor_species[ocor_species$species %in% carn_species,]
common_species=ocor_carn[ocor_carn$occ >4,]$species

#Changing names of the variables
names(carnivora_degree)[9]="degree"
names(carnivora_degree_relat)[9]="degree_relat"

#Merging both datasets to create a single object
carnivora_degree=merge(carnivora_degree,carnivora_degree_relat)
carnivora_degree=merge(carnivora_degree,func,by="species")

#Creating a dataset with the minimum age for each species
min_ages=data.frame(unique(carnivora_degree$species))
names(min_ages)=c("species")
for (i in 1:nrow(min_ages)){
min_ages$min[i]=max(carnivora_degree[carnivora_degree$species == min_ages$species[i],]$MIN.age)
}

###creating the longevity tables beggining at 0 (oldest occurence) to the species extinction
carnivora_degree$t_inf="NA"
carnivora_degree$t_sup="NA"
for (i in 1:nrow(carnivora_degree)){
  carnivora_degree$t_inf[i]=min_ages[min_ages$species == carnivora_degree$species[i],]$min-carnivora_degree$MIN.age[i]
  carnivora_degree$t_sup[i]=min_ages[min_ages$species == carnivora_degree$species[i],]$min-carnivora_degree$MAX.age[i]
}


##Determining if each species was extant or extinct at each time point
for (i in 1:nrow(carnivora_degree)){
  data=carnivora_degree[carnivora_degree$species == carnivora_degree$species[i],]
  if (carnivora_degree$MAX.age[i]==min(data$MAX.age)){carnivora_degree$state[i] = 2}else{carnivora_degree$state[i]=1}
  if (carnivora_degree$MAX.age[i] == 0){carnivora_degree$state[i]=1}
}

##Ordering per time
carnivora_degree=carnivora_degree[order(carnivora_degree$species,carnivora_degree$t_sup),]


#Restricting to the 3 FF identified by Blanco et al 2021
fes=c(3,4,8)

carnivora_degree=carnivora_degree[carnivora_degree$module %in% fes,]
carnivora_degree_common=carnivora_degree[carnivora_degree$species %in% common_species,]
carnivora_degree_common$species <- factor(carnivora_degree_common$species)

table(carnivora_degree_common$species)

require(dplyr)

sp.common <- levels(carnivora_degree_common$species)


plot(1, type = "n", ylim = c(0,2.5),xlim = c(-20,0))
abline(h = 1, lty = 2)

#Condensing the values per geological unit 
for(k in 1:length(sp.common)){
  
  temp <- filter(carnivora_degree_common, species == sp.common[k])
  k.temp <- tapply(temp$degree,INDEX = temp$Geological.unit,mean)
  kr.temp <- tapply(temp$degree_relat,INDEX = temp$Geological.unit,mean)
  
  age.temp <- tapply(temp$age,INDEX = temp$Geological.unit,mean)
  aux <- sum(!is.na(age.temp))
  
  temp.df <- data.frame(species = rep(unique(temp$species),aux),
                        family = rep(unique(temp$family),aux),
                        body.size = rep(unique(temp$body.size.x),aux),
                        unit = names(age.temp[!is.na(age.temp)]),
                        age = age.temp[!is.na(age.temp)],
                        degree = k.temp[!is.na(k.temp)],
                        degree_relat = kr.temp[!is.na(k.temp)],
                        t_sup = max(age.temp[!is.na(age.temp)])-age.temp[!is.na(age.temp)])
  temp.df <- temp.df[order(temp.df$age, decreasing = T),]
  temp.df <- cbind(temp.df, state = c(rep(1,nrow(temp.df)-1),2))
  temp.df$degree.standard = temp.df$degree/mean(temp.df$degree) #temp.df$degree/temp.df$degree[1]
  
  
  
  if(k==1){ 
    df.all <- temp.df
  }else{
    df.all = rbind(df.all,temp.df)
  }
  col = c(rep("black",nrow(temp.df)-1), "tomato")
  points(-temp.df$age,temp.df$degree.standard,type = "o", pch = 16, col = col)
}

#Removing extant species from the dataset
df.all <- df.all[-which(df.all$species == "Canis_lupus"),]
df.all <- df.all[-which(df.all$species == "Meles_meles"),]
df.all <- df.all[-which(df.all$species == "Lynx_pardinus"),]
#Removing species with k = 0
df.all <- df.all[-which(df.all$species == "Protictitherium_crassum"),]

#Calculating the Multi-State Markov models described in methodology, in order: Constant, time variable, degree, mean degree, body size and all biological covariates

q=rbind(c(0,0.5),c(0,1))
Q.crude <- crudeinits.msm(state ~ t_sup, species, data=df.all,qmatrix=q)

ext.msm.constant<- msm(state ~ t_sup, subject=species, data = df.all,
                           qmatrix = Q.crude, deathexact = 2)
ext.msm.time <- msm(state ~ t_sup, subject=species, data = df.all,
                      qmatrix = Q.crude, deathexact = 2, covariates =~t_sup)
ext.msm.degree <- msm(state ~ t_sup, subject=species, data = df.all,
                      qmatrix = Q.crude, deathexact = 2, covariates =~degree)
ext.msm.degree.relat <- msm(state ~ t_sup, subject=species, data = df.all,
                      qmatrix = Q.crude, deathexact = 2, covariates =~degree_relat)
ext.msm.body.size <- msm(state ~ t_sup, subject=species, data = df.all,
                      qmatrix = Q.crude, deathexact = 2, covariates =~body.size)
ext.msm.cov <- msm(state ~ t_sup, subject=species, data = df.all,
                         qmatrix = Q.crude, deathexact = 2, covariates =~body.size+degree)

#Calculating the AIC for each model
aics.common=data.frame(model=c("ext.msm.constant","ext.msm.time","ext.msm.degree","ext.msm.degree.relat","ext.msm.body.size","ext.msm.cov"),
                       aics=c(AIC(ext.msm.constant),AIC(ext.msm.time),AIC(ext.msm.degree),AIC(ext.msm.degree.relat),AIC(ext.msm.body.size),AIC(ext.msm.cov)))
aics.common

#Summary for the best model identified based on AIC
summary(ext.msm.degree)
qmatrix.msm(ext.msm.degree)

#Creating dataset with information on each model
degree_plot=rbind(hazard.msm(ext.msm.degree)[[1]],
                  hazard.msm(ext.msm.degree.relat)[[1]],
                  hazard.msm(ext.msm.time)[[1]],
                  hazard.msm(ext.msm.body.size)[[1]],
                  hazard.msm(ext.msm.cov)[[1]])
degree_plot=as.data.frame(degree_plot)
degree_plot$model=c("degree","relative degree","time","body size", "body size+degree")

#Saving results of model fit
write.table(degree_plot,"results/HR_models.txt",quote=F,row.names=F,sep="\t")

#Plotting the results of the MSM model (Figure 3c, d)
require(survival)
require(ggfortify)
s=read.table("results/HR_models.txt",h=T,sep="\t")
s=s[-2,]

m1=ggplot(s, aes(x=reorder(model,U), y=HR))+theme_classic()+
  ylim(c(0,2))+geom_errorbar(aes(ymin=L,ymax=U, width=0.05),size=0.8)+
  geom_hline(yintercept=1, linetype="dashed")+
  ylab("Hazard Ratio")+xlab("Model")+
  geom_point(fill="white",shape=21,size=2)

#Calculating and plotting the survival curve (Figure 3d)
data_survival=plot.survfit.msm(ext.msm.degree, survdata=T)
fit=survfit(Surv(data_survival$survtime,data_survival$died)~1, data=data_survival)
m2=autoplot(fit,surv.linetype = 'dashed',censor=F)+theme_classic()+xlab("Time (Ma)")+ylab("Survival (%)")
