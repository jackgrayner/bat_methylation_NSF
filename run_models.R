#run models

library(viridis)
library(car)
library(reshape2)
library(methylclock)
library(MammalMethylClock)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)

SexPalette1<-c("#fa8490","#5090bf","purple")
setwd("~/Desktop/umd/nsf_methylation/")

methmeta.final<-read.csv('~/Desktop/umd/nsf_methylation/methmeta_final.csv',row.names = 1)
batdat<-readRDS("~/Desktop/umd/nsf_methylation/batdat_sexspp.rds")

print(methmeta.final %>% group_by(Species,Sex) %>% dplyr::summarize(n=n()),n=Inf)
methmeta.final$Species<-factor(methmeta.final$Species,levels=c("Desmodus rotundus","Phyllostomus discolor","Phyllostomus hastatus","Carollia perspicillata",
                                                               "Rousettus aegyptiacus","Pteropus pumilus","Pteropus vampyrus","Pteropus hypomelanus",
                                                               "Eptesicus fuscus","Tadarida brasiliensis",NA,"Saccopteryx bilineata"))
final.species<-levels(factor(methmeta.final$Species))

#######  
#### FIRST RUN MODELS OF SEX + AGE (NON-SEX-SPECIFIC)
######
#test
#Anova(lm(m ~ Sex * scale(Age,scale = FALSE),data=batdat[batdat$Species=="Tadarida brasiliensis" & batdat$site=="cg00000165",] ),type="III")
#add lab info for Eptesicus models
batdat<-merge(batdat,methmeta.final[,c("SID","Lab")],by="SID")

final.species<-"Eptesicus fuscus"
for (species_full in final.species) {
 species1<-gsub(" ","_",species_full)
 if (species_full=="Eptesicus fuscus"){
  temp<-batdat[batdat$Species==species_full,] %>%
   nest_by(site) %>% 
   mutate(spp=species1,
          model.coefs=list(summary(lm(m ~ Sex + scale(Age,scale = FALSE) + factor(Lab), data = data))$coefficients),
          model.test=list(Anova(lm(m ~ Sex + scale(Age,scale = FALSE) + factor(Lab), data = data))),
          intercept = model.coefs[1,1],
          sex_coef = model.coefs[2,1],age_coef = model.coefs[3,1],
          sex_test_f = model.test[1,3],sex_test_p = model.test[1,4],         
          age_test_f = model.test[2,3],age_test_p = model.test[2,4],
          int.coef=(summary(lm(m ~ Sex * scale(Age,scale = FALSE) + factor(Lab), data = data))$coefficients[4,1]),
          int.test=list(Anova(lm(m ~ Sex * scale(Age,scale = FALSE) + factor(Lab), data = data),type="III")),
          int_test_f=int.test[4,3],int_test_p=int.test[4,4])
 }else{
  temp<-batdat[batdat$Species==species_full,] %>%
   nest_by(site) %>% 
   mutate(spp=species1,
          model.coefs=list(summary(lm(m ~ Sex + scale(Age,scale = FALSE), data = data))$coefficients),
          model.test=list(Anova(lm(m ~ Sex + scale(Age,scale = FALSE), data = data))),
          intercept = model.coefs[1,1],
          sex_coef = model.coefs[2,1],age_coef = model.coefs[3,1],
          sex_test_f = model.test[1,3],sex_test_p = model.test[1,4],         
          age_test_f = model.test[2,3],age_test_p = model.test[2,4],
          int.coef=(summary(lm(m ~ Sex * scale(Age,scale = FALSE), data = data))$coefficients[4,1]),
          int.test=list(Anova(lm(m ~ Sex * scale(Age,scale = FALSE), data = data),type="III")),
          int_test_f=int.test[4,3],int_test_p=int.test[4,4])
 }
 write.csv(data.frame(temp[,c(1,3,6:13,15:16)]),file=paste("~/Desktop/umd/nsf_methylation/",species1,"_","meth_sex_age.csv",sep=""),row.names=FALSE,quote=FALSE)
 rm(temp)
}

#######  
#### NEXT RUN SEX-SPECIFIC AGE MODELS
######

for (species_full in final.species) {
 species1<-gsub(" ","_",species_full)
 if (species_full=="Eptesicus fuscus"){
  temp<-batdat[batdat$Species==species_full,] %>%
   nest_by(site) %>% 
   mutate(spp=species1,
          Fmodel.coefs=list(summary(lm(m ~ scale(Age,scale = FALSE)+factor(Lab), data = data %>% filter(Sex=="F")))$coefficients),
          Fmodel.test=list(Anova(lm(m ~ scale(Age,scale = FALSE)+factor(Lab), data = data %>% filter(Sex=="F")))),
          Fintercept = Fmodel.coefs[1,1],
          Fage_coef = Fmodel.coefs[2,1],
          Fage_test_f = Fmodel.test[1,3],Fage_test_p = Fmodel.test[1,4],
          Mmodel.coefs=list(summary(lm(m ~ scale(Age,scale = FALSE)+factor(Lab), data = data %>% filter(Sex=="M")))$coefficients),
          Mmodel.test=list(Anova(lm(m ~ scale(Age,scale = FALSE)+factor(Lab), data = data %>% filter(Sex=="M")))),
          Mintercept = Mmodel.coefs[1,1],
          Mage_coef = Mmodel.coefs[2,1],
          Mage_test_f = Mmodel.test[1,3],Mage_test_p = Mmodel.test[1,4])
 }else{
  temp<-batdat[batdat$Species==species_full,] %>%
   nest_by(site) %>% 
   mutate(spp=species1,
          Fmodel.coefs=list(summary(lm(m ~ scale(Age,scale = FALSE), data = data %>% filter(Sex=="F")))$coefficients),
          Fmodel.test=list(Anova(lm(m ~ scale(Age,scale = FALSE), data = data %>% filter(Sex=="F")))),
          Fintercept = Fmodel.coefs[1,1],
          Fage_coef = Fmodel.coefs[2,1],
          Fage_test_f = Fmodel.test[1,3],Fage_test_p = Fmodel.test[1,4],
          Mmodel.coefs=list(summary(lm(m ~ scale(Age,scale = FALSE), data = data %>% filter(Sex=="M")))$coefficients),
          Mmodel.test=list(Anova(lm(m ~ scale(Age,scale = FALSE), data = data %>% filter(Sex=="M")))),
          Mintercept = Mmodel.coefs[1,1],
          Mage_coef = Mmodel.coefs[2,1],
          Mage_test_f = Mmodel.test[1,3],Mage_test_p = Mmodel.test[1,4])
 }
 write.csv(data.frame(temp[,c(1,3,6:9,12:15)]),file=paste("~/Desktop/umd/nsf_methylation/",species1,"_","meth_sex-specific_age.csv",sep=""),row.names=FALSE,quote=FALSE)
 rm(temp)
}
