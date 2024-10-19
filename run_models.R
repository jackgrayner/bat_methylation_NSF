library(viridis)
library(reshape2)
library(methylclock)
library(MammalMethylClock)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(sesame)


SexPalette1<-c("#fa8490","#5090bf","purple")
setwd("~/Documents/UMD_new/DNAm-sex_diffs/")

#read original 11 spp, get rid of juves, NWO samples, and old PhHa females
methmeta.final<-read.csv('~/Documents/UMD_new/DNAm-sex_diffs/methmeta_final.csv',row.names = 1)
methmeta.final<-methmeta.final[!methmeta.final$Lab=="Northwest Ohio Med",]#get rid of NWO as strongly skewed sample
methmeta.final<-methmeta.final[!(methmeta.final$Species=="Phyllostomus hastatus" & methmeta.final$Age>10),]#get rid of PhHa older females
print(methmeta.final %>% group_by(Species,Sex) %>% dplyr::summarize(n=n()),n=Inf)

#add more species, get rid of juves, combine
# meta1<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/methmeta_geo.csv")
# meta1<-meta1[meta1$Species %in% c("Pteropus poliocephalus","Pteropus rodricensis") &
#                meta1$Age>0.99,]
# print(meta1 %>% group_by(Species,Sex) %>% dplyr::summarize(n=n()),n=Inf)
# meta.all<-data.frame(SID=c(methmeta.final$SID,meta1$SID),
#                      Species=c(methmeta.final$Species,meta1$Species),
#                      Sex=c(methmeta.final$Sex,meta1$Sex),
#                      Basename=c(methmeta.final$Basename,meta1$Basename),
#                      Age=c(methmeta.final$Age,meta1$Age))
# methmeta.final<-meta.all


#run models in SeSAMe or whatever it's called
manifest_sesame <- read.csv("~/Documents/UMD_new/DNAm-sex_diffs/geo_rawdat/HorvathMammal40.CanonicalManifest.3.2019.sesame.csv",)

normalized_betas_sesame<-
 openSesame("~/Documents/UMD_new/DNAm-sex_diffs/geo_rawdat/GSE164127_RAW/", 
            prep="SHCDPB",func = getBetas, manifest=manifest_sesame,
            mask=FALSE)

colnames(normalized_betas_sesame)<-gsub("GSM.*\\_20","20",colnames(normalized_betas_sesame))
normalized_betas_sesame<-normalized_betas_sesame[,colnames(normalized_betas_sesame) %in% methmeta.final$Basename]

rownames(methmeta.final)<-methmeta.final$Basename
summary(rownames(methmeta.final) == colnames(normalized_betas_sesame))#check order
methmeta.final<-methmeta.final[colnames(normalized_betas_sesame),]#reorder
summary(rownames(methmeta.final) == colnames(normalized_betas_sesame))#check order

final.species<-levels(factor(methmeta.final$Species))

#get rid of non-cg sites
normalized_betas_sesame<-normalized_betas_sesame[grep("cg",rownames(normalized_betas_sesame)),]

# ADDITIVE BOTH SEX MODELS
for (species in final.species){
 species1<-gsub(" ","_",species)
 methmeta.1<-methmeta.final[methmeta.final$Species==species,]
 methmeta.1$Age<-scale(methmeta.1$Age,scale=FALSE)
 normalized_betas_sesame1<-normalized_betas_sesame[,colnames(normalized_betas_sesame) %in% methmeta.1$Basename]
 smry<-DML(betas=normalized_betas_sesame1, 
           meta=methmeta.1,fm=~ Sex + Age)
 smry<-data.frame(summaryExtractTest(smry))
 smry$Species<-species
 write.csv(smry,
           file=paste("~/Documents/UMD_new/DNAm-sex_diffs/results_files/",species1,"_","meth_sex_age_sesame.csv",sep=""),
           row.names=FALSE,quote=FALSE)
 
}

# SINGLE-SEX MODELS FOR SLOPE COMPARISON
for (species in final.species){
  species1<-gsub(" ","_",species) 
 #females
 methmeta.1<-methmeta.final[methmeta.final$Species==species & methmeta.final$Sex=="F",]
 methmeta.1$Age<-scale(methmeta.1$Age,scale=FALSE)
 smry.f<-DML(betas=normalized_betas_sesame[,colnames(normalized_betas_sesame) %in% methmeta.1$Basename], 
             meta=methmeta.1[,c("Sex","Age")],fm=~ Age)
 smry.f<-data.frame(summaryExtractTest(smry.f))
 colnames(smry.f)<-paste(colnames(smry.f),".F")
 #males
 methmeta.1<-methmeta.final[methmeta.final$Species==species & methmeta.final$Sex=="M",]
 methmeta.1$Age<-scale(methmeta.1$Age,scale=FALSE)
 normalized_betas_sesame1<-normalized_betas_sesame[,methmeta.1$Sex=="M"]
 smry.m<-DML(betas=normalized_betas_sesame[,colnames(normalized_betas_sesame) %in% methmeta.1$Basename],
             meta=methmeta.1[,c("Sex","Age")],fm=~ Age)
 smry.m<-data.frame(summaryExtractTest(smry.m))
 smry.m$Species<-species
 colnames(smry.m)<-paste(colnames(smry.f),".M")
 
 smry<-data.frame(cbind(smry.f,smry.m))
 write.csv(smry,
           file=paste("~/Documents/UMD_new/DNAm-sex_diffs/results_files/",species1,"_","meth_sex-specific_age_sesame.csv",sep=""),
           row.names=FALSE,quote=FALSE)
 
}



rhino.ann<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/Rhinolophus_ferrumequinum.hlrhifer5.HorvathMammalMethylChip40.v1.csv",h=T)
xlinked<-rhino.ann[rhino.ann$geneChr==1,]$CGid
autosomal<-rhino.ann[rhino.ann$geneChr>1,]$CGid

#combine species additive model results
CaPe<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/results_files/Carollia_perspicillata_meth_sex_age_sesame.csv",h=T) %>% 
  filter(Probe_ID %in% rhino.ann[!is.na(rhino.ann$SYMBOL),]$CGid) %>%
  mutate(Pval_SexMadj=p.adjust(Pval_SexM),Pval_Ageadj=p.adjust(Pval_Age),SD="M-biased",spp1="CaPe",spp2="Carollia perspicillata",
         mating="harem") 
DeRo<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/results_files/Desmodus_rotundus_meth_sex_age_sesame.csv",h=T) %>%
  filter(Probe_ID %in% rhino.ann[!is.na(rhino.ann$SYMBOL),]$CGid) %>% 
  mutate(Pval_SexMadj=p.adjust(Pval_SexM),Pval_Ageadj=p.adjust(Pval_Age),SD="F-biased",spp1="DeRo",spp2="Desmodus rotundus",
         mating="polygynous")
EpFu<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/results_files/Eptesicus_fuscus_meth_sex_age_sesame.csv",h=T) %>%
  filter(Probe_ID %in% rhino.ann[!is.na(rhino.ann$SYMBOL),]$CGid) %>% 
  mutate(Pval_SexMadj=p.adjust(Pval_SexM),Pval_Ageadj=p.adjust(Pval_Age),SD="F-biased",spp1="EpFu",spp2="Eptesicus fuscus",
         mating="?")
PhDi<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/results_files/Phyllostomus_discolor_meth_sex_age_sesame.csv",h=T) %>%
  filter(Probe_ID %in% rhino.ann[!is.na(rhino.ann$SYMBOL),]$CGid) %>% 
  mutate(Pval_SexMadj=p.adjust(Pval_SexM),Pval_Ageadj=p.adjust(Pval_Age),SD="F-biased",spp1="PhDi",spp2="Phyllostomus discolor",
         mating="harem")
PhHa<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/results_files/Phyllostomus_hastatus_meth_sex_age_sesame.csv",h=T) %>%
  filter(Probe_ID %in% rhino.ann[!is.na(rhino.ann$SYMBOL),]$CGid) %>% 
  mutate(Pval_SexMadj=p.adjust(Pval_SexM),Pval_Ageadj=p.adjust(Pval_Age),SD="M-biased",spp1="PhHa",spp2="Phyllostomus hastatus",
         mating="harem")
PtHy<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/results_files/Pteropus_hypomelanus_meth_sex_age_sesame.csv",h=T) %>%
  filter(Probe_ID %in% rhino.ann[!is.na(rhino.ann$SYMBOL),]$CGid) %>% 
  mutate(Pval_SexMadj=p.adjust(Pval_SexM),Pval_Ageadj=p.adjust(Pval_Age),SD="M-biased",spp1="PtPy",spp2="Pteropus hypomelanus",
         mating="polygynous")
PtPu<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/results_files/Pteropus_pumilus_meth_sex_age_sesame.csv",h=T)  %>%
  filter(Probe_ID %in% rhino.ann[!is.na(rhino.ann$SYMBOL),]$CGid) %>% 
  mutate(Pval_SexMadj=p.adjust(Pval_SexM),Pval_Ageadj=p.adjust(Pval_Age),SD="M-biased",spp1="PtPu",spp2="Pteropus pumilus",
         mating="polygynous")
PtVa<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/results_files/Pteropus_vampyrus_meth_sex_age_sesame.csv",h=T) %>%
  filter(Probe_ID %in% rhino.ann[!is.na(rhino.ann$SYMBOL),]$CGid) %>% 
  mutate(Pval_SexMadj=p.adjust(Pval_SexM),Pval_Ageadj=p.adjust(Pval_Age),SD="M-biased",spp1="PtVa",spp2="Pteropus vampyrus",
         mating="polygynous")
RoAe<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/results_files/Rousettus_aegyptiacus_meth_sex_age_sesame.csv",h=T) %>%
  filter(Probe_ID %in% rhino.ann[!is.na(rhino.ann$SYMBOL),]$CGid) %>% 
  mutate(Pval_SexMadj=p.adjust(Pval_SexM),Pval_Ageadj=p.adjust(Pval_Age),SD="M-biased",spp1="RoAe",spp2="Rousettus aegyptiacus",
         mating="polygynandrous")
SaBi<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/results_files/Saccopteryx_bilineata_meth_sex_age_sesame.csv",h=T) %>%
  filter(Probe_ID %in% rhino.ann[!is.na(rhino.ann$SYMBOL),]$CGid) %>% 
  mutate(Pval_SexMadj=p.adjust(Pval_SexM),Pval_Ageadj=p.adjust(Pval_Age),SD="F-biased",spp1="SaBi",spp2="Saccopteryx bilineata",
         mating="harem")
TaBr<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/results_files/Tadarida_brasiliensis_meth_sex_age_sesame.csv",h=T) %>%
  filter(Probe_ID %in% rhino.ann[!is.na(rhino.ann$SYMBOL),]$CGid) %>% 
  mutate(Pval_SexMadj=p.adjust(Pval_SexM),Pval_Ageadj=p.adjust(Pval_Age),SD="F-biased",spp1="TaBr",spp2="Tadarida brasiliensis",
         mating="polygynandrous")

allspp<-rbind(CaPe,DeRo,EpFu,PhDi,PhHa,PtHy,PtPu,PtVa,RoAe,SaBi,TaBr)
allspp$XA<-NA
allspp[allspp$Probe_ID %in% xlinked,]$XA<-"X-linked"
allspp[allspp$Probe_ID %in% autosomal,]$XA<-"autosomal"

write.csv(allspp,"~/Documents/UMD_new/DNAm-sex_diffs/results_files/AllSpp_sex_age.csv")


#combine species sex-specific results
CaPe<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/results_files/Carollia_perspicillata_meth_sex-specific_age_sesame.csv",h=T) %>% 
 filter(Probe_ID..F %in% rhino.ann[!is.na(rhino.ann$SYMBOL),]$CGid) %>%
 mutate(Fage_test_padj=p.adjust(Pval_Age..F),Mage_test_padj=p.adjust(Pval_Age..F..M),SD="M-biased",spp1="CaPe",spp2="Carollia perspicillata",
        mating="harem") 
DeRo<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/results_files/Desmodus_rotundus_meth_sex-specific_age_sesame.csv",h=T) %>%
 filter(Probe_ID..F %in% rhino.ann[!is.na(rhino.ann$SYMBOL),]$CGid) %>% 
 mutate(Fage_test_padj=p.adjust(Pval_Age..F),Mage_test_padj=p.adjust(Pval_Age..F..M),SD="F-biased",spp1="DeRo",spp2="Desmodus rotundus",
        mating="polygynous")
EpFu<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/results_files/Eptesicus_fuscus_meth_sex-specific_age_sesame.csv",h=T) %>%
 filter(Probe_ID..F %in% rhino.ann[!is.na(rhino.ann$SYMBOL),]$CGid) %>% 
 mutate(Fage_test_padj=p.adjust(Pval_Age..F),Mage_test_padj=p.adjust(Pval_Age..F..M),SD="F-biased",spp1="EpFu",spp2="Eptesicus fuscus",
        mating="?")
PhDi<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/results_files/Phyllostomus_discolor_meth_sex-specific_age_sesame.csv",h=T) %>%
 filter(Probe_ID..F %in% rhino.ann[!is.na(rhino.ann$SYMBOL),]$CGid) %>% 
 mutate(Fage_test_padj=p.adjust(Pval_Age..F),Mage_test_padj=p.adjust(Pval_Age..F..M),SD="F-biased",spp1="PhDi",spp2="Phyllostomus discolor",
        mating="harem")
PhHa<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/results_files/Phyllostomus_hastatus_meth_sex-specific_age_sesame.csv",h=T) %>%
 filter(Probe_ID..F %in% rhino.ann[!is.na(rhino.ann$SYMBOL),]$CGid) %>% 
 mutate(Fage_test_padj=p.adjust(Pval_Age..F),Mage_test_padj=p.adjust(Pval_Age..F..M),SD="M-biased",spp1="PhHa",spp2="Phyllostomus hastatus",
        mating="harem")
PtHy<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/results_files/Pteropus_hypomelanus_meth_sex-specific_age_sesame.csv",h=T) %>%
 filter(Probe_ID..F %in% rhino.ann[!is.na(rhino.ann$SYMBOL),]$CGid) %>% 
 mutate(Fage_test_padj=p.adjust(Pval_Age..F),Mage_test_padj=p.adjust(Pval_Age..F..M),SD="M-biased",spp1="PtPy",spp2="Pteropus hypomelanus",
        mating="polygynous")
PtPu<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/results_files/Pteropus_pumilus_meth_sex-specific_age_sesame.csv",h=T)  %>%
 filter(Probe_ID..F %in% rhino.ann[!is.na(rhino.ann$SYMBOL),]$CGid) %>% 
 mutate(Fage_test_padj=p.adjust(Pval_Age..F),Mage_test_padj=p.adjust(Pval_Age..F..M),SD="M-biased",spp1="PtPu",spp2="Pteropus pumilus",
        mating="polygynous")
PtVa<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/results_files/Pteropus_vampyrus_meth_sex-specific_age_sesame.csv",h=T) %>%
 filter(Probe_ID..F %in% rhino.ann[!is.na(rhino.ann$SYMBOL),]$CGid) %>% 
 mutate(Fage_test_padj=p.adjust(Pval_Age..F),Mage_test_padj=p.adjust(Pval_Age..F..M),SD="M-biased",spp1="PtVa",spp2="Pteropus vampyrus",
        mating="polygynous")
RoAe<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/results_files/Rousettus_aegyptiacus_meth_sex-specific_age_sesame.csv",h=T) %>%
 filter(Probe_ID..F %in% rhino.ann[!is.na(rhino.ann$SYMBOL),]$CGid) %>% 
 mutate(Fage_test_padj=p.adjust(Pval_Age..F),Mage_test_padj=p.adjust(Pval_Age..F..M),SD="M-biased",spp1="RoAe",spp2="Rousettus aegyptiacus",
        mating="polygynandrous")
SaBi<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/results_files/Saccopteryx_bilineata_meth_sex-specific_age_sesame.csv",h=T) %>%
 filter(Probe_ID..F %in% rhino.ann[!is.na(rhino.ann$SYMBOL),]$CGid) %>% 
 mutate(Fage_test_padj=p.adjust(Pval_Age..F),Mage_test_padj=p.adjust(Pval_Age..F..M),SD="F-biased",spp1="SaBi",spp2="Saccopteryx bilineata",
        mating="harem")
TaBr<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/results_files/Tadarida_brasiliensis_meth_sex-specific_age_sesame.csv",h=T) %>%
 filter(Probe_ID..F %in% rhino.ann[!is.na(rhino.ann$SYMBOL),]$CGid) %>% 
 mutate(Fage_test_padj=p.adjust(Pval_Age..F),Mage_test_padj=p.adjust(Pval_Age..F..M),SD="F-biased",spp1="TaBr",spp2="Tadarida brasiliensis",
        mating="polygynandrous")

allspp<-rbind(CaPe,DeRo,EpFu,PhDi,PhHa,PtHy,PtPu,PtVa,RoAe,SaBi,TaBr)
allspp$sig<-"NS"
allspp[allspp$Pval_Age..F<0.05,]$sig<-"F"
allspp[allspp$Pval_Age..F..M<0.05,]$sig<-"M"
allspp[allspp$Pval_Age..F..M<0.05 &
        allspp$Pval_Age..F<0.05,]$sig<-"F, M"
allspp$sig<-factor(allspp$sig,levels=c("F","M","F, M"))

allspp$XA<-NA
allspp[allspp$Probe_ID..F %in% xlinked,]$XA<-"X-linked"
allspp[allspp$Probe_ID..F %in% autosomal,]$XA<-"autosomal"

colnames(allspp)[c(1:11)]<-c("Fem_Probe_ID","Fem_intercept","Fem_age_est","Fem_intercept_P","Fem_age_P",
                    "Male_Probe_ID","Male_intercept","Male_age_est","Male_intercept_P","Male_age_P",
                    "Spp")

write.csv(allspp,"~/Documents/UMD_new/DNAm-sex_diffs/results_files/AllSpp_sex_specific_age.csv")
