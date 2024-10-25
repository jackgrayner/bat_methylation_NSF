#process metadata and filter methylation data

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

methmeta<-read.csv("~/Desktop/umd/nsf_methylation/all_meth_sampinfo.csv",h=T)
methmeta<-methmeta[!duplicated(methmeta$SID),]
#keep species with >= 8 of each sex, skin, and non-hibernating

methmeta.m<-methmeta[methmeta$Sex=="M",]
methmeta.f<-methmeta[methmeta$Sex=="F",]
methsum.m<-(data.frame(N.M=summary(factor(methmeta.m$Species))))
methsum.f<-(data.frame(N.F=summary(factor(methmeta.f$Species))))

methsum<-merge(methsum.m,methsum.f,by="row.names")
methsum<-methsum[methsum$N.M>7 & methsum$N.F>7,]
final.species<-methsum$Row.names
methmeta.final<-methmeta[methmeta$Species %in% final.species,]
methmeta.final<-methmeta.final[!is.na(methmeta.final$Age),]

# plot known age distributions
# ggplot(methmeta.final,aes(x=Age,fill=Sex))+geom_density()+
#   facet_grid(Sex~Species,space='free_x',scales='free_x')+
#   theme_minimal()+theme(panel.grid=element_blank(),
#                         plot.background=element_rect(fill='white'),
#                         strip.text=element_text(size=9.5,face='italic',family='sans'),
#                         axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),
#                         axis.ticks.x=element_line(),axis.ticks.y=element_line())+
#   scale_fill_manual(values=SexPalette1)+theme(panel.grid=element_blank())

#do some more filtering
#get rid of non-skin samples
methmeta.final<-methmeta.final[methmeta.final$Tissue=="Skin",]
#remove unaged individuals
methmeta.final<-methmeta.final[!is.na(methmeta.final$Age),]

#save files
write.csv(methmeta.final,file='~/Desktop/umd/nsf_methylation/methmeta_final.csv',row.names=TRUE,quote=FALSE)
rm(list=ls())
