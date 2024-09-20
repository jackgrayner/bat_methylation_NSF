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
#remove JHU samples of Carollia due to small sample
methmeta.final<-methmeta.final[!(methmeta.final$Lab %in% "JHU" & methmeta.final$Species=="Carollia perspicillata"),]
#get rid of non-skin samples
methmeta.final<-methmeta.final[methmeta.final$Tissue=="Skin",]
#remove unaged individuals
methmeta.final<-methmeta.final[!is.na(methmeta.final$Age),]

#read in methylation data
batdat<-readRDS("~/Desktop/umd/nsf_methylation/BatTib_25Oct23.rds")
#keep samples in final metadata file
batdat<-batdat[batdat$SID %in% methmeta.final$SID,]
#remove duplicates
batdat.dupes<-batdat.dupes[duplicated(batdat.dupes$SID),]$Basename
#in all cases, it looks like the first record in batdat df is the record in methmeta.final
#therefore, we can just remove the second copy
batdat<-batdat[!batdat$Basename %in% batdat.dupes,]

batclocks<-read.csv("~/Desktop/umd/nsf_methylation/Bats_Coef.CombAnn.csv",h=T)

#remove individuals with low (<80) confidence known ages
lowconf<-read.csv("low_confidence_samples.csv",h=T)
batdat<-batdat[!batdat$Basename %in% lowconf$Basename,]
methmeta.final<-methmeta.final[methmeta.final$SID %in% batdat$SID,]

#output summary of remaining samples
print(methmeta.final %>% group_by(Species,Sex) %>% dplyr::summarize(n=n()),n=Inf)

#build a function to estimate ages for each species
fun_sqrt.inv <- function(y, ...) y^2-1
ageEst <- function(species, clockname){
 #retrieve clock coefficients
 species.clock<-data.frame(batclocks[,c("var",clockname)])
 species.clock<-data.frame(species.clock[!is.na(species.clock[,2]),])
 
 #rearrange meth data for species
 batdat.temp<-batdat[batdat$Species==species & batdat$tissue=="Skin",]
 batdat.temp<-batdat.temp[batdat.temp$site %in% species.clock$var,c("Basename","SID","tissue","site","m")]
 batdat.temp<-spread(batdat.temp, key = site, value = m)
 
 #predict ages
 species1<-gsub(" ","_",species)
 return(data.frame(Basename=batdat.temp$Basename,SID=batdat.temp$SID,
                   est.age=as.numeric(predictClockSimple(species.clock,batdat.temp[,-c(1,2,3)],fun_sqrt.inv))))
 rm(batdat.temp)
}

colnames(batclocks)
final.species
all.est.ages<-rbind(
 ageEst("Carollia perspicillata","Coef.BatSkin.Sqrt"),#all bat looks ok
 ageEst("Desmodus rotundus","Coef.BatSkin.Sqrt"),#all bat looks ok
 ageEst("Eptesicus fuscus","Coef.BatSkin.Sqrt"),#not ideal
 ageEst("Phyllostomus discolor","Coef.Phyllostomus.discolor_Skin.Sqrt"),#all bat not ideal, try species
 ageEst("Phyllostomus hastatus","Coef.BatSkin.Sqrt"),#all bat looks ok
 ageEst("Pteropus hypomelanus","Coef.BatSkin.Sqrt"),#all bat looks ok
 ageEst("Pteropus pumilus","Coef.BatSkin.Sqrt"),#all bat looks ok
 ageEst("Pteropus vampyrus","Coef.BatSkin.Sqrt"),#all bat looks ok
 ageEst("Rousettus aegyptiacus","Coef.BatSkin.Sqrt"),#all bat looks ok
 ageEst("Saccopteryx bilineata","Coef.Saccopteryx.bilineata_Skin.Sqrt"),#all bat not ideal
 ageEst("Tadarida brasiliensis","Coef.Tadarida.brasiliensis_Skin.Sqrt"))#all bat terrible, try species

methmeta.final<-merge(methmeta.final,all.est.ages,by="SID")

#relevel
methmeta.final$Species<-factor(methmeta.final$Species,levels=c("Desmodus rotundus","Phyllostomus discolor","Phyllostomus hastatus","Carollia perspicillata",
                                                               "Rousettus aegyptiacus","Pteropus pumilus","Pteropus vampyrus","Pteropus hypomelanus",
                                                               "Eptesicus fuscus","Tadarida brasiliensis",NA,"Saccopteryx bilineata"))

#plot age vs estimated age
g.agecor<-ggplot(methmeta.final,aes(x=Age,y=est.age,colour=Sex))+
 geom_abline(intercept=0,slope=1,colour='black',linewidth=0.5)+
 geom_point(size=1.5,aes(colour=Sex),alpha=0.5)+
 geom_smooth(method='lm',se=TRUE,aes(fill=Sex),alpha=0.25)+
 facet_wrap(.~Species,scales='free')+
 scale_colour_manual(values=SexPalette1[c(1,2)])+
 scale_fill_manual(values=SexPalette1[c(1,2)])+
 theme_minimal()+theme(panel.grid=element_blank(),
                       panel.background=element_rect(fill='#fcf7f2',colour='#cccccc'),
                       plot.background=element_rect(fill='white',colour='white'),
                       strip.text=element_text(size=9.5,face='italic',family='sans'),
                       axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),
                       axis.ticks.x=element_line(),axis.ticks.y=element_line())+
 ylab("Estimated age")+xlab("Recorded age")+
 scale_x_continuous(breaks=seq(0,25,3))+scale_y_continuous(breaks=seq(0,25,3))+
 labs(caption="Species-specific clocks: PhDi, SaBi, TaBr\n
       All-bat clock: CaPe, DeRo, PhHa, PtHy, PtPu, PtVa, RoAe")

ggsave('~/Desktop/umd/nsf_methylation/Age_correlations.png',plot=g.agecor,dpi=600,height=6.5,width=7.75)

#plot age residuals
resid.df<-data.frame(spp=NA,sex=NA,resid=NA,age=NA)
for (spp.temp in levels(factor(methmeta.final$Species))){
 spp2<-gsub("\\s","_",spp.temp)
 assign('temp.df',data.frame(spp=spp.temp,
                             sex=methmeta.final[methmeta.final$Species==spp.temp & !is.na(methmeta.final$Age),]$Sex,
                             age=methmeta.final[methmeta.final$Species==spp.temp & !is.na(methmeta.final$Age),]$Age,
                             resid=resid(lm(est.age ~ Age,data=methmeta.final[methmeta.final$Species==spp.temp,]))))
 resid.df<-rbind(resid.df,temp.df)
}

resid.df$spp<-factor(resid.df$spp,levels=c("Desmodus rotundus","Phyllostomus discolor","Phyllostomus hastatus","Carollia perspicillata",
                                           "Rousettus aegyptiacus","Pteropus pumilus","Pteropus vampyrus","Pteropus hypomelanus",
                                           "Eptesicus fuscus","Tadarida brasiliensis",NA,"Saccopteryx bilineata"))

library(ggbeeswarm)
#plot residual from age vs est.age regression
g.resid<-ggplot(resid.df[!is.na(resid.df$spp),],aes(x=sex,y=resid,colour=sex,fill=sex))+
 geom_hline(yintercept=0,linetype='dashed')+geom_quasirandom(size=0.75)+
 stat_summary(fun.data='mean_sdl',fun.args = list(mult = 1),fill='white',colour='black',size=0.5,shape=21)+
 facet_wrap(.~spp,scales='free')+theme_minimal()+
 scale_colour_manual(values=SexPalette1)+
 scale_fill_manual(values=SexPalette1)+
 xlab("Sex")+ylab('residuals (est.age ~ known_age)')+
 theme(panel.grid=element_blank(),panel.background=element_rect(fill='#fcf7f2',colour='#cccccc'),
       plot.background=element_rect(fill='white',colour='white'),
       strip.text=element_text(size=9.5,face='italic',family='sans'))

ggsave('~/Desktop/umd/nsf_methylation/Age_residuals.png',plot=g.resid,dpi=600,height=6.5,width=7.75)

#save files
write.csv(methmeta.final,file='~/Desktop/umd/nsf_methylation/methmeta_final.csv',row.names=TRUE,quote=FALSE)
batdat<-batdat[batdat$Basename %in% methmeta.final$Basename,]
batclocks<-read.csv("~/Desktop/umd/nsf_methylation/Bats_Coef.CombAnn.csv",h=T)
batdat<-batdat[!batdat$site %in% batclocks$var,]
saveRDS(batdat,file="batdat_sexspp.rds")
rm(list=ls())
