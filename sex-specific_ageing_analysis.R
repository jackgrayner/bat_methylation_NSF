library(tidyverse)
library(RColorBrewer)
library(nlme)
library(dplyr)
library(ggbeeswarm)

#theme for plots
SexPalette<-c("#C74955","#5090bf")
cust.theme<-function(){
  theme_minimal() %+replace%
  theme(panel.grid=element_blank(),
        axis.ticks=element_line(),
        axis.text.y=element_text(face='italic',size=10),
        plot.background=element_rect(fill='white',colour='white'),
        panel.background=element_rect(fill="#fcfbfa",colour='#aaaaaa'))
}

#read in Rhino annotation and results files (allspp)
setwd("~/Documents/UMD_new/DNAm-sex_diffs/")
rhino.ann<-read.csv("Rhinolophus_ferrumequinum.hlrhifer5.HorvathMammalMethylChip40.v1.csv",row.names = 1)
allspp<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/results_files/AllSpp_sex_specific_age.csv",h=T,row.names = 1)#removed PhHa>10yr and NWO EpFu samples, juveniles retained

#merge
colnames(allspp)[1]<-"CGid"
allspp<-merge(allspp,rhino.ann,by="CGid")

#subset sites showing age-associated change in both sexes, in same direction
allspp.sig<-allspp[allspp$Fem_age_P<0.05 & allspp$Male_age_P<0.05 & 
                     (allspp$Fem_age_est * allspp$Male_age_est) > 0,]

#read in clock sites (a posteriori designation of age-associated sites)
batclocks<-read.csv("~/Documents/UMD_new/DNAm-sex_diffs/Bats_Coef.CombAnn.csv",h=T)
#allspp.sig<-allspp[allspp$CGid %in% batclocks[!is.na(batclocks$Coef.BatSkin.Sqrt),]$var,]#for clock sites isntead

#calculate relative male slopes
allspp.sig$rel.age.m<-(abs(allspp.sig$Male_age_est)/abs(allspp.sig$Fem_age_est))

#define male/female-biased slopes
allspp.sig$dir<-NA
allspp.sig[abs(allspp.sig$Male_age_est) < abs(allspp.sig$Fem_age_est),]$dir<-"Female-biased change"
allspp.sig[abs(allspp.sig$Male_age_est) > abs(allspp.sig$Fem_age_est),]$dir<-"Male-biased change"

#relevel consistent with phylogeny
allspp.sig$Spp<-factor(allspp.sig$spp2,levels=c("Eptesicus fuscus","Tadarida brasiliensis","Desmodus rotundus",
                                              "Phyllostomus discolor","Phyllostomus hastatus","Carollia perspicillata",
                                              "Saccopteryx bilineata",
                                              "Rousettus aegyptiacus","Pteropus pumilus", "Pteropus vampyrus",
                                              "Pteropus hypomelanus"))

#plot relative slopes
g.relageing<-ggplot(allspp.sig,aes(y=log10(rel.age.m),x=Spp))+theme_minimal()+
  ylab("Log10 relative rate of male aging")+#facet_grid(.~Fem_age_est>0)+
  scale_x_discrete(limits = rev(levels(allspp.sig$Spp)))+
  theme(axis.text.x=element_text(face='italic',hjust=0.95,vjust=0.2))+
  scale_colour_manual(values=SexPalette)+
  scale_fill_manual(values=SexPalette)+
  geom_quasirandom(size=1,aes(colour=dir))+
  geom_hline(yintercept=0,colour='black',linetype='dashed')+
  geom_boxplot(width=0.25,outlier.shape=NA,coef = 0)+
  theme_minimal()+ 
  theme(panel.grid=element_blank(),
        axis.ticks=element_line(),
        plot.background=element_rect(fill='white',colour='white'),
        panel.background=element_rect(fill="#fcfbfa",colour='#aaaaaa'))+
  theme(legend.position='none',axis.title.y=element_blank(),axis.text.y=element_text(face='italic',hjust=0.95,vjust=0.2))+
  coord_flip()

#ggsave('MF_age_slopes1.png',plot=g.relageing,dpi=600,height=6,width=5)

#plot EWAS
g.relage.ewas<-ggplot(allspp.sig,aes(x=probeStart,y=log10(abs(rel.age.m))))+
  theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank())+
  geom_point(size=0.5,aes(colour=Spp))+
  geom_hline(yintercept=0,linetype='solid',colour='white')+
  facet_grid(.~geneChr,space='free',scales='free',switch='x')+
  ylab("Log10 relative rate of male aging")+
  xlab("")+
  scale_color_brewer(palette = "RdYlBu")+
  cust.theme()+
  theme(panel.background=element_rect(fill="#fcfbfa",colour='white'),
        axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        legend.position='top')+guides(shape = guide_legend(override.aes = list(size = 1)),
                  color = guide_legend(override.aes = list(size = 1)))+
    theme(legend.title = element_blank(), 
          legend.text  = element_text(size = 8),
          legend.key.size = unit(0.5, "lines"))

#ggsave('MF_relage_ewas.png',dpi=600,width=9,height=5,
#       plot=(g.relage.ewas))

#subset sites showing male-biased slope in >5 spp
mbias.age.sites<-allspp.sig[(abs(allspp.sig$Male_age_est)>abs(allspp.sig$Fem_age_est)),] %>% 
  group_by(CGid) %>% summarise(n=n()) %>% filter(n>5)
mbias.age.sites$CGid
mbias.age.sites<-rhino.ann[rhino.ann$CGid %in% mbias.age.sites$CGid,]$SYMBOL
#write.table(file='Mbiased_age_genes.txt',mbias.age.sites,row.names = FALSE,col.names = FALSE,quote=FALSE)
all.list<-unique(rhino.ann$SYMBOL)


library("KEGGREST")
library("clusterProfiler")
library("tidyverse")
library(AnnotationHub)
hub <- AnnotationHub()
organism<-'org.Cf.eg.db'
library(organism, character.only = TRUE)
go.all.mb<-enrichGO(gene = mbias.age.sites, 
                    universe = all.list,
                    keyType = "SYMBOL",
                    OrgDb = organism,
                    ont = "BP",
                    pAdjustMethod = "none",
                    qvalueCutoff = 0.05,
                    pvalueCutoff = 0.05,
                    readable = TRUE)
go.all.mb.plot<-dotplot(go.all.mb)+theme(panel.grid=element_blank(),
                         axis.ticks=element_line(),
                         plot.background=element_rect(fill='white',colour='white'),
                         panel.background=element_rect(fill="#fcfbfa",colour='#aaaaaa'))

ggsave('MBaging_GOs.png',dpi=600,width=7,height=5,
       plot=go.all.mb.plot)


#Run PGLS
library(phytools)
library(dplyr)
library(ggrepel)
phylo<-read.newick(file="bats_species.nwk")
species<-gsub(" ","_",levels(factor(allspp$Spp)))
species[3]<-"Eptesicus_fuscus_dutertreus"
pruned.tree<-drop.tip(phylo,phylo$tip.label[-match(species, phylo$tip.label)])
pruned.tree$tip.label<-gsub("Eptesicus_fuscus_dutertreus","Eptesicus_fuscus",pruned.tree$tip.label)

#png(filename="bat_phylo.png",width=3,height=6,units="in",res=600)
plotTree(pruned.tree,fsize=1,lwd=2,plot=TRUE,ftype="i")
#dev.off()

summary.allspp<-allspp.sig %>% group_by(Spp=gsub(" ","_",Spp)) %>% 
  summarize(med.m.bias=median(rel.age.m))

#summary.allspp<-allspp.sig[!allspp.sig$geneChr==1,] %>% 
#  group_by(Spp=gsub(" ","_",Spp)) %>% 
#  summarize(med.m.bias=median(rel.age.m))

pglsModel.sexbias<-gls(med.m.bias ~ 1, summary.allspp[!is.na(summary.allspp$med.m.bias),], 
               correlation=corBrownian(1, pruned.tree, form = ~Spp))
anova(pglsModel.sexbias)

summary.allspp$Spp
summary.allspp$rel.testes.size<-c(0.00616,0.00359,0.02076,0.0133,0.00285,
                                  0.00408,0.00819,0.00447,0.0493,0.00142,0.00813)
summary.allspp$Slope<-"M-biased"
summary.allspp[log10(summary.allspp$med.m.bias)<0,]$Slope<-"F-biased"
pglsModel.testes<-gls(log10(med.m.bias) ~ log10(rel.testes.size), summary.allspp[!is.na(summary.allspp$rel.testes.size),], 
               correlation=corBrownian(1, pruned.tree, form = ~Spp))
anova(pglsModel.testes)

# what if we just classify spp as polygynous/promiscuous
# summary.allspp$matingsystem<-"polygynous"
# summary.allspp[summary.allspp$Spp %in% c("Rousettus_aegyptiacus","Tadarida_brasiliensis","Eptesicus_fuscus"),]$matingsystem<-"promiscuous"
# pglsModel.msystem<-gls(log10(med.m.bias) ~ matingsystem, summary.allspp[!is.na(summary.allspp$rel.testes.size),], 
#                       correlation=corBrownian(1, pruned.tree, form = ~Spp))
# anova(pglsModel.msystem)

library(ggeffects)
g <- predict_response(pglsModel.testes, terms=c("rel.testes.size"),back_transform = FALSE,ci_level = 0.95) 
g1<-data.frame(g)
g1$med.m.bias<-summary.allspp$med.m.bias

#plot association btwn testes and m.biased slope
g.cor<-ggplot(summary.allspp,aes(x=log10(rel.testes.size),y=log10(med.m.bias)))+
  cust.theme()+
  geom_hline(yintercept=0,colour='red',linetype='dotted')+
  geom_point(aes(colour=Slope),size=2)+
  scale_colour_manual(values=SexPalette)+
  #geom_ribbon(data=g1,aes(x=log10(x),y=predicted,ymin=conf.low,ymax=conf.high,group=1),fill='#cccccc',alpha=0.5)+
  geom_line(data=g1,aes(x=log10(x),y=predicted),colour="#1f4f70")+
  geom_text_repel(aes(label=gsub("_"," ",Spp),colour=Slope),fontface = "italic",
                  min.segment.length = 0.0001,alpha=1,show.legend=FALSE,size=3.5)+
  scale_x_continuous(breaks=seq(-3,-1.2,by=0.2))+
  scale_y_continuous(breaks=seq(-0.2,0.3,by=0.1),limits=c(-0.15,0.3))+
  xlab("Log10 relative testes size")+
  ylab("Log10 relate male slope")+
  annotate(geom='text',label="y = -0.39 - 0.22x",x = -1.5,y=0.30,
           size=4,colour="#1f4f70",alpha=0.5)

ggsave("./final_figs/rel_Mageing_testes.png",plot=g.cor,dpi=600,height=5,width=6.3)


