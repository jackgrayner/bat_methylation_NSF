library(pheatmap)
library(patchwork)
library(ggbeeswarm)
library(viridis)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(dplyr)

SexPalette<-c("#C74955","#5090bf","purple")
cust.theme<-function(){
  theme_minimal() %+replace%
    theme(panel.grid=element_blank(),
          axis.ticks=element_line(),
          axis.text.y=element_text(face='italic',size=10),
          plot.background=element_rect(fill='white',colour='white'),
          panel.background=element_rect(fill="#fcfbfa",colour='#aaaaaa'))
}
addSmallLegend <- function(myPlot, pointSize = 1, textSize = 8, spaceLegend = 0.5) {
  myPlot + guides(shape = guide_legend(override.aes = list(size = pointSize)),
                  color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_blank(), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

setwd("~/Documents/UMD_new/DNAm-sex_diffs/")
rhino.ann<-read.csv("Rhinolophus_ferrumequinum.hlrhifer5.HorvathMammalMethylChip40.v1.csv",row.names = 1) %>%
  filter(!is.na(SYMBOL))
colnames(rhino.ann)[1]<-"Probe_ID"
allspp<-read.csv("./results_files/AllSpp_sex_age.csv",row.names=1)

allspp$spp2<-factor(allspp$spp2,levels=c("Eptesicus fuscus","Tadarida brasiliensis","Desmodus rotundus",
                                         "Phyllostomus discolor","Phyllostomus hastatus","Carollia perspicillata",
                                         "Saccopteryx bilineata",
                                         "Rousettus aegyptiacus","Pteropus pumilus","Pteropus vampyrus","Pteropus hypomelanus"))

allspp<-merge(allspp,rhino.ann,by="Probe_ID")

#add sig details
allspp$sigsb<-"NSB"
allspp[allspp$Pval_SexMadj < 0.05,]$sigsb<-"SB"

allspp$sigab<-"NAB"
allspp[allspp$Pval_Ageadj < 0.05,]$sigab<-"AB"

allspp$sig<-NA
allspp[allspp$Pval_SexMadj<0.05 & allspp$XA=="autosomal" & allspp$Est_SexM<0,]$sig<-"F-biased, autosomal"
allspp[allspp$Pval_SexMadj<0.05 & allspp$XA=="autosomal" & allspp$Est_SexM>0,]$sig<-"M-biased, autosomal"
allspp[allspp$Pval_SexMadj<0.05 & allspp$XA=="X-linked" & allspp$Est_SexM<0,]$sig<-"F-biased, X-linked"
allspp[allspp$Pval_SexMadj<0.05 & allspp$XA=="X-linked" & allspp$Est_SexM>0,]$sig<-"M-biased, X-linked"


#test overrepresentation of X among sbDMPs
allspp.sexsig.sum<-allspp[allspp$sigsb=="SB",] %>% group_by(Probe_ID,XA,sigsb) %>% summarize(n=n())
allspp.all.sum<-allspp[allspp$sigsb=="NSB",] %>% group_by(Probe_ID,XA,sigsb) %>% summarize(n=n()) %>% filter(n==11)
sum1<-rbind(allspp.sexsig.sum,allspp.all.sum)
fisher.test(table(sum1$XA, sum1$sigsb))
nrow(allspp[allspp$sigsb=="SB",] %>% group_by(Probe_ID,XA,sigab) %>% summarize(n=n()) %>% filter(n==11))

#test overrepresentation of X among ageDMPs
allspp$sigab<-factor(allspp$sigab)
allspp$sigab<-relevel(allspp$sigab,ref="NAB")
allspp.agesig.sum<-allspp[allspp$sigab=="AB",] %>% group_by(Probe_ID,XA,sigab) %>% summarize(n=n())
allspp.nonsig.sum<-allspp[allspp$sigab=="NAB",] %>% group_by(Probe_ID,XA,sigab) %>% summarize(n=n()) %>% filter(n==11)
sum2<-rbind(allspp.agesig.sum,allspp.nonsig.sum)
fisher.test(table(sum2$XA, sum2$sigab))

#plot gene feature proportions
SexPalette2<-c("#944149","#C74955","#195582","#5090bf")
allspp$feature<-"Other feature"
allspp[allspp$main_Categories=="Promoter",]$feature<-"Promoter"
g.num.sb<-ggplot(allspp[!is.na(allspp$sig),],aes(x=spp2,fill=sig))+
  #geom_bar(colour='black')+
  scale_x_discrete(expand=c(0,0))+scale_x_discrete(limits=rev)+
  geom_bar(colour='black',position='fill',alpha=0.75)+
  theme_minimal()+scale_fill_manual(values=SexPalette2,name="Sex-bias")+
  theme(panel.grid=element_blank(),axis.ticks.y = element_line(),axis.text.y=element_text(face='italic'),
        panel.background = element_rect(fill='#fcfbfa',colour='#aaaaaa'),legend.title=element_blank())+
  labs(y="Proportion")+labs(tag="B")+xlab(" ")+coord_flip()+
  ylab("Proportion")+facet_grid(.~feature)

g.feature.sb<-ggplot(allspp[!is.na(allspp$sig),],aes(x=spp2,fill=main_Categories))+
  geom_bar(colour='black',position='fill',alpha=0.75)+
  scale_x_discrete(expand=c(0,0))+scale_x_discrete(limits=rev)+
  theme_minimal()+scale_fill_viridis(discrete=TRUE,name="Feature")+
  theme(axis.title.x=element_blank(),axis.text.y=element_text(face='italic'),
        panel.grid=element_blank(),legend.position='right',axis.ticks.y = element_line(),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill='#fcfbfa',colour='#aaaaaa'),legend.title=element_blank())+
  labs(y="Proportion")+labs(tag="C")+coord_flip()+xlab(" ")

g.feature.all<-ggplot(allspp[allspp$spp1=="CaPe",],
                      aes(x=factor("Reference"),fill=main_Categories))+
  geom_bar(colour='black',position='fill',alpha=0.75)+scale_x_discrete(limits=rev)+
  theme_minimal()+scale_fill_viridis(discrete=TRUE)+
  theme(axis.title.y=element_blank(),axis.ticks.x = element_line(),legend.title=element_blank(),
        panel.grid=element_blank(),legend.position='none',
        plot.background = element_rect(fill='white',colour='white'))+
  labs(y="Proportion")+coord_flip()+xlab("")+ylab("Proportion")

ggsave('SB_DMPs_XA_features.png',dpi=600,width=7,height=5,
       plot=(g.num.sb/g.feature.sb/g.feature.all)+plot_layout(heights=c(1,1,0.1)))

#test enrichment of gene features --SB
fish.results<-data.frame(spp=NA,group=NA,p=NA,odds=NA)
for (spp1 in levels(factor(allspp$spp1))){
  for (i in c(1:7)){
    t<-table(allspp[allspp$spp1==spp1,]$main_Categories,
             allspp[allspp$spp1==spp1,]$sigsb)
    fish.results<-rbind(fish.results,
                        data.frame(spp=spp1,group=rownames(t)[i],
                                   p=fisher.test(matrix(c(t[i,2],sum(t[-i,2]), t[i,1], sum(t[-i,1])), nrow = 2, byrow = TRUE))$p.value,
                                   odds=fisher.test(matrix(c(t[i,2],sum(t[-i,2]), t[i,1], sum(t[-i,1])), nrow = 2, byrow = TRUE))$estimate)
    )
  }
}
fish.results$padj<-p.adjust(fish.results$p)
fish.results$sig<-"N"
fish.results[fish.results$padj<0.05 & !is.na(fish.results$padj) & fish.results$odds>1,]$sig<-"Y_up"
fish.results[fish.results$padj<0.05 & !is.na(fish.results$padj) & fish.results$odds<1,]$sig<-"Y_dn"
fish.results[fish.results$sig=="Y_up",] %>% group_by(group) %>% summarize(n=n())
fish.results[fish.results$sig=="Y_dn",] %>% group_by(group) %>% summarize(n=n())

#test enrichment of gene features --AB
fish.results<-data.frame(spp=NA,group=NA,p=NA,odds=NA)
allspp$sigab<-factor(allspp$sigab)
allspp$sigab<-relevel(allspp$sigab,ref="NAB")
for (spp1 in levels(factor(allspp$spp1))){
  for (i in c(1:7)){
    t<-table(allspp[allspp$spp1==spp1,]$main_Categories,
             allspp[allspp$spp1==spp1,]$sigab)
    fish.results<-rbind(fish.results,
                        data.frame(spp=spp1,group=rownames(t)[i],
                                   p=fisher.test(matrix(c(t[i,2],sum(t[-i,2]), t[i,1], sum(t[-i,1])), nrow = 2, byrow = TRUE))$p.value,
                                   odds=fisher.test(matrix(c(t[i,2],sum(t[-i,2]), t[i,1], sum(t[-i,1])), nrow = 2, byrow = TRUE))$estimate)
    )
  }
}

fish.results$padj<-p.adjust(fish.results$p)
fish.results$sig<-"N"
fish.results[fish.results$padj<0.05 & !is.na(fish.results$padj) & fish.results$odds>1,]$sig<-"Y_up"
fish.results[fish.results$padj<0.05 & !is.na(fish.results$padj) & fish.results$odds<1,]$sig<-"Y_dn"
fish.results[fish.results$sig=="Y_up",] %>% group_by(group) %>% summarize(n=n())
fish.results[fish.results$sig=="Y_dn",] %>% group_by(group) %>% summarize(n=n())

#plot enrichment of gene features
library(ggrepel)
ggplot(fish.results,aes(x=odds,y=-log10(padj),colour=group))+
  geom_point(alpha=1)+
  geom_label_repel(data=fish.results[fish.results$padj<0.05,],aes(label=spp),show.legend = FALSE)+
  geom_vline(xintercept=1,linetype='dashed',alpha=0.5)+
  geom_hline(yintercept=-log10(0.05),linetype='dashed',alpha=0.5)+
  cust.theme()


#M or F biased in 5 or more
allspp$sig<-NA
allspp[allspp$Pval_SexMadj<0.05 & allspp$Est_SexM<0,]$sig<-"F-biased"
allspp[allspp$Pval_SexMadj<0.05 & allspp$Est_SexM>0,]$sig<-"M-biased"

allspp.f.biased<-allspp[allspp$sig=="F-biased",] %>% group_by(Probe_ID,XA,sig) %>%
  summarize(n=n()) %>% filter(n>5)
allspp.m.biased<-allspp[allspp$sig=="M-biased",] %>% group_by(Probe_ID,XA,sig) %>%
  summarize(n=n()) %>% filter(n>5)
sb.list<-unique(allspp[allspp$Pval_SexMadj<0.05,]
                $SYMBOL)

write.table(file='f_biased_sites.txt',allspp.f.biased$Probe_ID,quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(file='m_biased_sites.txt',allspp.m.biased$Probe_ID,quote=FALSE,row.names=FALSE,col.names=FALSE)


library("clusterProfiler")
organism<-'org.Cf.eg.db'
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
fbiased.go<-enrichGO(gene = rhino.ann[rhino.ann$Probe_ID %in% allspp.f.biased$Probe_ID,]$SYMBOL, 
                    universe = unique(rhino.ann$SYMBOL),
                    keyType = "SYMBOL",
                    OrgDb = organism,
                    ont = "BP",
                    pAdjustMethod = "none",
                    qvalueCutoff = 0.2,
                    pvalueCutoff = 0.5,
                    readable = TRUE)
dotplot(fbiased.go)
mbiased.go<-enrichGO(gene = rhino.ann[rhino.ann$Probe_ID %in% allspp.m.biased$Probe_ID,]$SYMBOL, 
                     universe = unique(rhino.ann$SYMBOL),
                     keyType = "SYMBOL",
                     OrgDb = organism,
                     ont = "BP",
                     pAdjustMethod = "none",
                     qvalueCutoff = 0.2,
                     pvalueCutoff = 0.5,
                     readable = TRUE)
dotplot(mbiased.go)
go.all.sb<-enrichGO(gene = sb.list, 
                       universe = unique(rhino.ann$SYMBOL),
                       keyType = "SYMBOL",
                       OrgDb = organism,
                       ont = "all",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       pvalueCutoff = 1,
                       readable = TRUE)
go.all.sb.plot<-dotplot(go.all.sb)+theme(panel.grid=element_blank(),
                         axis.ticks=element_line(),
                         plot.background=element_rect(fill='white',colour='white'),
                         panel.background=element_rect(fill="#fcfbfa",colour='#aaaaaa'))

ggsave('SB_GOs.png',dpi=600,width=7,height=5,
       plot=go.all.sb.plot)


g.ewas<-ggplot()+
  theme_bw()+
  geom_point(data=allspp,
             aes(x=probeStart,y=-log10(Pval_SexM)),
             size=0.35,alpha=1,colour='#aaaaaa')+
  geom_point(data=allspp[(allspp$sigsb=="SB"),],
             aes(x=probeStart,y=-log10(Pval_SexM),colour=spp2),
             size=0.35,alpha=0.75)+
  facet_grid(.~geneChr,scales='free',space='free',switch='x')+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        panel.grid=element_blank(),legend.position='right',axis.title.x=element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(fill='#fcfbfa',colour='#fcfbfa'),
        panel.background = element_rect(fill='#fcfbfa',colour='#fcfbfa'),
        panel.margin.x=unit(0.02, "lines") , panel.margin.y=unit(0,"lines"),
        strip.text.x=element_text(size=6,angle = 90),legend.title=element_blank(),
        axis.text.y=element_text(size=8),legend.text=element_text(face='italic',size=8))+
  guides(colour = guide_legend(override.aes = list(size=1)))+
  labs(tag="A")+ylab("-log10(P_sex)")+scale_color_brewer(palette = "RdYlBu")
ggsave('sex_ewas.png',dpi=600,width=7,height=2.5,
       plot=addSmallLegend(g.ewas))


#plot gene feature proportions
SexPalette2<-c("#944149","#C74955","#195582","#5090bf")
allspp$feature<-"Other feature"
allspp[allspp$main_Categories=="Promoter",]$feature<-"Promoter"
g.num.sb<-ggplot(allspp[!is.na(allspp$sig),],aes(x=spp2,fill=sig))+
  #geom_bar(colour='black')+
  scale_x_discrete(expand=c(0,0))+scale_x_discrete(limits=rev)+
  geom_bar(colour='black',position='fill',alpha=0.75)+
  theme_minimal()+scale_fill_manual(values=SexPalette2,name="Sex-bias")+
  theme(panel.grid=element_blank(),axis.ticks.y = element_line(),axis.text.y=element_text(face='italic'),
        panel.background = element_rect(fill='#fcfbfa',colour='#aaaaaa'),legend.title=element_blank())+
  labs(y="Proportion")+labs(tag="B")+xlab(" ")+coord_flip()+
  ylab("Proportion")+facet_grid(.~feature)

g.feature.sb<-ggplot(allspp[!is.na(allspp$sig),],aes(x=spp2,fill=main_Categories))+
  geom_bar(colour='black',position='fill',alpha=0.75)+
  scale_x_discrete(expand=c(0,0))+scale_x_discrete(limits=rev)+
  theme_minimal()+scale_fill_viridis(discrete=TRUE,name="Feature")+
  theme(axis.title.x=element_blank(),axis.text.y=element_text(face='italic'),
        panel.grid=element_blank(),legend.position='right',axis.ticks.y = element_line(),
        axis.text.x=element_blank(),
        panel.background = element_rect(fill='#fcfbfa',colour='#aaaaaa'),legend.title=element_blank())+
  labs(y="Proportion")+labs(tag="C")+coord_flip()+xlab(" ")

g.feature.all<-ggplot(allspp[allspp$spp1=="CaPe",],
                      aes(x=factor("Reference"),fill=main_Categories))+
  geom_bar(colour='black',position='fill',alpha=0.75)+scale_x_discrete(limits=rev)+
  theme_minimal()+scale_fill_viridis(discrete=TRUE)+
  theme(axis.title.y=element_blank(),axis.ticks.x = element_line(),legend.title=element_blank(),
        panel.grid=element_blank(),legend.position='none',
        plot.background = element_rect(fill='white',colour='white'))+
  labs(y="Proportion")+coord_flip()+xlab("")+ylab("Proportion")

ggsave('SB_DMPs_XA_features.png',dpi=600,width=7,height=5,
       plot=(g.num.sb/g.feature.sb/g.feature.all)+plot_layout(heights=c(1,1,0.1)))


#test representation of androgen-sensitive probes
andros<-read.csv("elife-64932-supp3-v2.csv",h=T)

#SEX
allspp$andros<-"NAS"
allspp[allspp$Probe_ID %in% andros$CpG.ID,]$andros<-"AS"
allspp.sexsig1.sum<-allspp[allspp$Pval_SexMadj<0.05,] %>%
  group_by(Probe_ID,XA,sigsb,andros) %>% summarize(n=n()) %>% filter(n>0)
allspp.nosexsig1.sum<-allspp[allspp$Pval_SexMadj>0.05,] %>%
  group_by(Probe_ID,XA,sigsb,andros) %>% summarize(n=n()) %>% filter(n==11)
allspp.sexsig1.sum<-rbind(allspp.sexsig1.sum,allspp.nosexsig1.sum)

table(allspp.sexsig1.sum$sigsb,allspp.sexsig1.sum$andros)
fisher.test(table(allspp.sexsig1.sum$sigsb,allspp.sexsig1.sum$andros))

#AGE
allspp$andros<-"NAS"
allspp[allspp$Probe_ID %in% andros$CpG.ID,]$andros<-"AS"
allspp.agesig1.sum<-allspp[allspp$Pval_Ageadj<0.05,] %>%
  group_by(Probe_ID,XA,sigab,andros) %>% summarize(n=n()) %>% filter(n>0)
allspp.noagesig1.sum<-allspp[allspp$Pval_Ageadj>0.05,] %>%
  group_by(Probe_ID,XA,sigab,andros) %>% summarize(n=n()) %>% filter(n==11)
allspp.agesig1.sum<-rbind(allspp.agesig1.sum,allspp.noagesig1.sum)

table(allspp.agesig1.sum$sigab,allspp.agesig1.sum$andros)
fisher.test(table(allspp.agesig1.sum$sigab,allspp.agesig1.sum$andros))



#plot heatmaps - this gets messy. sorry

#sex

spp.list<-c("Eptesicus fuscus","Tadarida brasiliensis","Desmodus rotundus",
            "Phyllostomus discolor","Phyllostomus hastatus","Carollia perspicillata",
            "Saccopteryx bilineata","Rousettus aegyptiacus","Pteropus pumilus","Pteropus vampyrus","Pteropus hypomelanus")
allspp$sig1<-"NS"
allspp[allspp$Pval_SexMadj<0.05 & allspp$Est_SexM>0,]$sig1<-"sig.MB"
allspp[allspp$Pval_SexMadj<0.05 & allspp$Est_SexM<0,]$sig1<-"sig.FB"

allspp.a<-allspp[allspp$XA=="autosomal",]
allspp.wide.sex<-data.frame(
  EpFu=allspp.a[allspp.a$Species==spp.list[1],"sig1"],
  TaBr=allspp.a[allspp.a$Species==spp.list[2],"sig1"],
  DeRo=allspp.a[allspp.a$Species==spp.list[3],"sig1"],
  PhDi=allspp.a[allspp.a$Species==spp.list[4],"sig1"],
  PhHa=allspp.a[allspp.a$Species==spp.list[5],"sig1"],
  CaPe=allspp.a[allspp.a$Species==spp.list[6],"sig1"],
  SaBi=allspp.a[allspp.a$Species==spp.list[7],"sig1"],
  RoAe=allspp.a[allspp.a$Species==spp.list[8],"sig1"],
  PtPu=allspp.a[allspp.a$Species==spp.list[9],"sig1"],
  PtVa=allspp.a[allspp.a$Species==spp.list[10],"sig1"],
  PtHy=allspp.a[allspp.a$Species==spp.list[11],"sig1"])
rownames(allspp.wide.sex)<-allspp.a[allspp.a$Species=="Eptesicus fuscus",]$Probe_ID
sexdiff.matrix<-matrix(nrow=length(spp.list),
                       ncol=length(spp.list))
rownames(sexdiff.matrix)<-c("E.fuscus","T.brasiliensis","D.rotundus",
  "P.discolor","P.hastatus","C.perspicillata",
  "S.bilineata","R.aegyptiacus","Pt.pumilus","Pt.vampyrus","Pt.hypomelanus")
colnames(sexdiff.matrix)<-rownames(sexdiff.matrix)
for (x in c(1:11)){
  for (y in c(1:11)){
    sexdiff.matrix[x,y]<-nrow(allspp.wide.sex[!allspp.wide.sex[,x]=="NS" &
                                                !allspp.wide.sex[,y]=="NS"&
                                                allspp.wide.sex[,x]==allspp.wide.sex[,y],])/
      min(c(nrow(allspp.wide.sex[!allspp.wide.sex[,x]=="NS",]),
            nrow(allspp.wide.sex[!allspp.wide.sex[,y]=="NS",])))
  }
}
ph.sex.a<-pheatmap((sexdiff.matrix),cluster_rows = FALSE,cluster_cols = FALSE,show_rownames=TRUE,
                   main="Sex-biased, autosomal",
                   border_color = '#444444',breaks=seq(0,1, length.out=10),scale='none',
                   color=colorRampPalette(c("#3a7ca6","#e0c653", "#d13434"))(10)#  begin = 0, end = 1, option = "viridis")#
)

allspp.x<-allspp[allspp$XA=="X-linked",]
allspp.wide.sex<-data.frame(
  EpFu=allspp.x[allspp.x$Species==spp.list[1],"sig1"],
  TaBr=allspp.x[allspp.x$Species==spp.list[2],"sig1"],
  DeRo=allspp.x[allspp.x$Species==spp.list[3],"sig1"],
  PhDi=allspp.x[allspp.x$Species==spp.list[4],"sig1"],
  PhHa=allspp.x[allspp.x$Species==spp.list[5],"sig1"],
  CaPe=allspp.x[allspp.x$Species==spp.list[6],"sig1"],
  SaBi=allspp.x[allspp.x$Species==spp.list[7],"sig1"],
  RoAe=allspp.x[allspp.x$Species==spp.list[8],"sig1"],
  PtPu=allspp.x[allspp.x$Species==spp.list[9],"sig1"],
  PtVa=allspp.x[allspp.x$Species==spp.list[10],"sig1"],
  PtHy=allspp.x[allspp.x$Species==spp.list[11],"sig1"])
rownames(allspp.wide.sex)<-allspp.x[allspp.x$Species=="Eptesicus fuscus",]$Probe_ID
sexdiff.matrix<-matrix(nrow=length(spp.list),
                       ncol=length(spp.list))
rownames(sexdiff.matrix)<-c("E.fuscus","T.brasiliensis","D.rotundus",
                            "P.discolor","P.hastatus","C.perspicillata",
                            "S.bilineata","R.aegyptiacus","Pt.pumilus","Pt.vampyrus","Pt.hypomelanus")
colnames(sexdiff.matrix)<-rownames(sexdiff.matrix)
for (x in c(1:11)){
  for (y in c(1:11)){
    sexdiff.matrix[x,y]<-nrow(allspp.wide.sex[!allspp.wide.sex[,x]=="NS" &
                                                !allspp.wide.sex[,y]=="NS"&
                                                allspp.wide.sex[,x]==allspp.wide.sex[,y],])/
      min(c(nrow(allspp.wide.sex[!allspp.wide.sex[,x]=="NS",]),
            nrow(allspp.wide.sex[!allspp.wide.sex[,y]=="NS",])))
  }
}
ph.sex.x<-pheatmap((sexdiff.matrix),cluster_rows = FALSE,cluster_cols = FALSE,show_rownames=TRUE,
                   main="Sex-biased, X-linked",
                   border_color = '#444444',breaks=seq(0,1, length.out=10),scale='none',
                   color=colorRampPalette(c("#3a7ca6","#e0c653", "#d13434"))(10)#  begin = 0, end = 1, option = "viridis")#
)


#age

allspp$sig1<-"NS"
allspp[allspp$Pval_Ageadj<0.05 & allspp$Est_Age>0,]$sig1<-"sig.hyper"
allspp[allspp$Pval_Ageadj<0.05 & allspp$Est_Age<0,]$sig1<-"sig.hypo"
allspp.a<-allspp[allspp$XA=="autosomal",]
allspp.wide.age<-data.frame(
  EpFu=allspp.a[allspp.a$Species==spp.list[1],"sig1"],
  TaBr=allspp.a[allspp.a$Species==spp.list[2],"sig1"],
  DeRo=allspp.a[allspp.a$Species==spp.list[3],"sig1"],
  PhDi=allspp.a[allspp.a$Species==spp.list[4],"sig1"],
  PhHa=allspp.a[allspp.a$Species==spp.list[5],"sig1"],
  CaPe=allspp.a[allspp.a$Species==spp.list[6],"sig1"],
  SaBi=allspp.a[allspp.a$Species==spp.list[7],"sig1"],
  RoAe=allspp.a[allspp.a$Species==spp.list[8],"sig1"],
  PtPu=allspp.a[allspp.a$Species==spp.list[9],"sig1"],
  PtVa=allspp.a[allspp.a$Species==spp.list[10],"sig1"],
  PtHy=allspp.a[allspp.a$Species==spp.list[11],"sig1"])
rownames(allspp.wide.age)<-allspp.a[allspp.a$Species=="Eptesicus fuscus",]$Probe_ID
agediff.matrix<-matrix(nrow=length(spp.list),
                       ncol=length(spp.list))
rownames(agediff.matrix)<-c("E.fuscus","T.brasiliensis","D.rotundus",
                            "P.discolor","P.hastatus","C.perspicillata",
                            "S.bilineata","R.aegyptiacus","Pt.pumilus","Pt.vampyrus","Pt.hypomelanus")
colnames(agediff.matrix)<-rownames(agediff.matrix)
for (x in c(1:11)){
  for (y in c(1:11)){
    agediff.matrix[x,y]<-nrow(allspp.wide.age[!allspp.wide.age[,x]=="NS" &
                                                !allspp.wide.age[,y]=="NS"&
                                                allspp.wide.age[,x]==allspp.wide.age[,y],])/
      min(c(nrow(allspp.wide.age[!allspp.wide.age[,x]=="NS",]),
            nrow(allspp.wide.age[!allspp.wide.age[,y]=="NS",])))
  }
}
ph.age.a<-pheatmap((agediff.matrix),cluster_rows = FALSE,cluster_cols = FALSE,show_rownames=TRUE,
                   main="Age-associated, autosomal",
                   border_color = '#444444',breaks=seq(0,1, length.out=10),scale='none',
                   color=colorRampPalette(c("#3a7ca6","#e0c653", "#d13434"))(10)#  begin = 0, end = 1, option = "viridis")#
)

allspp.x<-allspp[allspp$XA=="X-linked",]
allspp.wide.age<-data.frame(
  EpFu=allspp.x[allspp.x$Species==spp.list[1],"sig1"],
  TaBr=allspp.x[allspp.x$Species==spp.list[2],"sig1"],
  DeRo=allspp.x[allspp.x$Species==spp.list[3],"sig1"],
  PhDi=allspp.x[allspp.x$Species==spp.list[4],"sig1"],
  PhHa=allspp.x[allspp.x$Species==spp.list[5],"sig1"],
  CaPe=allspp.x[allspp.x$Species==spp.list[6],"sig1"],
  SaBi=allspp.x[allspp.x$Species==spp.list[7],"sig1"],
  RoAe=allspp.x[allspp.x$Species==spp.list[8],"sig1"],
  PtPu=allspp.x[allspp.x$Species==spp.list[9],"sig1"],
  PtVa=allspp.x[allspp.x$Species==spp.list[10],"sig1"],
  PtHy=allspp.x[allspp.x$Species==spp.list[11],"sig1"])
rownames(allspp.wide.age)<-allspp.x[allspp.x$Species=="Eptesicus fuscus",]$Probe_ID
agediff.matrix<-matrix(nrow=length(spp.list),
                       ncol=length(spp.list))
rownames(agediff.matrix)<-c("E.fuscus","T.brasiliensis","D.rotundus",
                            "P.discolor","P.hastatus","C.perspicillata",
                            "S.bilineata","R.aegyptiacus","Pt.pumilus","Pt.vampyrus","Pt.hypomelanus")
colnames(agediff.matrix)<-rownames(agediff.matrix)
for (x in c(1:11)){
  for (y in c(1:11)){
    agediff.matrix[x,y]<-nrow(allspp.wide.age[!allspp.wide.age[,x]=="NS" &
                                                !allspp.wide.age[,y]=="NS"&
                                                allspp.wide.age[,x]==allspp.wide.age[,y],])/
      min(c(nrow(allspp.wide.age[!allspp.wide.age[,x]=="NS",]),
            nrow(allspp.wide.age[!allspp.wide.age[,y]=="NS",])))
  }
}
ph.age.x<-pheatmap((agediff.matrix),cluster_rows = FALSE,cluster_cols = FALSE,show_rownames=TRUE,
                   main="Age-associated, X-linked",
                   border_color = '#444444',breaks=seq(0,1, length.out=10),scale='none',
                   color=colorRampPalette(c("#3a7ca6","#e0c653", "#d13434"))(10)#  begin = 0, end = 1, option = "viridis")#
)

library(gridExtra)
ggsave('allheatmap.svg',grid.arrange(ph.sex.x[[4]],ph.sex.a[[4]],
                                     ph.age.x[[4]],ph.age.a[[4]],nrow=2),height=8,width=8.5)

