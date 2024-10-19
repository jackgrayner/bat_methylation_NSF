#The first part of this script is adapted from https://github.com/shorvath/MammalianMethylationConsortium/blob/main/Annotations%2C%20Amin%20Haghani/Alignment%20codes%2C%20Amin%20Haghani.pdf

library(easypackages)
libraries("tidyr", "dplyr","BiocManager","parallel", "QuasR", "Rsamtools", "ChIPseeker","GenomicFeatures","data.table")

setwd("~/Documents/UMD_new/DNAm-sex_diffs/probe alignments/")
path<-("~/Documents/UMD_new/DNAm-sex_diffs/probe alignments/PhHa_ncbi_dataset/ncbi_dataset/data/GCF_019186645.2/")

f <-paste0(path,list.files(path, pattern="\\.fna$")) 
alignment <- qAlign("~/Documents/UMD_new/DNAm-sex_diffs/probe alignments/sample.txt", genome = f, 
                      bisulfite = "undir", alignmentParameter = "-k 2 --strata --best -v 3")

#aln <- BamFile("GPL28271-57075_b1ac3e046a0.bam")#desro
#aln <- BamFile("GPL28271-57075_b1ac4ef66b7f.bam")#epfu
#aln <- BamFile("GPL28271-57075_b1aca7bfa42.bam")#phdi
#aln <- BamFile("GPL28271-57075_b1ac52f6ae3f.bam")#sabi
aln <- BamFile("GPL28271-57075_b1ac3a15caef.bam")#phha
aln <- scanBam(aln)
aln <- as.data.frame(aln[[1]])

manifest<-read.csv("Manifest, HorvathMammalMethylChip40.csv")

#create targetCG variable. I think this needs to refer to where in the SourceSeq the CG dimer is (1:2 or 49:50)
manifest$targetCG<-NA
manifest[substr(manifest$SourceSeq,0,2)=="CG",]$targetCG<-"1:2"
manifest[substr(manifest$SourceSeq,49,50)=="CG",]$targetCG<-"49:50"
summary(factor(manifest$targetCG))
head(manifest[is.na(manifest$targetCG),]$SourceSeq)

# Determination of CG location based on the probe design. The probe is designed by either top or bottom strand. 
aln1 <- manifest %>% dplyr::select(IlmnID, SourceSeq, targetCG) %>% dplyr::rename(qname = IlmnID) %>% right_join(aln, by="qname")%>% 
  mutate(targetCG = as.character(targetCG))

CGcount <- rbindlist(lapply(1:nrow(aln1), function(i){
  pattern <- DNAString(as.character(aln1$SourceSeq[i]))
  subject <- DNAString(aln1$seq[i])
  matches <- matchPattern(pattern, subject, max.mismatch = 0, algorithm = "naive-inexact")
  locations = paste(start(matches), end(matches), sep=":")
  pattern2 <-reverseComplement(DNAString(as.character(aln1$SourceSeq[i])))
  matches2 <- matchPattern(pattern2, subject, max.mismatch = 0, algorithm = "naive-inexact")
  locations2 = paste(start(matches2), end(matches2), sep=":")
  hits <- data.frame(qname=aln1$qname[i],
                     CGcount = length(start(matches))+length(start(matches2)), 
                     forward = paste(locations, collapse = " ; "),
                     reverse = paste(locations2, collapse = " ; "))
}))


aln1$alignedStand <- ifelse(CGcount$forward!="", "forward", "complementReverse")
aln1$targetCG <- ifelse(aln1$alignedStand=="forward", aln1$targetCG, 
                       ifelse(aln1$alignedStand=="complementReverse"&aln1$targetCG=="1:2", "49:50",
                              ifelse(aln1$alignedStand=="complementReverse"&aln1$targetCG=="49:50", "1:2",NA)))
aln1$targetCG <- as.numeric(as.character(factor(aln1$targetCG, levels = c("1:2", "49:50"), labels = c(0,48))))

# convert to GRange for annotation
input <- aln1 %>% dplyr::select(qname, rname, strand, pos) %>% dplyr::filter(complete.cases(.)) %>%
  mutate(start = pos) %>% mutate(end = pos+49)
input <- input[,c(2,5,6,1, 3)]
names(input) <- c("chr","start", "end", "CGid", "strand")
target <- with(input,
               GRanges( seqnames = Rle(chr),
                        ranges   = IRanges(start, end=end, names=CGid),
                        strand   = Rle(strand(strand)) ))

# create TxDB 
library(GenomicFeatures)
txdb <- makeTxDbFromGFF(paste0(path,"genomic.gff"), format = "gff3")

# annotating the probes and estimating the CG location
peakAnno <- annotatePeak(target, tssRegion=c(-10000, 1000),
                         TxDb=txdb,
                         sameStrand = FALSE, overlap = "all", addFlankGeneInfo=T)
genomeAnnotation <- data.frame(CGid = peakAnno@anno@ranges@NAMES, peakAnno@anno, 
                             peakAnno@detailGenomicAnnotation)
genomeAnnotation <- genomeAnnotation %>% dplyr::rename(probeStart = start, probeEnd = end)
genomeAnnotation <- aln1 %>% dplyr::select(qname, targetCG, seq) %>% 
  dplyr::rename(CGid = qname) %>% 
  right_join(genomeAnnotation, by="CGid") %>% 
  mutate(CGstart = probeStart+targetCG, CGend =probeStart+targetCG+1) %>%
  dplyr::select(-targetCG)


#I did not run these steps
# Confirming if the CG is real. This step is done by extracting the sequence from the original FASTA file
#BEDfile <- genomeAnnotation %>% dplyr::select(seqnames, CGstart, CGend, 
#                                     CGid, strand) %>% 
#  setnames(new = c("chrom", 'chromStart', 'chromEnd', 'name', "strand")) %>%
#  filter(!is.na(chromStart)) %>% mutate(chromStart = chromStart-1) 
#write.table(BEDfile, "BEDfile.bed", 
#            sep = "\t", row.names=F, col.names=F, quote = F)

#bedtools getfasta -fi GCF_022682495.1_HLdesRot8A_genomic.fna -bed BEDfile.bed -fo BEDfile.fasta

#CGs <- readDNAStringSet("BEDfile.fasta")
#seq_name = names(CGs)
#I added the substr bit -- assume this was what was intended (checking it starts with CG)
#sequence = paste(CGs)
#df <- data.frame(seq_name, sequence) %>% dplyr::rename(CG = sequence) %>% 
#  mutate(CG = ifelse(CG %in% c("CG", "GC"), TRUE, FALSE))

#for some reason, have to replace ", " with " "
genomeAnnotation$annotation<-gsub(", "," ",genomeAnnotation$annotation)
f<-list.files(path, pattern="\\.fna$")
write.csv(x=genomeAnnotation,file=gsub("\\.fna",".HorvathMammalMethylChip40.csv",f),quote=FALSE)

#check for conserved rhinolophus/yangochiroptera probes
rhino<-read.csv("../Rhinolophus_ferrumequinum.hlrhifer5.HorvathMammalMethylChip40.v1.csv",h=T) %>% dplyr::select('CGid','SYMBOL')
colnames(rhino)

phha<-read.csv("PhHa_GCF_019186645.2_TTU_PhHast_1.1_genomic.HorvathMammalMethylChip40.csv",h=T) %>% 
  mutate(SYMBOL=geneId) %>% dplyr::select('CGid','SYMBOL')
phdi<-read.csv("PhDi_GCF_004126475.2_mPhyDis1.pri.v3_genomic.HorvathMammalMethylChip40.csv",h=T) %>% 
  mutate(SYMBOL=geneId) %>% dplyr::select("CGid","SYMBOL")
sabi<-read.csv("SaBi_GCF_036850765.1_mSacBil1_pri_phased_curated_genomic.HorvathMammalMethylChip40.csv",h=T) %>% 
  mutate(SYMBOL=geneId) %>% dplyr::select("CGid","SYMBOL")
dero<-read.csv("DeRo_GCF_022682495.1_HLdesRot8A_genomic.HorvathMammalMethylChip40.csv",h=T) %>% 
  mutate(SYMBOL=geneId) %>% dplyr::select("CGid","SYMBOL")
epfu<-read.csv("EpFu_GCF_027574615.1_DD_ASM_mEF_20220401_genomic.HorvathMammalMethylChip40.csv",h=T) %>% 
  mutate(SYMBOL=geneId) %>% dplyr::select("CGid","SYMBOL")

df_list <- list(rhino,dero,epfu,phdi,phha,sabi)
all_annos<-df_list %>% reduce(full_join, by='CGid')
colnames(all_annos)<-c("CGid","rhino", "dero", "epfu","phdi","phha","sabi")
all_annos<-all_annos[!is.na(all_annos$rhino),]

all_annos.conserved<-subset(all_annos,
                  rhino==dero | rhino==epfu | rhino==phdi | rhino==phha | rhino==sabi)

all_annos$conserved<-"N"
all_annos[all_annos$CGid %in% all_annos.conserved$CGid,]$conserved<-"Y"

summary(factor(all_annos$conserved))
all.annos.notconserved<-all_annos[all_annos$conserved=="N",]
nrow(all.annos.notconserved[grep("LOC",all.annos.notconserved$rhino),])

write.csv(all_annos,'rhino_yango_probe_annotations.csv',quote=FALSE)
