rm(list=ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library(data.table)
setwd("/data/zhang/pancan_cin/step6_compounds/")
load('data/ccle_cin_score.RData')
scores$Sample=toupper(scores$Sample)
#scores=subset(scores,scores$Genome.doublings==0)
library(ggplot2)
library(ggrepel)
format_name=function(x){
  x=gsub(')','',x)
  x_split = strsplit(x," \\(")[[1]]
  return(paste0(c(x_split),collapse = '_'))
}
achiles=fread('/storage/zhang/crc_mfgae/data/ccle/dwl/Achilles_gene_effect.csv') 
ind=achiles$V1
achiles=as.data.frame(achiles)
achiles$V1=NULL
rownames(achiles)=ind
colnames(achiles)=unlist(lapply(colnames(achiles), function(x){format_name(x)}))
head(achiles[,1:5])
d2achiles=fread('/storage/zhang/crc_mfgae/data/ccle/dwl/D2_DRIVE_gene_dep_scores.csv')
ind=d2achiles$V1
d2achiles$V1=NULL
head(d2achiles[,1:5])
d2achiles=t(d2achiles)
d2achiles=as.data.frame(d2achiles)
colnames(d2achiles)=unlist(lapply(ind, function(x){format_name(x)}))

d2drive=fread('/storage/zhang/crc_mfgae/data/ccle/dwl/D2_DRIVE_gene_dep_scores.csv')
ind=d2drive$V1
d2drive$V1=NULL
head(d2drive[,1:5])
d2drive=t(d2drive)
d2drive=as.data.frame(d2drive)
colnames(d2drive)=unlist(lapply(ind, function(x){format_name(x)}))


cell_line_info = fread("/storage/zhang/crc_mfgae/data/ccle/dwl/sample_info.csv")
ccle_achilles_map = as.list(cell_line_info$CCLE_Name)
names(ccle_achilles_map)=cell_line_info$DepMap_ID
colnames(d2drive)=unlist(lapply(colnames(d2drive), function(x){format_name(x)}))
achiles$Sample=unlist(lapply(rownames(achiles),function(x){
  if(!(x %in% names(ccle_achilles_map)))
  {return(x)}
  else{return(ccle_achilles_map[[x]])}
}))
d2achiles$Sample=rownames(d2achiles)
d2drive$Sample=rownames(d2drive)


ginsgene=c("GINS1_9837", "GINS2_51659", "GINS3_64785", "GINS4_84296")
ginsgene=c("CKS1B_1163", "CKS2_1164")
tmpmat=merge(scores, achiles, by='Sample')
scatter1<- ggplot(tmpmat, aes(x=scin, y=CKS1B_1163))+ #scale_x_log10() +
  #scale_y_log10() +
  xlab("scin score") + 
  ylab(paste("CKS1B dependency score")) + 
  ggtitle(paste0("")) +
  geom_point(colour="grey70", size=1, alpha=1) +
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(low="white",high="grey20") +
  scale_alpha(range = c(0.1,0.3)) +
  guides(alpha="none", fill="none") +
  theme(panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), axis.text.y=element_text(size=13, colour="black"), axis.text.x=element_text(size=13, colour="black"), axis.title=element_text(size=16, colour="black"), plot.title = element_text(size=17, hjust = 0.5, colour="black"), legend.key=element_rect(fill=NA), legend.text = element_text(size=13, colour="black"), legend.title = element_text(size=16, colour="black"), axis.title.x=element_text(margin=margin(10,0,0,0)))

scatter2<- ggplot(tmpmat, aes(x=scin, y=CKS2_1164))+ #scale_x_log10() +
  #scale_y_log10() +
  xlab("scin score") + 
  ylab(paste("CKS2 dependency score")) + 
  ggtitle(paste0("")) +
  geom_point(colour="grey70", size=1, alpha=1) +
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(low="white",high="grey20") +
  scale_alpha(range = c(0.1,0.3)) +
  guides(alpha="none", fill="none") +
  theme(panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), axis.text.y=element_text(size=13, colour="black"), axis.text.x=element_text(size=13, colour="black"), axis.title=element_text(size=16, colour="black"), plot.title = element_text(size=17, hjust = 0.5, colour="black"), legend.key=element_rect(fill=NA), legend.text = element_text(size=13, colour="black"), legend.title = element_text(size=16, colour="black"), axis.title.x=element_text(margin=margin(10,0,0,0)))
pdf("plot/scin_vs_CKS_dependency.pdf", height  = 4,width=8)
plot_grid(scatter1,scatter2,labels = c("A", "B"),ncol=2)
dev.off()



ginsgene=c("GINS1_9837", "GINS2_51659", "GINS3_64785", "GINS4_84296")
tmpmat=merge(scores, achiles, by='Sample')
scatter1<- ggplot(tmpmat, aes(x=wgii, y=GINS1_9837))+ #scale_x_log10() +
  #scale_y_log10() +
  xlab("wgii score") + 
  ylab(paste("GINS1 dependency score")) + 
  ggtitle(paste0("")) +
  geom_point(colour="grey70", size=1, alpha=1) +
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(low="white",high="grey20") +
  scale_alpha(range = c(0.1,0.3)) +
  guides(alpha="none", fill="none") +
  theme(panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), axis.text.y=element_text(size=13, colour="black"), axis.text.x=element_text(size=13, colour="black"), axis.title=element_text(size=16, colour="black"), plot.title = element_text(size=17, hjust = 0.5, colour="black"), legend.key=element_rect(fill=NA), legend.text = element_text(size=13, colour="black"), legend.title = element_text(size=16, colour="black"), axis.title.x=element_text(margin=margin(10,0,0,0)))

scatter2<- ggplot(tmpmat, aes(x=wgii, y=GINS2_51659))+ #scale_x_log10() +
  #scale_y_log10() +
  xlab("wgii score") + 
  ylab(paste("GINS2 dependency score")) + 
  ggtitle(paste0("")) +
  geom_point(colour="grey70", size=1, alpha=1) +
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(low="white",high="grey20") +
  scale_alpha(range = c(0.1,0.3)) +
  guides(alpha="none", fill="none") +
  theme(panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), axis.text.y=element_text(size=13, colour="black"), axis.text.x=element_text(size=13, colour="black"), axis.title=element_text(size=16, colour="black"), plot.title = element_text(size=17, hjust = 0.5, colour="black"), legend.key=element_rect(fill=NA), legend.text = element_text(size=13, colour="black"), legend.title = element_text(size=16, colour="black"), axis.title.x=element_text(margin=margin(10,0,0,0)))

scatter3<- ggplot(tmpmat, aes(x=wgii, y=GINS3_64785))+ #scale_x_log10() +
  #scale_y_log10() +
  xlab("wgii score") + 
  ylab(paste("GINS3 dependency score")) + 
  ggtitle(paste0("")) +
  geom_point(colour="grey70", size=1, alpha=1) +
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(low="white",high="grey20") +
  scale_alpha(range = c(0.1,0.3)) +
  guides(alpha="none", fill="none") +
  theme(panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), axis.text.y=element_text(size=13, colour="black"), axis.text.x=element_text(size=13, colour="black"), axis.title=element_text(size=16, colour="black"), plot.title = element_text(size=17, hjust = 0.5, colour="black"), legend.key=element_rect(fill=NA), legend.text = element_text(size=13, colour="black"), legend.title = element_text(size=16, colour="black"), axis.title.x=element_text(margin=margin(10,0,0,0)))

scatter4<- ggplot(tmpmat, aes(x=wgii, y=GINS4_84296))+ #scale_x_log10() +
  #scale_y_log10() +
  xlab("wgii score") + 
  ylab(paste("GINS4 dependency score")) + 
  ggtitle(paste0("")) +
  geom_point(colour="grey70", size=1, alpha=1) +
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(low="white",high="grey20") +
  scale_alpha(range = c(0.1,0.3)) +
  guides(alpha="none", fill="none") +
  theme(panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), axis.text.y=element_text(size=13, colour="black"), axis.text.x=element_text(size=13, colour="black"), axis.title=element_text(size=16, colour="black"), plot.title = element_text(size=17, hjust = 0.5, colour="black"), legend.key=element_rect(fill=NA), legend.text = element_text(size=13, colour="black"), legend.title = element_text(size=16, colour="black"), axis.title.x=element_text(margin=margin(10,0,0,0)))

pdf("plot/wgii_vs_GINS_dependency.pdf", height  = 4,width=16)
plot_grid(scatter1,scatter2,scatter3,scatter4,labels = c("A", "B","C", "D"),ncol=4)
dev.off()