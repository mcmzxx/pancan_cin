rm(list=ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library(data.table)
library(RColorBrewer)
setwd("/data/zhang/pancan_cin/step3_mutation/")
fns = c(
  "DDR gene mutations.csv",
  "DDR gene alterations ONCOPRINT.csv",
  "DDR gene alterations.csv",
  "DDR footprints.csv",
  "DDR epigenetic silencing.csv",
  "DDR deep deletions.csv"
)
pathdir = "../step8_tme/data/TCGA_DDR_Data_Resources_xls"
datalist = list()
for (i in fns) {
  if (i != 'DDR footprints.csv') {
    tmp = fread(file.path(pathdir, i))
    tmp = data.frame(tmp, stringsAsFactors = F)
    coltmp = colnames(tmp)[-c(1:2)]
    rowtmp = as.character(tmp[, 2][-c(1:2)])
    mattmp = tmp[-c(1:2), -c(1:2)]
    rownames(mattmp) = rowtmp
    colnames(mattmp) = coltmp
    datalist[[i]] = mattmp
  }
}


i = "DDR footprints.csv"
tmp = fread(file.path(pathdir, i))
tmp = data.frame(tmp, stringsAsFactors = F)
coltmp = c('Cohort', colnames(tmp)[-c(1:4)])
rowtmp = as.character(tmp[, 2][-c(1:4)])
mattmp = tmp[-c(1:4), -c(1:2, 4)]
rownames(mattmp) = rowtmp
colnames(mattmp) = coltmp
mattmp = mattmp[,c("LOH_n_seg", "LOH_frac_altered", "purity", "ploidy", "genome_doublings","tp53_score")]
mattmp$sample_id=substr(gsub('-','.',rownames(mattmp)),1,12)

setwd("/data/zhang/pancan_cin/step1_calculate_cin//")
#Load Table S2 from data from Taylor et al, Cancer Cell 2018 (https://www.sciencedirect.com/science/article/pii/S1535610818301119)
table <- read.delim("../step2_genome_instability/data/TaylorCancerCell_TableS2.txt", na.strings = c("", " ", "na", "NA", "n.a.", "#N/A"))
rownames(table) <- as.character(table$Sample)
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer','Primary Solid Tumor'),]
Final_table$Sample <- substr(Final_table$Sample, 1, 15)
rownames(Final_table)=Final_table$Sample
sels=intersect(rownames(table),rownames(Final_table))
Final_table=Final_table[sels,]

# Proliferation rates from paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5186797/pdf/fphys-07-00644.pdf
# The table used is in https://github.com/cdiener/proliferation/blob/master/results/pred_rates.csv
Prolif_score <- read.csv("../step1_calculate_cin/data/pred_rates.csv")
Prolif_score$sample_id <- gsub("-",".",Prolif_score$patient_barcode)
Prolif_score <- Prolif_score[Prolif_score$tumor %in% "TRUE",]
# summarize by homen
Prolif_score_2 <- aggregate(Prolif_score[,3],by=list(sample_id=Prolif_score$sample_id),data=Prolif_score,FUN=mean)
names(Prolif_score_2)[2]<- "rates"

Final_table$sample_id <-gsub('-','.',substr(Final_table$Sample,1,12))
Final_table<- join(Final_table, Prolif_score_2,type='left',by='sample_id')
Final_table<- join(Final_table, mattmp,type='left',by='sample_id')

#Final_table=Final_table[which(Final_table$genome_doublings==0),]
res1=cor.test(Final_table$wgii, Final_table$rates, method = "spearman")
tmp=names(res1)
res1=lapply(tmp, function(x){res1[[x]]})
names(res1)=tmp
res2=cor.test(Final_table$ncin, Final_table$rates, method = "spearman")
tmp=names(res2)
res2=lapply(tmp, function(x){res2[[x]]})
names(res2)=tmp
res3=cor.test(Final_table$scin, Final_table$rates, method = "spearman")
tmp=names(res3)
res3=lapply(tmp, function(x){res3[[x]]})
names(res3)=tmp
sink('table/wgii_vs_proliferation_cor.txt')
print(res1)
sink()

sink('table/ncin_vs_proliferation_cor.txt')
print(res2)
sink()

sink('table/scin_vs_proliferation_cor.txt')
print(res3)
sink()

scatter1 <- ggplot(Final_table, aes(x=wgii, y=rates)) + #scale_x_log10() +
  #scale_y_log10() +
  xlab("wgii score") + 
  ylab(paste("proliferation rates [1/h]")) + 
  ggtitle(paste0("")) +
  geom_point(colour="grey70", size=1, alpha=1) +
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(low="white",high="grey20") +
  scale_alpha(range = c(0.1,0.3)) +
  guides(alpha="none", fill="none") +
  theme(panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), axis.text.y=element_text(size=13, colour="black"), axis.text.x=element_text(size=13, colour="black"), axis.title=element_text(size=16, colour="black"), plot.title = element_text(size=17, hjust = 0.5, colour="black"), legend.key=element_rect(fill=NA), legend.text = element_text(size=13, colour="black"), legend.title = element_text(size=16, colour="black"), axis.title.x=element_text(margin=margin(10,0,0,0)))


scatter2 <- ggplot(Final_table, aes(y=ncin, x=rates))+ #scale_x_log10() +
  #scale_y_log10() +
  ylab("ncin score") + 
  xlab(paste("proliferation rates [1/h]")) + 
  ggtitle(paste0("")) +
  geom_point(colour="grey70", size=1, alpha=1) +
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(low="white",high="grey20") +
  scale_alpha(range = c(0.1,0.3)) +
  guides(alpha="none", fill="none") +
  theme(panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), axis.text.y=element_text(size=13, colour="black"), axis.text.x=element_text(size=13, colour="black"), axis.title=element_text(size=16, colour="black"), plot.title = element_text(size=17, hjust = 0.5, colour="black"), legend.key=element_rect(fill=NA), legend.text = element_text(size=13, colour="black"), legend.title = element_text(size=16, colour="black"), axis.title.x=element_text(margin=margin(10,0,0,0)))

scatter3 <- ggplot(Final_table, aes(y=scin, x=rates))+ #scale_x_log10() +
  #scale_y_log10() +
  ylab("scin score") + 
  xlab(paste("proliferation rates [1/h]")) + 
  ggtitle(paste0("")) +
  geom_point(colour="grey70", size=1, alpha=1) +
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(low="white",high="grey20") +
  scale_alpha(range = c(0.1,0.3)) +
  guides(alpha="none", fill="none") +
  theme(panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), axis.text.y=element_text(size=13, colour="black"), axis.text.x=element_text(size=13, colour="black"), axis.title=element_text(size=16, colour="black"), plot.title = element_text(size=17, hjust = 0.5, colour="black"), legend.key=element_rect(fill=NA), legend.text = element_text(size=13, colour="black"), legend.title = element_text(size=16, colour="black"), axis.title.x=element_text(margin=margin(10,0,0,0)))

pdf("plot/cin_score_vs_proliferation.pdf", height  = 4,width=12)
plot_grid(scatter1,scatter2,scatter3, labels = c("A", "B","C"),ncol=3)
dev.off()



pdf("plot/cin_score_vs_proliferation.pdf", height  = 4,width=9)
plot_grid(scatter2,scatter3, labels = c("N-CIN vs proliferation", "S-CIN vs proliferation"),ncol=2)
dev.off()
