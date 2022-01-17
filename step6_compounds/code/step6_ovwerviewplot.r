rm(list=ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library(RColorBrewer)
setwd("/data/zhang/pancan_cin/step1_calculate_cin/")
#Load Table S2 from data from Taylor et al, Cancer Cell 2018 (https://www.sciencedirect.com/science/article/pii/S1535610818301119)
table <- read.delim("../step2_genome_instability/data/TaylorCancerCell_TableS2.txt", na.strings = c("", " ", "na", "NA", "n.a.", "#N/A"))
rownames(table) <- as.character(table$Sample)
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c("Primary Blood Derived Cancer", "Metastatic", "Solid Tissue Normal", 
                                                    "Blood Derived Normal", "Primary Solid Tumor"),]
Final_table$Sample <- substr(Final_table$Sample, 1, 15)
rownames(Final_table)=Final_table$Sample

res=wilcox.test(Final_table$wgii[Final_table$sample_type_detail %in% "Primary Blood Derived Cancer"], Final_table$wgii[Final_table$sample_type_detail %in% "Blood Derived Normal"])
tmp=names(res)
res=lapply(tmp, function(x){res[[x]]})
names(res)=tmp
sink('table/cmp_wgii_blood_tum_vs_blood_norm.txt')
print(res)
sink()
res=wilcox.test(Final_table$wgii[Final_table$sample_type_detail %in% "Metastatic"], Final_table$wgii[Final_table$sample_type_detail %in% "Blood Derived Normal"])
tmp=names(res)
res=lapply(tmp, function(x){res[[x]]})
names(res)=tmp
sink('table/cmp_wgii_metastatic_tum_vs_blood_norm.txt')
print(res)
sink()
res=wilcox.test(Final_table$wgii[Final_table$sample_type_detail %in% "Primary Solid Tumor"], Final_table$wgii[Final_table$sample_type_detail %in% "Blood Derived Normal"])
tmp=names(res)
res=lapply(tmp, function(x){res[[x]]})
names(res)=tmp
sink('table/cmp_wgii_solid_tum_vs_blood_norm.txt')
print(res)
sink()
res=wilcox.test(Final_table$wgii[Final_table$sample_type_detail %in% "Primary Blood Derived Cancer"], Final_table$wgii[Final_table$sample_type_detail %in% "Solid Tissue Normal"])
tmp=names(res)
res=lapply(tmp, function(x){res[[x]]})
names(res)=tmp
sink('table/cmp_wgii_blood_tum_vs_solid_norm.txt')
print(res)
sink()
res=wilcox.test(Final_table$wgii[Final_table$sample_type_detail %in% "Metastatic"], Final_table$wgii[Final_table$sample_type_detail %in% "Solid Tissue Normal"])
tmp=names(res)
res=lapply(tmp, function(x){res[[x]]})
names(res)=tmp
sink('table/cmp_wgii_metstatic_tum_vs_solid_norm.txt')
print(res)
sink()
res=wilcox.test(Final_table$wgii[Final_table$sample_type_detail %in% "Primary Solid Tumor"], Final_table$wgii[Final_table$sample_type_detail %in% "Solid Tissue Normal"])
tmp=names(res)
res=lapply(tmp, function(x){res[[x]]})
names(res)=tmp
sink('table/cmp_wgii_solid_tum_vs_solid_norm.txt')
print(res)
sink()

res=wilcox.test(Final_table$ncin[Final_table$sample_type_detail %in% "Primary Blood Derived Cancer"], Final_table$ncin[Final_table$sample_type_detail %in% "Blood Derived Normal"])
tmp=names(res)
res=lapply(tmp, function(x){res[[x]]})
names(res)=tmp
sink('table/cmp_ncin_blood_tum_vs_blood_norm.txt')
print(res)
sink()
res=wilcox.test(Final_table$ncin[Final_table$sample_type_detail %in% "Metastatic"], Final_table$ncin[Final_table$sample_type_detail %in% "Blood Derived Normal"])
tmp=names(res)
res=lapply(tmp, function(x){res[[x]]})
names(res)=tmp
sink('table/cmp_ncin_metastatic_tum_vs_blood_norm.txt')
print(res)
sink()
res=wilcox.test(Final_table$ncin[Final_table$sample_type_detail %in% "Primary Solid Tumor"], Final_table$ncin[Final_table$sample_type_detail %in% "Blood Derived Normal"])
tmp=names(res)
res=lapply(tmp, function(x){res[[x]]})
names(res)=tmp
sink('table/cmp_ncin_solid_tum_vs_blood_norm.txt')
print(res)
sink()
res=wilcox.test(Final_table$ncin[Final_table$sample_type_detail %in% "Primary Blood Derived Cancer"], Final_table$ncin[Final_table$sample_type_detail %in% "Solid Tissue Normal"])
tmp=names(res)
res=lapply(tmp, function(x){res[[x]]})
names(res)=tmp
sink('table/cmp_ncin_blood_tum_vs_solid_norm.txt')
print(res)
sink()
res=wilcox.test(Final_table$ncin[Final_table$sample_type_detail %in% "Metastatic"], Final_table$ncin[Final_table$sample_type_detail %in% "Solid Tissue Normal"])
tmp=names(res)
res=lapply(tmp, function(x){res[[x]]})
names(res)=tmp
sink('table/cmp_ncin_metstatic_tum_vs_solid_norm.txt')
print(res)
sink()
res=wilcox.test(Final_table$ncin[Final_table$sample_type_detail %in% "Primary Solid Tumor"], Final_table$ncin[Final_table$sample_type_detail %in% "Solid Tissue Normal"])
tmp=names(res)
res=lapply(tmp, function(x){res[[x]]})
names(res)=tmp
sink('table/cmp_ncin_solid_tum_vs_solid_norm.txt')
print(res)
sink()

res=wilcox.test(Final_table$scin[Final_table$sample_type_detail %in% "Primary Blood Derived Cancer"], Final_table$scin[Final_table$sample_type_detail %in% "Blood Derived Normal"])
tmp=names(res)
res=lapply(tmp, function(x){res[[x]]})
names(res)=tmp
sink('table/cmp_scin_blood_tum_vs_blood_norm.txt')
print(res)
sink()
res=wilcox.test(Final_table$scin[Final_table$sample_type_detail %in% "Metastatic"], Final_table$scin[Final_table$sample_type_detail %in% "Blood Derived Normal"])
tmp=names(res)
res=lapply(tmp, function(x){res[[x]]})
names(res)=tmp
sink('table/cmp_scin_metastatic_tum_vs_blood_norm.txt')
print(res)
sink()
res=wilcox.test(Final_table$scin[Final_table$sample_type_detail %in% "Primary Solid Tumor"], Final_table$scin[Final_table$sample_type_detail %in% "Blood Derived Normal"])
tmp=names(res)
res=lapply(tmp, function(x){res[[x]]})
names(res)=tmp
sink('table/cmp_scin_solid_tum_vs_blood_norm.txt')
print(res)
sink()
res=wilcox.test(Final_table$scin[Final_table$sample_type_detail %in% "Primary Blood Derived Cancer"], Final_table$scin[Final_table$sample_type_detail %in% "Solid Tissue Normal"])
tmp=names(res)
res=lapply(tmp, function(x){res[[x]]})
names(res)=tmp
sink('table/cmp_scin_blood_tum_vs_solid_norm.txt')
print(res)
sink()
res=wilcox.test(Final_table$scin[Final_table$sample_type_detail %in% "Metastatic"], Final_table$scin[Final_table$sample_type_detail %in% "Solid Tissue Normal"])
tmp=names(res)
res=lapply(tmp, function(x){res[[x]]})
names(res)=tmp
sink('table/cmp_scin_metstatic_tum_vs_solid_norm.txt')
print(res)
sink()
res=wilcox.test(Final_table$scin[Final_table$sample_type_detail %in% "Primary Solid Tumor"], Final_table$scin[Final_table$sample_type_detail %in% "Solid Tissue Normal"])
tmp=names(res)
res=lapply(tmp, function(x){res[[x]]})
names(res)=tmp
sink('table/cmp_scin_solid_tum_vs_solid_norm.txt')
print(res)
sink()


data_summary <- function(x) {
  m <- median(x)
  ymin <- as.numeric(quantile(x)[2])
  ymax <- as.numeric(quantile(x)[4])
  return(c(y=m,ymin=ymin,ymax=ymax))
}

dataf=ddply(Final_table, .(sample_type_detail), summarise,py=median(wgii))
dataf$sample_type_detail=as.character(dataf$sample_type_detail)
ii=dataf$sample_type_detail[order(dataf$py)]
Final_table$sample_type_detail=factor(Final_table$sample_type_detail,levels = c("Blood Derived Normal", "Solid Tissue Normal", "Primary Blood Derived Cancer", 
                                                                    "Primary Solid Tumor", "Metastatic")
)


box1<- ggplot(Final_table, aes(x=factor(sample_type_detail), y=wgii)) + 
  scale_y_continuous(limits = quantile(Final_table$wgii, c(0.01, 0.99))) +
  geom_boxplot(aes(fill=factor(sample_type_detail)),width=0.5)+  xlab('sample type')+
  ylab('wgii score')+
  theme(panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), axis.ticks=element_line(colour="black"),
        axis.text.x=element_text(size=5, colour="black",angle = 60, hjust = 1), axis.text.y=element_text(size=8, colour="black"), axis.title=element_text(size=10, colour="black"),
        axis.title.x = element_text(margin = margin(5,0,0,0)), legend.position = "none") + 
  scale_fill_manual('sample type',values=brewer.pal(n = 5, name = "Greens")) +
  geom_signif(comparisons = list(c("Blood Derived Normal", "Primary Blood Derived Cancer")), annotations="****", y_position =0.25, tip_length = 0.01, size = 0.5, textsize = 3, vjust=0.05) +
  geom_signif(comparisons = list(c("Blood Derived Normal", "Primary Solid Tumor")), annotations="****", y_position =0.9, tip_length = 0.01, size = 0.5, textsize = 3, vjust=0.05) +
  geom_signif(comparisons = list(c("Blood Derived Normal", "Metastatic")), annotations="****", y_position =0.93, tip_length = 0.01, size = 0.5, textsize = 3, vjust=0.05) +
  geom_signif(comparisons = list(c("Solid Tissue Normal", "Primary Blood Derived Cancer")), annotations="****", y_position =0.35, tip_length = 0.01, size = 0.5, textsize = 3, vjust=0.05) +
  geom_signif(comparisons = list(c("Solid Tissue Normal", "Primary Solid Tumor")), annotations="****", y_position =0.96, tip_length = 0.01, size = 0.5, textsize = 3, vjust=0.05) +
  geom_signif(comparisons = list(c("Solid Tissue Normal", "Metastatic")), annotations="****", y_position =0.99, tip_length = 0.01, size = 0.5, textsize = 3, vjust=0.05)+
  stat_summary(fun.data=data_summary, color="black", size=0.4)
box2<- ggplot(Final_table, aes(x=factor(sample_type_detail), y=ncin)) + 
  #scale_y_continuous(limits = quantile(Final_table$ncin, c(0.01, 0.99))) +
  geom_boxplot(aes(fill=factor(sample_type_detail)),width=0.5)+  xlab('sample type')+
  ylab('ncin score')+
  theme(panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), axis.ticks=element_line(colour="black"),
        axis.text.x=element_text(size=5, colour="black",angle = 60, hjust = 1), axis.text.y=element_text(size=8, colour="black"), axis.title=element_text(size=10, colour="black"),
        axis.title.x = element_text(margin = margin(5,0,0,0)), legend.position = "none") + 
  scale_fill_manual('sample type',values=brewer.pal(n = 5, name = "Greens")) +
  geom_signif(comparisons = list(c("Blood Derived Normal", "Primary Blood Derived Cancer")), annotations="****", y_position =3, tip_length = 0.03, size = 0.5, textsize = 3, vjust=0.05) +
  geom_signif(comparisons = list(c("Blood Derived Normal", "Primary Solid Tumor")), annotations="****", y_position =19, tip_length = 0.03, size = 0.5, textsize = 3, vjust=0.05) +
  geom_signif(comparisons = list(c("Blood Derived Normal", "Metastatic")), annotations="****", y_position =20, tip_length = 0.03, size = 0.5, textsize = 3, vjust=0.05) +
  geom_signif(comparisons = list(c("Solid Tissue Normal", "Primary Blood Derived Cancer")), annotations="****", y_position =4, tip_length = 0.03, size = 0.5, textsize = 3, vjust=0.05) +
  geom_signif(comparisons = list(c("Solid Tissue Normal", "Primary Solid Tumor")), annotations="****", y_position =21, tip_length = 0.03, size = 0.5, textsize = 3, vjust=0.05) +
  geom_signif(comparisons = list(c("Solid Tissue Normal", "Metastatic")), annotations="****", y_position =22, tip_length = 0.03, size = 0.5, textsize = 3, vjust=0.05)+
  stat_summary(fun.data=data_summary, color="black", size=0.4)

box3<- ggplot(Final_table, aes(x=factor(sample_type_detail), y=scin)) + 
  scale_y_continuous(limits = quantile(Final_table$scin, c(0.01, 0.99))) +
  geom_boxplot(aes(fill=factor(sample_type_detail)),width=0.5)+  xlab('sample type')+
  ylab('scin score')+
  theme(panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), axis.ticks=element_line(colour="black"),
        axis.text.x=element_text(size=5, colour="black",angle = 60, hjust = 1), axis.text.y=element_text(size=8, colour="black"), axis.title=element_text(size=10, colour="black"),
        axis.title.x = element_text(margin = margin(5,0,0,0)), legend.position = "none") + 
  scale_fill_manual('sample type',values=brewer.pal(n = 5, name = "Greens")) +
  geom_signif(comparisons = list(c("Blood Derived Normal", "Primary Blood Derived Cancer")), annotations="****", y_position =50, tip_length = 0.003, size = 0.5, textsize = 3, vjust=0.05) +
  geom_signif(comparisons = list(c("Blood Derived Normal", "Primary Solid Tumor")), annotations="****", y_position =54, tip_length = 0.003, size = 0.5, textsize = 3, vjust=0.05) +
  geom_signif(comparisons = list(c("Blood Derived Normal", "Metastatic")), annotations="****", y_position =58, tip_length = 0.003, size = 0.5, textsize = 3, vjust=0.05) +
  geom_signif(comparisons = list(c("Solid Tissue Normal", "Primary Blood Derived Cancer")), annotations="****", y_position =65, tip_length = 0.003, size = 0.5, textsize = 3, vjust=0.05) +
  geom_signif(comparisons = list(c("Solid Tissue Normal", "Primary Solid Tumor")), annotations="****", y_position =69, tip_length = 0.003, size = 0.5, textsize = 3, vjust=0.05) +
  geom_signif(comparisons = list(c("Solid Tissue Normal", "Metastatic")), annotations="****", y_position =73, tip_length = 0.003, size = 0.5, textsize = 3, vjust=0.05)+
  stat_summary(fun.data=data_summary, color="black", size=0.4)

pdf("plot/cin_vs_sample_type.pdf", width =15, height = 4.5)
plot_grid(box1,box2,box3, labels = c("A", "B","C"),ncol=3)
dev.off()
pdf("plot/wgii_vs_sample_type_cohort.pdf", width =15, height = 4.5)
plot_sample_type_violin <- ggplot(Final_table, aes(x=Cohort, y=wgii, fill=factor(sample_type_detail)))+geom_boxplot(outlier.size = 0.1)+scale_y_log10()+
  labs(x="", y = "wgii score")+
  theme(panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), axis.text.x=element_text(size=14, angle=45, hjust=1, colour="black"), axis.text.y=element_text(size=12, colour="black"), axis.title.y=element_text(size=16, colour="black", margin=margin(0,10,0,0)), legend.key = element_rect(fill="White"), legend.text = element_text(size=14), plot.title = element_text(size=17, hjust=0.5, colour="black")) +
  #guides(fill=F, colour=F) +
  scale_fill_manual('sample type',values=brewer.pal(n = 5, name = "Set3")) 
print(plot_sample_type_violin)
dev.off()
pdf("plot/ncin_vs_sample_type_cohort.pdf", width =15, height = 4.5)
plot_sample_type_violin <- ggplot(Final_table, aes(x=Cohort, y=ncin+1, fill=factor(sample_type_detail)))+geom_boxplot(outlier.size = 0.1)+scale_y_log10()+
  labs(x="", y = "ncin score")+
  theme(panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), axis.text.x=element_text(size=14, angle=45, hjust=1, colour="black"), axis.text.y=element_text(size=12, colour="black"), axis.title.y=element_text(size=16, colour="black", margin=margin(0,10,0,0)), legend.key = element_rect(fill="White"), legend.text = element_text(size=14), plot.title = element_text(size=17, hjust=0.5, colour="black")) +
  #guides(fill=F, colour=F) +
  scale_fill_manual('sample type',values=brewer.pal(n = 5, name = "Set3")) 
print(plot_sample_type_violin)
dev.off()
pdf("plot/scin_vs_sample_type_cohort.pdf", width =15, height = 4.5)
plot_sample_type_violin <- ggplot(Final_table, aes(x=Cohort, y=scin+1, fill=factor(sample_type_detail)))+geom_boxplot(outlier.size = 0.1)+scale_y_log10()+
  labs(x="", y = "scin score")+
  theme(panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), axis.text.x=element_text(size=14, angle=45, hjust=1, colour="black"), axis.text.y=element_text(size=12, colour="black"), axis.title.y=element_text(size=16, colour="black", margin=margin(0,10,0,0)), legend.key = element_rect(fill="White"), legend.text = element_text(size=14), plot.title = element_text(size=17, hjust=0.5, colour="black")) +
  #guides(fill=F, colour=F) +
  scale_fill_manual('sample type',values=brewer.pal(n = 5, name = "Set3")) 
print(plot_sample_type_violin)
dev.off()

