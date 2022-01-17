rm(list=ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library(RColorBrewer)
setwd("/data/zhang/pancan_cin/step1_calculate_cin/")
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
#write.table(scores,file='/storage/zhang/experiments/DrugCell/data/raw/bulk_tumor/bulktumor_cin_scores.txt',sep='\t',quote = F,row.names = F,col.names = T)

Final_table=scores[scores$sample_type_detail %in% c("Primary Blood Derived Cancer", "Metastatic", "Solid Tissue Normal", 
                                                    "Blood Derived Normal", "Primary Solid Tumor"),]
Final_table$sample_type_detail[Final_table$sample_type_detail %in% c("Solid Tissue Normal","Blood Derived Normal")]='Normal'
Final_table$sample_type_detail[Final_table$sample_type_detail %in% c("Primary Blood Derived Cancer","Primary Solid Tumor")]='Tumor'
Final_table$Sample <- substr(Final_table$Sample, 1, 15)
rownames(Final_table)=Final_table$Sample


Final_table$sample_type_detail=factor(Final_table$sample_type_detail,levels =c("Normal", "Tumor", "Metastatic"))
gcolor=unname(yarrr::piratepal("pony"))[1:nlevels(Final_table$sample_type_detail)]
my_comparisons <- list(c("Normal","Tumor"),c("Normal","Metastatic"),c("Tumor","Metastatic"))
pl <- ggboxplot(Final_table, x = "sample_type_detail", y = "wgii",fill = "sample_type_detail",width=0.3, size=0.1, outlier.colour = 'black')+ylab('WGII score')+
  scale_fill_manual("Sample Type", labels = c("Normal Tissue", "Tumor Tisue", "Metastatic Tumor"),values=gcolor)+
  scale_x_discrete(labels=c("Normal"="normal", "Tumor"="tumor", "Metastatic"="metastatic"))+
  stat_compare_means(aes(label = paste('p',..p.format..)),label.y= 1.4,label.x.npc = 'left')+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(legend.position = 'none',axis.text.y=element_text(size=10, colour="black"), 
        axis.text.x=element_text(size=10, colour="black",angle = 90, vjust=0.5),
        axis.title.x=element_blank())
pl
pl1 <-ggboxplot(Final_table, x = "sample_type_detail", y = "scin",fill = "sample_type_detail",width=0.3, size=0.1, outlier.colour = 'black')+ylab('SCIN score')+
  scale_fill_manual("Sample Type", labels = c("Normal Tissue", "Tumor Tisue", "Metastatic Tumor"),values=gcolor)+
  scale_x_discrete(labels=c("Normal"="normal", "Tumor"="tumor", "Metastatic"="metastatic"))+
  stat_compare_means(aes(label = paste('p',..p.format..)),label.y= 580,label.x.npc = 'left')+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(legend.position = 'none',axis.text.y=element_text(size=10, colour="black"), 
        axis.text.x=element_text(size=10, colour="black",angle = 90, vjust=0.5),
        axis.title.x=element_blank())
pl1=pl1+theme(plot.margin=margin(l=-0.1,t=0.1,unit="cm"))
pl=pl+theme(plot.margin=margin(r=-0.1,t=0.1,l=0.1,unit="cm"))
plot_grid(pl,pl1)
pdf("plot/cin_vs_sample_type.pdf",width =4,height =4)
plot_grid(pl,pl1)
dev.off()
