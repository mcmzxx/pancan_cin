rm(list=ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library(data.table)
library("colorspace")
tropiccolor=c("deepskyblue","grey69","deeppink")
mycolor=c('darkblue','white','darkred','black')
setwd("/data/zhang/pancan_cin/step3_mutation/")
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer','Primary Solid Tumor'),]
Final_table$Sample = substr(Final_table$Sample, 1, 15)
rownames(Final_table)=Final_table$Sample
Final_table$Genome_doublings = as.factor(Final_table$`Genome doublings`)

  Results = read.delim("table/scin_allGeneExp_metacor_results.txt")
Results = Results[order(Results$Pvalue),]
Results = Results[!Results$Gene %in% ".",]

Results$significant = ifelse(Results$FDR < 0.05 & Results$Coef > 0, "significant", "not significant")
Results$significant[Results$FDR < 0.05 & Results$Coef < 0] = "significant"
Results$significant2 = ifelse(Results$FDR < 0.05 & Results$Coef > 0, "postive", "not significant")
Results$significant2[Results$FDR < 0.05 & Results$Coef < 0] = "negtive"
Results$text = "no"
selt=unique(c(which(Results$Gene %in% c("MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "CDC7", "DBF4", 
                                        "GINS1", "POLD1", "POLD2", "POLD3", "POLE", "PCNA", "GINS2", 
                                        "GINS3", "GINS4", "CDC45", "CDK1", "CCNE1", "CCNE2", "CDK2", 
                                        "CCNA1", "CCNA2", "WDHD1", "RECQL4", "C15orf42", "TOPBP1","AURKB","BUB1",'BUB3','AURKA'))))
#"TPX2", "E2F1", "UBE2C",'PIK3CA','KRAS','BRAF','NRAS', 'HRAS', "AURKB","BUB1",'MYC','RB1','BUB3','AURKA','PTEN'
selt=intersect(selt,which(Results$significant2!='not significant'))
mitotic=which(Results$Gene %in% c("RAE1","TPX2","UBE2C","AURKA","AURB","BUB1","CDK1"))
mitotic=intersect(mitotic,which(Results$significant2!='not significant'))
replication=which(Results$Gene %in% c("GINS1","CDC45","MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7","DBF4","CDC7","GINS3","GINS2"))
replication=intersect(replication,which(Results$significant2!='not significant'))
Results$text[selt] = "others"
Results$text[replication] = "replication factor"
Results$text[mitotic] = "mitotic factor"
Results$text = as.factor(Results$text)
represults=subset(Results, text %in% c("replication factor"))
mitresults=subset(Results, text %in% c("mitotic factor"))
otherresults=subset(Results, text %in% c("others"))
highresults=subset(Results, Gene %in% c("GINS1","CDC45"))
plot1 =ggplot(Results, aes(x = Coef, y = -log10(FDR))) +
  geom_point(aes(colour = significant),shape=19,size=0.2) +
  xlab("mean difference high WGII vs low WGII (meta-analysis)") +
  ylab("-log10(FDR adjusted p-value)") +
  theme_cowplot(12)+ theme(axis.text.x=element_text(size=15, colour="black"),  legend.position = 'bottom',
                           axis.text.y=element_text(size=15, colour="black"), 
                           axis.title=element_text(size=15, colour="black"), 
                           legend.key = element_rect(fill="White"), 
                           legend.text = element_text(size=15), 
                           legend.title = element_text(size=15)) +
  geom_text_repel(
    data = represults,
    aes(label = Gene,colour=text),
    size =3,max.overlaps = Inf,segment.size=0.2,
    nudge_x= 0.02, direction = "y",angle=0,
    arrow = arrow (length = unit (0.02, "npc"), ends = "last", angle = 0),
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines"),
    hjust=0
  )+
  geom_text_repel(
    data = mitresults,
    aes(label = Gene,colour=text),
    size =3,max.overlaps = Inf,segment.size=0.2,
    nudge_y=55, direction = "both",angle=90,
    arrow = arrow (length = unit (0.02, "npc"), ends = "last", angle = 0),
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines"),
    hjust=0
  )+
  geom_text_repel(
    data = otherresults,
    aes(label = Gene,colour=text),
    size =3,max.overlaps = Inf,segment.size=0.2,
    nudge_x = -0.02, direction = "y",
    arrow = arrow (length = unit (0.02, "npc"), ends = "last", angle = 0),
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines"),
    hjust=0

  )+
  geom_label_repel (    data = highresults,
                        aes(label = Gene,colour=text),
                        size =3,max.overlaps = Inf,segment.size=0.4,
                        nudge_x= 0.02, direction = "y",
                        arrow = arrow (length = unit (0.02, "npc"), ends = "last", angle = 0),
                        box.padding = unit(0.1, "lines"),
                        point.padding = unit(0.1, "lines"),
                        hjust=0)+scale_color_manual(values = c('mitotic factor'='darkblue','replication factor'='darkred',"others"="black","not significant"=tropiccolor[2], "significant"=tropiccolor[1]))
plot1=plot1+guides(color=guide_legend(title="significant (FDR < 0.05)",nrow=2,byrow=T,override.aes = list(shape =19)))
plot1
prow=plot_grid(plot1, align = "h", axis = "b",ncol=1,rel_widths = c(1,1), labels = c("WGII score vs gene expressions"),label_fontface = 'plain')
print(prow)
pdf('../holger_bastian/plot/volcano_wgii_geneexp_metattest.pdf',width =8,height = 6)
print(prow)
dev.off()