rm(list=ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library(data.table)
library("colorspace")
tropiccolor=c("#0374be","grey69","#cc0f0f")
mycolor=c('darkblue','white','darkred','black')
setwd("/data/zhang/pancan_cin/step3_mutation/")
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer','Primary Solid Tumor'),]
Final_table$Sample = substr(Final_table$Sample, 1, 15)
rownames(Final_table)=Final_table$Sample
Final_table$Genome_doublings = as.factor(Final_table$`Genome doublings`)
Results = read.delim("table/wgii_allGeneExp_metacor_results.txt")
Results = Results[order(Results$Pvalue),]
Results = Results[!Results$Gene %in% ".",]

Results$significant = ifelse(Results$FDR < 0.05 & Results$Coef > 0, "significant", "insignificant")
Results$significant[Results$FDR < 0.05 & Results$Coef < 0] = "significant"
Results$significant2 = ifelse(Results$FDR < 0.05 & Results$Coef > 0, "postive", "insignificant")
Results$significant2[Results$FDR < 0.05 & Results$Coef < 0] = "negtive"
Results$text = "no"
selt=unique(c(which(Results$Gene %in% c("CKS1B", "TPX2", "NEK2", "AURKA", "UBE2C", "CENPA", "CSE1L", 
                                        "ZNF334", "GINS1", "PKMYT1", "CABLES2", "STX1A.3", "STX1A", "DSN1", 
                                        "RAE1"))))
#selt=intersect(selt,which(Results$significant2!='insignificant'))
# mitotic=which(Results$Gene %in% c("RAE1","TPX2","UBE2C","AURKA","AURB","BUB1","CDK1"))
# mitotic=intersect(mitotic,which(Results$significant2!='insignificant'))
replication=which(Results$Gene %in% c("GINS1","CDC45","MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7","DBF4","CDC7","GINS3","GINS2"))
replication=intersect(replication,which(Results$significant2!='insignificant'))
Results$text[selt] = "PARADIGM feature"
Results$text[replication] = "replication factor"
Results$text = as.factor(Results$text)
represults=subset(Results, text %in% c("replication factor"))
represults=subset(represults,!(Gene %in% c("GINS1","CDC45")))
otherresults=subset(Results, text %in% c("PARADIGM feature"))
highresults=subset(Results, Gene %in% c("GINS1","CDC45"))
pwgii =ggplot(Results, aes(x = Coef, y = -log10(FDR))) +
  geom_point(aes(colour = significant),shape=19,size=0.2) +
  xlab("correlation coefficient(meta-analysis)") +
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
    nudge_x= 0.2, direction = "y",angle=0,
    arrow = arrow (length = unit (0.02, "npc"), ends = "last", angle = 0),
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines"),
    hjust=0
  )+
  geom_text_repel(
    data = otherresults,
    aes(label = Gene,colour=text),
    size =3,max.overlaps = Inf,segment.size=0.2,
    nudge_x = -0.2, direction = "y",
    arrow = arrow (length = unit (0.02, "npc"), ends = "last", angle = 0),
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines"),
    hjust=0
  )+
  geom_label_repel (    data = highresults,
                        aes(label = Gene,colour=text),
                        size =3,max.overlaps = Inf,segment.size=0.4,
                        nudge_x= -0.01, direction = "x",
                        arrow = arrow (length = unit (0.02, "npc"), ends = "last", angle = 0),
                        box.padding = unit(0.1, "lines"),
                        point.padding = unit(0.1, "lines"),
                        hjust=0)+scale_color_manual(values = c('PARADIGM feature'='#006747FF','replication factor'='#4F2C1DFF',"insignificant"=tropiccolor[2], "significant"=tropiccolor[1]))
pwgii=pwgii+guides(color=guide_legend(title="significant (FDR < 0.05)",nrow=2,byrow=T,override.aes = list(shape =19)))
pwgii


Results = read.delim("table/scin_allGeneExp_metacor_results.txt")
Results = Results[order(Results$Pvalue),]
Results = Results[!Results$Gene %in% ".",]

Results$significant = ifelse(Results$FDR < 0.05 & Results$Coef > 0, "significant", "insignificant")
Results$significant[Results$FDR < 0.05 & Results$Coef < 0] = "significant"
Results$significant2 = ifelse(Results$FDR < 0.05 & Results$Coef > 0, "postive", "insignificant")
Results$significant2[Results$FDR < 0.05 & Results$Coef < 0] = "negtive"
Results$text = "no"
selt=unique(c(which(Results$Gene %in% c("CKS1B", "MYBL2", "DNAJB11", "SHC1.3", "SNRPG", "E2F1", "NEK2", 
                                        "AURKA", "VIF", "TAF4", "TBCE", "LIN9", "SKP2", "CDK4", "UBE2C"))))
#"TPX2", "E2F1", "UBE2C",'PIK3CA','KRAS','BRAF','NRAS', 'HRAS', "AURKB","BUB1",'MYC','RB1','BUB3','AURKA','PTEN'
#selt=intersect(selt,which(Results$significant2!='insignificant'))
replication=which(Results$Gene %in% c("GINS1","CDC45","MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7","DBF4","CDC7","GINS3","GINS2"))
replication=intersect(replication,which(Results$significant2!='insignificant'))
Results$text[selt] = "PARADIGM feature"
Results$text[replication] = "replication factor"
Results$text = as.factor(Results$text)
represults=subset(Results, text %in% c("replication factor"))
represults=subset(represults,!(Gene %in% c("GINS1","CDC45")))
otherresults=subset(Results, text %in% c("PARADIGM feature"))
highresults=subset(Results, Gene %in% c("GINS1","CDC45"))
pscin =ggplot(Results, aes(x = Coef, y = -log10(FDR))) +
  geom_point(aes(colour = significant),shape=19,size=0.2) +
  xlab("correlation coefficient(meta-analysis)") +
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
    nudge_x= 0.2, direction = "y",angle=0,
    arrow = arrow (length = unit (0.02, "npc"), ends = "last", angle = 0),
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines"),
    hjust=0
  )+
  geom_text_repel(
    data = otherresults,
    aes(label = Gene,colour=text),
    size =3,max.overlaps = Inf,segment.size=0.2,
    nudge_x = -0.2, direction = "y",
    arrow = arrow (length = unit (0.02, "npc"), ends = "last", angle = 0),
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines"),
    hjust=0
  )+
  geom_label_repel (    data = highresults,
                        aes(label = Gene,colour=text),
                        size =3,max.overlaps = Inf,segment.size=0.4,
                        nudge_x= -0.3, direction = "x",
                        arrow = arrow (length = unit (0.02, "npc"), ends = "last", angle = 0),
                        box.padding = unit(0.1, "lines"),
                        point.padding = unit(0.1, "lines"),
                        hjust=0)+scale_color_manual(values = c('PARADIGM feature'='#006747FF','replication factor'='#4F2C1DFF',"insignificant"=tropiccolor[2], "significant"=tropiccolor[1]))
pscin=pscin+guides(color=guide_legend(title="significant (FDR < 0.05)",nrow=2,byrow=T,override.aes = list(shape =19)))
pscin


legend1=get_legend(pscin+ theme(legend.key = element_rect(colour = NA, fill = NA),legend.box.margin = margin(-1, -1, -1, 0)))
pwgii=pwgii+theme(legend.position = 'none')
pscin=pscin+theme(legend.position = 'none')
prow=plot_grid(pwgii,pscin, align = "v", axis = "l",ncol=1,rel_heights= c(1,1), labels = c("B","D"),label_x = -0.005)
pdf('plot/volcano_cin_geneexp_metacor_of.pdf',width=5,height = 9.35)
#plot_grid(prow,legend, align = "v",axis='c',ncol=1,rel_heights =c(10,0.7))
print(prow)
dev.off()
dumdata=data.frame(x=1:10,y=1:10)
dumplot=ggplot(dumdata,aes(x,y))+geom_point()
legendplot=plot_grid(dumplot,legend1, align = "v",axis='c',ncol=1,rel_heights =c(0.1,0.7))
pdf('plot/volcano_cin_geneexp_metacor_legend.pdf',width=4.39,height = 0.66)
print(legendplot)
dev.off()

