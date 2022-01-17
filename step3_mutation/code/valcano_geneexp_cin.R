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
Results$Significant2 = ifelse(Results$FDR < 0.05 & Results$Coef > 0, "positive", "insignificant")
Results$Significant2[Results$FDR < 0.05 & Results$Coef < 0] = "negative"
Results$text = "no"
selt=unique(c(which(Results$Gene %in% c("CKS1B", "TPX2", "NEK2", "AURKA", "UBE2C", "CENPA", "CSE1L", 
                                            "ZNF334", "GINS1", "PKMYT1", "CABLES2", "STX1A.3", "STX1A", "DSN1", 
                                            "RAE1"))))
selt=intersect(selt,which(Results$Significant2!="insignificant"))
Results$text[selt] = "yes"
Results$text = as.factor(Results$text)
pwgii =  ggplot(Results, aes(x = Coef, y = -log10(FDR))) +
  geom_point(aes(color = Significant2),size=0.2) +
  xlab("correlation coefficient(meta-analysis)") +
  ylab("-log10(FDR adjusted p-value)") +
  scale_color_manual(values = c("negative"=tropiccolor[1], "insignificant"=tropiccolor[2], "positive"=tropiccolor[3])) +
  theme_cowplot(12)+ theme(axis.text.x=element_text(size=15, colour="black"),  legend.position = 'bottom',
                           axis.text.y=element_text(size=15, colour="black"), 
                           axis.title=element_text(size=15, colour="black"), 
                           legend.key = element_rect(fill="White"), 
                           legend.text = element_text(size=15), 
                           legend.title = element_text(size=15)) +
  geom_text_repel(
    data = subset(Results, text == "yes"),
    aes(label = Gene),
    size =3,
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines")
  ) +
  guides(color=guide_legend(title="significant(FDR<0.05)",ncol=2,override.aes = list(size=2)))


Results = read.delim("table/scin_allGeneExp_metacor_results.txt")
Results$Significant2 = ifelse(Results$FDR < 0.05 & Results$Coef > 0, "positive", "insignificant")
Results$Significant2[Results$FDR < 0.05 & Results$Coef < 0] = "negative"
Results$text = "no"
selt=unique(c(which(Results$Gene %in% c("CKS1B", "MYBL2", "DNAJB11", "SHC1.3", "SNRPG", "E2F1", "NEK2", 
                                            "AURKA", "VIF", "TAF4", "TBCE", "LIN9", "SKP2", "CDK4", "UBE2C"))))
selt=intersect(selt,which(Results$Significant2!="insignificant"))
Results$text[selt] = "yes"
Results$text = as.factor(Results$text)
pscin =  ggplot(Results, aes(x = Coef, y = -log10(FDR))) +
  geom_point(aes(color = Significant2),size=0.2) +
  xlab("correlation coefficient(meta-analysis)") +
  ylab("-log10(FDR adjusted p-value)") +
  scale_color_manual(values = c("negative"=tropiccolor[1], "insignificant"=tropiccolor[2], "positive"=tropiccolor[3])) +
  theme_cowplot(12)+ theme(axis.text.x=element_text(size=15, colour="black"),  legend.position = 'bottom',
                           axis.text.y=element_text(size=15, colour="black"), 
                           axis.title=element_text(size=15, colour="black"), 
                           legend.key = element_rect(fill="White"), 
                           legend.text = element_text(size=15), 
                           legend.title = element_text(size=15)) +
  geom_text_repel(
    data = subset(Results, text == "yes"),
    aes(label = Gene),
    size =3,
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines")
  ) +
  guides(color=guide_legend(title="significant(FDR<0.05)",ncol=2,override.aes = list(size=2)))
legend=get_legend(pscin+ theme(legend.key = element_rect(colour = NA, fill = NA),legend.box.margin = margin(-1, -1, -1, 0)))
pwgii=pwgii+theme(legend.position = 'none')
pscin=pscin+theme(legend.position = 'none')
prow=plot_grid(pwgii,pscin, align = "v", axis = "l",ncol=1,rel_heights= c(1,1), labels = c("B","D"),label_x = -0.005)
pdf('plot/volcano_cin_geneexp_metacor.pdf',width=5,height = 10)
plot_grid(prow,legend, align = "v",axis='c',ncol=1,rel_heights =c(10,0.7))
dev.off()
