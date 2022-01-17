rm(list=ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library(data.table)
library("colorspace")
tropiccolor=c("#0374be","grey69","#cc0f0f")
setwd("/data/zhang/pancan_cin/step9_driver_interaction/")
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer','Primary Solid Tumor'),]
Final_table$Sample <- substr(Final_table$Sample, 1, 15)
rownames(Final_table)=Final_table$Sample
Results <- read.delim("table/wgii_allAmplifications_linear_model_results.txt")
Results <- Results[order(Results$Pvalue),]
Results <- Results[!Results$Gene %in% ".",]
Results$Significant = ifelse(Results$FDR < 0.05 & Results$MUTminusWT_diff > 0, "positive", "insignificant")
Results$Significant[Results$FDR < 0.05 & Results$MUTminusWT_diff < 0] = "negative"
Results$Significant2 = ifelse(Results$FDR < 0.05 & Results$Coef > 0, "positive", "insignificant")
Results$Significant2[Results$FDR < 0.05 & Results$Coef < 0] = "negative"
Results$text = "no"
selt=unique(c(1:5,which(Results$Gene %in% c('PIK3CA','KRAS','BRAF','NRAS', 'HRAS', "AURKB","BUB1",'MYC','RB1','BUB3','AURKA','PTEN')),which(Results$Significant2 %in% c("positive"))[1:5]))
selt=intersect(selt,which(Results$Significant2!="insignificant"))
Results$text[selt] = "yes"
Results$text = as.factor(Results$text)
plot1 =  ggplot(Results, aes(x = Coef, y = -log10(Pvalue))) +
  geom_point(aes(color = Significant2),size=0.2) +
  xlab("coefficient (linear model)") +
  ylab("-log10(p-value)") +
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
  guides(color=guide_legend(title="significant(FDR<0.05)",override.aes = list(size=2)))

Results <- read.delim("table/scin_allAmplifications_linear_model_results.txt")
Results <- Results[order(Results$Pvalue),]
Results <- Results[!Results$Gene %in% ".",]
Results$Significant = ifelse(Results$FDR < 0.05 & Results$MUTminusWT_diff > 0, "positive", "insignificant")
Results$Significant[Results$FDR < 0.05 & Results$MUTminusWT_diff < 0] = "negative"
Results$Significant2 = ifelse(Results$FDR < 0.05 & Results$Coef > 0, "positive", "insignificant")
Results$Significant2[Results$FDR < 0.05 & Results$Coef < 0] = "negative"
Results$text = "no"
selt=unique(c(1:5,which(Results$Gene %in% c('PIK3CA','KRAS','BRAF','NRAS', 'HRAS', "AURKB","BUB1",'MYC','RB1','BUB3','AURKA','PTEN')),which(Results$Significant2 %in% c("positive"))[1:5]))
selt=intersect(selt,which(Results$Significant2!="insignificant"))
Results$text[selt] = "yes"
Results$text = as.factor(Results$text)
plot2 = ggplot(Results, aes(x = Coef, y = -log10(Pvalue))) +
  geom_point(aes(color = Significant2),size=0.2) +
  xlab("coefficient (linear model)") +
  ylab("-log10(p-value)") +
  scale_color_manual(values = c("negative"=tropiccolor[1], "insignificant"=tropiccolor[2], "positive"=tropiccolor[3])) +
  theme_cowplot(12)+ theme(axis.text.x=element_text(size=15, colour="black"), legend.position = 'bottom',
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
  guides(color=guide_legend(title="significant(FDR<0.05)",override.aes = list(size=2)))
legend <- get_legend(
  plot2 + theme(legend.box.margin = margin(-1, -1, -1, 0))
)

plot1=plot1+theme(legend.position = "none",plot.margin=margin(r=-0.1,t=0.4,l=0.1,unit="cm"))
plot2=plot2+theme(legend.position = "none",plot.margin=margin(t=0.4,r=0.2,unit="cm"))
prow=plot_grid(plot1,plot2, align = "h", axis = "b",ncol=2,rel_widths = c(1,1), labels = c("WGII score vs amplifications","SCIN score vs amplifications"),label_fontface = 'plain',label_x = -0.1)
pdf('plot/lm_cin_amplification.pdf',width =7,height = 6)
plot_grid(prow,legend, align = "v",ncol=1,rel_heights=c(1,0.1))
dev.off()