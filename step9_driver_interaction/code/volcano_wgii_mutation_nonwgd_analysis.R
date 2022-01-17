rm(list=ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library(data.table)
setwd("/data/zhang/pancan_cin/step9_driver_interaction/")
#Load Table S2 from data from Taylor et al, Cancer Cell 2018 (https://www.sciencedirect.com/science/article/pii/S1535610818301119)
table <- read.delim("../step2_genome_instability/data/TaylorCancerCell_TableS2.txt", na.strings = c("", " ", "na", "NA", "n.a.", "#N/A"))
rownames(table) <- as.character(table$Sample)
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer','Primary Solid Tumor'),]
Final_table$Sample <- substr(Final_table$Sample, 1, 15)
rownames(Final_table)=Final_table$Sample
sels=intersect(rownames(table),rownames(Final_table))
# Mut_table <- fread("../step8_tme/data/mc3.v0.2.8.PUBLIC.maf")
# Mut_table$Sample=substr(Mut_table$Tumor_Sample_Barcode,1,15)
Final_expression_table_tum=Final_table[sels,]
Results <- read.delim("table/scin_allAmplifications_linear_model_results.txt")
Results <- Results[order(Results$Pvalue),]
Results <- Results[!Results$Gene %in% ".",]
Results$text <- "no"
Results$text[c(1:10,grep('PIK3CA',Results$Gene))] <- "yes"
Results$text <- as.factor(Results$text)
Results$Significant <- ifelse(Results$FDR < 0.05 & Results$MUTminusWT_diff > 0, "Pos", "Not Sig")
Results$Significant[Results$FDR < 0.05 & Results$MUTminusWT_diff < 0] <- "Neg"

Results$Significant2 <- ifelse(Results$FDR < 0.05 & Results$Coef > 0, "Pos", "Not Sig")
Results$Significant2[Results$FDR < 0.05 & Results$Coef < 0] <- "Neg"
Results[grep('PTEN',Results$Gene),]
pdf("plot/lm_scin_amplification.pdf")
plot <- ggplot(Results, aes(x = Coef, y = -log10(Pvalue))) +
  geom_point(aes(color = Significant2)) +
  xlab("Coefficient (linear model)") +
  ylab("-log10(p-value)") +
  scale_color_manual(values = c("Neg"=tropiccolor[1], "Not Sig"=tropiccolor[1], "Pos"=tropiccolor[3])) +
  theme_bw(base_size = 12) + theme(axis.text.x=element_text(size=13, colour="black"), axis.text.y=element_text(size=13, colour="black"), axis.title=element_text(size=16, colour="black"), legend.key = element_rect(fill="White"), legend.text = element_text(size=13), legend.title = element_text(size=14), legend.position = "bottom") +
  geom_text_repel(
    data = subset(Results, text == "yes"),
    aes(label = Gene),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) +
  guides(color=guide_legend(title="Significant (FDR < 0.05)"))
print(plot)
dev.off()