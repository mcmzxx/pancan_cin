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

Final_expression_table_tum=Final_table[sels,]
Results <- read.delim("table/wgii_allMutations_linear_model_results.txt")
Results1 <- read.delim("../step3_mutation/table/wgii_allGeneExp_linear_model_results.txt")
Results <- Results[order(Results$Pvalue),]
Results <- Results[!Results$Gene %in% ".",]
Results$Significant2 <- ifelse(Results$FDR < 0.05 & Results$Coef > 0, "Pos", "Not Sig")
Results$Significant2[Results$FDR < 0.05 & Results$Coef < 0] <- "Neg"

Results1 <- Results1[order(Results1$Pvalue),]
Results1 <- Results1[!Results1$Gene %in% ".",]
Results1$Significant2 <- ifelse(Results1$FDR < 0.05 & Results1$Coef > 0, "Pos", "Not Sig")
Results1$Significant2[Results1$FDR < 0.05 & Results1$Coef < 0] <- "Neg"
colnames(Results1)=paste0(colnames(Results1),'_geneExp')
Results1$Gene=Results1$Gene_geneExp
Results=merge(Results,Results1,by='Gene')
Results=Results[which(Results$Significant2==Results$Significant2_geneExp),]
Results=Results[which(Results$Significant2!='Not Sig'),]
write.table(Results, "table//mutation_dependent_mechnisim_wgii.txt", quote=F, sep="\t", row.names = F)


