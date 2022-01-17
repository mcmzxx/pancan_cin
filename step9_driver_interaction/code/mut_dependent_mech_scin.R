rm(list=ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library(data.table)
threshold=0.1
setwd("/data/zhang/pancan_cin/step9_driver_interaction/")
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer','Primary Solid Tumor'),]
Final_table$Sample <- substr(Final_table$Sample, 1, 15)
rownames(Final_table)=Final_table$Sample
Results <- read.delim("table/scin_allMutations_linear_model_results.txt")
Results1 <- read.delim("../step3_mutation/table/scin_allGeneExp_linear_model_results.txt")
Results <- Results[order(Results$Pvalue),]
Results <- Results[!Results$Gene %in% ".",]
Results$Significant2 <- ifelse(Results$Pvalue < 0.05 & Results$Coef > threshold, "Pos", "Not Sig")
Results$Significant2[Results$Pvalue < 0.05 & Results$Coef < -threshold] <- "Neg"

Results1 <- Results1[order(Results1$Pvalue),]
Results1 <- Results1[!Results1$Gene %in% ".",]
Results1$Significant2 <- ifelse(Results1$Pvalue < 0.05 & Results1$Coef > threshold, "Pos", "Not Sig")
Results1$Significant2[Results1$Pvalue < 0.05 & Results1$Coef < -threshold] <- "Neg"
colnames(Results1)=paste0(colnames(Results1),'_geneExp')
Results1$Gene=Results1$Gene_geneExp
Results=merge(Results,Results1,by='Gene')
Results=Results[which(Results$Significant2==Results$Significant2_geneExp),]
Results=Results[which(Results$Significant2!='Not Sig'),]
write.table(Results, "table//mutation_dependent_mechnisim_scin.txt", quote=F, sep="\t", row.names = F)
