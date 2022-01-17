rm(list=ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library(data.table)
setwd("/data/zhang/pancan_cin/step6_compounds/")
load('data/ccle_cin_score.RData')
scores$Genome.doublings=as.factor(scores$Genome.doublings)
data=scores[which(!is.na(scores$tcga_code)),]
data$Cohort=as.character(data$tcga_code)
data=data[which(complete.cases(data$Genome.doublings)),]
Final_expression_table_tum=data[which(data$tcga_code!="UNABLE TO CLASSIFY"),]
Mut_table <- fread("data/CCLE_DepMap_18q3_maf_20180718.txt")
Mut_table$Sample=Mut_table$Tumor_Sample_Barcode
Results <- data.frame()
for(i in 1:length(unique(Mut_table$Hugo_Symbol))){
  gene <- paste(unique(Mut_table$Hugo_Symbol)[i])
  if(length(which(Final_expression_table_tum$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene])) >=20 & length(which(!Final_expression_table_tum$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene])) >=20){
    Final_expression_table_tum$Mutation_selected <- Final_expression_table_tum$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene]
    fit <- lm(scin ~ Mutation_selected + Cohort, Final_expression_table_tum)
    pvalue = summary(fit)$coefficients[2,4]  
    coef = summary(fit)$coefficients[2,1]
    wt <- mean(Final_expression_table_tum$scin[!Final_expression_table_tum$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene]])
    mut <- mean(Final_expression_table_tum$scin[Final_expression_table_tum$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene]])
    WT_samples <- length(which(!Final_expression_table_tum$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene]))
    MUT_samples <- length(which(Final_expression_table_tum$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene]))
    Results=rbind(Results, data.frame(gene, coef, pvalue, wt, mut, mut-wt, WT_samples, MUT_samples))
  }
}
Results$FDR <- p.adjust(Results$pvalue, method = "fdr")
Results <- Results[order(Results$pvalue),c(1:3,9,4:8)]
names(Results) <- c("Gene", "Coef", "Pvalue", "FDR", "WTmean", "MUTmean", "MUTminusWT_diff", "WT_samples", "MUT_samples")
write.table(Results, "table/scin_allMutations_linear_model_results.txt", quote=F, sep="\t", row.names = F)