rm(list=ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library(data.table)
setwd("/data/zhang/pancan_cin/step3_mutation/")
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer','Primary Solid Tumor'),]
Final_table$Sample <- substr(Final_table$Sample, 1, 15)
rownames(Final_table)=Final_table$Sample
msitable <- readxl::read_xlsx("../step2_genome_instability/data/1-s2.0-S009286741830237X-mmc1.xlsx",skip=8 ,sheet = 'Table S5')
msitable$Sample=msitable$TCGA_Barcode
msitable=as.data.frame(msitable)
ii=c( "Tota_Mutation_cnt", "MSIscore", 
      "POLE", "MLH1", "MLH3", "MGMT", "MSH6", "MSH3", "MSH2", "PMS1", 
      "PMS2")
msitable[,ii]=apply(msitable[,ii],2,as.numeric)

Final_expression_table_tum=Final_table[which(!(Final_table$Patient %in% msitable$TCGA_Barcode)),]
Mut_table <- fread("../step8_tme/data/mc3.v0.2.8.PUBLIC.maf.gz")
Mut_table$Sample=substr(Mut_table$Tumor_Sample_Barcode,1,15)
Results <- data.frame()
for(i in 1:length(unique(Mut_table$Hugo_Symbol))){
  gene <- paste(unique(Mut_table$Hugo_Symbol)[i])
  # at least 20 samples per group
  if(length(which(Final_expression_table_tum$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene])) >=20 & length(which(!Final_expression_table_tum$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene])) >=20){
    
    Final_expression_table_tum$Mutation_selected <- Final_expression_table_tum$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene]
    
    fit <- lm(wgii ~ Mutation_selected + Cohort, Final_expression_table_tum)
    
    pvalue = summary(fit)$coefficients[2,4]  
    # delta mean
    coef = summary(fit)$coefficients[2,1]
    
    wt <- mean(Final_expression_table_tum$wgii[!Final_expression_table_tum$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene]])
    mut <- mean(Final_expression_table_tum$wgii[Final_expression_table_tum$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene]])
    
    WT_samples <- length(which(!Final_expression_table_tum$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene]))
    MUT_samples <- length(which(Final_expression_table_tum$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene]))
    
    Results=rbind(Results, data.frame(gene, coef, pvalue, wt, mut, mut-wt, WT_samples, MUT_samples))
    
  }
}
Results$FDR <- p.adjust(Results$pvalue, method = "fdr")
Results <- Results[order(Results$pvalue),c(1:3,9,4:8)]
names(Results) <- c("Gene", "Coef", "Pvalue", "FDR", "WTmean", "MUTmean", "MUTminusWT_diff", "WT_samples", "MUT_samples")
write.table(Results, "table/wgii_allMutations_linear_model_results_nmsi.txt", quote=F, sep="\t", row.names = F)