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
Final_expression_table_tum=merge(Final_table, table, by='Sample')
# Final_expression_table_tum=Final_expression_table_tum[which(Final_expression_table_tum$Genome_doublings==0),]
Mut_table <- fread("data/data_CNA.oncokb.txt")
Mut_table=subset(Mut_table,Mut_table$ALTERATION=='Deletion',)
Mut_table$Sample=substr(Mut_table$SAMPLE_ID,1,15)
Results <- data.frame()
unigene=unique(Mut_table$HUGO_SYMBOL)
for(i in 1:length(unigene)){
  gene <-unigene[i]
  # at least 20 samples per group
  print(i)
  sampleif=Mut_table$Sample[Mut_table$HUGO_SYMBOL %in% gene]
  gmut=which(Final_expression_table_tum$Sample %in% sampleif)
  gwt=which(!Final_expression_table_tum$Sample %in% sampleif)
  if(length(gmut) >=20 & length(gwt) >=20){
    
    Final_expression_table_tum$Mutation_selected <- factor(Final_expression_table_tum$Sample %in% sampleif,labels = c('wt','mut'))
    
    fit <- lm(wgii ~ Mutation_selected + Cohort, Final_expression_table_tum)
    tmpcoef=summary(fit)$coefficients
    pvalue = tmpcoef[2,4]  
    # delta mean
    coef = tmpcoef[2,1]
    
    wt <- mean(Final_expression_table_tum$wgii[gwt])
    mut <- mean(Final_expression_table_tum$wgii[gmut])
    
    WT_samples <- length(gwt)
    MUT_samples <- length(gmut)
    Results=rbind(Results, data.frame(gene, coef, pvalue, wt, mut, mut-wt, WT_samples, MUT_samples))
  }
}
Results$FDR <- p.adjust(Results$pvalue, method = "fdr")
Results <- Results[order(Results$pvalue),c(1:3,9,4:8)]
names(Results) <- c("Gene", "Coef", "Pvalue", "FDR", "WTmean", "MUTmean", "MUTminusWT_diff", "WT_samples", "MUT_samples")
write.table(Results, "table/allsamples_wgii_allDeletions_linear_model_results.txt", quote=F, sep="\t", row.names = F)