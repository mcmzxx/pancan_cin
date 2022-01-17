rm(list=ls())
library(ggplot2)
library(ggsignif)
library(matrixStats)
library(ggrepel)
library(cowplot)
library(data.table)
setwd("/data/zhang/pancan_cin/step5_pathway/")
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer','Primary Solid Tumor'),]
Final_table$Sample <- substr(Final_table$Sample, 1, 12)
rownames(Final_table)=Final_table$Sample
imm=fread('data/merge_merged_reals.txt')
ii=imm$Gene
rownames(imm)=imm$Gene
imm$Gene=NULL
imm=t(as.data.frame(imm))
colnames(imm)=ii
imm=data.frame(imm)
imm$Sample=rownames(imm)
# Final_table=merge(Final_table,til,by='Sample')
Final_table=merge(Final_table,imm,by='Sample')
tmpc=unique(Final_table$Cohort)
ii=colnames(Final_table)[51:length(colnames(Final_table))]

reslist=list()
for(j in tmpc[21:23]){
  Results <- data.frame()
  Results1 <- data.frame()
  tmpmat=subset(Final_table,Final_table$Cohort==j)
  for(i in ii){
    if(length(unique(tmpmat[,i]))>5){
      cell_type=i
      cor_p <- cor.test(tmpmat$wgii, tmpmat[,i], method = "pearson")
      cor_s <- cor.test(tmpmat$wgii, tmpmat[,i], method = "spearman")
      Results=rbind(Results, data.frame(j,cell_type, as.numeric(cor_s$estimate), as.numeric(cor_s$p.value), as.numeric(cor_p$estimate), as.numeric(cor_p$p.value)))
      cor_p <- cor.test(tmpmat$scin, tmpmat[,i], method = "pearson")
      cor_s <- cor.test(tmpmat$scin, tmpmat[,i], method = "spearman")
      Results1=rbind(Results1, data.frame(j,cell_type, as.numeric(cor_s$estimate), as.numeric(cor_s$p.value), as.numeric(cor_p$estimate), as.numeric(cor_p$p.value)))
    }}
  # Results$FDR_s <- p.adjust(Results$as.numeric.cor_s.p.value., method = "fdr")
  # Results$FDR_p <- p.adjust(Results$as.numeric.cor_p.p.value., method = "fdr")
  names(Results) <- c("Cohort","Pathway","Spearman_coef", "Spearman_pvalue", "Pearson_coef", "Pearson_pvalue")#,"Spearman_FDR", "Pearson_FDR")
  names(Results1) <- c("Cohort","Pathway","Spearman_coef", "Spearman_pvalue", "Pearson_coef", "Pearson_pvalue")#,"Spearman_FDR", "Pearson_FDR")
  save(Results,file = file.path('table/paradigm_pathway',paste0(c('cor_paradigm_wgii_',j,'.RData'),collapse = '')))
  save(Results1,file = file.path('table/paradigm_pathway',paste0(c('cor_paradigm_scin_',j,'.RData'),collapse = '')))
}