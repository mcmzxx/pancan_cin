rm(list=ls())
library(ggplot2)
library(ggsignif)
library(matrixStats)
library(ggrepel)
library(cowplot)
library(data.table)
library(readxl)
setwd("/data/zhang/pancan_cin/step8_tme/")
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer','Primary Solid Tumor'),]
rownames(Final_table)=Final_table$Sample
imm=fread('data/Scores_160_Signatures.tsv')
ii=imm$SetName
imm=imm[,-c(1:3)]
imm=t(imm)
jj=rownames(imm)
imm=apply(imm,2,function(x){as.numeric(trimws(x))})
imm=data.frame(imm)
jj=substr(jj,1,15)
idx=which(!duplicated(jj))
imm=imm[idx,]
rownames(imm)=jj[idx]
colnames(imm)=ii
imm$Sample=rownames(imm)
Final_table=merge(Final_table,imm,by='Sample')
tmpc=unique(Final_table$Cohort)
reslist=list()
reslist1=list()
for(j in tmpc){
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
  Results$FDR_s <- p.adjust(Results$as.numeric.cor_s.p.value., method = "fdr")
  Results$FDR_p <- p.adjust(Results$as.numeric.cor_p.p.value., method = "fdr")
  Results <- Results[order(Results$as.numeric.cor_s.p.value.),]
  names(Results) <- c("Cohort","Cell_type","Spearman_coef", "Spearman_pvalue", "Pearson_coef", "Pearson_pvalue", "Spearman_FDR", "Pearson_FDR")
  reslist[[j]]=Results
  
  Results1$FDR_s <- p.adjust(Results1$as.numeric.cor_s.p.value., method = "fdr")
  Results1$FDR_p <- p.adjust(Results1$as.numeric.cor_p.p.value., method = "fdr")
  Results1 <- Results1[order(Results1$as.numeric.cor_s.p.value.),]
  names(Results1) <- c("Cohort","Cell_type","Spearman_coef", "Spearman_pvalue", "Pearson_coef", "Pearson_pvalue", "Spearman_FDR", "Pearson_FDR")
  reslist1[[j]]=Results1
}
Results=data.frame(do.call(rbind,reslist))
Results2 <-Results[order(Results$Spearman_pvalue),]
Results2$text <- "no"
Results2$text[ abs(Results2$Spearman_coef)>0.3] <- "yes"
Results2$text <- as.factor(Results2$text)
Results2$Significant <- ifelse(Results2$Spearman_FDR < 0.05 & Results2$Spearman_coef > 0, "Pos", "Not Sig")
Results2$Significant[Results2$Spearman_FDR < 0.05 & Results2$Spearman_coef < 0] <- "Neg"
write.table(Results2, "table/immune_corr_results_wgii.txt", quote=F, sep="\t", row.names = F)
Results=data.frame(do.call(rbind,reslist1))
Results2 <-Results[order(Results$Spearman_pvalue),]
Results2$text <- "no"
Results2$text[ abs(Results2$Spearman_coef)>0.3] <- "yes"
Results2$text <- as.factor(Results2$text)
Results2$Significant <- ifelse(Results2$Spearman_FDR < 0.05 & Results2$Spearman_coef > 0, "Pos", "Not Sig")
Results2$Significant[Results2$Spearman_FDR < 0.05 & Results2$Spearman_coef < 0] <- "Neg"
write.table(Results2, "table/immune_corr_results_scin.txt", quote=F, sep="\t", row.names = F)