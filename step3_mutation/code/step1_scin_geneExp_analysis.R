rm(list = ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library(data.table)
setwd("/data/zhang/pancan_cin/step3_mutation/")
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
Final_table = scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer', 'Primary Solid Tumor'), ]
Final_table$Sample <- substr(Final_table$Sample, 1, 15)
rownames(Final_table) = Final_table$Sample
rna_table=fread("../step8_tme/data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv")
ii = unlist(lapply(rna_table$gene_id, function(x) {
  strsplit(as.character(x), '\\|')[[1]][1]
}))
iidx = which(!duplicated(ii))
ii = ii[iidx]
rna_table = rna_table[iidx,]
rownames(rna_table) = ii
rna_table = rna_table[which(rownames(rna_table) != "?"),]
rna_table$gene_id = NULL
rna_table = t(as.data.frame(rna_table))
colnames(rna_table) = ii[-1]
rna_table = data.frame(rna_table)
rna_table$Sample =  substr(rownames(rna_table), 1, 15)
rna_table = rna_table[!duplicated(rna_table$Sample), ]
rownames(rna_table) = rna_table$Sample
Final_expression_table_tum = merge(Final_table, rna_table, by = 'Sample')
Results <- data.frame()
for (i in 1:(length(colnames(rna_table))-1)) {
  gene=colnames(rna_table)[i]
  Final_expression_table_tum$Mutation_selected=Final_expression_table_tum[, gene]
  fit=lm(scin ~ Mutation_selected + Cohort, Final_expression_table_tum)
  pvalue = summary(fit)$coefficients[2, 4]
  coef = summary(fit)$coefficients[2, 1]
  Results = rbind(Results, data.frame(gene, coef, pvalue))
}
Results$FDR <- p.adjust(Results$pvalue, method = "fdr")
Results <- Results[order(Results$pvalue),]
names(Results) <- c("Gene", "Coef", "Pvalue", "FDR")
write.table(Results, "table/scin_allGeneExp_linear_model_results.txt", quote=F, sep="\t", row.names = F)