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
Final_expression_table_bak = merge(Final_table, rna_table, by = 'Sample')

organs=c("ACC", 
         "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", 
         "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", 
         "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", 
         "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

finalresults=list()
for(organ in organs){
  Final_expression_table_tum=subset(Final_expression_table_bak,Final_expression_table_bak$Cohort==organ)
  Results <- data.frame()
  for (i in 1:(length(colnames(rna_table))-1)) {
    print(i)
    gene=colnames(rna_table)[i]
    if(!all(is.na(Final_expression_table_tum[, gene])))
    {if(var(Final_expression_table_tum[, gene],na.rm = T)!=0)
    {
      Final_expression_table_tum$Mutation_selected=ifelse(Final_expression_table_tum[, gene]<=median(Final_expression_table_tum[, gene],na.rm = T), "Lower", "Higher")
      Final_expression_table_tum$Mutation_selected <- factor(Final_expression_table_tum$Mutation_selected, levels = c("Lower", "Higher"))
      fit=lm(wgii ~ Mutation_selected, Final_expression_table_tum)
      pvalue = summary(fit)$coefficients[2, 4]
      coef = summary(fit)$coefficients[2, 1]
      Results = rbind(Results, data.frame(gene, coef, pvalue))
    }
    }
  }
  Results$FDR <- p.adjust(Results$pvalue, method = "fdr")
  Results <- Results[order(Results$pvalue),]
  names(Results) <- c("Gene", "Coef", "Pvalue", "FDR")
  Results$Organ=organ
  finalresults[[organ]]=Results
}
save(finalresults,file ='table/wgii_gexp.RData')