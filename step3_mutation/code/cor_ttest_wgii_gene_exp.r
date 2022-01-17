rm(list = ls())
library(arules)
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
print('start calculating')
corfn=function(organ){
  print(organ)
  usedmat=subset(Final_expression_table_tum,Final_expression_table_tum$Cohort==organ)
  usedmat$cin_type=discretize(usedmat$wgii,method = 'cluster',categories = 2)
  usedmat$cin_type=factor(usedmat$cin_type)
  levels(usedmat$cin_type) = c('low','high')
  Results <- data.frame()
  for (i in 1:(length(colnames(rna_table))-1)) {
    
    gene=colnames(rna_table)[i]
    print(gene)
    if(length(which(is.na(usedmat[, gene])))<=10)
    {
      if(var(usedmat[, gene],na.rm = T)!=0)
      {
        usedmat$selected_term=(usedmat[, gene]-min(usedmat[, gene]))/(max(usedmat[, gene])-min(usedmat[, gene]))
        hind=which(usedmat$cin_type=='high')
        lind=which(usedmat$cin_type=='low')
     
        hi=usedmat$selected_term[hind]
        li=usedmat$selected_term[lind]
        m1i=mean(hi, na.rm=TRUE)
        m2i=mean(li, na.rm=TRUE)
        sd1i=sd(hi, na.rm=TRUE)
        sd2i=sd(li, na.rm=TRUE)
        n1i=length(hi)
        n2i=length(li)
        wt=t.test(hi,li)
        pvalue =wt$p.value
        coef = mean(usedmat$selected_term[hind])-mean(usedmat$selected_term[lind])
        Results = rbind(Results, data.frame(gene, coef, pvalue,m1i, sd1i, n1i, m2i, sd2i, n2i))
      }
    }
  }
  Results$FDR <- p.adjust(Results$pvalue, method = "fdr")
  Results <- Results[order(Results$pvalue),]
  names(Results) <- c("Gene", "Coef", "Pvalue", "m1i","sd1i","n1i","m2i","sd2i","n2i","FDR")
  write.table(Results, paste0(c("table/wgii_gene_exp_ttest_",organ,".txt"),collapse = ''), quote=F, sep="\t", row.names = F)
}
args = commandArgs(trailingOnly=TRUE)
corfn(organ=args[1])