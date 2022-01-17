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
cnv_table=fread("../step8_tme/data/all_data_by_genes_whitelisted.tsv")
cnv_table=cnv_table[which(!duplicated(cnv_table$`Gene Symbol`)),]
ii= cnv_table$`Gene Symbol`
cnv_table[,1:3] = NULL
cnv_table=data.frame(cnv_table)
rownames(cnv_table)=ii
cnv_table=cnv_table[,which(substr(colnames(cnv_table), 14, 15)!=11)]
cnv_table = data.frame(t(cnv_table))
rownames(cnv_table)=unlist(lapply(rownames(cnv_table), function(x){gsub('\\.','-',x)}))
cnv_table$Sample =  substr(rownames(cnv_table), 1, 15)
cnv_table = cnv_table[!duplicated(cnv_table$Sample), ]
rownames(cnv_table) = cnv_table$Sample
Final_expression_table_tum = merge(Final_table, cnv_table, by = 'Sample')
print('start calculating')
corfn=function(organ){
  print(organ)
  usedmat=subset(Final_expression_table_tum,Final_expression_table_tum$Cohort==organ)
  usedmat$cin_type=discretize(usedmat$wgii,method = 'cluster',categories = 2)
  usedmat$cin_type=factor(usedmat$cin_type)
  levels(usedmat$cin_type) = c('low','high')
  Results <- data.frame()
  for (i in 1:(length(colnames(cnv_table))-1)) {
    
    gene=colnames(cnv_table)[i]
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
  write.table(Results, paste0(c("table/wgii_cnv_ttest_",organ,".txt"),collapse = ''), quote=F, sep="\t", row.names = F)
}
args = commandArgs(trailingOnly=TRUE)
corfn(organ=args[1])
