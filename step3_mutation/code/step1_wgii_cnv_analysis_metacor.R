rm(list = ls())
library(arules)
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library(data.table)
library("robumeta")
library("metafor")
library("dplyr")
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
clist=split(Final_expression_table_tum,Final_expression_table_tum$Cohort)
cname=c("GBM", "OV", "LUAD", "LUSC", "PRAD", "UCEC", "BLCA", "TGCT", 
        "ESCA", "PAAD", "KIRP", "LIHC", "CESC", "SARC", "BRCA", "THYM", 
        "MESO", "COAD", "STAD", "CHOL", "KIRC", "THCA", "HNSC", "LAML", 
        "READ", "SKCM", "LGG", "DLBC", "KICH", "UCS", "ACC", "PCPG", 
        "UVM")
dflist=list()
for(i in cname){
  fn=paste0(c('table/wgii_cnv_ttest_',i,'.txt'),collapse = '')
  df=fread(fn)
  df$ni=dim(clist[[i]])[1]
  df=data.frame(df)
  rownames(df)=df$Gene
  dflist[[i]]=df
}
allgenes=Reduce(intersect,lapply(dflist, function(x){return(x$Gene)}))
Results <- data.frame()
for(gene in allgenes){
  print(gene)
  dat=data.frame(do.call(rbind,lapply(dflist, function(x){return(x[gene,])})))
  dat$organ=rownames(dat)
  #dat=escalc(measure="COR", ri=Coef, ni=ni, data=dat, slab=organ,vtype="HO") 
  dat=escalc(measure="MD", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat, slab=organ) 
  res=rma(yi, vi, data=dat,method="FE")
  resest=predict(res, digits=3)
  Results=rbind(Results,data.frame(Gene=gene,Coef=resest$pred,pvalue=res$pval))
}
Results$FDR <- p.adjust(Results$pvalue, method = "fdr")
Results <- Results[order(Results$pvalue),]
names(Results) <- c("Gene", "Coef", "Pvalue", "FDR")
write.table(Results, "table/wgii_allCnv_metattest_results.txt", quote=F, sep="\t", row.names = F)
