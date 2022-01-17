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
clist=split(Final_expression_table_tum,Final_expression_table_tum$Cohort)
cname=c("GBM", "OV", "LUAD", "LUSC", "PRAD", "UCEC", "BLCA", "TGCT", 
        "ESCA", "PAAD", "KIRP", "LIHC", "CESC", "SARC", "BRCA", "THYM", 
        "MESO", "COAD", "STAD", "CHOL", "KIRC", "THCA", "HNSC", "LAML", 
        "READ", "SKCM", "LGG", "DLBC", "KICH", "UCS", "ACC", "PCPG", 
        "UVM")
dflist=list()
for(i in cname){
fn=paste0(c('table/scin_gene_exp_cor_',i,'.txt'),collapse = '')
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
dat=escalc(measure="COR", ri=Coef, ni=ni, data=dat, slab=organ) 
#dat=escalc(measure="MD", m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i, data=dat, slab=organ) 
res=rma(yi, vi, data=dat,method="FE")
resest=predict(res, digits=3)
Results=rbind(Results,data.frame(Gene=gene,Coef=resest$pred,pvalue=res$pval))
}
Results$FDR <- p.adjust(Results$pvalue, method = "fdr")
Results <- Results[order(Results$pvalue),]
names(Results) <- c("Gene", "Coef", "Pvalue", "FDR")
write.table(Results, "table/scin_allGeneExp_metacor_results.txt", quote=F, sep="\t", row.names = F)
