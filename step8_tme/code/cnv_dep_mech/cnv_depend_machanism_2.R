rm(list=ls())
library(ggplot2)
library(ggsignif)
library(matrixStats)
library(ggrepel)
library(cowplot)
library(data.table)
setwd("/data/zhang/pancan_cin/step8_tme/")
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer','Primary Solid Tumor'),]
#Final_table=scores[scores$sample_type%in% c('Tumor'),]
Final_table$Sample <- substr(Final_table$Sample, 1, 15)
rownames(Final_table)=Final_table$Sample
rownames(Final_table)=gsub('-','.',rownames(Final_table))
mrna_tcga = fread('data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv')
colnames(mrna_tcga) = substr(colnames(mrna_tcga), 1, 15)
mrna_tcga$gene_id = unlist(lapply(as.character(mrna_tcga$gene_id), function(x) {
  strsplit(x, '\\|')[[1]][1]
}))
mrna_tcga = mrna_tcga[!duplicated(mrna_tcga$gene_id),]
mrna_tcga_tmp = data.frame(mrna_tcga[,-1])
rownames(mrna_tcga_tmp) = mrna_tcga$gene_id
mrna_tcga = mrna_tcga_tmp

cnv_tcga = fread('data/all_data_by_genes_whitelisted.tsv')
colnames(cnv_tcga) = substr(colnames(cnv_tcga), 1, 15)
cnv_tcga = cnv_tcga[!duplicated(cnv_tcga$`Gene Symbol`),]
cnv_tcga_tmp = data.frame(cnv_tcga[,-c(1:3)])
rownames(cnv_tcga_tmp) = cnv_tcga$`Gene Symbol`
cnv_tcga = cnv_tcga_tmp
selg=intersect(rownames(cnv_tcga),rownames(mrna_tcga))
sels=Reduce(intersect,list(rownames(Final_table),colnames(mrna_tcga),colnames(cnv_tcga)))
cnv_tcga=cnv_tcga[selg,sels]
mrna_tcga=mrna_tcga[selg,sels]
Final_table=Final_table[sels,]
tmpc=unique(Final_table$Cohort)
reslist=list()
for(j in tmpc[4:6]){
  Results <- data.frame()
  tmpmat=subset(Final_table,Final_table$Cohort==j)
  tmpcnv=cnv_tcga[,rownames(tmpmat)]
  tmpmrna=mrna_tcga[,rownames(tmpmat)]
  res_mol=unlist(lapply(1:nrow(tmpcnv), function(x){cor(as.numeric(tmpcnv[x,]),as.numeric(tmpmrna[x,]),method = 'spearman')}))
  res_wgii=unlist(lapply(1:nrow(tmpcnv), function(x){cor(as.numeric(tmpcnv[x,]),as.numeric(tmpmat$wgii),method = 'spearman')}))
  res_ncin=unlist(lapply(1:nrow(tmpcnv), function(x){cor(as.numeric(tmpcnv[x,]),as.numeric(tmpmat$ncin),method = 'spearman')}))
  res_scin=unlist(lapply(1:nrow(tmpcnv), function(x){cor(as.numeric(tmpcnv[x,]),as.numeric(tmpmat$scin),method = 'spearman')}))
  restmp=data.frame(cor_cnv_mrna=res_mol,cor_cnv_wgii=res_wgii,cor_cnv_ncin=res_ncin,cor_cnv_scin=res_scin)
  rownames(restmp)=selg
  save(restmp,file = file.path('table/cnv_dep_mech',paste0(c('cor_cv_rna_cin_',j,'.RData'),collapse = '')))
}