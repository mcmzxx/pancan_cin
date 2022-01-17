rm(list=ls())
library("tidyverse") # loads ggplot2, tibble, tidyr, readr, purr, dplyr
library("genefilter")
library("edgeR")
library("ggrepel")
library("gplots")
library("reshape2")
library("data.table")
library("magrittr")
library("GSVA")
library("GSA")
library("limma")
library("edgeR")
library("rio")
library(gdata)
library(gplots)
library(RColorBrewer)

setwd("/data/zhang/pancan_cin/step9_driver_interaction")
table <-read.delim("../step2_genome_instability/data/TaylorCancerCell_TableS2.txt",na.strings = c("", " ", "na", "NA", "n.a.", "#N/A"))
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
score_table = scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer', 'Primary Solid Tumor'), ]
Final_table <- merge(table, score_table, by = 'Sample')
Final_table$Sample=substr(Final_table$Sample,1,15)
rownames(Final_table) = Final_table$Sample

rna_idx=c("GINS1", "GINS2", "GINS3", "GINS4","CKS1B", "CKS2")
pi3k=c("ATR", "CHK1", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "PCNA", "CDKN1A", "FEN1", "EXO1", 
       "MRE11", "BRCA2", "RAD50", "RECQL4", "WRN", "CENPJ", "CEP152", 
       "PCNT", "POLD2", "FANCD2", "FANCM", "FANCI", "FANCJ", "FANCN", 
       "FANCD1", "BLM", "RMI1", "RMI2", "TOPOIII", "PICH", "SMC1", "POLE", 
       "BRCA1", "MUS81", "SLX4", "WEE1", "SLX1", "EME1", "EME2", "RAD51", 
       "RAD52", "ATRIP", "RPA", "RIF1")
pro_idx=c("PI3KP110ALPHA","PI3KP85","PTEN")
cnv_idx=c('PIK3CA','PTEN')
cnv_tcga = fread('../step8_tme/data/all_data_by_genes_whitelisted.tsv')
colnames(cnv_tcga) = substr(colnames(cnv_tcga), 1, 15)
cnv_tcga = cnv_tcga[!duplicated(cnv_tcga$`Gene Symbol`),]
cnv_tcga_tmp = cnv_tcga[,-c(1:3)]
rownames(cnv_tcga_tmp) = cnv_tcga$`Gene Symbol`
cnv_tcga_tmp = data.frame(t(cnv_tcga_tmp))
colnames(cnv_tcga_tmp)=cnv_tcga$`Gene Symbol`
cnv_tcga=cnv_tcga_tmp
cnv_tcga$Sample=rownames(cnv_tcga)
rna_table <-fread("../step8_tme/data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv")
ii = unlist(lapply(rna_table$gene_id, function(x) {
  strsplit(as.character(x), '\\|')[[1]][1]
}))
iidx = which(!duplicated(ii))
ii = ii[iidx]
rna_table = rna_table[iidx, ]
rownames(rna_table) = ii
rna_table = rna_table[which(rownames(rna_table) != "?"), ]
rna_table$gene_id = NULL
rna_table = t(as.data.frame(rna_table))
colnames(rna_table) = ii[-1]
rna_table = data.frame(rna_table)
rna_table$Sample =  substr(rownames(rna_table), 1, 15)
rna_table=rna_table[!duplicated(rna_table$Sample),]
rownames(rna_table)=rna_table$Sample

onkotablebak=fread('/data/zhang/pancan_cin/step9_driver_interaction/data/onkokb.csv')
onkotable=data.frame(apply(onkotablebak[,-1], 2, function(x){ifelse(x==FALSE,0,1)}),stringsAsFactors = F)
onkotable[is.na(onkotable)]=0
rownames(onkotable)=onkotablebak$SAMPLE_BARCODE
Mut_table <- fread("../step8_tme/data/mc3.v0.2.8.PUBLIC.maf.gz")
Mut_table$Sample=substr(Mut_table$Tumor_Sample_Barcode,1,15)
protein_table <-fread("../step8_tme/data/TCGA-RPPA-pancan-clean.txt")
protein_table=data.frame(protein_table)
rownames(protein_table)=substr(protein_table$SampleID,1,15)
cnv_table <- fread("data/data_CNA.oncokb.txt")
cnv_table=subset(cnv_table,cnv_table$ALTERATION=='Deletion',)
cnv_table=cnv_table[which(cnv_table$HUGO_SYMBOL %in% cnv_idx),]
cnv_table$Sample=substr(cnv_table$SAMPLE_ID,1,15)

Final_expression_table_tum=merge(Final_table,cnv_tcga[,c('Sample',cnv_idx)],by='Sample')
colnames(Final_expression_table_tum)[which(colnames(Final_expression_table_tum) %in% cnv_idx)]=paste(cnv_idx,'_cnv',sep='')
Final_expression_table_tum=merge(Final_expression_table_tum,rna_table[,c('Sample',rna_idx)],by='Sample')
colnames(Final_expression_table_tum)[which(colnames(Final_expression_table_tum) %in% rna_idx)]=paste(rna_idx,'_gene',sep='')
Final_expression_table_tum$TP53_mut <- "wt"
Final_expression_table_tum$TP53_mut[Final_expression_table_tum$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% "TP53"]] <- "mut"
Final_expression_table_tum$TP53_mut<- factor(Final_expression_table_tum$TP53, levels = c("wt", "mut"))
study_var=c("PIK3CA_cnv", "PTEN_cnv", "GINS1_gene", 
           # "GINS2_gene", "GINS3_gene", "GINS4_gene",
           "CKS1B_gene")
           #, "CKS2_gene")
Results_TP53 <- data.frame()

for(i in 1:length(unique(Final_expression_table_tum$Cohort))){
  
  cohort <- paste(unique(Final_expression_table_tum$Cohort)[i])
  
  if(length(which(Final_expression_table_tum$TP53 %in% "mut" & Final_expression_table_tum$Cohort %in% cohort)) >=20 & length(which(Final_expression_table_tum$TP53 %in% "mut" & Final_expression_table_tum$Cohort %in% cohort)) >=20){
    
    Final_expression_table_tum_cohort <- Final_expression_table_tum[Final_expression_table_tum$Cohort %in% cohort,]
    
   for(j in study_var){
     fit <- lm(as.formula(paste0(j," ~ TP53_mut")), Final_expression_table_tum_cohort)
     
     pvalue = summary(fit)$coefficients[2,4]  
     coef = summary(fit)$coefficients[2,1]
     
     wt <- mean(Final_expression_table_tum_cohort$ncin[Final_expression_table_tum_cohort$TP53 %in% "wt"])
     mut <- mean(Final_expression_table_tum_cohort$ncin[Final_expression_table_tum_cohort$TP53 %in% "mut"])
     
     WT_samples <- length(which(Final_expression_table_tum$TP53 %in% "wt" & Final_expression_table_tum$Cohort %in% cohort))
     MUT_samples <- length(which(Final_expression_table_tum$TP53 %in% "mut" & Final_expression_table_tum$Cohort %in% cohort))
     
     Results_TP53=rbind(Results_TP53, data.frame(cohort, coef, pvalue, wt, mut, mut-wt, WT_samples, MUT_samples,j))
   }
    
  }
}

Results_TP53$FDR <- p.adjust(Results_TP53$pvalue, method = "fdr")
Results_TP53 <- Results_TP53[,c(1:3,10,9,4:8)]
names(Results_TP53) <- c("Cohort", "Coef", "Pvalue", "FDR","Gene", "WTmean", "MUTmean", "MUTminusWT_diff", "WT_samples", "MUT_samples")


### plot for all cohorts

Results_TP53 <- Results_TP53[order(Results_TP53$Coef, decreasing = T),]

order_names=as.character(Results_TP53$Cohort)
Results_TP53$Cohort <- factor(Results_TP53$Cohort, levels=order_names)

Results_TP53$Significant <- ifelse(Results_TP53$FDR < 0.05 & Results_TP53$Coef > 0, "Pos", "Not Sig")
Results_TP53$Significant[Results_TP53$FDR < 0.05 & Results_TP53$Coef < 0] <- "Neg"
Results_TP53$Significant=factor(Results_TP53$Significant,levels=c('Not Sig','Pos'))