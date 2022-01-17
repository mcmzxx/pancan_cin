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
tablex <-readxl::read_excel("../step2_genome_instability/data/1-s2.0-S1535610818301119-mmc2.xlsx",skip = 1) %>% as.data.frame()
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
score_table = scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer', 'Primary Solid Tumor'), ]
Final_table <- merge(tablex, score_table, by = 'Sample')
Final_table$Sample=substr(Final_table$Sample,1,15)
rownames(Final_table) = Final_table$Sample

pi3k=c("GINS1", "GINS2", "GINS3", "GINS4","CKS1B", "CKS2")
pi3k=c("ATR", "CHK1", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "PCNA", "CDKN1A", "FEN1", "EXO1", 
       "MRE11", "BRCA2", "RAD50", "RECQL4", "WRN", "CENPJ", "CEP152", 
       "PCNT", "POLD2", "FANCD2", "FANCM", "FANCI", "FANCJ", "FANCN", 
       "FANCD1", "BLM", "RMI1", "RMI2", "TOPOIII", "PICH", "SMC1", "POLE", 
       "BRCA1", "MUS81", "SLX4", "WEE1", "SLX1", "EME1", "EME2", "RAD51", 
       "RAD52", "ATRIP", "RPA", "RIF1")

cnv_tcga = fread('../step8_tme/data/all_data_by_genes_whitelisted.tsv')
colnames(cnv_tcga) = substr(colnames(cnv_tcga), 1, 15)
cnv_tcga = cnv_tcga[!duplicated(cnv_tcga$`Gene Symbol`),]
cnv_tcga_tmp = cnv_tcga[,-c(1:3)]
rownames(cnv_tcga_tmp) = cnv_tcga$`Gene Symbol`
cnv_tcga = t(cnv_tcga_tmp)
cnv_tcga$Sample=rownames(cnv_tcga)
cnv_idx=colnames(cnv_tcga)[which(colnames(cnv_tcga) %in%  pi3k)]
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
rna_idx=colnames(rna_table)[which(colnames(rna_table) %in%  pi3k)]
onkotablebak=fread('/data/zhang/pancan_cin/step9_driver_interaction/data/onkokb.csv')
onkotable=data.frame(apply(onkotablebak[,-1], 2, function(x){ifelse(x==FALSE,0,1)}),stringsAsFactors = F)
onkotable[is.na(onkotable)]=0
rownames(onkotable)=onkotablebak$SAMPLE_BARCODE

protein_table <-fread("../step8_tme/data/TCGA-RPPA-pancan-clean.txt")
protein_table=data.frame(protein_table)
rownames(protein_table)=substr(protein_table$SampleID,1,15)
pro_idx=colnames(protein_table)[which(unlist(lapply(colnames(protein_table), function(x){strsplit(x,'_')[[1]][1]})) %in% c(pi3k,"PI3KP110ALPHA","PI3KP85","AKT"))]
pro_idx=c("PI3KP110ALPHA","PI3KP85","PTEN")
rna_idx=c('PIK3CA','PTEN')
Mut_table=onkotable[,onco_idx]
Mut_table$Sample=rownames(Mut_table)
Final_expression_table_tum=merge(Final_table, Mut_table, by='Sample')
Results <- data.frame()
for(i in colnames(Final_expression_table_tum)[109:dim(Final_expression_table_tum)[2]]){
  gmut=which(Final_expression_table_tum[,i]==1)
  gwt=which(Final_expression_table_tum[,i]==0)
  if(length(gmut) >=20 & length(gwt) >=20){
    
    Final_expression_table_tum$Mutation_selected <- factor(Final_expression_table_tum[,i],labels = c('wt','mut'))
    
    fit <- lm(scin ~ Mutation_selected + Cohort, Final_expression_table_tum)
    tmpcoef=summary(fit)$coefficients
    pvalue = tmpcoef[2,4]  
    # delta mean
    coef = tmpcoef[2,1]
    
    wt <- mean(Final_expression_table_tum$ncin[gwt])
    mut <- mean(Final_expression_table_tum$ncin[gmut])
    
    WT_samples <- length(gwt)
    MUT_samples <- length(gmut)
    Results=rbind(Results, data.frame(i, coef, pvalue, wt, mut, mut-wt, WT_samples, MUT_samples))
  }
}
Results$FDR <- p.adjust(Results$pvalue, method = "fdr")
Results <- Results[order(Results$pvalue),c(1:3,9,4:8)]
names(Results) <- c("Gene", "Coef", "Pvalue", "FDR", "WTmean", "MUTmean", "MUTminusWT_diff", "WT_samples", "MUT_samples")
Results=Results[order(Results$Pvalue),]
write.table(Results, "table/allsamples_ncin_onkokb_pathway_linear_model_results.txt", quote=F, sep="\t", row.names = F)



ii=colnames(Mut_table)[-grep('Sample',colnames(Mut_table))]
z <-merge(Final_table,Mut_table,by='Sample')
rownames(z)=z$Sample
z1=z[,c('wgii','scin')]
# z1$wgii=sapply(z1$wgii,function(x){ifelse(x>median(z1$wgii),1,0)})
# z1$scin=sapply(z1$scin,function(x){ifelse(x>median(z1$scin),1,0)})
z=z[,ii]
#z=z[,which(colSums(z)>70)]
z=z[which(rowSums(z)>0),]
z1=z1[rownames(z),]
#z=cbind(z,z1)
ii=colnames(z)
col_order = order(colSums(z), decreasing = TRUE)
scoreCol = function(x) {
  score = 0
  for(i in 1:length(x)) {
    if(x[i]) {
      score = score + 2^(length(x)-i*1/x[i])
    }
  }
  return(score)
}
scores = apply(z[,col_order], 1, scoreCol)
row_order=order(scores, decreasing=TRUE)

colMutations = colorRampPalette(yarrr::piratepal(palette = "pony",trans = 0)[-c(6:7)])(length(ii))
names(colMutations) <- ii
par(bty="n", mgp = c(2,.33,0), mar=rep(3,4), las=1, tcl=-.25, xpd=NA)
plot(NA,NA, xlim=c(0,ncol(z)+2), ylim=c(0,nrow(z)), xaxt="n", yaxt="n", xlab="",ylab="", xaxs="i", yaxs="i")
xtick<-seq(0.5, length(ii)+2+length(c(pro_idx,rna_idx))-0.5, by=1)
axis(side=2, at=xtick, labels = FALSE)
score_mat=z1[row_order,]
score_mat$scin[which(score_mat$scin>40)]=40
score_mat=apply(score_mat,2,function(x){rev(brewer.pal(11,"PiYG"))[cut(x, quantile(x+runif(length(x), -0.01,0.01), seq(0,1,l=11), na.rm=TRUE), right=F)]})
zmat=sapply(1:ncol(z), function(i) ifelse(z[,i]>0, colMutations[i], "#FFFFFF"))[row_order,col_order]

pdf('plot/onco_pi3k.pdf',height = 10,width = 18)
par(bty="n", mgp = c(2,.33,0), mar=c(0,6,0,0), las=1, tcl=-.25, xpd=NA)
imagemat=t(cbind(score_mat,zmat))
rownames(imagemat)=c('wgii','scin',colnames(z)[col_order])
colnames(imagemat)=rownames(z)[row_order]
tmpmat=rna_table[colnames(imagemat),rna_idx]
tmpmat=apply(tmpmat,2,function(x){rev(brewer.pal(11,"RdBu"))[cut(x, quantile(x+runif(length(x), -0.01,0.01), seq(0,1,l=11), na.rm=TRUE), right=F)]})
imagemat=rbind(t(tmpmat),imagemat)
tmpmat=protein_table[colnames(imagemat),pro_idx]
tmpmat=apply(tmpmat,2,function(x){rev(brewer.pal(11,"PuOr"))[cut(x, quantile(x+runif(length(x), -0.01,0.01), seq(0,1,l=11), na.rm=TRUE), right=F)]})
imagemat=rbind(t(tmpmat),imagemat)
plot(NA,NA, xlim=c(0,ncol(imagemat)), ylim=c(0,nrow(imagemat)), xaxt="n", yaxt="n", xlab="",ylab="", xaxs="i", yaxs="i")
rasterImage(imagemat, 0, 0, ncol(imagemat), nrow(imagemat), interpolate=FALSE)
text(y=xtick,  par("usr")[4], labels = rev(c(paste0(pro_idx,'_pro'),paste0(rna_idx,'_gene'),'wgii','scin',colnames(z)[col_order])), srt = 0, pos = 2, xpd = T,offset = 1.5,cex = 0.5)
dev.off()