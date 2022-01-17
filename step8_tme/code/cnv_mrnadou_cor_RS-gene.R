rm(list=ls())
library(ggplot2)
library(ggsignif)
library(matrixStats)
library(ggrepel)
library(cowplot)
library(data.table)
setwd("/data/zhang/pancan_cin/step8_tme/")
table <- read.delim("../step2_genome_instability/data/TaylorCancerCell_TableS2.txt", na.strings = c("", " ", "na", "NA", "n.a.", "#N/A"))
rownames(table) <- as.character(table$Sample)
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer','Primary Solid Tumor'),]
rownames(Final_table)=Final_table$Sample
sels=intersect(rownames(table),rownames(Final_table))
Final_table=Final_table[sels,]
Final_table$Sample=substr(Final_table$Sample,1,15)
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
pi3k=c("GINS1", "GINS2", "GINS3", "GINS4","CKS1B", "CKS2","ATR", "CHK1", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "PCNA", "CDKN1A", "FEN1", "EXO1", 
       "MRE11", "BRCA2", "RAD50", "RECQL4", "WRN", "CENPJ", "CEP152", 
       "PCNT", "POLD2", "FANCD2", "FANCM", "FANCI", "FANCJ", "FANCN", 
       "FANCD1", "BLM", "RMI1", "RMI2", "TOPOIII", "PICH", "SMC1", "POLE", 
       "BRCA1", "MUS81", "SLX4", "WEE1", "SLX1", "EME1", "EME2", "RAD51", 
       "RAD52", "ATRIP", "RPA", "RIF1")
selg=intersect(rownames(cnv_tcga),rownames(mrna_tcga))
selg=intersect(pi3k,selg)
sels=Reduce(intersect,list(rownames(Final_table),colnames(mrna_tcga),colnames(cnv_tcga)))
cnv_tcga=cnv_tcga[selg,sels]
mrna_tcga=mrna_tcga[selg,sels]
Final_table=Final_table[sels,]
tmpc=unique(Final_table$Cohort)
reslist=list()
threshold_val=0.3
for(j in tmpc){
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
  restmp=data.frame(apply(restmp,2,abs))

    wgii=restmp[which((restmp$cor_cnv_wgii>threshold_val & restmp$cor_cnv_mrna>threshold_val)),]
    #wgii=rbind(wgii,restmp[which((restmp$cor_cnv_wgii< -threshold_val & restmp$cor_cnv_mrna< -threshold_val)),])
    
    ncin=restmp[which((restmp$cor_cnv_ncin>threshold_val & restmp$cor_cnv_mrna>threshold_val)),]
    #ncin=rbind(ncin,restmp[which((restmp$cor_cnv_ncin< -threshold_val & restmp$cor_cnv_mrna< -threshold_val)),])
    
    scin=restmp[which((restmp$cor_cnv_scin>threshold_val & restmp$cor_cnv_mrna>threshold_val)),]
    #scin=rbind(scin,restmp[which((restmp$cor_cnv_scin< -threshold_val & restmp$cor_cnv_mrna< -threshold_val)),])
    
    restmp$sig_wgii=rep('',dim(restmp)[1])
    restmp[which((restmp$cor_cnv_wgii>threshold_val & restmp$cor_cnv_mrna>threshold_val)),'sig_wgii']='yes'
    restmp$sig_ncin=rep('',dim(restmp)[1])
    restmp[which((restmp$cor_cnv_ncin>threshold_val & restmp$cor_cnv_mrna>threshold_val)),'sig_ncin']='yes'
    restmp$sig_scin=rep('',dim(restmp)[1])
    restmp[which((restmp$cor_cnv_scin>threshold_val & restmp$cor_cnv_mrna>threshold_val)),'sig_scin']='yes'
    reslist[[j]]=restmp
    rm(restmp,wgii,ncin,scin)
  
}
ACC=reslist[['ACC']]
BLCA=reslist[['BLCA']]
BRCA=reslist[['BRCA']]
CESC=reslist[['CESC']]
CHOL=reslist[['CHOL']]
COAD=reslist[['COAD']]
DLBC=reslist[['DLBC']]
ESCA=reslist[['ESCA']]
GBM=reslist[['GBM']]
HNSC=reslist[['HNSC']]
KICH=reslist[['KICH']]
KIRC=reslist[['KIRC']]
KIRP=reslist[['KIRP']]
LAML=reslist[['LAML']]
LGG=reslist[['LGG']]
LIHC=reslist[['LIHC']]
LUAD=reslist[['LUAD']]
LUSC=reslist[['LUSC']]
MESO=reslist[['MESO']]
OV=reslist[['OV']]
PAAD=reslist[['PAAD']]
PCPG=reslist[['PCPG']]
PRAD=reslist[['PRAD']]
READ=reslist[['READ']]
SARC=reslist[['SARC']]
STAD=reslist[['STAD']]
TGCT=reslist[['TGCT']]
THCA=reslist[['THCA']]
THYM=reslist[['THYM']]
UCEC=reslist[['UCEC']]
UCS=reslist[['UCS']]
UVM=reslist[['UVM']]
xls=names(reslist)
WriteXLS::WriteXLS(xls,ExcelFileName = "table/cnv_mrna_dou_cor_RSgene.xls",AdjWidth = TRUE, BoldHeaderRow = TRUE, row.names = T)
for(col_name in colnames(reslist[[1]])[1:4])
{
  heatdat=data.frame(do.call(cbind,lapply(reslist, function(x){x[,col_name]})))
  rownames(heatdat)=rownames(reslist[[1]])
  
  library(pheatmap)
  library(RColorBrewer)
  pheatmap::pheatmap(heatdat,scale = "none",
                     #clustering_method = "ward.D",
                     clustering_distance_cols = "euclidean",
                     clustering_distance_rows = "euclidean",
                     border_color = NA,main =col_name,cluster_cols = T,cluster_rows = F,
                     show_rownames =T,show_colnames = T,treeheight_row = 0,fontsize=3,fontsize_row =3,cellwidth =3,cellheight =3,treeheight_col = 0,
                     file =paste0(c("/data/zhang/pancan_cin/step5_pathway/plot/duo_correlation_",col_name,'.pdf'),collapse = ''),
                     width =5.5, height =5,color = rev(colorRampPalette(brewer.pal(11,"RdBu"))(2000)))
  
}
