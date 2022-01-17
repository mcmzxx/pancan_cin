rm(list = ls())
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
library(BiocStyle)
library(knitr)
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)
library(gdata)
library(snow)
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
selg=intersect(rownames(cnv_tcga),rownames(mrna_tcga))
sels=Reduce(intersect,list(rownames(Final_table),colnames(mrna_tcga),colnames(cnv_tcga)))
cnv_tcga=cnv_tcga[selg,sels]
mrna_tcga=mrna_tcga[selg,sels]
Final_table=Final_table[sels,]
tmpc=unique(Final_table$Cohort)
setwd("/data/zhang/pancan_cin/step8_tme/")
flists=list.files('data',pattern = 'cor_cv_rna_vin',full.names = T)
barlist=list()
threshold_val=0.4
for(fl in flists){
  #fl=flists[6]
  tmpname=gsub('.RData','',gsub('cor_cv_rna_vin_','',strsplit(fl,'\\/')[[1]][2]))
  load(fl)
  restmp=data.frame(apply(restmp,2,abs))
  rownames(restmp)=rownames(cnv_tcga)
  wgii=restmp[which((restmp$cor_cnv_wgii>threshold_val & restmp$cor_cnv_mrna>threshold_val)),]
  #wgii=rbind(wgii,restmp[which((restmp$cor_cnv_wgii< -threshold_val & restmp$cor_cnv_mrna< -threshold_val)),])
  
  ncin=restmp[which((restmp$cor_cnv_ncin>threshold_val & restmp$cor_cnv_mrna>threshold_val)),]
  #ncin=rbind(ncin,restmp[which((restmp$cor_cnv_ncin< -threshold_val & restmp$cor_cnv_mrna< -threshold_val)),])
  
  scin=restmp[which((restmp$cor_cnv_scin>threshold_val & restmp$cor_cnv_mrna>threshold_val)),]
  #scin=rbind(scin,restmp[which((restmp$cor_cnv_scin< -threshold_val & restmp$cor_cnv_mrna< -threshold_val)),])
  
  barlist[[tmpname]]=list(wgii=rownames(wgii),ncin=rownames(ncin),scin=rownames(scin))
  rm(restmp,wgii,ncin,scin)
}
hallmark=GSA.read.gmt('/data/zhang/pancan_cin/step5_pathway/data/h.all.v7.0.symbols.gmt')
names(hallmark$genesets)=hallmark$geneset.names
hallmark=hallmark$genesets
get.enrichment.list = function(KEGG.PathwayGeneSets2, Gene.Modules, genes.entrezId){
  table = list()
  for(kk in 1:length(Gene.Modules)){
    module.gene.set = Gene.Modules[[kk]]
    pvalues = NULL
    for(i in 1:length(KEGG.PathwayGeneSets2)){
      term.gene.set = KEGG.PathwayGeneSets2[[i]]
      k = sum(module.gene.set%in%term.gene.set)
      M = length(term.gene.set)
      N = length(genes.entrezId)
      n = length(module.gene.set)
      pvalues[i] = phyper(k-1,M, N-M, n, lower.tail=FALSE)
      # e.g. phyper(29,165,423-165,46,lower.tail=FALSE)
    }
    
    pvalues.adjust = p.adjust(pvalues, "BH")
    pvalues = pvalues
    
    enrichment.score = data.frame(pvalues = pvalues, pvalues.adjust = pvalues.adjust)
    row.names(enrichment.score) = names(KEGG.PathwayGeneSets2)
    enrichment.score=enrichment.score[order(enrichment.score$pvalues),]
    table[[kk]] = enrichment.score
  }
  names(table)=c('wgii','ncin','scin')
  return(table)
}
corlist=lapply(barlist, function(x){get.enrichment.list(hallmark,x, rownames(cnv_tcga))})
ACC=corlist[['ACC']][['ncin']]
BLCA=corlist[['BLCA']][['ncin']]
BRCA=corlist[['BRCA']][['ncin']]
CESC=corlist[['CESC']][['ncin']]
CHOL=corlist[['CHOL']][['ncin']]
COAD=corlist[['COAD']][['ncin']]
DLBC=corlist[['DLBC']][['ncin']]
ESCA=corlist[['ESCA']][['ncin']]
GBM=corlist[['GBM']][['ncin']]
HNSC=corlist[['HNSC']][['ncin']]
KICH=corlist[['KICH']][['ncin']]
KIRC=corlist[['KIRC']][['ncin']]
KIRP=corlist[['KIRP']][['ncin']]
LAML=corlist[['LAML']][['ncin']]
LGG=corlist[['LGG']][['ncin']]
LIHC=corlist[['LIHC']][['ncin']]
LUAD=corlist[['LUAD']][['ncin']]
LUSC=corlist[['LUSC']][['ncin']]
MESO=corlist[['MESO']][['ncin']]
OV=corlist[['OV']][['ncin']]
PAAD=corlist[['PAAD']][['ncin']]
PCPG=corlist[['PCPG']][['ncin']]
PRAD=corlist[['PRAD']][['ncin']]
READ=corlist[['READ']][['ncin']]
SARC=corlist[['SARC']][['ncin']]
STAD=corlist[['STAD']][['ncin']]
TGCT=corlist[['TGCT']][['ncin']]
THCA=corlist[['THCA']][['ncin']]
THYM=corlist[['THYM']][['ncin']]
UCEC=corlist[['UCEC']][['ncin']]
UCS=corlist[['UCS']][['ncin']]
UVM=corlist[['UVM']][['ncin']]
xls=names(corlist)
WriteXLS::WriteXLS(xls,"table/cnv_denpend_mechanism_ncin.xls",AdjWidth = TRUE, BoldHeaderRow = TRUE, row.names = T)




