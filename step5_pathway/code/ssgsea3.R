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
library(readxl)
library(gplots)
library(RColorBrewer)
setwd("/data/zhang/pancan_cin/step5_pathway")
onco_template=readxl::excel_sheets('../step9_driver_interaction/data/NIHMS957693-supplement-7.xlsx')[1:10]
oncosignature=list()
for(i in onco_template[-1])
{
  oncosignature[[i]]=readxl::read_xlsx('../step9_driver_interaction/data/NIHMS957693-supplement-7.xlsx',sheet = i)
}

oncosignature$`Cell Cycle`=readxl::read_xlsx('../step9_driver_interaction/data/NIHMS957693-supplement-7.xlsx',sheet = onco_template[1],skip = 1)
oncosignature=lapply(oncosignature, function(x){as.character(x$Gene)})
names(oncosignature)=gsub(' ','_',names(oncosignature))

ddr_v1=readxl::read_xlsx('../step8_tme/data/TCGA_DDR_Data_Resources.xlsx',sheet='DDR genes and pathways',skip = 1)
colnames(ddr_v1)[11:20]=paste0('ddr_',colnames(ddr_v1)[11:20])
colnames(ddr_v1)[21:29]=paste0('ddr_core_',colnames(ddr_v1)[21:29])
colnames(ddr_v1)=unlist(lapply(colnames(ddr_v1), function(x){strsplit(x,'\\.\\.\\.')[[1]][1]}))
names(ddr_v1)=gsub(' ','_',names(ddr_v1))
names(ddr_v1)=gsub('\\(','',names(ddr_v1))
names(ddr_v1)=gsub('\\)','',names(ddr_v1))
pp=c(  "ddr_Base_Excision_Repair_BER", "ddr_Nucleotide_Excision_Repair_NER_-_includes_TC-NER_and_GC-NER", 
       "ddr_Mismatch_Repair_MMR", "ddr_Fanconi_Anemia_FA", "ddr_Homology-dependent_recombination_HDR", 
       "ddr_Non-homologous_End_Joining_NHEJ", "ddr_Direct_Repair_DR", 
       "ddr_Translesion_Synthesis_TLS", "ddr_Nucleotide_pools_NP", "ddr_Others", 
       "ddr_core_Base_Excision_Repair_BER", "ddr_core_Nucleotide_Excision_Repair_NER,_including_TC-NER_and_GC-NER", 
       "ddr_core_Mismatch_Repair_MMR", "ddr_core_Fanconi_Anemia_FA", 
       "ddr_core_Homologous_Recomination_HR", "ddr_core_Non-homologous_End_Joining_NHEJ", 
       "ddr_core_Direct_Repair_DR", "ddr_core_Translesion_Synthesis_TLS", 
       "ddr_core_Damage_Sensor_etc.")
ddr_v1=lapply(pp,function(x){ddr_v1[complete.cases(ddr_v1[,x]),x][[1]]})
names(ddr_v1)= pp

ddr_v2=readxl::read_xls('../step5_pathway/data/ddr_nrc3891.xls')
pp=c("DSR", "CS", "TLS", "MMR", "UR", "TM", "FA", "NHEJ", "CR", 
     "BER", "HR", "CF", "NER", "PDDR")
ddr_v2=lapply(pp,function(x){ddr_v2[which(ddr_v2$Abbreviation==x),'Gene.ID'][[1]]})
names(ddr_v2)=pp

biocart <- GSA.read.gmt("../step5_pathway/data/c2.cp.biocarta.v7.0.symbols.gmt")
names(biocart$genesets) <- biocart$geneset.names
biocart=biocart$genesets
react <- GSA.read.gmt("../step5_pathway/data/c2.cp.reactome.v7.0.symbols.gmt")
names(react$genesets) <- react$geneset.names
react=react$genesets
kegg <- GSA.read.gmt("../step5_pathway/data/c2.cp.kegg.v7.0.symbols.gmt")
names(kegg$genesets) <- kegg$geneset.names
kegg=kegg$genesets
cgp <- GSA.read.gmt("../step5_pathway/data/c2.cgp.v7.0.symbols.gmt")
names(cgp$genesets) <- cgp$geneset.names
cgp=cgp$genesets
hallmark <- GSA.read.gmt("../step5_pathway/data/h.all.v7.0.symbols.gmt")
names(hallmark$genesets) <- hallmark$geneset.names
hallmark=hallmark$genesets
go <- GSA.read.gmt("../step5_pathway/data/c5.go.v7.4.symbols.gmt")
names(go$genesets) <- go$geneset.names
go=go$genesets
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer','Primary Solid Tumor'),]
Final_table$Sample <- substr(Final_table$Sample, 1, 15)
rownames(Final_table)=Final_table$Sample

rna_table =fread("../step8_tme/data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv")
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
sels=intersect(rna_table$Sample,Final_table$Sample)
Final_table <- Final_table[sels,]
rna_table=rna_table[sels,]
cname=unique(Final_table$Cohort)

gsva_oncokb=list()
gsva_biocart=list()
gsva_react=list()
gsva_kegg=list()
gsva_cgp=list()
gsva_hallmark=list()
gsva_ddr_v1=list()
gsva_ddr_v2=list()
gsva_oncosignature=list()
gsva_go=list()
for(i in cname[21:22]){
  print(i)
  gene_exp=rna_table[Final_table$Sample[which(Final_table$Cohort==i)],]
  gene_exp=gene_exp[,which(apply(gene_exp, 2, function(x){length(unique(x))>2}))]
  gsva_go[[i]] <- gsva(t(gene_exp),go,method='ssgsea',kcdf='Gaussian', mx.diff = TRUE,parallel.sz=32)
  gsva_biocart[[i]] <- gsva(t(gene_exp),biocart,method='ssgsea',kcdf='Gaussian', mx.diff = TRUE,parallel.sz=32)
  gsva_react[[i]] <- gsva(t(gene_exp),react,method='ssgsea',kcdf='Gaussian', mx.diff = TRUE,parallel.sz=32)
  gsva_kegg[[i]] <- gsva(t(gene_exp),kegg,method='ssgsea',kcdf='Gaussian', mx.diff = TRUE,parallel.sz=32)
  gsva_cgp[[i]] <- gsva(t(gene_exp),cgp,method='ssgsea',kcdf='Gaussian', mx.diff = TRUE,parallel.sz=32)
  gsva_hallmark[[i]] <- gsva(t(gene_exp),hallmark,method='ssgsea',kcdf='Gaussian', mx.diff = TRUE,parallel.sz=32)
  gsva_oncokb[[i]] <- gsva(t(gene_exp),oncosignature,method='ssgsea',kcdf='Gaussian', mx.diff = TRUE,parallel.sz=32)
  gsva_ddr_v1[[i]] <- gsva(t(gene_exp),ddr_v1,method='ssgsea',kcdf='Gaussian', mx.diff = TRUE,parallel.sz=32)
  gsva_ddr_v2[[i]] <- gsva(t(gene_exp),ddr_v2,method='ssgsea',kcdf='Gaussian', mx.diff = TRUE,parallel.sz=32)
}
save(gsva_biocart,file ='table/gsva_biocart_pathway6.RData')
save(gsva_react,file ='table/gsva_react_pathway6.RData')
save(gsva_kegg,file ='table/gsva_kegg_pathway6.RData')
save(gsva_cgp,file ='table/gsva_cgp_pathway6.RData')
save(gsva_hallmark,file ='table/gsva_hallmark_pathway6.RData')
save(gsva_oncokb,file ='table/gsva_oncokb_pathway6.RData')
save(gsva_ddr_v1,file ='table/gsva_ddrv1_pathway6.RData')
save(gsva_ddr_v2,file ='table/gsva_ddrv2_pathway6.RData')
save(gsva_go,file ='table/gsva_go_pathway6.RData')
