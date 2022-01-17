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
library(gdata)
library(snow)
library(cetcolor)
library(WriteXLS)
setwd("/data/zhang/pancan_cin/step8_tme/")
flists=list.files('table/cnv_dep_mech',full.names = T)
barlist=list()
threshold_val=0.4
for(fl in flists){
  #fl=flists[6]
  tmpname=gsub('.RData','',gsub('cor_cv_rna_cin_','',strsplit(fl,'\\/')[[1]][3]))
  load(fl)
  restmp=data.frame(apply(restmp,2,abs))
  wgii=restmp[which((restmp$cor_cnv_wgii>threshold_val & restmp$cor_cnv_mrna>threshold_val)),]
  #wgii=rbind(wgii,restmp[which((restmp$cor_cnv_wgii< -threshold_val & restmp$cor_cnv_mrna< -threshold_val)),])
  
  scin=restmp[which((restmp$cor_cnv_scin>threshold_val & restmp$cor_cnv_mrna>threshold_val)),]
  #scin=rbind(scin,restmp[which((restmp$cor_cnv_scin< -threshold_val & restmp$cor_cnv_mrna< -threshold_val)),])
  
  barlist[[tmpname]]=list(wgii=rownames(wgii),scin=rownames(scin))
  rm(restmp,wgii,scin)
}
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
cpg=cgp$genesets
hallmark <- GSA.read.gmt("../step5_pathway/data/h.all.v7.0.symbols.gmt")
names(hallmark$genesets) <- hallmark$geneset.names
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
    enrichment.score$pathway = names(KEGG.PathwayGeneSets2)
    enrichment.score=enrichment.score[order(enrichment.score$pvalues),]
    table[[kk]] = enrichment.score
  }
  names(table)=c('wgii','scin')
  return(table)
}
load(fl)
geneset=c(oncosignature,ddr_v1,ddr_v2,hallmark)
corlist=lapply(barlist, function(x){get.enrichment.list(geneset,x,rownames(restmp))})
# ACC=corlist[['ACC']][['scin']]
# BLCA=corlist[['BLCA']][['scin']]
# BRCA=corlist[['BRCA']][['scin']]
# CESC=corlist[['CESC']][['scin']]
# CHOL=corlist[['CHOL']][['scin']]
# COAD=corlist[['COAD']][['scin']]
# DLBC=corlist[['DLBC']][['scin']]
# ESCA=corlist[['ESCA']][['scin']]
# GBM=corlist[['GBM']][['scin']]
# HNSC=corlist[['HNSC']][['scin']]
# KICH=corlist[['KICH']][['scin']]
# KIRC=corlist[['KIRC']][['scin']]
# KIRP=corlist[['KIRP']][['scin']]
# LAML=corlist[['LAML']][['scin']]
# LGG=corlist[['LGG']][['scin']]
# LIHC=corlist[['LIHC']][['scin']]
# LUAD=corlist[['LUAD']][['scin']]
# LUSC=corlist[['LUSC']][['scin']]
# MESO=corlist[['MESO']][['scin']]
# OV=corlist[['OV']][['scin']]
# PAAD=corlist[['PAAD']][['scin']]
# PCPG=corlist[['PCPG']][['scin']]
# PRAD=corlist[['PRAD']][['scin']]
# READ=corlist[['READ']][['scin']]
# SARC=corlist[['SARC']][['scin']]
# STAD=corlist[['STAD']][['scin']]
# TGCT=corlist[['TGCT']][['scin']]
# THCA=corlist[['THCA']][['scin']]
# THYM=corlist[['THYM']][['scin']]
# UCEC=corlist[['UCEC']][['scin']]
# UCS=corlist[['UCS']][['scin']]
# UVM=corlist[['UVM']][['scin']]
 xls=names(corlist)
# WriteXLS::WriteXLS(xls,"table/cnv_denpend_mechanism_scin.xls",AdjWidth = TRUE, BoldHeaderRow = TRUE, row.names = T)
Results=do.call(rbind,lapply(xls, function(x){
  tmp=corlist[[x]][['scin']];
  tmp$Cohort=rep(x,dim(tmp)[1]);
  return(tmp)
  }))
Results=data.frame(Results)
length(which(duplicated(with(Results,paste(pathway,Cohort)))))
Results <- dcast(Results, pathway~ Cohort, value.var="pvalues")
rownames(Results)=Results$pathway
Results$pathway=NULL
Results=Results[order(apply(Results, 1, function(x){length(which(x<=0.05))}),decreasing = T)[1:25],]
Results=Results[hclust(dist(Results))$order,hclust(dist(t(Results)))$order]
ii=rownames(Results)
uu=melt(Results)
colnames(uu)=c('Cohort','pvalue')
uu$pathway=rep(ii,32)
pheatmap::pheatmap(Results)

#uu$pvalue=ifelse(uu$pvalue>0.1,0.05,1)
pdf('plot/cn_dep_mech_scin.pdf',width = 16)
ggplot(uu, aes(Cohort, pathway)) + 
  geom_tile(aes(fill = pvalue)) + 
  geom_text(aes(fill = uu$pvalue, label = '*',color=ifelse(uu$pvalue<=0.05,'Sig','Not-sig'))) +scale_color_manual(name='Significance',values=c('Sig'='black','Not-sig'=NA))+
  scale_fill_gradientn(colours  = cet_pal(200, name = "d9")) +
  theme(panel.grid.major.x=element_blank(), 
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"),
        axis.text.x = element_text(angle=90, hjust = 0.2,vjust=0.4,size = 10),
        axis.text.y = element_text(size = 10)) + 
  ggtitle("CN dependent mechanisms driving SCIN") + 
  theme(legend.title=element_text(size=10)) + 
  scale_x_discrete(name="") +
  scale_y_discrete(name="") 
dev.off()

Results=do.call(rbind,lapply(xls, function(x){
  tmp=corlist[[x]][['wgii']];
  tmp$Cohort=rep(x,dim(tmp)[1]);
  return(tmp)
}))
Results=data.frame(Results)
length(which(duplicated(with(Results,paste(pathway,Cohort)))))
Results <- dcast(Results, pathway~ Cohort, value.var="pvalues")
rownames(Results)=Results$pathway
Results$pathway=NULL
Results=Results[order(apply(Results, 1, function(x){length(which(x<=0.05))}),decreasing = T)[1:25],]
Results=Results[hclust(dist(Results))$order,hclust(dist(t(Results)))$order]
ii=rownames(Results)
uu=melt(Results)
colnames(uu)=c('Cohort','pvalue')
uu$pathway=rep(ii,32)
pheatmap::pheatmap(Results)

pdf('plot/cn_dep_mech_wgii.pdf',width = 16)
ggplot(uu, aes(Cohort, pathway)) + 
  geom_tile(aes(fill = pvalue)) + 
  geom_text(aes(fill = uu$pvalue, label = '*',color=ifelse(uu$pvalue<=0.05,'Sig','Not-sig'))) +scale_color_manual(name='Significance',values=c('Sig'='black','Not-sig'=NA))+
  scale_fill_gradientn(colours  = cet_pal(200, name = "d9")) +
  theme(panel.grid.major.x=element_blank(), 
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"),
        axis.text.x = element_text(angle=90, hjust = 0.2,vjust=0.4,size = 10),
        axis.text.y = element_text(size = 10)) + 
  ggtitle("CN dependent mechanisms driving NCIN") + 
  theme(legend.title=element_text(size=10)) + 
  scale_x_discrete(name="") +
  scale_y_discrete(name="") 
dev.off()

