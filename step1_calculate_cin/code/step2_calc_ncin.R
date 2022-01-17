rm(list=ls())
library(data.table)
library(plyr)
setwd("/data/zhang/pancan_cin/step1_calculate_cin/")
numComplexity <- function( sample_name, seg){
  require(Hmisc)
  seg <- seg[which(!seg$Chromosome %in% 23),]
  tmpn=unique(seg$Chromosome)
  vect.whole.chrom <- rep(0, length( unique(seg$Chromosome) ))
  names(vect.whole.chrom) <- tmpn
  seg.sample <- subset(seg, seg$Sample == sample_name)
  #seg.sample$ploidy=rep(3,dim(seg.sample)[1])
  seg.sample$Copy_Number <- seg.sample$Copy_Number -as.integer(round(seg.sample$ploidy))
  for (chr in tmpn){
    seg.chr <- subset(seg.sample, seg.sample$Chromosome == chr)
    idx=c(which(seg.chr$Copy_Number>0),which(seg.chr$Copy_Number<0))
    if(sum(seg.chr$End[idx]-seg.chr$Start[idx])/sum(seg.chr$End-seg.chr$Start)>=0.75)
      vect.whole.chrom[chr] <- 1
  }
  return(vect.whole.chrom)
}
cn.segs=fread("/data/zhang/pancan_cin/step8_tme/data/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg")
cn.segs$Sample=substr(cn.segs$Sample,1,15)
ploidy_abs=fread('/data/zhang/pancan_cin/step8_tme/data/TCGA_mastercalls.abs_tables_JSedit.fixed.txt')
ploidy_abs$Sample=ploidy_abs$array
clinic=gdata::read.xls('/data/zhang/pancan_cin/step8_tme/data/TCGA-CDR-SupplementalTableS1.xlsx')
clinic=data.frame(apply(clinic, 2, as.character))
cn.segs=plyr::join (cn.segs,ploidy_abs, by='Sample',type='left')
cn.segs$ploidy[which(is.na(cn.segs$ploidy))]=2
cn.segs$sample_type <- NA
cn.segs$sample_type[grep("-01$", cn.segs$Sample)] <- "Tumor"
cn.segs$sample_type[grep("-02$", cn.segs$Sample)] <- "Tumor"
cn.segs$sample_type[grep("-03$", cn.segs$Sample)] <- "Tumor"
cn.segs$sample_type[grep("-04$", cn.segs$Sample)] <- "Tumor"
cn.segs$sample_type[grep("-05$", cn.segs$Sample)] <- "Tumor"
cn.segs$sample_type[grep("-06$", cn.segs$Sample)] <- "Tumor"
cn.segs$sample_type[grep("-07$", cn.segs$Sample)] <- "Tumor"
cn.segs$sample_type[grep("-08$", cn.segs$Sample)] <- "Tumor"
cn.segs$sample_type[grep("-09$", cn.segs$Sample)] <- "Tumor"
cn.segs$sample_type[grep("-10$", cn.segs$Sample)] <- "Normal"
cn.segs$sample_type[grep("-11$", cn.segs$Sample)] <- "Normal"
cn.segs$sample_type[grep("-12$", cn.segs$Sample)] <- "Normal"
cn.segs$sample_type[grep("-13$", cn.segs$Sample)] <- "Normal"
cn.segs$sample_type[grep("-14$", cn.segs$Sample)] <- "Normal"
cn.segs$sample_type_detail <- NA
cn.segs$sample_type_detail[grep("-01$", cn.segs$Sample)] <- "Primary Solid Tumor"
cn.segs$sample_type_detail[grep("-02$", cn.segs$Sample)] <- "Recurrent Solid Tumor"
cn.segs$sample_type_detail[grep("-03$", cn.segs$Sample)] <- "Primary Blood Derived Cancer"
cn.segs$sample_type_detail[grep("-04$", cn.segs$Sample)] <- "Recurrent Blood Derived Cancer"
cn.segs$sample_type_detail[grep("-05$", cn.segs$Sample)] <- "Additional - New Primary"
cn.segs$sample_type_detail[grep("-06$", cn.segs$Sample)] <- "Metastatic"
cn.segs$sample_type_detail[grep("-07$", cn.segs$Sample)] <- "Additional Metastatic"
cn.segs$sample_type_detail[grep("-08$", cn.segs$Sample)] <- "Human Tumor Original Cells"
cn.segs$sample_type_detail[grep("-09$", cn.segs$Sample)] <- "Primary Blood Derived Cancer - Bone Marrow"
cn.segs$sample_type_detail[grep("-10$", cn.segs$Sample)] <- "Blood Derived Normal"
cn.segs$sample_type_detail[grep("-11$", cn.segs$Sample)] <- "Solid Tissue Normal"
cn.segs$sample_type_detail[grep("-12$", cn.segs$Sample)] <- "Buccal Cell Normal"
cn.segs$sample_type_detail[grep("-13$", cn.segs$Sample)] <- "EBV Immortalized Normal"
cn.segs$sample_type_detail[grep("-14$", cn.segs$Sample)] <- "Bone Marrow Normal"
tmpclinic=clinic[,2:3]
colnames(tmpclinic)=c('Patient','Cohort')
cn.segs$Patient=substr(cn.segs$Sample,1,12)
cn.segs=plyr::join (cn.segs,tmpclinic, by='Patient',type='left')
cn.segs$Copy_Number <- as.integer(round(2*(2^cn.segs$Segment_Mean)))
sample.names <- unique(cn.segs$Sample)
sample.names <- sample.names
num.complexity.chr <- sapply(sample.names, numComplexity, seg=cn.segs)
num.complexity.score <- colSums(num.complexity.chr)
colnames(num.complexity.chr) <- sample.names
names(num.complexity.score) <- sample.names
save(num.complexity.score,file ="pan_can_ncin_score.RData")
