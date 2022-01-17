rm(list=ls())
library(data.table)
library(plyr)
setwd("/data/zhang/pancan_cin/step6_compounds/")
structComplexity <- function( sample_name, seg, thresh = 1e6, autosomal=TRUE){
  seg <- seg[which(!seg$Chromosome %in% c("X","Y")),]
  tmpn=unique(seg$Chromosome)
  vect.num.arm <- rep(0, length(tmpn))
  names(vect.num.arm) <- tmpn
  seg.sample <- subset(seg, seg$Sample== sample_name)
  for (chr in tmpn){
    seg.chr <- subset(seg.sample, seg.sample$Chromosome == chr)
    modal.cn <- seg.chr$Copy_Number[which.max(seg.chr$Num_Probes)]
    a <- ((seg.chr$Copy_Number > modal.cn) | (seg.chr$Copy_Number < modal.cn))
    a <- a & ((seg.chr$End - seg.chr$Start) > thresh)
    vect.num.arm[chr] <- sum(a)
  }
  return(vect.num.arm)
}
cn.segs=fread("data/CCLE_copynumber_2013-12-03.seg.txt")
cn.segs$Sample=cn.segs$CCLE_name
ploidy_abs=gdata::read.xls('data/CCLE_ABSOLUTE_combined_20181227.xlsx',sheet = 'ABSOLUTE_combined.table')
ploidy_abs$Sample=as.character(ploidy_abs$CCLE_ID)
clinic=fread('data/Cell_lines_annotations_20181226.txt')
clinic$Sample=clinic$CCLE_ID
cn.segs=plyr::join (cn.segs,ploidy_abs, by='Sample',type='left')
cn.segs=plyr::join (cn.segs,clinic, by='Sample',type='left')
cn.segs$ploidy[which(is.na(cn.segs$ploidy))]=2
cn.segs$Copy_Number <- as.integer(round(2*(2^cn.segs$Segment_Mean)))
sample.names <- unique(cn.segs$Sample)
sample.names <- sample.names
struct.complexity.chr <- sapply(sample.names, structComplexity, seg=cn.segs)
struct.complexity.score <- colSums(struct.complexity.chr)
colnames(struct.complexity.chr) <- sample.names
names(struct.complexity.score) <- sample.names
save(struct.complexity.score,file ="ccle_scin_score.RData")
