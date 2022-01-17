rm(list=ls())
library(data.table)
library(plyr)
setwd("/data/zhang/pancan_cin/step6_compounds/")
numComplexity <- function( sample_name, seg){
  require(Hmisc)
  seg <- seg[which(!seg$Chromosome %in% c("X","Y")),]
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
cn.segs=fread("data/CCLE_copynumber_2013-12-03.seg.txt")
cn.segs$Sample=cn.segs$CCLE_name
ploidy_abs=gdata::read.xls('data/CCLE_ABSOLUTE_combined_20181227.xlsx',sheet = 'ABSOLUTE_combined.table')
ploidy_abs$Sample=ploidy_abs$CCLE_ID
ploidy_abs$Sample=as.character(ploidy_abs$CCLE_ID)
clinic=fread('data/Cell_lines_annotations_20181226.txt')
clinic$Sample=clinic$CCLE_ID
cn.segs=plyr::join (cn.segs,ploidy_abs, by='Sample',type='left')
cn.segs=plyr::join (cn.segs,clinic, by='Sample',type='left')
cn.segs$ploidy[which(is.na(cn.segs$ploidy))]=2
cn.segs$Copy_Number <- as.integer(round(2*(2^cn.segs$Segment_Mean)))
sample.names <- unique(cn.segs$Sample)
sample.names <- sample.names
num.complexity.chr <- sapply(sample.names, numComplexity, seg=cn.segs)
num.complexity.score <- colSums(num.complexity.chr)
colnames(num.complexity.chr) <- sample.names
names(num.complexity.score) <- sample.names
save(num.complexity.score,file ="ccle_ncin_score.RData")
