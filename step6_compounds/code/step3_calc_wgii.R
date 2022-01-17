rm(list=ls())
library(data.table)
library(plyr)
setwd("/data/zhang/pancan_cin/step6_compounds/")
getWgii <- function(x,cnsegs){
  # data should have columns chr, start, end, cn
  # For each chromosome (except X,Y) find the percentage for which CN !=2
  # Then take the mean across all chromosomes
  data=subset(cnsegs,cnsegs$Sample==x)
  chrdata <- c()
  ploidytmp=as.integer(round(unique(data$ploidy)))
  for(chr in unique(data$Chromosome)){
    if(!(chr %in% 1:22)){ next }
    segs <- data[which(data$Chromosome ==chr),]
    seglengths <- as.numeric(segs$End) - as.numeric(segs$Start)
    x <- sum(seglengths[which(segs$Copy_Number !=ploidytmp)]) / sum(seglengths)
    chrdata <- c(chrdata, x)
  }
  return(mean(chrdata))
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
clinic=cn.segs[which(!duplicated(cn.segs$Sample)),]
cn.segs$Copy_Number <- as.integer(round(2*(2^cn.segs$Segment_Mean)))
sample.names <- unique(cn.segs$Sample)
sample.names <- sample.names
WGII.chr <- sapply(sample.names, getWgii, cnsegs=cn.segs)
WGII.score <- WGII.chr
names(WGII.score) <- sample.names
save(WGII.score,clinic,file ="ccle_wgii_score.RData")

