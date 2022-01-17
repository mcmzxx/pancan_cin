rm(list = ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library(data.table)
setwd("/data/zhang/pancan_cin/step3_mutation/")
load("table/wgii_gexp.RData")
#load("table/scin_gexp.RData")
gpos=lapply(finalresults,function(x){
  fn=unique(x$Organ)
  cind=x$Coef[1:50]
  pind=x$FDR[1:50]
  posg=data.frame(x[which(cind>0&pind<=0.1),'Gene'])
  negg=data.frame(x[which(cind<0&pind<=0.1),'Gene'])
  write.table(posg,file.path('table/cmap_query/wgii',paste0(fn,'_pos.txt')),quote = F,col.names = F,row.names = F)
write.table(negg,file.path('table/cmap_query/wgii',paste0(fn,'_neg.txt')),quote = F,col.names = F,row.names = F)})
