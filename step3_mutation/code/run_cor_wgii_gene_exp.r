rm(list = ls())
library(arules)
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library(data.table)
setwd("/data/zhang/pancan_cin/step3_mutation/")
cname=c("GBM", "OV", "LUAD", "LUSC", "PRAD", "UCEC", "BLCA", "TGCT", 
        "ESCA", "PAAD", "KIRP", "LIHC", "CESC", "SARC", "BRCA", "THYM", 
        "MESO", "COAD", "STAD", "CHOL", "KIRC", "THCA", "HNSC", "LAML", 
        "READ", "SKCM", "LGG", "DLBC", "KICH", "UCS", "ACC", "PCPG", 
        "UVM")
for(i in cname[25:30]){
  #shell_cmd=paste0(c('screen -dmS "session" sh -c "Rscript code/cor_ttest_wgii_gene_exp.r ',i,'; exec bash"'),collapse = '')
  shell_cmd=paste0(c('screen -dmS "session" sh -c "Rscript code/cor_ttest_wgii_cnv.r ',i,'; exec bash"'),collapse = '')
  print(shell_cmd)
  system(shell_cmd)
  }
