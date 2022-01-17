rm(list=ls())
library(doParallel)
library("foreach")
library(ABSOLUTE)
setwd("/data/zhang/pancan_cin/PanCancer_CentrosomeAmplification/Genomic_instability/")
file_dir="/media/xiaoxiao/tiger1/cnv_absolute/single_files_pancancer"
tmpc=c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COADREAD", "DLBC", 
       "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", 
       "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", 
       "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", 
       "UVM")
data.dir='/media/xiaoxiao/tiger3/cnv_absolute_pancan'
for(j in tmpc){
  results.dir <-file.path(data.dir,j)
  absolute_rdata <- list.files(results.dir, pattern=".ABSOLUTE.RData",full.names = T)
  
  CreateReviewObject(obj.name =j, 
                     absolute.files = absolute_rdata,
                     indv.results.dir = paste0(results.dir, "/merged"), plot.modes = T,
                     copy_num_type = "total", verbose = TRUE)
  calls.path = file.path(results.dir,"merged", paste0(j,".PP-calls_tab.txt"))
  
  modes.path =  file.path(results.dir,"merged", paste0(j,".PP-modes.data.RData"))
  
  output.path =  file.path(results.dir,"merged", paste0(j,".abs_extract"))
  
  ExtractReviewedResults(reviewed.pp.calls.fn = calls.path, analyst.id = j, 
                         modes.fn = modes.path, out.dir.base = output.path, 
                         obj.name = j, copy_num_type="total")
  
  
}
