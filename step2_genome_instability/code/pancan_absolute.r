rm(list=ls())
library(doParallel)
library("foreach")
library(ABSOLUTE)
setwd("/data/zhang/pancan_cin/PanCancer_CentrosomeAmplification/Genomic_instability/")
# Table_final_results<- read.delim("Firebrowse_download/CNV_per_sample_all_TCGA.txt")
# tmpc=unique(as.character(Table_final_results$Cohort))
# resdir='/media/xiaoxiao/tiger1/cnv_absolute/single_files_pancancer'
# for(j in tmpc){
# cgh=Table_final_results[which(Table_final_results$Cohort==j),]
# dir.create(paste0(file.path(resdir, j)))
# tmpn=unique(as.character(cgh$Sample))
# for(i in tmpn){
#   tmp=cgh[which(cgh$Sample==i),]
#   write.table(tmp,file.path(file.path(resdir,j),paste0(i,'.hg19.cnv.txt')),col.names = T,sep='\t',quote = F,row.names=F)
# }}
file_dir="/media/xiaoxiao/tiger1/cnv_absolute/single_files_pancancer"
tmpc=c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COADREAD", "DLBC", 
       "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", 
       "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", 
       "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", 
       "UVM")
results.dir='/media/xiaoxiao/tiger3/cnv_absolute_pancan'
for(j in tmpc){

  dir.create(paste0(file.path(results.dir, j)))
  
  seg_files <- dir(file.path(file_dir,j), pattern = "\\.hg19.cnv.txt",recursive = T)
  files=file.path(file.path(file_dir,j), seg_files)
  doParallel::registerDoParallel(parallel::makePSOCKcluster(6))
  foreach(i = 1:length(files), .errorhandling = "pass",
          .packages = c("ABSOLUTE", "numDeriv"), .verbose = TRUE) %dopar%
    {
      # redefine function in each loop iteration to avoid error; run with default settings
      absolute_run_complete  <- function(seg.dat.fn, sample.name, results.dir) {
        sigma.p = 0
        max.sigma.h = 0.02
        min.ploidy <- 0.95
        max.ploidy <- 10
        primary.disease = j
        platform = "SNP_6.0"
        max.as.seg.count <- 1500
        max.non.clonal <- 0
        max.neg.genome <- 0
        copy_num_type <- "total"
        
        # creat log and results directory 
        log.dir <- paste0(results.dir, "/logs")
        # if (!file.exists(log.dir)){
        #   dir.create(log.dir, recursive = TRUE)
        # } else {
        #   message(paste0("log directory already made"))
        # }
        # if (!file.exists(results.dir)){
        #   for(i in seq_along(results.dir)){
        #     dir.create(results.dir[i], recursive = TRUE)}
        # } else {
        #   message(paste0("results directories already made"))
        # }
        # 
        # open log file and run aboslute - close sink after finishing file
        sink(file = file.path(log.dir, paste(sample.name, ".abs.log", sep="")), type = "output", append = T)
        RunAbsolute(seg.dat.fn, sigma.p, max.sigma.h, min.ploidy, max.ploidy,
                    primary.disease, platform, sample.name, results.dir,
                    max.as.seg.count, max.non.clonal, max.neg.genome,
                    copy_num_type, maf.fn =NULL, min.mut.af =NULL, verbose = TRUE)
       sink()
      }
      
      absolute_run_complete(
        seg.dat.fn = files[i],
        results.dir = file.path(results.dir,j),
        sample.name =substr(strsplit(files[i],'/')[[1]][8],1,27)
      )
    }
  }
results.dir <-paste0("/media/xiaoxiao/tiger1/cnv_absolute/results/absolute")

# register a cluster - will find how many nodes you have local



# combine purity estimates with edited CreateReviewObject from Andrew
source(file = "~/Projects/src/createreviewobject.R")
library(ABSOLUTE)
results.dir <-paste0("/media/xiaoxiao/cnv/cnv_absolute/results/absolute")
absolute_rdata <- list.files(results.dir, pattern=".ABSOLUTE.RData",full.names = T)

CreateReviewObject(obj.name = "COADREAD", 
                   absolute.files = absolute_rdata,
                   indv.results.dir = paste0(results.dir, "/merged"), plot.modes = T,
                   copy_num_type = "total", verbose = TRUE)
calls.path = file.path(results.dir,"merged", "COADREAD.PP-calls_tab.txt")

modes.path =  file.path(results.dir,"merged", "COADREAD.PP-modes.data.RData")

output.path =  file.path(results.dir,"merged", "COADREAD.abs_extract")

ExtractReviewedResults(reviewed.pp.calls.fn = calls.path, analyst.id = "COADREAD", 
                       modes.fn = modes.path, out.dir.base = output.path, 
                       obj.name = "COADREAD", copy_num_type="total")

