# CNA analyses

```{r load TCGA copy number data}
setwd("/data/zhang/pancan_cin/PanCancer_CentrosomeAmplification/Genomic_instability/")
system("mkdir Firebrowse_download")
setwd("Firebrowse_download")

TCGA_cohorts_table <- read.delim("../../TCGA_cohorts.txt")
TCGA_cohorts <- TCGA_cohorts_table$Cohort

Table_final_results <- data.frame()

for(i in 1:length(TCGA_cohorts)) {
  
  cohort_TCGA <- TCGA_cohorts[i]
  
  ##### Download from firebrowse
  system(paste("wget http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/", cohort_TCGA,"/20160128/gdac.broadinstitute.org_", cohort_TCGA,".Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0.tar.gz", sep=""))
  system(paste("tar -zxvf gdac.broadinstitute.org_", cohort_TCGA,".Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0.tar.gz", sep=""))
  system(paste("rm gdac.broadinstitute.org_", cohort_TCGA,".Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0.tar.gz", sep=""))
  
  CNV_table <- read.delim(paste("gdac.broadinstitute.org_", cohort_TCGA,".Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0/", cohort_TCGA,".snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt", sep=""), na.strings = c("null", "", " ", "NA", "na"))
  
  CNV_table$sample_name=substr(CNV_table$Sample,1,15)
  
  CNV_table$sample_type <- NA
  CNV_table$sample_type[grep(".01$", CNV_table$sample_name)] <- "Tumor"
  CNV_table$sample_type[grep(".02$", CNV_table$sample_name)] <- "Tumor"
  CNV_table$sample_type[grep(".03$", CNV_table$sample_name)] <- "Tumor"
  CNV_table$sample_type[grep(".04$", CNV_table$sample_name)] <- "Tumor"
  CNV_table$sample_type[grep(".05$", CNV_table$sample_name)] <- "Tumor"
  CNV_table$sample_type[grep(".06$", CNV_table$sample_name)] <- "Tumor"
  CNV_table$sample_type[grep(".07$", CNV_table$sample_name)] <- "Tumor"
  CNV_table$sample_type[grep(".08$", CNV_table$sample_name)] <- "Tumor"
  CNV_table$sample_type[grep(".09$", CNV_table$sample_name)] <- "Tumor"
  CNV_table$sample_type[grep(".10$", CNV_table$sample_name)] <- "Normal"
  CNV_table$sample_type[grep(".11$", CNV_table$sample_name)] <- "Normal"
  CNV_table$sample_type[grep(".12$", CNV_table$sample_name)] <- "Normal"
  CNV_table$sample_type[grep(".13$", CNV_table$sample_name)] <- "Normal"
  CNV_table$sample_type[grep(".14$", CNV_table$sample_name)] <- "Normal"
  
  
  CNV_table$sample_type_detail <- NA
  CNV_table$sample_type_detail[grep(".01$", CNV_table$sample_name)] <- "Primary Solid Tumor"
  CNV_table$sample_type_detail[grep(".02$", CNV_table$sample_name)] <- "Recurrent Solid Tumor"
  CNV_table$sample_type_detail[grep(".03$", CNV_table$sample_name)] <- "Primary Blood Derived Cancer"
  CNV_table$sample_type_detail[grep(".04$", CNV_table$sample_name)] <- "Recurrent Blood Derived Cancer"
  CNV_table$sample_type_detail[grep(".05$", CNV_table$sample_name)] <- "Additional - New Primary"
  CNV_table$sample_type_detail[grep(".06$", CNV_table$sample_name)] <- "Metastatic"
  CNV_table$sample_type_detail[grep(".07$", CNV_table$sample_name)] <- "Additional Metastatic"
  CNV_table$sample_type_detail[grep(".08$", CNV_table$sample_name)] <- "Human Tumor Original Cells"
  CNV_table$sample_type_detail[grep(".09$", CNV_table$sample_name)] <- "Primary Blood Derived Cancer - Bone Marrow"
  CNV_table$sample_type_detail[grep(".10$", CNV_table$sample_name)] <- "Blood Derived Normal"
  CNV_table$sample_type_detail[grep(".11$", CNV_table$sample_name)] <- "Solid Tissue Normal"
  CNV_table$sample_type_detail[grep(".12$", CNV_table$sample_name)] <- "Buccal Cell Normal"
  CNV_table$sample_type_detail[grep(".13$", CNV_table$sample_name)] <- "EBV Immortalized Normal"
  CNV_table$sample_type_detail[grep(".14$", CNV_table$sample_name)] <- "Bone Marrow Normal"
  
  
  CNV_table$Cohort <- cohort_TCGA
  
  Table_final_results=rbind(Table_final_results, CNV_table)
  
  print(cohort_TCGA)
  
} # cohorts

# Merge all data
write.table(Table_final_results, "CNV_per_sample_all_TCGA.txt", quote=F, sep="\t", row.names = F)

# return to directory of analyses
setwd("../")