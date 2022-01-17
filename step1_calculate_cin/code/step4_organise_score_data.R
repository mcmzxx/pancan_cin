rm(list=ls())
setwd('/data/zhang/pancan_cin/step1_calculate_cin/')
data1='pan_can_wgii_score.RData'
data2='pan_can_scin_score.RData'
data3='pan_can_ncin_score.RData'
load(data1)
load(data2)
load(data3)
scores = data.frame(scin=struct.complexity.score,ncin=num.complexity.score,wgii=WGII.score)
scores$Sample=rownames(scores)
clinic=gdata::read.xls('/data/zhang/pancan_cin/step8_tme/data/TCGA-CDR-SupplementalTableS1.xlsx')
clinic=data.frame(apply(clinic, 2, as.character))
scores$sample_type <- NA
scores$sample_type[grep("-01$", scores$Sample)] <- "Tumor"
scores$sample_type[grep("-02$", scores$Sample)] <- "Tumor"
scores$sample_type[grep("-03$", scores$Sample)] <- "Tumor"
scores$sample_type[grep("-04$", scores$Sample)] <- "Tumor"
scores$sample_type[grep("-05$", scores$Sample)] <- "Tumor"
scores$sample_type[grep("-06$", scores$Sample)] <- "Tumor"
scores$sample_type[grep("-07$", scores$Sample)] <- "Tumor"
scores$sample_type[grep("-08$", scores$Sample)] <- "Tumor"
scores$sample_type[grep("-09$", scores$Sample)] <- "Tumor"
scores$sample_type[grep("-10$", scores$Sample)] <- "Normal"
scores$sample_type[grep("-11$", scores$Sample)] <- "Normal"
scores$sample_type[grep("-12$", scores$Sample)] <- "Normal"
scores$sample_type[grep("-13$", scores$Sample)] <- "Normal"
scores$sample_type[grep("-14$", scores$Sample)] <- "Normal"
scores$sample_type_detail <- NA
scores$sample_type_detail[grep("-01$", scores$Sample)] <- "Primary Solid Tumor"
scores$sample_type_detail[grep("-02$", scores$Sample)] <- "Recurrent Solid Tumor"
scores$sample_type_detail[grep("-03$", scores$Sample)] <- "Primary Blood Derived Cancer"
scores$sample_type_detail[grep("-04$", scores$Sample)] <- "Recurrent Blood Derived Cancer"
scores$sample_type_detail[grep("-05$", scores$Sample)] <- "Additional - New Primary"
scores$sample_type_detail[grep("-06$", scores$Sample)] <- "Metastatic"
scores$sample_type_detail[grep("-07$", scores$Sample)] <- "Additional Metastatic"
scores$sample_type_detail[grep("-08$", scores$Sample)] <- "Human Tumor Original Cells"
scores$sample_type_detail[grep("-09$", scores$Sample)] <- "Primary Blood Derived Cancer - Bone Marrow"
scores$sample_type_detail[grep("-10$", scores$Sample)] <- "Blood Derived Normal"
scores$sample_type_detail[grep("-11$", scores$Sample)] <- "Solid Tissue Normal"
scores$sample_type_detail[grep("-12$", scores$Sample)] <- "Buccal Cell Normal"
scores$sample_type_detail[grep("-13$", scores$Sample)] <- "EBV Immortalized Normal"
scores$sample_type_detail[grep("-14$", scores$Sample)] <- "Bone Marrow Normal"
scores$Patient=substr(scores$Sample,1,12)
colnames(clinic)[2:3]=c('Patient','Cohort')
scores$Type=clinic$Cohort[match(scores$Patient,clinic$Patient)]
scores=plyr::join (scores,clinic[,-1], by='Patient',type='left')
scores=scores[which(!is.na(scores$Cohort)),]
#scores=scores[which(scores$scin<135),]
ploidy_abs=fread('../step8_tme/data/TCGA_mastercalls.abs_tables_JSedit.fixed.txt')
ploidy_abs$Sample=ploidy_abs$array
scores=plyr::join (scores,ploidy_abs, by='Sample',type='left')
scores$ploidy[which(is.na(scores$ploidy))]=2
save(scores,file ="data/pan_can_cin_score.RData")