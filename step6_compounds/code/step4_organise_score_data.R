rm(list=ls())
setwd('/data/zhang/pancan_cin/step6_compounds/')
#lohtab=gdata::read.xls('data/CCLE_ABSOLUTE_combined_20181227.xlsx',sheet = 'ABSOLUTE_combined.segtab')
data1='ccle_wgii_score.RData'
data2='ccle_scin_score.RData'
data3='ccle_ncin_score.RData'
load(data1)
load(data2)
load(data3)
scores = data.frame(scin=struct.complexity.score,ncin=num.complexity.score,wgii=WGII.score)
scores$Sample=rownames(scores)

scores=plyr::join (scores,clinic, by='Sample',type='left')
ii=c("scin", "ncin", "wgii", "Sample",  "CCLE_ID", "depMapID", 
     "call.status", "purity", "ploidy", "Genome.doublings", "delta", 
     "Coverage.for.80..power", "Cancer.DNA.fraction", "Subclonal.genome.fraction", 
     "tau", "E_CR", "Name", "Pathology", "Site_Primary", "Site_Subtype1", 
     "Site_Subtype2", "Site_Subtype3", "Histology", "Hist_Subtype1", 
     "Hist_Subtype2", "Hist_Subtype3", "Gender", "Life_Stage", "Age", 
     "Race", "Geo_Loc", "inferred_ethnicity", "Site_Of_Finding", "Disease", 
     "Annotation_Source", "Original.Source.of.Cell.Line", "Characteristics", 
     "Growth.Medium", "Supplements", "Freezing.Medium", "Doubling.Time.from.Vendor", 
     "Doubling.Time.Calculated.hrs", "type", "type_refined", "PATHOLOGIST_ANNOTATION", 
     "mutRate", "tcga_code")
scores=scores[,ii]
scores=scores[which(!is.na(scores$Site_Primary)),]
#scores=scores[which(scores$scin<165),]
save(scores,file ="data/ccle_cin_score.RData")