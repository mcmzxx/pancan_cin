rm(list=ls())
library(data.table)
library(gridExtra)
library(gr)
setwd("/data/zhang/pancan_cin/step4_survival_analyses/")
TCGA_cohorts_table <- read.delim("../TCGA_cohorts.txt")
TCGA_cohorts <- TCGA_cohorts_table$Cohort
# setwd("Firebrowse_download")
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer','Primary Solid Tumor'),]

rownames(Final_table)=Final_table$Sample
clinic_list=list()
for(i in 1:length(TCGA_cohorts)) {
  cohort_TCGA <- as.character(TCGA_cohorts[i])
  # if(length(grep(paste("wget http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/", cohort_TCGA, "/20160128/gdac.broadinstitute.org_", cohort_TCGA,".Merge_Clinical.Level_1.2016012800.0.0.tar.gz", sep=""), list.files())) ==0 ){
  #   system(paste("wget http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/", cohort_TCGA, "/20160128/gdac.broadinstitute.org_", cohort_TCGA,".Merge_Clinical.Level_1.2016012800.0.0.tar.gz", sep=""))
  #   system(paste("tar -zxvf gdac.broadinstitute.org_", cohort_TCGA,".Merge_Clinical.Level_1.2016012800.0.0.tar.gz", sep=""))
  #   system(paste("rm gdac.broadinstitute.org_", cohort_TCGA,".Merge_Clinical.Level_1.2016012800.0.0.tar.gz", sep=""))
  # }
  file_str=paste("gdac.broadinstitute.org_", cohort_TCGA,".Merge_Clinical.Level_1.2016012800.0.0/", cohort_TCGA,".merged_only_clinical_clin_format.txt", sep="")
  d=read.table(file.path('Firebrowse_download',file_str),header=T,sep='\t',row.names=1,fill=T,stringsAsFactors=F,quote = '$')
  colnames(d) = toupper(unlist(lapply(as.character(d["patient.bcr_patient_barcode", ]), function(x) {substr(x, 1, 12)})))
  mat = subset(Final_table, Final_table$Patient %in% colnames(d))
  mat$radicr = as.character(d["patient.radiations.radiation.measure_of_response", mat$Patient])
  mat$chemcr = as.character(d["patient.drugs.drug.measure_of_response", mat$Patient])
  mat$chemcr[which(mat$therapy != "chemotherapy")] = NA
  mat$chemcr[which(mat$chemcr == "complete response")] = "responder"
  mat$chemcr[which(mat$chemcr == "partial response")] = "responder"
  mat$chemcr[which(mat$chemcr != "responder" & (!is.na(mat$chemcr)))] = "non-responder"
  mat$radicr[which(mat$radicr == "complete response")] = "responder"
  mat$radicr[which(mat$radicr == "partial response")] = "responder"
  mat$radicr[which(mat$radicr != "responder" & (!is.na(mat$radicr)))] = "non-responder"
  clinic_list[[cohort_TCGA]] = mat
}
Final_table=data.frame(do.call(rbind,clinic_list))
Final_table$chemcr=as.factor(Final_table$chemcr)
Final_table$radicr=as.factor(Final_table$radicr)
radimat=Final_table[complete.cases(Final_table$radicr),]
chemmat=Final_table[complete.cases(Final_table$chemcr),]
gcolor=unname(yarrr::piratepal("info2"))[1:2]
my_comparisons <- list(c("responder","non-responder"))
pl <- ggboxplot(radimat, x = "radicr", y = "wgii",fill = "radicr",width=0.3,size=0.1, outlier.colour = NA)+ylab('WGII score')+
  scale_fill_manual("radicr",values=gcolor)+
  stat_compare_means(comparisons = my_comparisons,aes(label = paste0("p = ", ..p.format..)))+
  theme(legend.position = 'none',axis.text.y=element_text(size=12, colour="black"), 
        axis.text.x=element_blank(),
        axis.title.x=element_blank())+ylim(-0.1,1.2)
pl
pl1 <- ggboxplot(radimat, x = "radicr", y = "scin",fill = "radicr",width=0.3,size=0.1, outlier.colour = NA)+ylab('SCIN score')+
  scale_fill_manual("radicr", values=gcolor)+
  stat_compare_means(comparisons = my_comparisons,aes(label = paste0("p = ", ..p.format..)),label.y = 80)+
  theme(legend.position = 'none',axis.text.y=element_text(size=12, colour="black"), 
        axis.text.x=element_blank(),
        axis.title.x=element_blank())
pl1
prow1=plot_grid(pl,pl1)
pdf("plot/cin_vs_radicr.pdf",width =4,height =2)
plot_grid(pl,pl1)
dev.off()


pl <- ggboxplot(chemmat, x = "chemcr", y = "wgii",fill = "chemcr",width=0.3,size=0.1, outlier.colour = NA)+ylab('WGII score')+
  scale_fill_manual("chemcr", values=gcolor)+
  stat_compare_means(comparisons = my_comparisons,aes(label = paste0("p = ", ..p.format..)))+
  theme(legend.position = 'none',axis.text.y=element_text(size=12, colour="black"), 
        axis.text.x=element_blank(),
        axis.title.x=element_blank())+ylim(-0.1,1.2)
pl
pl1 <- ggboxplot(chemmat, x = "chemcr", y = "scin",fill = "chemcr",width=0.3,size=0.1, outlier.colour = NA)+ylab('SCIN score')+
  scale_fill_manual("chemcr", values=gcolor)+
  stat_compare_means(comparisons = my_comparisons,aes(label = paste0("p = ", ..p.format..)),label.y = 80)+
  theme(legend.position = 'right',axis.text.y=element_text(size=12, colour="black"), 
        axis.text.x=element_blank(),
        axis.title.x=element_blank())+guides(fill=guide_legend(title = 'therapy\nresponse'))
legend=get_legend(pl1)
pl1=pl1+theme(legend.position = 'none')
prow2=plot_grid(pl,pl1)
prow=plot_grid(prow1,prow2,nrow=2,rel_widths = c(1,1),labels = c('C','D'))
pdf("plot/cin_vs_chemcr.pdf",width =4,height =2)
plot_grid(prow2,legend,rel_widths = c(1,0.1),nrow = 1)
dev.off()
pdf("plot/cin_vs_therapy.pdf",width =6,height =4)
plot_grid(prow,legend,rel_widths = c(1,0.3),nrow = 1)
dev.off()
