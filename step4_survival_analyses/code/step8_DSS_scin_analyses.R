rm(list=ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
require(ggplot2)
library(survival)
library(survminer)
library("colorspace")
tropiccolor=diverging_hcl(401, "Tropic")[c(1,201,401)]
setwd("/data/zhang/pancan_cin/step4_survival_analyses/")

load('../step1_calculate_cin/data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer','Primary Solid Tumor'),]
rownames(Final_table)=Final_table$Sample

ii=c("age_at_initial_pathologic_diagnosis", 
     "initial_pathologic_dx_year",  "birth_days_to",
     "last_contact_days_to", "death_days_to","new_tumor_event_dx_days_to", 
     "OS", "OS.time", "DSS", "DSS.time", "DFI", "DFI.time", "PFI", 
     "PFI.time")
Final_table[,ii]=apply(Final_table[,ii], 2, as.numeric)
TCGA_cohorts=unique(Final_table$Cohort)


Table_survival_results <- data.frame()
plist=list()
for(c in 1:length(TCGA_cohorts)) {
  cohort_TCGA <- TCGA_cohorts[c]
  Final_table_merged_tumor <-subset(Final_table,Final_table$Cohort==cohort_TCGA)
  if(!all(is.na(Final_table_merged_tumor$DSS.time))){
    Final_table_merged_tumor$DSS.time <- Final_table_merged_tumor$DSS.time/365
    
    Final_table_merged_tumor$scin_group <- ifelse(Final_table_merged_tumor$scin<=median(Final_table_merged_tumor$scin), "Lower", "Higher")
    Final_table_merged_tumor$scin_group <- factor(Final_table_merged_tumor$scin_group, levels = c("Lower", "Higher"))
    
    Final_table_merged_tumor$SurvObj <- with(Final_table_merged_tumor, Surv(Final_table_merged_tumor$DSS.time  , Final_table_merged_tumor$DSS))
    km_cluster <- survfit(SurvObj ~ scin_group, data = Final_table_merged_tumor)
    surv_results=survdiff(SurvObj ~ scin_group, data = Final_table_merged_tumor)
    surv_results_p.value <- 1 - pchisq(surv_results$chisq, length(surv_results$n) - 1)
    
    if(surv_results$obs[1] > surv_results$exp[1]) worst_median <- "Higher scin"
    if(surv_results$obs[2] > surv_results$exp[2]) worst_median <- "Lower scin"
    
    higher_median <- as.numeric(length(which(complete.cases(Final_table_merged_tumor$DSS.time) & complete.cases(Final_table_merged_tumor$DSS) & Final_table_merged_tumor$scin_group %in% "Higher")))
    lower_median <- as.numeric(length(which(complete.cases(Final_table_merged_tumor$DSS.time) & complete.cases(Final_table_merged_tumor$DSS) & Final_table_merged_tumor$scin_group %in% "Lower")))
    
    surv_summ <- summary(km_cluster,times=c(1,5))
    Lower_surv5 <- ifelse(length(surv_summ$surv) > 2, surv_summ$surv[2], surv_summ$surv[1])
    Higher_surv5 <- ifelse(length(surv_summ$surv) > 2, surv_summ$surv[4], surv_summ$surv[2])
    plot <- ggsurvplot(km_cluster, legend = "none", xlab="Years",ylab='Disease Specific Survival',font.main = 18,size=0.2,
                       risk.table = TRUE, # absolute number and percentage at risk.
                       risk.table.col = "strata", # Change risk table color by groups
                       risk.table.fontsize = 4.5,censor.size=1,
                       pval = TRUE,
                       font.x =  16,
                       font.y = 16,
                       font.tickslab = 14,
                       xlim=c(0,10),
                       break.time.by=2.5,
                       palette = c(tropiccolor[1], tropiccolor[3]))
    plot$plot <- plot$plot +labs(title = cohort_TCGA) +
      theme(plot.title = element_text(hjust = 0.5))
    plot$plot <- plot$plot +
      geom_segment(x=5,xend=5,y=0, yend=Lower_surv5, linetype=2, size=0.2) +
      geom_segment(x=0,xend=5,y=Lower_surv5, yend=Lower_surv5, linetype=2, size=0.2) +
      geom_segment(x=0,xend=5,y=Higher_surv5, yend=Higher_surv5, linetype=2, size=0.2) +
      geom_point(x = 5, y = Lower_surv5, fill=tropiccolor[1], colour="black", shape=21, size=2, stroke=0.1) +
      geom_point(x = 5, y = Higher_surv5, fill=tropiccolor[3], colour="black", shape=21, size=2, stroke=0.1)
    plot$table <- plot$table +
      theme(panel.border = element_rect(colour="black", fill=NA),
            plot.title = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            legend.position = 'none'
      )
    names(plot$table$theme$axis.text.y$colour)=gsub('_group=','_',names(plot$table$theme$axis.text.y$colour))
    
    plist[[as.character(cohort_TCGA)]]=plot
    Table_survival_results_cohort <- data.frame(TCGA_cohort=cohort_TCGA,
                                                Sample_number= nrow(Final_table_merged_tumor[complete.cases(Final_table_merged_tumor$OS.time) & complete.cases(Final_table_merged_tumor$OS),]),
                                                Median_cutoff_pvalue=surv_results_p.value,Lower_surv5=Lower_surv5,Higher_surv5=Higher_surv5
    )
    
    Table_survival_results=rbind(Table_survival_results, Table_survival_results_cohort)
    rm(plot,Lower_surv5,Higher_surv5)
  }
}
Table_survival_results$FDR <- p.adjust(Table_survival_results$Median_cutoff_pvalue, method = "fdr")
out <- Table_survival_results[order(Table_survival_results$Median_cutoff_pvalue),]
write.table(out, "table/TCGA_DSS_analyses_scin_results.txt", sep="\t", row.names = F, quote=F)

ACC=plist[['ACC']]
BLCA=plist[['BLCA']]
BRCA=plist[['BRCA']]
CESC=plist[['CESC']]
CHOL=plist[['CHOL']]
COAD=plist[['COAD']]
DLBC=plist[['DLBC']]
ESCA=plist[['ESCA']]
GBM=plist[['GBM']]
HNSC=plist[['HNSC']]
KICH=plist[['KICH']]
KIRC=plist[['KIRC']]
KIRP=plist[['KIRP']]
LAML=plist[['LAML']]
LGG=plist[['LGG']]
LIHC=plist[['LIHC']]
LUAD=plist[['LUAD']]
LUSC=plist[['LUSC']]
MESO=plist[['MESO']]
OV=plist[['OV']]
PAAD=plist[['PAAD']]
PCPG=plist[['PCPG']]
PRAD=plist[['PRAD']]
READ=plist[['READ']]
SARC=plist[['SARC']]
SKCM=plist[['SKCM']]
STAD=plist[['STAD']]
TGCT=plist[['TGCT']]
THCA=plist[['THCA']]
THYM=plist[['THYM']]
UCEC=plist[['UCEC']]
UCS=plist[['UCS']]
UVM=plist[['UVM']]
xls=as.character(out$TCGA_cohort[which(out$Median_cutoff_pvalue<=0.05)])
xls
length(xls)
pr1=plot_grid(UCEC$plot,THCA$plot,KIRP$plot,ACC$plot,SARC$plot,ncol = 5)
pr2=plot_grid(LGG$plot,OV$plot,KICH$plot,ESCA$plot,READ$plot,ncol = 5)
pdf('plot/scin_dss_sig.pdf',height=6,width=12.5)
plot_grid(pr1,pr2,ncol = 1,rel_heights = c(1,1))
dev.off()