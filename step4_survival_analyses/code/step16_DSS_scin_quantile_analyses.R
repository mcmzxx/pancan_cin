rm(list=ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
require(ggplot2)
library(survival)
library(survminer)
library(RColorBrewer)
setwd("/data/zhang/pancan_cin/step4_survival_analyses/")
table <- read.delim("../step2_genome_instability/data/TaylorCancerCell_TableS2.txt", na.strings = c("", " ", "na", "NA", "n.a.", "#N/A"))
rownames(table) <- as.character(table$Sample)
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer','Primary Solid Tumor'),]
rownames(Final_table)=Final_table$Sample
sels=intersect(rownames(table),rownames(Final_table))
Final_table=Final_table[sels,]
ii=c("age_at_initial_pathologic_diagnosis", 
     "initial_pathologic_dx_year",  "birth_days_to",
     "last_contact_days_to", "death_days_to","new_tumor_event_dx_days_to", 
     "OS", "OS.time", "DSS", "DSS.time", "DFI", "DFI.time", "PFI", 
     "PFI.time")
Final_table[,ii]=apply(Final_table[,ii], 2, as.numeric)
TCGA_cohorts=unique(Final_table$Cohort)
Table_survival_results <- data.frame()
pdf("plot/OS_wgii_quantile_plots.pdf", width = 8, height = 8)
for(c in 1:length(TCGA_cohorts)) {
  
  cohort_TCGA <- TCGA_cohorts[c]
  Final_table_merged_tumor <-subset(Final_table,Final_table$Cohort==cohort_TCGA)
  if(!all(is.na(Final_table_merged_tumor$DSS.time))){
    Final_table_merged_tumor$OS.time <- Final_table_merged_tumor$OS.time/365
    Final_table_merged_tumor$wgii_group=cut(Final_table_merged_tumor$wgii, quantile(Final_table_merged_tumor$wgii+runif(length(Final_table_merged_tumor$wgii), -0.01,0.01), seq(0,1,l=5), na.rm=TRUE), right=F,labels=c("0_25%","25_50%","50_75%","75_100%"))

Final_table_merged_tumor$SurvObj <- with(Final_table_merged_tumor, Surv(Final_table_merged_tumor$OS.time  , Final_table_merged_tumor$OS))
km_cluster <- survfit(SurvObj ~ wgii_group, data = Final_table_merged_tumor)
surv_results=survdiff(SurvObj ~ wgii_group, data = Final_table_merged_tumor)
surv_results_p.value <- 1 - pchisq(surv_results$chisq, length(surv_results$n) - 1)

msurv=Surv(Final_table_merged_tumor$OS.time, Final_table_merged_tumor$OS)
f=survfit(msurv~wgii_group,data= Final_table_merged_tumor);
km=coxph(msurv~wgii_group,data= Final_table_merged_tumor);
print(summary(km)$sctest['pvalue'])
color=c(c(brewer.pal(n = 9, name = 'Set1')),brewer.pal(n = 8, name = 'Dark2'))

#f=survfit(msurv~ cut( Final_table_merged_tumor$mut, quantile( Final_table_merged_tumor$mut+runif(length( Final_table_merged_tumor$mut), -0.01,0.01), seq(0,1,l=4), na.rm=TRUE), right=F))
sp=ggsurvplot(f, pval = TRUE, conf.int = F,
              risk.table = TRUE, risk.table.y.text.col = TRUE,xlim=c(0,5),palette ="Set1",main= cohort_TCGA,break.time.by=1,xlab="Time (years)")+ggtitle(cohort_TCGA)
print(sp)


group_n <- table(Final_table_merged_tumor$wgii_group)
surv_summ <- summary(km_cluster,times=c(1,3))
surv3=surv_summ$surv[seq(2,8,by=2)]

Table_survival_results_cohort <-cbind(data.frame(as.character(cohort_TCGA),nrow(Final_table_merged_tumor[complete.cases(Final_table_merged_tumor$OS.time) & complete.cases(Final_table_merged_tumor$OS),])),t(c(surv_results_p.value,group_n,surv3)))

Table_survival_results=rbind(Table_survival_results, Table_survival_results_cohort)
  }
} 
dev.off()
## multiple testing survival
ii=unlist(lapply(c("0_25%","25_50%","50_75%","75_100%"), function(x){paste0('sample_n_',x)}))
jj=unlist(lapply(c("0_25%","25_50%","50_75%","75_100%"), function(x){paste0('surv3_',x)}))
colnames(Table_survival_results)=c('Cohort','Non_na_sample_size','Log_rank_p',ii,jj)
Table_survival_results$FDR <- p.adjust(Table_survival_results$Log_rank_p, method = "fdr")
out <- Table_survival_results[order(Table_survival_results$Log_rank_p),]
write.table(out, "table/TCGA_OS_analyses_wgii_quantile_results.txt", sep="\t", row.names = F, quote=F)

