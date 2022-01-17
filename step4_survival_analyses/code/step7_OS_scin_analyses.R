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
tropiccolor=c("#0374be","#cc0f0f")
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
  cohort_TCGA <- as.character(TCGA_cohorts[c])
  Final_table_merged_tumor <-subset(Final_table,Final_table$Cohort==cohort_TCGA)
  if(!all(is.na(Final_table_merged_tumor$DSS.time))){
    Final_table_merged_tumor$OS.time <- Final_table_merged_tumor$OS.time/365
    Final_table_merged_tumor$scin_group <- ifelse(Final_table_merged_tumor$scin<=median(Final_table_merged_tumor$scin), "Lower", "Higher")
    Final_table_merged_tumor$scin_group <- factor(Final_table_merged_tumor$scin_group, levels = c("Lower", "Higher"))
    Final_table_merged_tumor$SurvObj <- with(Final_table_merged_tumor, Surv(Final_table_merged_tumor$OS.time  , Final_table_merged_tumor$OS))
    km_cluster <- survfit(SurvObj ~ scin_group, data = Final_table_merged_tumor)
    surv_results=survdiff(SurvObj ~ scin_group, data = Final_table_merged_tumor)
    surv_results_p.value <- 1 - pchisq(surv_results$chisq, length(surv_results$n) - 1)
    if(surv_results$obs[1] > surv_results$exp[1]) worst_median <- "Higher scin"
    if(surv_results$obs[2] > surv_results$exp[2]) worst_median <- "Lower scin"
    higher_median <- as.numeric(length(which(complete.cases(Final_table_merged_tumor$OS.time) & complete.cases(Final_table_merged_tumor$OS) & Final_table_merged_tumor$scin_group %in% "Higher")))
    lower_median <- as.numeric(length(which(complete.cases(Final_table_merged_tumor$OS.time) & complete.cases(Final_table_merged_tumor$OS) & Final_table_merged_tumor$scin_group %in% "Lower")))
    surv_summ <- summary(km_cluster,times=c(1,5))
    Lower_surv5 <- ifelse(length(surv_summ$surv) > 2, surv_summ$surv[2], surv_summ$surv[1])
    Higher_surv5 <- ifelse(length(surv_summ$surv) > 2, surv_summ$surv[4], surv_summ$surv[2])
    n_low5=ifelse(length(surv_summ$n.risk)>2,surv_summ$n.risk[2],surv_summ$n.risk[1])
    n_high5=ifelse(length(surv_summ$n.risk)>2,surv_summ$n.risk[4],surv_summ$n.risk[2])
    Table_survival_results_cohort <- data.frame(Cohort=cohort_TCGA,
                                                sample_number= nrow(Final_table_merged_tumor[complete.cases(Final_table_merged_tumor$OS.time) & complete.cases(Final_table_merged_tumor$OS),]),
                                                pvalue=surv_results_p.value,lower_surv5=Lower_surv5,higher_surv5=Higher_surv5,n_low5=n_low5,n_high5=n_high5
    )
    
    Table_survival_results=rbind(Table_survival_results, Table_survival_results_cohort)
    rm(plot,Lower_surv5,Higher_surv5)
  }
}
Table_survival_results$FDR <- p.adjust(Table_survival_results$pvalue, method = "fdr")
out <- Table_survival_results[order(Table_survival_results$pvalue),]
write.table(out, "table/TCGA_OS_analyses_scin_results.txt", sep="\t", row.names = F, quote=F)


out=subset(out,out$pvalue<=0.05)
sigcan=as.character(out$Cohort)

plist=split(Final_table,Final_table$Cohort)
plist=lapply(sigcan, function(x){return(plist[[x]])})
names(plist)=sigcan
plist=lapply(plist,function(x){
  if(!all(is.na(x$OS.time)))
  {x$scin_group <- ifelse(x$scin<=median(x$scin), "low", "high")
  x$scin_group <- factor(x$scin_group, levels = c("low", "high"))
  return(x)}
})
Final_table_merged_tumor = data.frame(do.call(rbind, plist))
Final_table_merged_tumor$Cohort=as.character(Final_table_merged_tumor$Cohort)
Final_table_merged_tumor$OS.time = Final_table_merged_tumor$OS.time / 365
Final_table_merged_tumor$SurvObj = with(
  Final_table_merged_tumor,
  Surv(
    Final_table_merged_tumor$OS.time,
    Final_table_merged_tumor$OS
  )
)
km_cluster=survfit(SurvObj ~ scin_group+Cohort, data = Final_table_merged_tumor)

pl=ggsurvplot_facet(fit = km_cluster,data = Final_table_merged_tumor,facet.by = "Cohort",ncol = 4,
                    xlab="years",ylab='overall survival',size=0.7,
                    censor.size=2.5,
                    pval = TRUE,
                    pval.coord = c(0,0.03),
                    xlim=c(0,10),
                    break.time.by=2,
                    palette = c(tropiccolor[1], tropiccolor[2]))+facet_wrap(.~Cohort,ncol = 4)+ 
  geom_text(aes(x=6,y=lower_surv5+0.1, label=paste0('n=',n_low5)),data=out,colour=tropiccolor[1])+
  geom_text(aes(x=4,y=higher_surv5-0.1, label=paste0('n=',n_high5)),data=out,colour=tropiccolor[2])+
  geom_segment(mapping = aes(x=5,xend=5,y=0, yend=lower_surv5),data = out,linetype=2, size=0.2) +
  geom_segment(mapping=aes(x=0,xend=5,y=lower_surv5, yend=out$lower_surv5),data = out, linetype=2, size=0.2) +
  geom_segment(mapping=aes(x=0,xend=5,y=higher_surv5, yend=higher_surv5),data=out, linetype=2, size=0.2) +
  geom_point(mapping = aes(x = 5, y = lower_surv5), fill=tropiccolor[1],data = out, colour="black", shape=21, size=2, stroke=0.1) +
  geom_point(mapping=aes(x = 5, y = higher_surv5), fill=tropiccolor[2],data=out, colour="black", shape=21, size=2, stroke=0.1)+
  theme(axis.line=element_line(),strip.background=element_blank(),legend.text = element_text(),legend.position = 'bottom')+guides(colour=guide_legend(title ='scin score'))
pl=pl+theme(text = element_text(size=11))

library(grid)
gt = ggplot_gtable(ggplot_build(pl))
modify_aix=function(str_aix='axis-l',gtname=gt){
  axl=grep(str_aix,gt$layout$name)
  groblist=lapply(axl, function(x){gt$grobs[[x]]})
  ind=axl[which(unlist(lapply(groblist,function(x){x$name}))=='NULL')]
  indm=axl[which(unlist(lapply(groblist,function(x){x$name}))!='NULL')][1]
  modgrob=gt$grobs[[indm]]
  for(i in ind)
  {
    gt$grobs[[i]]=modgrob
    gt$grobs[[i]]$children$axis$grobs$`1`$children[[1]]$label=NULL
  }
  return(gt)
}
gt=modify_aix('axis-l',gt)
gt=modify_aix('axis-b',gt)
grid.draw(gt)
pdf('plot/scin_os_sig.pdf',height=6,width=6.5)
grid.draw(gt)
dev.off()
