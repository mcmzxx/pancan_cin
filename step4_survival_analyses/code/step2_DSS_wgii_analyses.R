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
library(arules)
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
TCGA_Cohorts=unique(Final_table$Cohort)
Table_survival_results <- data.frame()
plist=list()
for(c in 1:length(TCGA_Cohorts)) {
  Cohort_TCGA <- as.character(TCGA_Cohorts[c])
  #Cohort_TCGA='SKCM'
  Final_table_merged_tumor <-subset(Final_table,Final_table$Cohort==Cohort_TCGA)
  if(!all(is.na(Final_table_merged_tumor$DSS.time))){
    Final_table_merged_tumor$DSS.time <- Final_table_merged_tumor$DSS.time/365
    # Final_table_merged_tumor$wgii_group=discretize(Final_table_merged_tumor$wgii,method='cluster',categories=2)
    # Final_table_merged_tumor$wgii_group=factor(Final_table_merged_tumor$wgii_group)
    # levels(Final_table_merged_tumor$wgii_group)=c('low','high')
    Final_table_merged_tumor$wgii_group <- ifelse(Final_table_merged_tumor$wgii<=median(Final_table_merged_tumor$wgii), "low", "high")
    Final_table_merged_tumor$wgii_group <- factor(Final_table_merged_tumor$wgii_group, levels=c("low", "high"))
    
    Final_table_merged_tumor$SurvObj <- with(Final_table_merged_tumor, Surv(Final_table_merged_tumor$DSS.time  , Final_table_merged_tumor$DSS))
    km_cluster <- survfit(SurvObj ~ wgii_group, data=Final_table_merged_tumor)
    surv_results=survdiff(SurvObj ~ wgii_group, data=Final_table_merged_tumor)
    surv_results_p.value <- 1 - pchisq(surv_results$chisq, length(surv_results$n) - 1)
    if(surv_results$obs[1] > surv_results$exp[1]) worst_median <- "high wgii"
    if(surv_results$obs[2] > surv_results$exp[2]) worst_median <- "low wgii"
    high_median <- as.numeric(length(which(complete.cases(Final_table_merged_tumor$DSS.time) & complete.cases(Final_table_merged_tumor$DSS) & Final_table_merged_tumor$wgii_group %in% "high")))
    low_median <- as.numeric(length(which(complete.cases(Final_table_merged_tumor$DSS.time) & complete.cases(Final_table_merged_tumor$DSS) & Final_table_merged_tumor$wgii_group %in% "low")))
    surv_summ1 <- summary(km_cluster,times=c(1))
    surv_summ5 <- summary(km_cluster,times=c(1,5))
    lensurv=length(surv_summ5$surv)
    if(lensurv==4 & (min(surv_summ5$n.risk)>2))
    {surv_summ=surv_summ5
    low_surv5=surv_summ$surv[2]
    high_surv5=surv_summ$surv[4]
    low_surv5_n=surv_summ$n.risk[2]
    high_surv5_n=surv_summ$n.risk[4]
    surv_end=5
    }else{
      surv_summ=surv_summ1
      low_surv5=surv_summ$surv[1]
      high_surv5=surv_summ$surv[2] 
      low_surv5_n=surv_summ$n.risk[1]
      high_surv5_n=surv_summ$n.risk[2] 
      surv_end=1
    }
    Table_survival_results_Cohort <- data.frame(Cohort=Cohort_TCGA,
                                                sample_number= nrow(Final_table_merged_tumor[complete.cases(Final_table_merged_tumor$DSS.time) & complete.cases(Final_table_merged_tumor$DSS),]),
                                                pvalue=surv_results_p.value,low_surv5=low_surv5,high_surv5=high_surv5,low_surv5_n=low_surv5_n,high_surv5_n=high_surv5_n,surv_end=surv_end
    )
    
    Table_survival_results=rbind(Table_survival_results, Table_survival_results_Cohort)
    
  }
}
Table_survival_results$FDR <- p.adjust(Table_survival_results$pvalue, method="fdr")
out <- Table_survival_results[order(Table_survival_results$pvalue),]
write.table(out, "table/TCGA_DSS_analyses_wgii_results.txt", sep="\t", row.names=F, quote=F)
print('wcin')
out$Cohort=as.character(out$Cohort)
outres=out[,c("Cohort", "sample_number", "pvalue", "low_surv5", "high_surv5","low_surv5_n", "high_surv5_n")]
outres$Cohort=unlist(lapply(1:length(out$Cohort),function(x){ifelse(out$surv_end[x]==1,paste0(out$Cohort[x],'$^*$'),out$Cohort[x])}))
print(xtable(outres, type="latex"),include.rownames=FALSE)
out=subset(out,out$pvalue<=0.05)
sigcan=as.character(out$Cohort)

plist=split(Final_table,Final_table$Cohort)
plist=lapply(sigcan, function(x){return(plist[[x]])})
names(plist)=sigcan
plist=lapply(plist,function(x){
  if(!all(is.na(x$DSS.time)))
  {
    x$wgii_group <- ifelse(x$wgii<=median(x$wgii), "low", "high")
    x$wgii_group <- factor(x$wgii_group, levels=c("low", "high"))
    # x$wgii_group=discretize(x$wgii,method='cluster',categories=2)
    # x$wgii_group=factor(x$wgii_group)
    # levels(x$wgii_group)=c('low','high')
    return(x)}
})
Final_table_merged_tumor=data.frame(do.call(rbind, plist))
Final_table_merged_tumor$Cohort=as.character(Final_table_merged_tumor$Cohort)
Final_table_merged_tumor$DSS.time=Final_table_merged_tumor$DSS.time / 365
Final_table_merged_tumor$SurvObj=with(
  Final_table_merged_tumor,
  Surv(
    Final_table_merged_tumor$DSS.time,
    Final_table_merged_tumor$DSS
  )
)
km_cluster=survfit(SurvObj ~ wgii_group+Cohort, data=Final_table_merged_tumor)

pl=ggsurvplot_facet(fit=km_cluster,data=Final_table_merged_tumor,facet.by="Cohort",ncol=4,
                    xlab="years",ylab='disease-specific survival',size=0.7,
                    censor.size=2.5,
                    pval=TRUE,
                    pval.coord=c(0,0.03),
                    xlim=c(0,10),
                    break.time.by=2,
                    palette=c(tropiccolor[1], tropiccolor[2]))+facet_wrap(.~Cohort,ncol=4)+ 
  geom_text(aes(x=surv_end+1,y=low_surv5+0.1, label=paste0('n=',low_surv5_n)),data=out,colour=tropiccolor[1])+
  geom_text(aes(x=surv_end-1,y=high_surv5-0.1, label=paste0('n=',high_surv5_n)),data=out,colour=tropiccolor[2])+
  geom_segment(mapping=aes(x=surv_end,xend=surv_end,y=0, yend=low_surv5),data=out,linetype=2, size=0.2) +
  geom_segment(mapping=aes(x=0,xend=surv_end,y=low_surv5, yend=out$low_surv5),data=out, linetype=2, size=0.2) +
  geom_segment(mapping=aes(x=0,xend=surv_end,y=high_surv5, yend=high_surv5),data=out, linetype=2, size=0.2) +
  geom_point(mapping=aes(x=surv_end, y=low_surv5), fill=tropiccolor[1],data=out, colour="black", shape=21, size=2, stroke=0.1) +
  geom_point(mapping=aes(x=surv_end, y=high_surv5), fill=tropiccolor[2],data=out, colour="black", shape=21, size=2, stroke=0.1)+
  theme(axis.line=element_line(),strip.background=element_blank(),legend.text=element_text(),legend.position='bottom',legend.margin=margin(l=-0.3,b=-0.2, t=-0.3,r=-0.3,unit='cm'))+guides(colour=guide_legend(title ='WGII score'))
pl=pl+theme(text=element_text(size=11))

library(grid)
gt=ggplot_gtable(ggplot_build(pl))
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
gtwgii=gt
Table_survival_results <- data.frame()
plist=list()
for(c in 1:length(TCGA_Cohorts)) {
  Cohort_TCGA <- as.character(TCGA_Cohorts[c])
  Final_table_merged_tumor <-subset(Final_table,Final_table$Cohort==Cohort_TCGA)
  if(!all(is.na(Final_table_merged_tumor$DSS.time))){
    Final_table_merged_tumor$DSS.time <- Final_table_merged_tumor$DSS.time/365
    # Final_table_merged_tumor$scin_group=discretize(Final_table_merged_tumor$scin,method='cluster',categories=2)
    # Final_table_merged_tumor$scin_group=factor(Final_table_merged_tumor$scin_group)
    # levels(Final_table_merged_tumor$scin_group)=c('low','high')
    Final_table_merged_tumor$scin_group <- ifelse(Final_table_merged_tumor$scin<=median(Final_table_merged_tumor$scin), "low", "high")
    Final_table_merged_tumor$scin_group <- factor(Final_table_merged_tumor$scin_group, levels=c("low", "high"))
    
    Final_table_merged_tumor$SurvObj <- with(Final_table_merged_tumor, Surv(Final_table_merged_tumor$DSS.time  , Final_table_merged_tumor$DSS))
    km_cluster <- survfit(SurvObj ~ scin_group, data=Final_table_merged_tumor)
    surv_results=survdiff(SurvObj ~ scin_group, data=Final_table_merged_tumor)
    surv_results_p.value <- 1 - pchisq(surv_results$chisq, length(surv_results$n) - 1)
    if(surv_results$obs[1] > surv_results$exp[1]) worst_median <- "high scin"
    if(surv_results$obs[2] > surv_results$exp[2]) worst_median <- "low scin"
    high_median <- as.numeric(length(which(complete.cases(Final_table_merged_tumor$DSS.time) & complete.cases(Final_table_merged_tumor$DSS) & Final_table_merged_tumor$scin_group %in% "high")))
    low_median <- as.numeric(length(which(complete.cases(Final_table_merged_tumor$DSS.time) & complete.cases(Final_table_merged_tumor$DSS) & Final_table_merged_tumor$scin_group %in% "low")))
    surv_summ1 <- summary(km_cluster,times=c(1))
    surv_summ5 <- summary(km_cluster,times=c(1,5))
    lensurv=length(surv_summ5$surv)
    if(lensurv==4 & (min(surv_summ5$n.risk)>2)){
      surv_summ=surv_summ5
      low_surv5=surv_summ$surv[2]
      high_surv5=surv_summ$surv[4]
      low_surv5_n=surv_summ$n.risk[2]
      high_surv5_n=surv_summ$n.risk[4]
      surv_end=5
    }else{
      surv_summ=surv_summ1
      low_surv5=surv_summ$surv[1]
      high_surv5=surv_summ$surv[2] 
      low_surv5_n=surv_summ$n.risk[1]
      high_surv5_n=surv_summ$n.risk[2] 
      surv_end=1
    }
    Table_survival_results_Cohort <- data.frame(Cohort=Cohort_TCGA,
                                                sample_number= nrow(Final_table_merged_tumor[complete.cases(Final_table_merged_tumor$DSS.time) & complete.cases(Final_table_merged_tumor$DSS),]),
                                                pvalue=surv_results_p.value,low_surv5=low_surv5,high_surv5=high_surv5,low_surv5_n=low_surv5_n,high_surv5_n=high_surv5_n,surv_end=surv_end
    )
    
    Table_survival_results=rbind(Table_survival_results, Table_survival_results_Cohort)
    rm(plot,low_surv5,high_surv5)
  }
}
Table_survival_results$FDR <- p.adjust(Table_survival_results$pvalue, method="fdr")
out <- Table_survival_results[order(Table_survival_results$pvalue),]
write.table(out, "table/TCGA_DSS_analyses_scin_results.txt", sep="\t", row.names=F, quote=F)
print('scin')
out$Cohort=as.character(out$Cohort)
outres=out[,c("Cohort", "sample_number", "pvalue", "low_surv5", "high_surv5","low_surv5_n", "high_surv5_n")]
outres$Cohort=unlist(lapply(1:length(out$Cohort),function(x){ifelse(out$surv_end[x]==1,paste0(out$Cohort[x],'$^*$'),out$Cohort[x])}))
print(xtable(outres, type="latex"),include.rownames=FALSE)
out=subset(out,out$pvalue<=0.05)
sigcan=as.character(out$Cohort)

plist=split(Final_table,Final_table$Cohort)
plist=lapply(sigcan, function(x){return(plist[[x]])})
names(plist)=sigcan
plist=lapply(plist,function(x){
  if(!all(is.na(x$DSS.time)))
  {
    x$scin_group <- ifelse(x$scin<=median(x$scin), "low", "high")
    x$scin_group <- factor(x$scin_group, levels=c("low", "high"))
    # x$scin_group=discretize(x$scin,method='cluster',categories=2)
    # x$scin_group=factor(x$scin_group)
    # levels(x$scin_group)=c('low','high')
    return(x)}
})
Final_table_merged_tumor=data.frame(do.call(rbind, plist))
Final_table_merged_tumor$Cohort=as.character(Final_table_merged_tumor$Cohort)
Final_table_merged_tumor$DSS.time=Final_table_merged_tumor$DSS.time / 365
Final_table_merged_tumor$SurvObj=with(
  Final_table_merged_tumor,
  Surv(
    Final_table_merged_tumor$DSS.time,
    Final_table_merged_tumor$DSS
  )
)
km_cluster=survfit(SurvObj ~ scin_group+Cohort, data=Final_table_merged_tumor)

pl=ggsurvplot_facet(fit=km_cluster,data=Final_table_merged_tumor,facet.by="Cohort",ncol=4,
                    xlab="years",ylab='disease-specific survival',size=0.7,
                    censor.size=2.5,
                    pval=TRUE,
                    pval.coord=c(0,0.03),
                    xlim=c(0,10),
                    break.time.by=2,
                    palette=c(tropiccolor[1], tropiccolor[2]))+facet_wrap(.~Cohort,ncol=4)+ 
  geom_text(aes(x=surv_end+1,y=low_surv5+0.1, label=paste0('n=',low_surv5_n)),data=out,colour=tropiccolor[1])+
  geom_text(aes(x=surv_end-1,y=high_surv5-0.1, label=paste0('n=',high_surv5_n)),data=out,colour=tropiccolor[2])+
  geom_segment(mapping=aes(x=surv_end,xend=surv_end,y=0, yend=low_surv5),data=out,linetype=2, size=0.2) +
  geom_segment(mapping=aes(x=0,xend=surv_end,y=low_surv5, yend=out$low_surv5),data=out, linetype=2, size=0.2) +
  geom_segment(mapping=aes(x=0,xend=surv_end,y=high_surv5, yend=high_surv5),data=out, linetype=2, size=0.2) +
  geom_point(mapping=aes(x=surv_end, y=low_surv5), fill=tropiccolor[1],data=out, colour="black", shape=21, size=2, stroke=0.1) +
  geom_point(mapping=aes(x=surv_end, y=high_surv5), fill=tropiccolor[2],data=out, colour="black", shape=21, size=2, stroke=0.1)+
  theme(axis.line=element_line(),strip.background=element_blank(),legend.text=element_text(),legend.position='bottom',legend.margin=margin(l=-0.3,b=-0.2, t=-0.3,r=-0.3,unit='cm'))+guides(colour=guide_legend(title ='scin score'))
pl=pl+theme(text=element_text(size=11))

library(grid)
gt=ggplot_gtable(ggplot_build(pl))
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
gtscin=gt
library(grid)
library(gridExtra)
pdf('plot/cin_dss_sig.pdf',height=9.7,width=6.5)
plot_grid(gtwgii,gtscin,ncol=1,labels=c('A','B'),rel_heights = c(2.13,3))
dev.off()
