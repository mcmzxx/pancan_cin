rm(list=ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library(data.table)
library("colorspace")
tropiccolor=c("grey69","#0374be","#cc0f0f")
setwd("/data/zhang/pancan_cin/step6_compounds/")
load('data/ccle_cin_score.RData')
ctrp_curve=fread('data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.data.curves_post_qc.txt')
ctrp_drug=fread('data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_compound.txt')
ctrp_cell=fread('data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_cell_line.txt')
ctrp_exp=fread('data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_experiment.txt')
ctrp_curve=plyr::join(ctrp_curve, ctrp_exp, by='experiment_id')
ctrp_curve=plyr::join(ctrp_curve, ctrp_cell, by='master_ccl_id')
ctrp_curve=plyr::join(ctrp_curve, ctrp_drug, by='master_cpd_id')
load('data/ccle_cin_score.RData')
scores$Sample=toupper(unlist(lapply(scores$Sample, function(x) {
  strsplit(x, '_')[[1]][1]
})))
ctrp_curve$Sample=ctrp_curve$ccl_name
library(ggplot2)
library(ggrepel)
normfn= function(x){(x-min(x))/(max(x)-min(x))}
ctrp_curve$area_under_curve=1-normfn(ctrp_curve$area_under_curve)
Results <- data.frame()
for(i in 1:length(unique(ctrp_curve$cpd_name))){
  cpd <- paste(unique(ctrp_curve$cpd_name)[i])
  tmpmat=ctrp_curve[which(ctrp_curve$cpd_name==cpd),]
  tmpmat=merge(scores,tmpmat,by='Sample')
  cor_p <- cor.test(tmpmat$wgii, tmpmat$area_under_curve, method = "pearson")
  cor_s <- cor.test(tmpmat$wgii, tmpmat$area_under_curve, method = "spearman")
  Results=rbind(Results, data.frame(cpd,max(tmpmat$area_under_curve),median(tmpmat$area_under_curve),as.numeric(cor_s$estimate), as.numeric(cor_s$p.value), as.numeric(cor_p$estimate), as.numeric(cor_p$p.value)))
}

Results$FDR_s <- p.adjust(Results$as.numeric.cor_s.p.value., method = "fdr")
Results$FDR_p <- p.adjust(Results$as.numeric.cor_p.p.value., method = "fdr")
Results <- Results[order(Results$as.numeric.cor_s.p.value.),]
head(Results)
names(Results) <- c("cpd_name","max_1minus_AUC","median_1minus_AUC","Spearman_coef", "Spearman_pvalue", "Pearson_coef", "Pearson_pvalue", "Spearman_FDR", "Pearson_FDR")
Results=join(Results,ctrp_drug,by='cpd_name')
ii=colnames(Results)
Results=Results[,c("cpd_name","max_1minus_AUC","median_1minus_AUC", "Spearman_coef", "Spearman_pvalue", "Pearson_coef", 
                   "Pearson_pvalue", "Spearman_FDR", "Pearson_FDR",  "cpd_status", "inclusion_rationale", "gene_symbol_of_protein_target", 
                   "target_or_activity_of_compound")]
Results=Results[with(Results,order(Spearman_coef,Spearman_pvalue)),]
WriteXLS::WriteXLS(
  Results,
  "table/drug_sensitivity_WCIN_all_sample.xlsx",
  AdjWidth=TRUE,
  BoldHeaderRow=TRUE,
  row.names=F
)
Results2=Results[order(Results$Spearman_pvalue), ]
Results2$Significant <-ifelse(Results2$Spearman_FDR < 0.05 & Results2$Spearman_coef > 0,"positive","insignificant")
Results2$Significant[Results2$Spearman_FDR  < 0.05 & Results2$Spearman_coef< 0]="negative"
Results2$text="no"
Results2$text[c(which(Results2$Significant == 'positive')[1:6],
                which(Results2$Significant == 'negative')[1:6])]="yes"
Results2$text=as.factor(Results2$text)
plot_drugs1=ggplot(Results2, aes(x = Spearman_coef, y = median_1minus_AUC)) +
  geom_point(aes(color=Significant), size=0.4) +
  xlab("correlation coeffecient") +
  ylab("drug sensitivity") +
  scale_color_manual(values=c(tropiccolor[1], tropiccolor[2], tropiccolor[3])) +
  theme_cowplot(12) + theme(
    axis.text.x=element_text(size=12, colour="black"),
    axis.text.y=element_text(size=12, colour="black"),
    axis.title=element_text(size=12, colour="black"),
    legend.key = element_rect(color = NA, fill = NA),
    legend.key.size = unit(0.05, "cm"),
    legend.text=element_text(size=12),
    legend.title=element_text(size=12),
    legend.position="right",legend.margin=margin(l=-0.3, unit='cm')
  ) +
  geom_text_repel(
    data=subset(Results2, Significant == 'positive' & median_1minus_AUC>0.5)[1:6,],
    aes(label=cpd_name),
    size =3,max.overlaps = Inf,segment.size=0.2,
    nudge_x= -0.03, direction = "y",angle=0,
    arrow = arrow (length = unit (0.02, "npc"), ends = "last", angle = 0),
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines"),
    hjust=0)+
  geom_text_repel(
    data=subset(Results2, Significant == 'negative' & median_1minus_AUC>0.5)[1:6,],
    aes(label=cpd_name),
    size =3,max.overlaps = Inf,segment.size=0.2,
    nudge_x= -0.03, direction = "y",angle=0,
    arrow = arrow (length = unit (0.02, "npc"), ends = "last", angle = 0),
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines"),
    hjust=0)+
  guides(color=guide_legend(title="significant",nrow=3,byrow=T,override.aes = list(size=2)))
plot_drugs1
Results <- data.frame()
for(i in 1:length(unique(ctrp_curve$cpd_name))){
  cpd <- paste(unique(ctrp_curve$cpd_name)[i])
  tmpmat=ctrp_curve[which(ctrp_curve$cpd_name==cpd),]
  tmpmat=merge(scores,tmpmat,by='Sample')
  cor_p <- cor.test(tmpmat$scin, tmpmat$area_under_curve, method = "pearson")
  cor_s <- cor.test(tmpmat$scin, tmpmat$area_under_curve, method = "spearman")
  Results=rbind(Results, data.frame(cpd,max(tmpmat$area_under_curve),median(tmpmat$area_under_curve), as.numeric(cor_s$estimate), as.numeric(cor_s$p.value), as.numeric(cor_p$estimate), as.numeric(cor_p$p.value)))
}

Results$FDR_s <- p.adjust(Results$as.numeric.cor_s.p.value., method = "fdr")
Results$FDR_p <- p.adjust(Results$as.numeric.cor_p.p.value., method = "fdr")
Results <- Results[order(Results$as.numeric.cor_s.p.value.),]
head(Results)
names(Results) <- c("cpd_name","max_1minus_AUC","median_1minus_AUC","Spearman_coef", "Spearman_pvalue", "Pearson_coef", "Pearson_pvalue", "Spearman_FDR", "Pearson_FDR")
Results=join(Results,ctrp_drug,by='cpd_name')
ii=colnames(Results)
Results=Results[,c("cpd_name","max_1minus_AUC","median_1minus_AUC", "Spearman_coef", "Spearman_pvalue", "Pearson_coef", 
                   "Pearson_pvalue", "Spearman_FDR", "Pearson_FDR",  "cpd_status", "inclusion_rationale", "gene_symbol_of_protein_target", 
                   "target_or_activity_of_compound")]
Results=Results[with(Results,order(Spearman_coef,Spearman_pvalue)),]
WriteXLS::WriteXLS(
  Results,
  "table/drug_sensitivity_SCIN_all_sample.xlsx",
  AdjWidth=TRUE,
  BoldHeaderRow=TRUE,
  row.names=F
)
Results2=Results[order(Results$Spearman_pvalue), ]
Results2$Significant <-ifelse(Results2$Spearman_FDR  < 0.05 & Results2$Spearman_coef > 0,"positive","insignificant")
Results2$Significant[Results2$Spearman_FDR  < 0.05 & Results2$Spearman_coef< 0]="negative"
Results2$text="no"
Results2$text[c(which(Results2$Significant == 'positive')[1:6],
                which(Results2$Significant == 'negative')[1:6])]="yes"
Results2$text=as.factor(Results2$text)
plot_drugs2=ggplot(Results2, aes(x = Spearman_coef, y = median_1minus_AUC)) +
  geom_point(aes(color=Significant), size=0.4) +
  xlab("correlation coeffecient") +
  ylab("drug sensitivity") +
  scale_color_manual(values=c(tropiccolor[1], tropiccolor[2], tropiccolor[3])) +
  theme_cowplot(12) + theme(
    axis.text.x=element_text(size=12, colour="black"),
    axis.text.y=element_text(size=12, colour="black"),
    axis.title=element_text(size=12, colour="black"),
    legend.key = element_rect(color = NA, fill = NA),
    legend.key.size = unit(0.05, "cm"),
    legend.text=element_text(size=12),
    legend.title=element_text(size=12),
    legend.position="right",legend.margin=margin(l=-0.3, unit='cm')
  ) +
  geom_text_repel(
    data=subset(Results2, Significant == 'positive' & median_1minus_AUC>0.5)[1:6,],
    aes(label=cpd_name),
    size =3,max.overlaps = Inf,segment.size=0.2,
    nudge_x= -0.03, direction = "y",angle=0,
    arrow = arrow (length = unit (0.02, "npc"), ends = "last", angle = 0),
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines"),
    hjust=0)+
  geom_text_repel(
    data=subset(Results2, Significant == 'negative' & median_1minus_AUC>0.5)[1:6,],
    aes(label=cpd_name),
    size =3,max.overlaps = Inf,segment.size=0.2,
    nudge_x= -0.03, direction = "y",angle=0,
    arrow = arrow (length = unit (0.02, "npc"), ends = "last", angle = 0),
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines"),
    hjust=0)+
  guides(color=guide_legend(title="significant",nrow=3,byrow=T,override.aes = list(size=2)))
plot_drugs2
drug_legend=get_legend(plot_drugs2)
plot_drugs1=plot_drugs1+theme(legend.position = 'none')
plot_drugs2=plot_drugs2+theme(legend.position = 'none')
prow=plot_grid(plot_drugs1,plot_drugs2,nrow=2,rel_heights = c(1,1),labels = c('E','F'))
pdf('plot/volcano_cin_vs_drug_sensitivity.pdf', width=6, height=5.7)
plot_grid(prow,drug_legend,nrow=1,rel_widths= c(1,0.23))+theme(plot.margin=margin(r=-0.5,l=-0.1,unit="cm"))
dev.off()