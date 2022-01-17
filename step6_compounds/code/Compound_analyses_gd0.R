rm(list=ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library(data.table)
setwd("/data/zhang/pancan_cin/step6_compounds/")
load('data/ccle_cin_score.RData')
ctrp_curve=fread('data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.data.curves_post_qc.txt')
ctrp_drug=fread('data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_compound.txt')
ctrp_cell=fread('data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_cell_line.txt')
ctrp_exp=fread('data/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_experiment.txt')
ctrp_curve=join(ctrp_curve,ctrp_exp,by='experiment_id')
ctrp_curve=join(ctrp_curve,ctrp_cell,by='master_ccl_id')
ctrp_curve=join(ctrp_curve,ctrp_drug,by='master_cpd_id')
load('data/ccle_cin_score.RData')
scores$Sample=toupper(unlist(lapply(scores$Sample, function(x){strsplit(x,'_')[[1]][1]})))
scores=subset(scores,scores$Genome.doublings==0)
ctrp_curve$Sample=ctrp_curve$ccl_name
library(ggplot2)
library(ggrepel)
Results <- data.frame()
for(i in 1:length(unique(ctrp_curve$cpd_name))){
  cpd <- paste(unique(ctrp_curve$cpd_name)[i])
  tmpmat=ctrp_curve[which(ctrp_curve$cpd_name==cpd),]
  tmpmat=merge(scores,tmpmat,by='Sample')
  cor_p <- cor.test(tmpmat$wgii, tmpmat$area_under_curve, method = "pearson")
  cor_s <- cor.test(tmpmat$wgii, tmpmat$area_under_curve, method = "spearman")
  Results=rbind(Results, data.frame(cpd, as.numeric(cor_s$estimate), as.numeric(cor_s$p.value), as.numeric(cor_p$estimate), as.numeric(cor_p$p.value)))
}

Results$FDR_s <- p.adjust(Results$as.numeric.cor_s.p.value., method = "fdr")
Results$FDR_p <- p.adjust(Results$as.numeric.cor_p.p.value., method = "fdr")
Results <- Results[order(Results$as.numeric.cor_s.p.value.),]
head(Results)
names(Results) <- c("cpd_name","Spearman_coef", "Spearman_pvalue", "Pearson_coef", "Pearson_pvalue", "Spearman_FDR", "Pearson_FDR")
Results=join(Results,ctrp_drug,by='cpd_name')
ii=colnames(Results)
Results=Results[,c("cpd_name", "Spearman_coef", "Spearman_pvalue", "Pearson_coef", 
                   "Pearson_pvalue", "Spearman_FDR", "Pearson_FDR",  "cpd_status", "inclusion_rationale", "gene_symbol_of_protein_target", 
                   "target_or_activity_of_compound")]
Results=Results[with(Results,order(Spearman_coef,Spearman_pvalue)),]
WriteXLS::WriteXLS(Results, "table/drug_sensitivity_NCIN_gd0.xlsx",AdjWidth = TRUE, BoldHeaderRow = TRUE, row.names = F)
Results2 <-Results[order(Results$Spearman_pvalue),]
Results2$Significant <- ifelse(Results2$Spearman_FDR < 0.05 & Results2$Spearman_coef > 0, "Pos", "Not Sig")
Results2$Significant[Results2$Spearman_FDR < 0.05 & Results2$Spearman_coef < 0] <- "Neg"
Results2$text <- "no"
Results2$text[c(which(Results2$Significant=='Pos')[1:6],which(Results2$Significant=='Neg')[1:6])] <- "yes"
Results2$text <- as.factor(Results2$text)
plot_drugs1 <- ggplot(Results2, aes(x = Spearman_coef, y = -log10(Spearman_pvalue))) +
  geom_point(aes(color = Significant)) +
  xlab("cor_spearman_wgii_auc") +
  ylab("-log10(p-value)") +
  scale_color_manual(values = c(tropiccolor[1], tropiccolor[1], tropiccolor[3])) +
  theme_cowplot(12) + theme(axis.text.x=element_text(size=13, colour="black"), axis.text.y=element_text(size=13, colour="black"), axis.title=element_text(size=16, colour="black"), legend.key = element_rect(fill="White"), legend.text = element_text(size=13), legend.title = element_text(size=14), legend.position = "bottom") +
  geom_text_repel(
    data = subset(Results2, text == "yes"),
    aes(label = target_or_activity_of_compound),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) +
  guides(color=guide_legend(title="Significant (FDR < 0.05)"))
plot_drugs1
Results <- data.frame()
for(i in 1:length(unique(ctrp_curve$cpd_name))){
  cpd <- paste(unique(ctrp_curve$cpd_name)[i])
  tmpmat=ctrp_curve[which(ctrp_curve$cpd_name==cpd),]
  tmpmat=merge(scores,tmpmat,by='Sample')
  cor_p <- cor.test(tmpmat$scin, tmpmat$area_under_curve, method = "pearson")
  cor_s <- cor.test(tmpmat$scin, tmpmat$area_under_curve, method = "spearman")
  Results=rbind(Results, data.frame(cpd, as.numeric(cor_s$estimate), as.numeric(cor_s$p.value), as.numeric(cor_p$estimate), as.numeric(cor_p$p.value)))
}

Results$FDR_s <- p.adjust(Results$as.numeric.cor_s.p.value., method = "fdr")
Results$FDR_p <- p.adjust(Results$as.numeric.cor_p.p.value., method = "fdr")
Results <- Results[order(Results$as.numeric.cor_s.p.value.),]
head(Results)
names(Results) <- c("cpd_name","Spearman_coef", "Spearman_pvalue", "Pearson_coef", "Pearson_pvalue", "Spearman_FDR", "Pearson_FDR")
Results=join(Results,ctrp_drug,by='cpd_name')
ii=colnames(Results)
Results=Results[,c("cpd_name", "Spearman_coef", "Spearman_pvalue", "Pearson_coef", 
                   "Pearson_pvalue", "Spearman_FDR", "Pearson_FDR",  "cpd_status", "inclusion_rationale", "gene_symbol_of_protein_target", 
                   "target_or_activity_of_compound")]
Results=Results[with(Results,order(Spearman_coef,Spearman_pvalue)),]
WriteXLS::WriteXLS(Results, "table/drug_sensitivity_SCIN_gd0.xlsx",AdjWidth = TRUE, BoldHeaderRow = TRUE, row.names = F)
Results2 <-Results[order(Results$Spearman_pvalue),]
Results2$Significant <- ifelse(Results2$Spearman_FDR < 0.05 & Results2$Spearman_coef > 0, "Pos", "Not Sig")
Results2$Significant[Results2$Spearman_FDR < 0.05 & Results2$Spearman_coef < 0] <- "Neg"
Results2$text <- "no"
Results2$text[c(which(Results2$Significant=='Pos')[1:6],which(Results2$Significant=='Neg')[1:6])] <- "yes"
Results2$text <- as.factor(Results2$text)
plot_drugs2 <- ggplot(Results2, aes(x = Spearman_coef, y = -log10(Spearman_pvalue))) +
  geom_point(aes(color = Significant)) +
  xlab("cor_spearman_scin_auc") +
  ylab("-log10(p-value)") +
  scale_color_manual(values = c(tropiccolor[1], tropiccolor[1], tropiccolor[3])) +
  theme_cowplot( 12) + theme(axis.text.x=element_text(size=13, colour="black"), axis.text.y=element_text(size=13, colour="black"), axis.title=element_text(size=16, colour="black"), legend.key = element_rect(fill="White"), legend.text = element_text(size=13), legend.title = element_text(size=14), legend.position = "bottom") +
  geom_text_repel(
    data = subset(Results2, text == "yes"),
    aes(label = target_or_activity_of_compound),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) +
  guides(color=guide_legend(title="Significant (FDR < 0.05)"))
plot_drugs2
plot2by2 <- plot_grid(plot_drugs1,plot_drugs2,labels=c("A", "B"), ncol = 2)
save_plot("plot/volcano_cin_vs_drug_sensitivity_gd0.pdf", plot2by2,ncol = 3, # we're saving a grid plot of 2 columns
          nrow = 1, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3)