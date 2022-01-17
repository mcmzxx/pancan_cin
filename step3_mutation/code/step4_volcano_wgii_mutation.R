rm(list=ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library(data.table)
library("colorspace")
tropiccolor=c("#0374be","grey69","#cc0f0f")
setwd("/data/zhang/pancan_cin/step3_mutation/")
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer','Primary Solid Tumor'),]
Final_table$Sample = substr(Final_table$Sample, 1, 15)
rownames(Final_table)=Final_table$Sample
Final_table$Genome_doublings = as.factor(Final_table$`Genome doublings`)

Results = read.delim("table/wgii_allMutations_linear_model_results.txt")
Results = Results[order(Results$Pvalue),]
Results = Results[!Results$Gene %in% ".",]

Results$Significant = ifelse(Results$FDR < 0.05 & Results$MUTminusWT_diff > 0, "positive", "insignificant")
Results$Significant[Results$FDR < 0.05 & Results$MUTminusWT_diff < 0] = "negative"
Results$Significant2 = ifelse(Results$FDR < 0.05 & Results$Coef > 0, "positive", "insignificant")
Results$Significant2[Results$FDR < 0.05 & Results$Coef < 0] = "negative"
Results$text = "no"
selt=unique(c(1:10,which(Results$Gene %in% c('PIK3CA','KRAS','BRAF','NRAS', 'HRAS', "AURKB","BUB1",'MYC','RB1','BUB3','AURKA','PTEN')),which(Results$Significant2 %in% c("positive"))[1:5]))
selt=intersect(selt,which(Results$Significant2!='insignificant'))
Results$text[selt] = "yes"
Results$text = as.factor(Results$text)
plot1 =  ggplot(Results, aes(x = Coef, y = -log10(Pvalue))) +
  geom_point(aes(color = Significant2),size=0.2) +
  xlab("coefficient (linear model)") +
  ylab("-log10(p-value)") +
  scale_color_manual(values = c("negative"=tropiccolor[1], "insignificant"=tropiccolor[2], "positive"=tropiccolor[3])) +
  theme_cowplot(12)+ theme(axis.text.x=element_text(size=15, colour="black"),  legend.position = 'bottom',
                           axis.text.y=element_text(size=15, colour="black"), 
                           axis.title=element_text(size=15, colour="black"), 
                           legend.key = element_rect(fill="White"), 
                           legend.text = element_text(size=15), 
                           legend.title = element_text(size=15)) +
  geom_text_repel(
    data = subset(Results, text == "yes"),
    aes(label = Gene),
    size =3,max.overlaps =Inf,
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines")
  ) +
  guides(color=guide_legend(title="significant(FDR<0.05)"))

Results = read.delim("table/scin_allMutations_linear_model_results.txt")
Results = Results[order(Results$Pvalue),]
Results = Results[!Results$Gene %in% ".",]

Results$Significant = ifelse(Results$FDR < 0.05 & Results$MUTminusWT_diff > 0, "positive", "insignificant")
Results$Significant[Results$FDR < 0.05 & Results$MUTminusWT_diff < 0] = "negative"
Results$Significant2 = ifelse(Results$FDR < 0.05 & Results$Coef > 0, "positive", "insignificant")
Results$Significant2[Results$FDR < 0.05 & Results$Coef < 0] = "negative"
Results$text = "no"
selt=unique(c(1:10,which(Results$Gene %in% c('PIK3CA','KRAS','BRAF','NRAS', 'HRAS', "AURKB","BUB1",'MYC','RB1','BUB3','AURKA','PTEN')),which(Results$Significant2 %in% c("positive"))[1:5]))
selt=intersect(selt,which(Results$Significant2!='insignificant'))
Results$text[selt] = "yes"
Results$text = as.factor(Results$text)
plot2 = ggplot(Results, aes(x = Coef, y = -log10(Pvalue))) +
  geom_point(aes(color = Significant2),size=0.2) +
  xlab("coefficient (linear model)") +
  ylab("-log10(p-value)") +
  scale_color_manual(values = c("negative"=tropiccolor[1], "insignificant"=tropiccolor[2], "positive"=tropiccolor[3])) +
  theme_cowplot(12)+ theme(axis.text.x=element_text(size=15, colour="black"), legend.position = 'bottom',
                           axis.text.y=element_text(size=15, colour="black"), 
                           axis.title=element_text(size=15, colour="black"), 
                           legend.key = element_rect(fill="White"), 
                           legend.text = element_text(size=15), 
                           legend.title = element_text(size=15)) +
  geom_text_repel(
    data = subset(Results, text == "yes"),
    aes(label = Gene),
    size =3,max.overlaps =Inf,
  box.padding = unit(0.1, "lines"),
   point.padding = unit(0.1, "lines")
  ) +
  guides(color=guide_legend(title="significant(FDR<0.05)",override.aes = list(size=2)))
legend <- get_legend(
  plot2 + theme(legend.box.margin = margin(-1, -1, -1, 0))
)


plot1=plot1+theme(legend.position = "none",plot.margin=margin(r=-0.1,t=0.1,l=0.1,unit="cm"))
plot2=plot2+theme(legend.position = "none",plot.margin=margin(t=0.1,unit="cm"))
prow=plot_grid(plot1,plot2, align = "h", axis = "b",ncol=2,rel_widths = c(1,1), labels = c("WGII score vs mutations","SCIN score vs mutations"),label_fontface = 'plain',label_x=-0.05)
plot_grid(prow,legend, align = "v",ncol=1,rel_heights=c(1,0.1))
pdf('plot/lm_cin_mutation.pdf',width =7,height = 6)
plot_grid(prow,legend, align = "v",ncol=1,rel_heights=c(1,0.1))
dev.off()



Mut_table = fread("../step8_tme/data/mc3.v0.2.8.PUBLIC.maf.gz")
Mut_table$Sample = substr(Mut_table$Tumor_Sample_Barcode, 1, 15)

Final_table$TP53 = "wt"
Final_table$TP53[Final_table$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% "TP53"]] ="mut"
Final_table$TP53=factor(Final_table$TP53, levels = c("wt", "mut"))

Results_TP531 = data.frame()
Results_TP532 = data.frame()
for (i in 1:length(unique(Final_table$Cohort))) {
  cohort = paste(unique(Final_table$Cohort)[i])
  if (length(which(Final_table$TP53 %in% "wt" &Final_table$Cohort %in% cohort)) >= 20 &
      length(which(Final_table$TP53 %in% "mut" & Final_table$Cohort %in% cohort)) >= 20) {
    Final_table_cohort = Final_table[Final_table$Cohort %in% cohort, ]
    fit = lm(wgii ~ TP53, Final_table_cohort)
    pvalue = summary(fit)$coefficients[2, 4]
    coef = summary(fit)$coefficients[2, 1]
    wt= mean(Final_table_cohort$wgii[Final_table_cohort$TP53 %in% "wt"])
    mut=mean(Final_table_cohort$wgii[Final_table_cohort$TP53 %in% "mut"])
    WT_samples =length(which(Final_table$TP53 %in% "wt" &Final_table$Cohort %in% cohort))
    MUT_samples =length(which(Final_table$TP53 %in% "mut" &Final_table$Cohort %in% cohort))
    Results_TP531 = rbind(Results_TP531,data.frame(cohort,coef,pvalue,wt,mut,mut - wt,WT_samples,MUT_samples))
    fit = lm(scin ~ TP53, Final_table_cohort)
    pvalue = summary(fit)$coefficients[2, 4]
    coef = summary(fit)$coefficients[2, 1]
    Results_TP532 = rbind(Results_TP532,data.frame(cohort,coef,pvalue,wt,mut,mut - wt,WT_samples,MUT_samples))
  }
}
Results_TP53=Results_TP531
Results_TP53$FDR = p.adjust(Results_TP53$pvalue, method = "fdr")
Results_TP53 = Results_TP53[,c(1:3,9,4:8)]
names(Results_TP53) = c("Cohort", "Coef", "Pvalue", "FDR", "WTmean", "MUTmean", "MUTminusWT_diff", "WT_samples", "MUT_samples")
Results_TP53 = Results_TP53[order(Results_TP53$Coef, decreasing = T),]
order_names=as.character(Results_TP53$Cohort)
Results_TP53$Cohort = factor(Results_TP53$Cohort, levels=order_names)
Results_TP53$Significant = ifelse(Results_TP53$FDR < 0.05 & Results_TP53$Coef > 0, "positive", "insignificant")
Results_TP53$Significant[Results_TP53$FDR < 0.05 & Results_TP53$Coef < 0] = "negative"
Results_TP53$Significant=factor(Results_TP53$Significant,levels=c('insignificant','positive'))
plot1=ggplot(Results_TP53, aes(Cohort, Coef,fill=Significant)) +
  geom_bar(width = 0.9, stat = "identity") +
  labs(x = "", y = "coefficient (linear model)") +
  theme_cowplot(12)+ theme(axis.text.x = element_text(size =10, angle = -270, vjust=0.5),
                           axis.text.y=element_text(size=15, colour="black"), legend.position = 'bottom',
                           axis.title=element_text(size=15, colour="black"), 
                           legend.key = element_rect(fill="White"), 
                           legend.text = element_text(size=15), 
                           legend.title = element_text(size=15))
plot1=plot1+scale_fill_manual(values = c("negative"=tropiccolor[1], "insignificant"=tropiccolor[2], "positive"=tropiccolor[3])) 
Results_TP53=Results_TP532
Results_TP53$FDR = p.adjust(Results_TP53$pvalue, method = "fdr")
Results_TP53 = Results_TP53[,c(1:3,9,4:8)]
names(Results_TP53) = c("Cohort", "Coef", "Pvalue", "FDR", "WTmean", "MUTmean", "MUTminusWT_diff", "WT_samples", "MUT_samples")
Results_TP53 = Results_TP53[order(Results_TP53$Coef, decreasing = T),]
order_names=as.character(Results_TP53$Cohort)
Results_TP53$Cohort = factor(Results_TP53$Cohort, levels=order_names)
Results_TP53$Significant = ifelse(Results_TP53$FDR < 0.05 & Results_TP53$Coef > 0, "positive", "insignificant")
Results_TP53$Significant[Results_TP53$FDR < 0.05 & Results_TP53$Coef < 0] = "negative"
Results_TP53$Significant=factor(Results_TP53$Significant,levels=c('insignificant','positive'))
plot2=ggplot(Results_TP53, aes(Cohort, Coef,fill=Significant)) +
  geom_bar(width = 0.9, stat = "identity") +
  labs(x = "", y = "coefficient (linear model)") +
  theme_cowplot(12)+ theme(axis.text.x = element_text(size =10, angle = -270, vjust=0.5),
                           axis.text.y=element_text(size=15, colour="black"), legend.position = 'bottom',
                           axis.title=element_text(size=15, colour="black"), 
                           legend.key = element_rect(fill="White"), 
                           legend.text = element_text(size=15), 
                           legend.title = element_text(size=15))+guides(fill=guide_legend(title="significant(FDR<0.05)"))
plot2=plot2+scale_fill_manual(values = c("negative"=tropiccolor[1], "insignificant"=tropiccolor[2], "positive"=tropiccolor[3])) 
legend <- get_legend(
  plot2 + theme(legend.box.margin = margin(-1, -1, -1, 0))
)

plot1=plot1+theme(legend.position = "none",plot.margin=margin(r=-0.1,t=0.1,l=0.1,unit="cm"))
plot2=plot2+theme(legend.position = "none",plot.margin=margin(t=0.1,unit="cm"))
prow=plot_grid(plot1,plot2, align = "h", axis = "b",ncol=2,rel_widths = c(1,1),labels = c("WGII score vs TP53 mutation","SCIN score vs TP53 mutation"),label_fontface = 'plain',label_x=-0.1)
pdf('plot/cin_tp53_mut.pdf',width =7,height = 6)
plot_grid(prow,legend, align = "v",ncol=1,rel_heights=c(1,0.1))
dev.off()
# driver mutations from Catalog of Validated Oncogenic Mutations (https://www.cancergenomeinterpreter.org/mutations)
Driver = read.delim("data/catalog_of_validated_oncogenic_mutations.tsv")
## which ones are driver?
Mut_tablebak=Mut_table
Mut_table$ID = paste0("chr", Mut_table$Chromosome, ":g.", Mut_table$Start_Position, Mut_table$Tumor_Seq_Allele1, ">", Mut_table$Tumor_Seq_Allele2)
Mut_table$Driver = Mut_table$ID %in% Driver$gdna
Mut_table=subset(Mut_table,Mut_table$Driver==T)

selg=as.character(unique(Mut_table$Hugo_Symbol[Mut_table$Driver %in% "TRUE"]))
Final_tablebak=Final_table
Results1 = data.frame()
Results2= data.frame()
for(gene in selg){
  # at least 10 samples per group
  Final_table=Final_tablebak
  if(length(which(Final_table$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene & Mut_table$Driver %in% "TRUE"])) >=10 & length(which(!Final_table$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene & Mut_table$Driver %in% "TRUE"])) >=10){
    Final_table$Mutation_selected = "WT"
    Final_table$Mutation_selected[Final_table$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene]] = "passenger"
    Final_table$Mutation_selected[Final_table$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene & Mut_table$Driver %in% "TRUE"]] = "driver"
    #only driver vs wt
    Final_table = Final_table[Final_table$Mutation_selected %in% c("WT", "driver"),]
    Final_table$Mutation_selected = droplevels(factor(Final_table$Mutation_selected))
    Final_table$Mutation_selected = factor(Final_table$Mutation_selected, levels=c("WT", "driver"))
    fit = lm(wgii ~ Mutation_selected + Cohort, Final_table)
    pvalue = summary(fit)$coefficients[2,4]  
    coef = summary(fit)$coefficients[2,1]
    wt = mean(Final_table$wgii[Final_table$Mutation_selected %in% "WT"])
    mut = mean(Final_table$wgii[Final_table$Mutation_selected %in% "driver"])
    WT_samples = length(which(!Final_table$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene]))
    MUT_samples = length(which(Final_table$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene]))
    Results1=rbind(Results1, data.frame(gene, coef, pvalue, wt, mut, mut-wt, WT_samples, MUT_samples))
    fit = lm(scin ~ Mutation_selected + Cohort, Final_table)
    pvalue = summary(fit)$coefficients[2, 4]
    coef = summary(fit)$coefficients[2, 1]
    Results2=rbind(Results2, data.frame(gene, coef, pvalue, wt, mut, mut-wt, WT_samples, MUT_samples))
    
  }
}
Results=Results1
Results$FDR = p.adjust(Results$pvalue, method = "fdr")
Results = Results[order(Results$pvalue),c(1:3,9,4:8)]
names(Results) = c("Gene", "Coef", "Pvalue", "FDR", "WTmean", "MUT_driver_mean", "MUTminusWT_diff", "WT_samples", "MUT_driver_samples")
Results1=Results
write.table(Results, "table/wgii_DriverMutations_linear_model_results.txt", quote=F, sep="\t", row.names = F)
Results=Results2
Results$FDR = p.adjust(Results$pvalue, method = "fdr")
Results = Results[order(Results$pvalue),c(1:3,9,4:8)]
names(Results) = c("Gene", "Coef", "Pvalue", "FDR", "WTmean", "MUT_driver_mean", "MUTminusWT_diff", "WT_samples", "MUT_driver_samples")
Results2=Results
write.table(Results, "table/scin_DriverMutations_linear_model_results.txt", quote=F, sep="\t", row.names = F)
Results = Results1
Results = Results[order(Results$Pvalue),]
Results = Results[!Results$Gene %in% ".",]
Results$text = "no"
Results$text[1:10] = "yes"
Results$text = as.factor(Results$text)
Results$Significant = ifelse(Results$FDR < 0.05 & Results$MUTminusWT_diff > 0, "positive", "insignificant")
Results$Significant[Results$FDR < 0.05 & Results$MUTminusWT_diff < 0] = "negative"
Results$Significant2 = ifelse(Results$FDR < 0.05 & Results$Coef > 0, "positive", "insignificant")
Results$Significant2[Results$FDR < 0.05 & Results$Coef < 0] = "negative"
plot1 = ggplot(Results, aes(x = Coef, y = -log10(Pvalue))) +
  geom_point(aes(color = Significant2),size=0.2) +scale_y_continuous( limits = c(0,43))+
  xlab("coefficient (linear model)") +
  ylab("-log10(p-value)") +
  scale_color_manual(values = c("negative"=tropiccolor[1], "insignificant"=tropiccolor[2], "positive"=tropiccolor[3])) +
  theme_cowplot(12)+ theme(axis.text.x=element_text(size=15, colour="black"),  legend.position = 'bottom',
                           axis.text.y=element_text(size=15, colour="black"), 
                           axis.title=element_text(size=15, colour="black"), 
                           legend.key = element_rect(fill="White"), 
                           legend.text = element_text(size=15), 
                           legend.title = element_text(size=15)) +
  geom_text_repel(
    data = subset(Results, text == "yes"),
    aes(label = Gene),
    size =3,max.overlaps =Inf,
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines")
  ) +
  guides(color=guide_legend(title="significant(FDR<0.05)"))
Results = Results2
Results = Results[order(Results$Pvalue),]
Results = Results[!Results$Gene %in% ".",]
Results$text = "no"
Results$text[1:10] = "yes"
Results$text = as.factor(Results$text)
Results$Significant = ifelse(Results$FDR < 0.05 & Results$MUTminusWT_diff > 0, "positive", "insignificant")
Results$Significant[Results$FDR < 0.05 & Results$MUTminusWT_diff < 0] = "negative"
Results$Significant2 = ifelse(Results$FDR < 0.05 & Results$Coef > 0, "positive", "insignificant")
Results$Significant2[Results$FDR < 0.05 & Results$Coef < 0] = "negative"
plot2 = ggplot(Results, aes(x = Coef, y = -log10(Pvalue))) +
  geom_point(aes(color = Significant2),size=0.2) +scale_y_continuous( limits = c(0,33))+
  xlab("coefficient (linear model)") +
  ylab("-log10(p-value)") +
  scale_color_manual(values = c("negative"=tropiccolor[1], "insignificant"=tropiccolor[2], "positive"=tropiccolor[3])) +
  theme_cowplot(12)+ theme(axis.text.x=element_text(size=15, colour="black"), legend.position = 'bottom',
                           axis.text.y=element_text(size=15, colour="black"), 
                           axis.title=element_text(size=15, colour="black"), 
                           legend.key = element_rect(fill="White"), 
                           legend.text = element_text(size=15), 
                           legend.title = element_text(size=15)) +
  geom_text_repel(
    data = subset(Results, text == "yes"),
    aes(label = Gene),
    size =3,max.overlaps =Inf,
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines")
  ) +
  guides(color=guide_legend(title="significant(FDR<0.05)",override.aes = list(size=2)))
legend <- get_legend(
  plot2 + theme(legend.box.margin = margin(-1, -1, -1, 0))
)

plot1=plot1+theme(legend.position = "none",plot.margin=margin(r=-0.1,t=0.1,l=0.1,unit="cm"))
plot2=plot2+theme(legend.position = "none",plot.margin=margin(r=-0.1,t=0.1,l=0.1,unit="cm"))
prow=plot_grid(plot1,plot2, align = "h", axis = "b",ncol=2,rel_widths = c(1,1),labels = c("WGII score vs driver mutations","SCIN score vs driver mutations"),label_fontface = 'plain',label_x=-0.18)

pdf('plot/cin_driver_mut.pdf',width =7,height = 6)
plot_grid(prow,legend, align = "v",ncol=1,rel_heights=c(1,0.1))
dev.off()
