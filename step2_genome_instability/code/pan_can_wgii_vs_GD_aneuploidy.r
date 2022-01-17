rm(list=ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library(sna)
library("ggbeeswarm")
library(cowplot)
library("colorspace")
library(readxl)
tropiccolor=c("grey69","#0374be","#cc0f0f")
gcolor=yarrr::piratepal("appletv") 
setwd("/data/zhang/pancan_cin/step2_genome_instability/")
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer','Primary Solid Tumor'),]
Final_table$Sample <- substr(Final_table$Sample, 1, 15)
rownames(Final_table)=Final_table$Sample
Final_table$Genome_doublings <- as.factor(Final_table$`Genome doublings`)
table <-read_xlsx("../step2_genome_instability/data/1-s2.0-S1535610818301119-mmc2.xlsx",skip = 1)
Final_table <- merge(Final_table, table, by = 'Sample')
tmpc=unique(Final_table$Cohort)
ii=c("1p", "1q", 
     "2p", "2q", "3p", "3q", "4p", "4q", "5p", "5q", "6p", "6q", "7p", 
     "7q", "8p", "8q", "9p", "9q", "10p", "10q", "11p", "11q", "12p", 
     "12q", "13q", "14q", "15q", "16p", "16q", "17p", "17q", "18p", 
     "18q", "19p", "19q", "20p", "20q", "21q", "22q", "1", "2", "3", 
     "4", "5", "6", "7", "8", "9", "10", "11", "12", "16", "17", "18", 
     "19", "20")
colnames(Final_table)[which(colnames(Final_table) %in% ii)]=paste0('X',ii)
ii=paste0('X',ii)
Results <- data.frame()
  tmpmat=Final_table
  for(i in ii){
    arm <- i
    tmpmat_arm <- tmpmat[tmpmat[,i] %in% c("-1", "1"), ]
    tmpmat_arm <- tmpmat_arm[complete.cases(tmpmat_arm[,i]),]
    tmpmat_arm[,i] <- droplevels(factor(tmpmat_arm[,i]))
    tmpmat_arm[,i] ='1'
    print(nlevels(tmpmat_arm[,i]))
    # tmpmat_arm=tmpmat
    # tmpmat_arm$Cohort=as.character(tmpmat_arm$Cohort)
    # tmpmat_arm[,'amp']=ifelse(tmpmat_arm[,i]==1 ,1,0)
    # tmpmat_arm[which(is.na(tmpmat_arm[,'amp'])),'amp']=0
    # tmpmat_arm[,'del']=ifelse(tmpmat_arm[,i]==-1 ,1,0)
    # tmpmat_arm[which(is.na(tmpmat_arm[,'del'])),'del']=0
    # tmpmat_arm$del=as.factor(tmpmat_arm$del)
    # tmpmat_arm$amp=as.factor(tmpmat_arm$amp)

    fit <- lm(wgii ~ factor(tmpmat_arm[,i])+Cohort,tmpmat_arm)
     pvalue = summary(fit)$coefficients[2,4]  
     coef = summary(fit)$coefficients[2,1]
     Results=rbind(Results, data.frame(arm, coef, pvalue))
   
}
Results$FDR <- p.adjust(Results$pvalue, method = "fdr")
Results$Significant2 <- ifelse(Results$FDR < 0.05 & Results$coef > 0, "positive", "insignificant")
Results$Significant2[Results$FDR < 0.05 & Results$coef < 0] <- "negative"
volcano_wgii <- ggplot(Results, aes(x = coef, y = -log10(pvalue)))+
  geom_point(aes(color = Significant2),size=0.5) +
  xlab("coefficient (linear model)") +
  ylab("-log10(p-value)") +
  scale_color_manual(values = c(tropiccolor[1], tropiccolor[2], tropiccolor[3])) +theme_cowplot(12)+
theme(axis.text.x=element_text(size=15, colour="black"), axis.text.y=element_text(size=15, colour="black"), 
      axis.title=element_text(size=15, colour="black"), legend.key = element_rect(fill="White"), 
      legend.text = element_text(size=15), legend.title = element_text(size=15)) +
  geom_text_repel(
    data = subset(Results, -log10(FDR) > 2),max.overlaps = Inf,
    aes(label = gsub("X", "", arm)),
    size = 5,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.1, "lines")
  ) +
  guides(color=guide_legend(title="significant(FDR<0.05)",override.aes = list(size=2)))

Results <- data.frame()
  for(i in ii){
    arm <- i
    tmpmat_arm <- tmpmat[tmpmat[,i] %in% c("-1", "1"), ]
    tmpmat_arm <- tmpmat_arm[complete.cases(tmpmat_arm[,i]),]
    tmpmat_arm[,i] <- droplevels(factor(tmpmat_arm[,i]))
    if(nlevels(tmpmat_arm[,i])>1){
      fit <- lm(scin ~ factor(tmpmat_arm[,i])+Cohort,tmpmat_arm)
      pvalue = summary(fit)$coefficients[2,4]  
      coef = summary(fit)$coefficients[2,1]
      Results=rbind(Results, data.frame(arm, coef, pvalue))
    }
  }
Results$FDR <- p.adjust(Results$pvalue, method = "fdr")
Results$Significant2 <- ifelse(Results$FDR < 0.05 & Results$coef > 0, "positive", "insignificant")
Results$Significant2[Results$FDR < 0.05 & Results$coef < 0] <- "negative"
volcano_scin <- ggplot(Results, aes(x = coef, y = -log10(pvalue))) +
  geom_point(aes(color = Significant2),size=0.5) +
  xlab("coefficient (linear model)") +
  ylab("-log10(p-value)") +
  scale_color_manual(values = c(tropiccolor[1], tropiccolor[2], tropiccolor[3])) +theme_cowplot(12)+
  theme(axis.text.x=element_text(size=15, colour="black"), axis.text.y=element_text(size=15, colour="black"), legend.position = 'bottom',
        axis.title=element_text(size=15, colour="black"), legend.key = element_rect(fill="White"), legend.text = element_text(size=15),
        legend.title = element_text(size=15)) +
  geom_text_repel(
    data = subset(Results, -log10(pvalue) > 2),max.overlaps = Inf,
    aes(label = gsub("X", "", arm)),
    size = 5,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.1, "lines")
  ) +
  guides(color=guide_legend(title="significant(FDR<0.05)",override.aes = list(size=2)))
legend <- get_legend(
  volcano_scin + theme(legend.box.margin = margin(-1, -1, -1, 0))
)

volcano_wgii=volcano_wgii+theme(legend.position = "none",plot.margin=margin(r=-0.1,t=0,l=0.1,unit="cm"))
volcano_scin=volcano_scin+theme(legend.position = "none",plot.margin=margin(t=0,unit="cm"))
volcano_plot=plot_grid(volcano_wgii,volcano_scin, align = "h", axis = "b",ncol=2,rel_widths = c(1,1), labels = c('WGII score','SCIN score'),label_x =0.3,label_fontface = 'plain')
pdf("plot/cin_arm_volcano.pdf",width =7,height =6)
plot_grid(volcano_plot,legend, align = "v",ncol=1,rel_heights=c(1,0.1))
dev.off()


library(ggpubr)
my_comparisons <- list(c("-1", "1"), c("-1", "0"), c("0", "1"))
df = Final_table[complete.cases(Final_table$X1q), ]# Add global p-value
box_X1q = ggboxplot(
  df,
  x = 'X1q',
  y = 'wgii',
  fill = "X1q",
  palette = c("-1" = tropiccolor[2], "0" = tropiccolor[1], "1" = tropiccolor[3]) ,
  order = c("0", "1", "-1"),
  width = 0.5,size=0.1,
  outlier.shape = NA
) + scale_x_discrete(name = "",labels = c("0" = "wt", "1" = "amp", "-1" = "del")) + theme(legend.position = "none",axis.title.y = element_blank())+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")
df = Final_table[complete.cases(Final_table$X7), ]# Add global p-value
box_X7_wgii = ggboxplot(
  df,
  x = 'X7',
  y = 'wgii',
  fill = "X7",
  palette = c("-1" = tropiccolor[2], "0" = tropiccolor[1], "1" = tropiccolor[3]) ,
  order = c("0", "1", "-1"),
  width = 0.5,size=0.1,
  outlier.shape = NA) + scale_x_discrete(name = "",labels = c("0" = "wt", "1" = "amp", "-1" = "del")) + theme(legend.position = "none",axis.title.y = element_blank())+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")

df = Final_table[complete.cases(Final_table$X19), ]# Add global p-value
box_X19 = ggboxplot(df,
  x = 'X19',
  y = 'wgii',
  fill = "X19",
  palette = c("-1" = tropiccolor[2], "0" = tropiccolor[1], "1" = tropiccolor[3]) ,
  order = c("0", "1", "-1"),
  width = 0.5,size=0.1,
  outlier.shape = NA) + scale_x_discrete(name = "",labels = c("0" = "wt", "1" = "amp", "-1" = "del")) + theme(legend.position = "none") +ylab('WGII score')+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")

df = Final_table[complete.cases(Final_table$X22q), ]# Add global p-value
box_X22q = ggboxplot(
  df,
  x = 'X22q',
  y = 'wgii',
  fill = "X22q",
  palette = c("-1" = tropiccolor[2], "0" = tropiccolor[1], "1" = tropiccolor[3]) ,
  order = c("0", "-1", "1"),
  width = 0.5,size=0.1,
  outlier.shape = NA) + scale_x_discrete(name = "",labels = c("0" = "wt", "1" = "amp", "-1" = "del")) + theme(legend.position = "none",axis.title.y = element_blank())+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")

df = Final_table[complete.cases(Final_table$X8p), ]# Add global p-value
box_X8p = ggboxplot(
  df,
  x = 'X8p',
  y = 'scin',
  fill = "X8p",
  palette = c("-1" = tropiccolor[2], "0" = tropiccolor[1], "1" = tropiccolor[3]) ,
  order = c("0", "1", "-1"),
  width = 0.5,size=0.1,
  outlier.shape = NA) + scale_x_discrete(name = "",labels = c("0" = "wt", "1" = "amp", "-1" = "del")) + theme(legend.position = "none",axis.title.y = element_blank())+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")

df = Final_table[complete.cases(Final_table$X18p), ]# Add global p-value
box_X18p = ggboxplot(
  df,
  x = 'X18p',
  y = 'scin',
  fill = "X18p",
  palette = c("-1" = tropiccolor[2], "0" = tropiccolor[1], "1" = tropiccolor[3]) ,
  order = c("0", "-1", "1"),
  width = 0.5,size=0.1,
  outlier.shape = NA) + scale_x_discrete(name = "",labels = c("0" = "wt", "1" = "amp", "-1" = "del")) + theme(legend.position = "none",axis.title.y = element_blank())+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")

df = Final_table[complete.cases(Final_table$X7), ]# Add global p-value
box_X7_scin = ggboxplot(
  df,
  x = 'X7',
  y = 'scin',
  fill = "X7",
  palette = c("-1" = tropiccolor[2], "0" = tropiccolor[1], "1" = tropiccolor[3]) ,
  order = c("0", "1", "-1"),
  width = 0.5,size=0.1,
  outlier.shape = NA) + scale_x_discrete(name = "",labels = c("0" = "wt", "1" = "amp", "-1" = "del")) + theme(legend.position = "none",axis.title.y = element_blank())+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")

df = Final_table[complete.cases(Final_table$X19p), ]
box_X19p = ggboxplot(
  df,
  x = 'X19p',
  y = 'scin',
  fill = "X19p",
  palette = c("-1" = tropiccolor[2], "0" = tropiccolor[1], "1" = tropiccolor[3]) ,
  order = c("0", "1", "-1"),
  width = 0.5,size=0.1,
  outlier.shape = NA) + scale_x_discrete(name = "",labels = c("0" = "wt", "1" = "amp", "-1" = "del")) + theme(legend.position = "none") +ylab('SCIN score')+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")
  
wgii_plot <- plot_grid(box_X19,
                      box_X7_wgii,
                      box_X1q,
                      box_X22q,
                      labels =c('Chr19','Chr7','Chr1q','Chr22q'),label_size =10,label_x = 0.6,label_fontface = 'plain',
                      ncol = 4)
scin_plot <- plot_grid(box_X19p,
                      box_X7_scin,
                      box_X8p,
                      box_X18p,
                      labels =c('Chr19','Chr7','Chr8p','Chr18p'),label_size =10,label_x = 0.6,label_fontface = 'plain',
                      ncol = 4)
wgii_plot=wgii_plot+theme(plot.margin=margin(l=-0.1,b=-0.3,t=0,unit="cm"))
scin_plot=scin_plot+theme(plot.margin=margin(l=-0.1,t=0,b=-0.3,unit="cm"))
ind_plot=plot_grid(wgii_plot,scin_plot, align = "v", axis = "l",ncol=1,rel_widths = c(1,1))
pdf("plot/cin_arm_individual.pdf",width =8.5,height =6.4)
print(ind_plot)
dev.off()
