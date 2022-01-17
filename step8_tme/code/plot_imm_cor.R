rm(list=ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library(data.table)
library(sna)
library("ggbeeswarm")
library(cowplot)
library(ggpmisc)
library(ggpubr)
library('unikn')  
library(colorspace)
tropiccolor=seecol(pal = c(rev(pal_seegruen), "white", pal_pinky),n=20)
setwd("/data/zhang/pancan_cin/step8_tme/")
Results2=fread('table/immune_corr_results_wgii.txt')
colnames(Results2)=c("Cohort", "Cell_type", "Coef", "Pvalue", "Pearson_coef", "Pearson_pvalue", "FDR", "Pearson_FDR","text", "Significant")
Results2$Significant <- ifelse(Results2$FDR <= 0.05 & Results2$Coef> 0.3 , "Pos", "Not Sig")
Results2$Significant[Results2$FDR < 0.05 & Results2$Coef< -0.3] <- "Neg"
library(reshape2)

data_wide= dcast(Results2,Cohort~Cell_type, value.var="Significant")
selc=apply(data_wide[,-1],2, function(x){length(which(x!='Not Sig'))})
selc=selc[which(selc>=6)]
print(length(selc))
Results2=Results2[which(Results2$Cell_type %in% names(selc)),]
pl1 <-ggplot(Results2, aes(x = Cohort,  y = Cell_type))
pl1=pl1+geom_point(aes(stroke=0.00001,size = -log(Pvalue),fill = Coef),shape=21)+
  scale_fill_gradientn(colours  = tropiccolor)+
  scale_size(range = c(1, 7)) +  
  theme_cowplot()+ theme(
    axis.text.x = element_text(size=13, angle =-270, vjust=0.5),
    axis.text.y= element_text(size=13),
    legend.text = element_text(size=13),
    legend.title = element_text(size=13),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    panel.background = element_rect(fill="white",colour="white")
  )
Results2=fread('table/immune_corr_results_scin.txt')
colnames(Results2)=c("Cohort", "Cell_type", "Coef", "Pvalue", "Pearson_coef", "Pearson_pvalue", "FDR", "Pearson_FDR","text", "Significant")
Results2$Significant <- ifelse(Results2$FDR <= 0.05 & Results2$Coef> 0.3 , "Pos", "Not Sig")
Results2$Significant[Results2$FDR < 0.05 & Results2$Coef< -0.3] <- "Neg"
data_wide= dcast(Results2,Cohort~Cell_type, value.var="Significant")
selc=apply(data_wide[,-1],2, function(x){length(which(x!='Not Sig'))})
selc=selc[which(selc>=11)]
print(length(selc))
Results2=Results2[which(Results2$Cell_type %in% names(selc)),]
pl2 <-ggplot(Results2, aes(x= Cohort,  y = Cell_type))
pl2=pl2+geom_point(aes(stroke=0,size = -log(Pvalue),fill = Coef),shape=21)+
  scale_fill_gradientn(colours  = tropiccolor,breaks=c(-0.5,0,0.5))+
  scale_size(range = c(1, 7)) +  
  theme_cowplot()+ theme(
    axis.text.x = element_text(size=13, angle =-270, vjust=0.5),
    axis.text.y= element_text(size=13),
    legend.text = element_text(size=13),
    legend.title = element_text(size=13),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),legend.position = 'bottom',
    panel.background = element_rect(fill="white",colour="white")
  )
legend <- get_legend(pl2 + theme(legend.box.margin = margin(0.1, 0.1, 0.1, 0,unit="cm")))
pl1=pl1+theme(legend.position = "none",plot.margin=margin(r=0.01,t=0.5,l=0,unit="cm"))
pl2=pl2+theme(legend.position = "none",plot.margin=margin(t=0.5,unit="cm"))
heatmap_plot=plot_grid(pl1,pl2, align = "h", axis = "b",ncol=2,rel_widths = c(1,1), labels = c('WGII score','SCIN score'),label_x = c(0.4,0),label_fontface = 'plain')
pdf("plot/imm_cor_cin_heatmap.pdf",width =15,height =13)
plot_grid(heatmap_plot,legend, align = "v",ncol=1,rel_heights=c(10,0.4),axis = 'l')
dev.off()

