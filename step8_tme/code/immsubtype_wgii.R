rm(list=ls())
library(ggplot2)
library(ggsignif)
library(matrixStats)
library(ggrepel)
library(cowplot)
library(data.table)
library(readxl)
library('unikn')  
library(colorspace)
tropiccolor=seecol(pal = c(rev(pal_seegruen), "white", pal_pinky),n=20)
setwd("/data/zhang/pancan_cin/step8_tme/")
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer','Primary Solid Tumor'),]
Final_table$Sample <- substr(Final_table$Sample, 1, 12)
rownames(Final_table)=Final_table$Sample
imm=read_xlsx('data/NIHMS958212-supplement-2.xlsx')
imm$Sample=imm$`TCGA Participant Barcode`
Final_table=merge(Final_table,imm,by='Sample')
Final_table$`Immune Subtype`[Final_table$`Immune Subtype`=='NA']=NA
Final_table=Final_table[complete.cases(Final_table$`Immune Subtype`),]
tmpc=unique(Final_table$`Immune Subtype`)
Results2 <- data.frame()
for(j in tmpc) {
  tmpmat = Final_table
  tmpmat$imm_bi = ifelse(Final_table$`Immune Subtype` == j, j, paste0('non_', j))
  fit <- lm(wgii ~ factor(imm_bi), tmpmat)
  pvalue = summary(fit)$coefficients[2, 4]
  coef = summary(fit)$coefficients[2, 1]
  Results2 = rbind(Results2, data.frame(j, coef, pvalue,cin_type='WGII Score'))
}

colnames(Results2)=c('imm_sub','coef','pvalue','cin_type')
pl <-ggplot(Results2, aes(y =cin_type ,  x = imm_sub))
pl=pl+geom_tile(aes(fill =coef)) + 
  geom_text(aes(label = '*',color=ifelse(Results2$pvalue<=0.05,'Sig','Not-sig'))) +scale_color_manual(name='Significance',values=c('Sig'='black','Not-sig'=NA))+
  scale_fill_gradientn(colours  = tropiccolor) +
  theme(panel.grid.major.x=element_blank(), 
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_blank(),
        legend.position = 'bottom',
        axis.title = element_blank(),axis.ticks.y =element_blank() )+ 
  theme(legend.title=element_text(size=18)) + 
  scale_x_discrete(name="") +
  scale_y_discrete(name="") 
legend <- get_legend(pl + theme(legend.box.margin = margin(-1, -1, -1, 0,unit = 'cm')))
pl=pl+theme(legend.position = 'none',plot.margin=margin(r=-0.1,t=-0.1,l=0,b=0.1,unit="cm"))
pdf('plot/imm_cor_wgii.pdf',width =5,height =2.3)
plot_grid(pl,legend, align = "v",axis = 'l',ncol=1,rel_heights=c(1,0.4))
dev.off()
