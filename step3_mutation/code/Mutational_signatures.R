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
testfn=function(data,ii,gi){
  Results <- data.frame()
  for(i in ii){
    for(j in unique(data$Cohort)){
      tmpdata=subset(data,data$Cohort==j)
      if(all(is.na(tmpdata[,i]))){
        tmpdf=data.frame(Cohort=j,Signature=i, Coef=NA, Pvalue=NA)
        Results=rbind(Results,tmpdf )
      }else{
        cor_s=cor.test(tmpdata[,gi], tmpdata[,i], method = "spearm")
        Results=rbind(Results, data.frame(Cohort=j,Signature=i, Coef=as.numeric(cor_s$estimate), Pvalue=as.numeric(cor_s$p.value)))
      }
    }
  }
  return(Results)
}
setwd("/data/zhang/pancan_cin/step3_mutation/")
mSignature_ann <-fread("data/signature_profile_project.txt")
mSignature_ann=data.frame(x=rep(1,30),y=1:30,text=unique(paste(mSignature_ann$Signature,mSignature_ann$Association,sep='_')))
mSignatureDB <-fread("data/signature_profile_sample.txt")
mSignatureDB_TCGA <- mSignatureDB[grep("TCGA", mSignatureDB$Tumor_Sample_Barcode),]
mSignatureDB_TCGA$project_code <- as.character(mSignatureDB_TCGA$project_code)
mSignatureDB_TCGA$Sample=gsub('\\.','-',as.character(mSignatureDB_TCGA$Tumor_Sample_Barcode))
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer','Primary Solid Tumor'),]
Final_table$Sample <- substr(Final_table$Sample, 1, 12)
rownames(Final_table)=Final_table$Sample
Final_table=merge(Final_table,mSignatureDB_TCGA,by='Sample')
data=dcast(Final_table, Sample+Cohort+scin+wgii~Signature, value.var="Contribution")
ii=colnames(data)[-c(1:4)]
res11=testfn(data,ii,'wgii')
res12=testfn(data,ii,'scin')
WriteXLS::WriteXLS(res11, "table/Mutational_signatures_wgii_cor_all_samples.xlsx",AdjWidth = TRUE, BoldHeaderRow = TRUE, row.names = F)
WriteXLS::WriteXLS(res12, "table/Mutational_signatures_scin_cor_all_samples.xlsx",AdjWidth = TRUE, BoldHeaderRow = TRUE, row.names = F)

Results2 <-res11
Results2$Significant <- ifelse(Results2$Pvalue <= 0.05 & Results2$Coef> 0.15, "Pos", "Not Sig")
Results2$Significant[Results2$Pvalue < 0.05 & Results2$Coef< -0.15] <- "Neg"
pl1 <-ggplot(Results2, aes(y = Cohort,  x = Signature))
pl1=pl1+geom_point(aes(stroke=0.01,size = -log(Pvalue),fill = Coef),shape=21)+
  scale_fill_gradientn(colours  = cet_pal(30, name = "d9"))+
  geom_point(mapping =aes(color=Significant),shape =22,size=5.5) +
  scale_color_manual(values=c('Pos'="#A89008" ,'Not Sig'="white",'Neg'="#3A90FE"))+
  scale_size(range = c(1, 4)) +  
  theme_cowplot()+ theme(
    axis.text.x = element_text(size =10, angle = -270, vjust=0.5),
    axis.text.y= element_text(size =10),
    legend.text = element_text(size=10),
    legend.title = element_text(size=10),
    #legend.position="none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    panel.background = element_rect(fill="white",colour="white")
  )+ggtitle("NCIN and Mut-sig Correlation")
pl1
Results2 <-res12
Results2$Significant <- ifelse(Results2$Pvalue <= 0.05 & Results2$Coef> 0.15, "Pos", "Not Sig")
Results2$Significant[Results2$Pvalue < 0.05 & Results2$Coef< -0.15] <- "Neg"

pl2 <-ggplot(Results2, aes(y = Cohort,  x = Signature))
pl2=pl2+geom_point(aes(stroke=0.01,size = -log(Pvalue),fill = Coef),shape=21)+
  scale_fill_gradientn(colours  = cet_pal(30, name = "d9"))+
  geom_point(mapping =aes(color=Significant),shape =22,size=5.5) +
  scale_color_manual(values=c('Pos'="#A89008" ,'Not Sig'="white",'Neg'="#3A90FE"))+
  scale_size(range = c(1, 4)) +  
  theme_cowplot()+ theme(
    axis.text.x = element_text(size =10, angle = -270, vjust=0.5),
    axis.text.y= element_text(size =10),
    legend.text = element_text(size=10),
    legend.title = element_text(size=10),
    #legend.position="none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    panel.background = element_rect(fill="white",colour="white")
  )+ggtitle("SCIN and Mut-sig Correlation")
pl2
legend=ggplot(mSignature_ann, aes(x, y, label = text))+geom_text(size=2)+ theme_cowplot()+ theme_void()+theme(plot.margin=margin(l=-0.5,unit="cm"))

pl1=pl1+theme(plot.margin=margin(r=-0.5,unit="cm"))
pl2=pl2+theme(plot.margin=margin(r=-0.5,unit="cm"))
prow=plot_grid(pl1,legend,ncol=2, align = "h", axis = "b",rel_widths  = c(3,0.5))
pdf('plot/NCIN_and_Mut-sig_Correlation_all_sample.pdf',width =12,height = 10)
print(prow)
dev.off()
prow=plot_grid(pl2,legend,ncol=2, align = "h", axis = "b",rel_widths  = c(3,0.5))
pdf('plot/SCIN_and_Mut-sig_Correlation_all_sample.pdf',width =12,height = 10)
print(prow)
dev.off()
