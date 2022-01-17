rm(list=ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library('unikn')  
library(RColorBrewer)
#library(data.table)
library(sna)
library("ggbeeswarm")
library(cowplot)
library(ggpmisc)
library(ggpubr)
library(cetcolor)
testfn=function(data,ii,gi){
  Results <- data.frame()
  for(i in ii){
    for(j in unique(data$Cohort)){
      tmpdata=subset(data,data$Cohort==j)
      if(all(is.na(tmpdata[,i]))){
        tmpdf=data.frame(Cohort=j,Genome_instability=i, correlation=NA, pvalue=NA)
        Results=rbind(Results,tmpdf )
      }else{
        cor_s=cor.test(tmpdata[,gi], tmpdata[,i], method = "spearm")
        Results=rbind(Results, data.frame(Cohort=j,Genome_instability=i, correlation=as.numeric(cor_s$estimate), pvalue=as.numeric(cor_s$p.value)))
      }
    }
  }
  return(Results)
}
gcolor=yarrr::piratepal("appletv") 
setwd("/data/zhang/pancan_cin/step1_calculate_cin/")
load('data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer','Primary Solid Tumor'),]
Final_table$Sample <- substr(Final_table$Sample, 1, 15)
rownames(Final_table)=Final_table$Sample
Final_table$Genome_doublings <- as.factor(Final_table$`Genome doublings`)
#Load Table S2 from data from Taylor et al, Cancer Cell 2018 (https://www.sciencedirect.com/science/article/pii/S1535610818301119)
astable <- readxl::read_xlsx("../step2_genome_instability/data/1-s2.0-S1535610818301119-mmc2.xlsx",skip = 1)
astable=as.data.frame(astable)
#rownames(astable) <- as.character(astable$Sample)
ii=c("AneuploidyScore(AS)", "AS_del", "AS_amp", 
     "Genome_doublings", "Leuk", "Purity", "Stroma", "Stroma_notLeukocyte", 
     "Stroma_notLeukocyte_Floor", "SilentMutationspeMb", "Non-silentMutationsperMb", 
     "1p", "1q", "2p", "2q", "3p", "3q", "4p", "4q", "5p", "5q", "6p", 
     "6q", "7p", "7q", "8p", "8q", "9p", "9q", "10p", "10q", "11p", 
     "11q", "12p", "12q", "13q", "14q", "15q", "16p", "16q", "17p", 
     "17q", "18p", "18q", "19p", "19q", "20p", "20q", "21q", "22q", 
     "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", 
     "16", "17", "18", "19", "20")
astable[,ii]=apply(astable[,ii],2,as.numeric)

data=merge(Final_table,astable,by='Sample')
#data$MutationBurden=data$SilentMutationspeMb+data$`Non-silentMutationsperMb`
ii=c("AneuploidyScore(AS)")
res0=testfn(data,ii,'wgii')
res02=testfn(data,ii,'scin')
msitable <- readxl::read_xlsx("../step2_genome_instability/data/1-s2.0-S009286741830237X-mmc1.xlsx",skip=8 ,sheet = 'Table S5')
msitable$Sample=msitable$TCGA_Barcode
msitable=as.data.frame(msitable)
ii=c( "Tota_Mutation_cnt", "MSIscore", 
      "POLE", "MLH1", "MLH3", "MGMT", "MSH6", "MSH3", "MSH2", "PMS1", 
      "PMS2")
msitable[,ii]=apply(msitable[,ii],2,as.numeric)
data=Final_table
data$Sample=substr(data$Sample,1,12)
library(data.table)
data=merge(data,msitable,by='Sample')
# data=subset(data,data$Cohort %in% c("LUSC", "LUAD", "UCEC", "COAD", "STAD"))
# data=subset(data,data$Cohort %in% c( "UCEC", "COAD", "STAD"))
data$Cohort=as.character(data$Cohort)
#data=subset(data,!(data$Patient %in% msitable$Sample))
#data$Cohort[which(data$Cohort %in% c("LUSC", "LUAD"))]='LUNG'
ii=c("MSIscore")
# res3=testfn(data,ii,'wgii')
# res32=testfn(data,ii,'scin')


require(tidyr)

formula =  y ~ x
Results2=data %>% tidyr::gather(key=CIN_type, value=CIN, c(scin,wgii))
Results2$Cohort=as.character(Results2$Cohort)
Results2$CIN_type=factor(Results2$CIN_type,labels = c('WGII score','SCIN score'))
Results2=Results2[complete.cases(Results2$Genome_doublings),]
pl <-ggplot(Results2, aes(x =MSIscore,  y = CIN))+facet_wrap(.~CIN_type,scales = 'free_y')
pl=pl+geom_point(aes(color=Genome_doublings),size=0.3)
pl=pl +ylab("CIN score") +
  xlab("MSI score") +
  scale_colour_manual("WGD", labels = c("0", "1","2"), values=unname(gcolor[c(1,5,6)]))+
  #stat_density2d(aes(fill = stat(level),alpha=0.5), geom = "polygon") + scale_fill_gradient(low="#1396DBFF", high="white", name="Distribution")+
  geom_smooth(method = "loess", formula = formula, se = F,span = 0.9, colour=unname(gcolor[1]),lwd=0.5) +
  stat_cor( aes(label = paste(..p.label.., sep = "~`,`~")),label.x = 0.5, label.y.npc=0.7,method = "spearman",label.sep = "\n", size = 4.5)+
  stat_cor( aes(label = paste(..rr.label.., sep = "~`,`~")),label.x = 0.5, label.y.npc=0.8,method = "spearman",label.sep = "\n", size = 4.5)+
  theme_cowplot(12)+
  theme(strip.background = element_rect(fill = 'white'),
        panel.background = element_rect(fill="white",colour="white"), 
        legend.position = 'bottom',axis.text.y=element_text(size=15, colour="black"), 
        axis.text.x=element_text(size=15, colour="black"),legend.box.margin = margin(-3, 0, -5, 0),
        axis.title=element_text(size=15, colour="black"),
        legend.key=element_rect(fill=NA), legend.text = element_text(size=15, colour="black"), 
        legend.title = element_text(size=15, colour="black"))+guides(color = guide_legend(override.aes = list(size = 3)))
pl

pdf('plot/pan_can_gi_msi.pdf',width =5,height = 4)
print(pl)
dev.off()
library(ggpubr)
data=Final_table
data$Sample=substr(data$Sample,1,12)
data=join(data,msitable,by='Sample',type='left')
data$Cohort=as.character(data$Cohort)
data$MSI=ifelse(is.na(data$MSIscore),'NMSI','MSI')
my_comparisons <- list(c('NMSI','MSI'))
pl <- ggboxplot(data, x = "MSI", y = "wgii",fill = "MSI",width=0.3, size=0.1,outlier.colour = NA)+ylab('WGII score')+
  scale_fill_manual("MSI", labels = c('NMSI','MSI'),values=colorRampPalette(brewer.pal(12, "Set3"), alpha=TRUE)(17)[10:11])+
  stat_compare_means(aes(label = paste0("p = ", ..p.format..)),label.y.npc ='top',label.x.npc = 'left')+
  theme(legend.position = 'none',axis.text.y=element_text(size=15, colour="black"), 
        axis.text.x=element_text(size=15, colour="black"),
        axis.title.x=element_blank())+ylim(-0.1,1.2)
pl1 <- ggboxplot(data, x = "MSI", y = "scin",fill = "MSI",width=0.3, size=0.1, outlier.colour = NA)+scale_y_log10()+ylab('SCIN score (log scaled)')+
  scale_fill_manual("MSI", labels = c('NMSI','MSI'),values=colorRampPalette(brewer.pal(12, "Set3"), alpha=TRUE)(17)[10:11])+
  stat_compare_means(aes(label = paste0("p = ", ..p.format..)),label.y.npc ='top',label.x.npc = 'left')+
  theme(legend.position = 'none',axis.text.y=element_text(size=15, colour="black"), 
        axis.text.x=element_text(size=15, colour="black"),
        axis.title.x=element_blank())
plot_grid(pl,pl1)
pdf('plot/pan_can_cin_msi.pdf',width =4,height = 4)
plot_grid(pl,pl1)
dev.off()





immtable=readxl::read_xlsx('../step8_tme/data/1-s2.0-S1074761318301213-mmc2.xlsx')
immtable$Sample=immtable$`TCGA Participant Barcode`
immtable=as.data.frame(immtable)
ii=c("Leukocyte Fraction", "Stromal Fraction", "Intratumor Heterogeneity", 
     "TIL Regional Fraction", "Proliferation", "Wound Healing", "Macrophage Regulation", 
     "Lymphocyte Infiltration Signature Score", "IFN-gamma Response", 
     "TGF-beta Response", "SNV Neoantigens", "Indel Neoantigens", 
     "Silent Mutation Rate", "Nonsilent Mutation Rate", "Number of Segments", 
     "Fraction Altered", "Aneuploidy Score", "Homologous Recombination Defects", 
     "BCR Evenness", "BCR Shannon", "BCR Richness", "TCR Shannon", 
     "TCR Richness", "TCR Evenness", "CTA Score", "Th1 Cells", "Th2 Cells", 
     "Th17 Cells", "OS", "OS Time", "PFI", "PFI Time", "B Cells Memory", 
     "B Cells Naive", "Dendritic Cells Activated", "Dendritic Cells Resting", 
     "Eosinophils...41", "Macrophages M0", "Macrophages M1", "Macrophages M2", 
     "Mast Cells Activated", "Mast Cells Resting", "Monocytes", "Neutrophils...48", 
     "NK Cells Activated", "NK Cells Resting", "Plasma Cells", "T Cells CD4 Memory Activated", 
     "T Cells CD4 Memory Resting", "T Cells CD4 Naive", "T Cells CD8", 
     "T Cells Follicular Helper", "T Cells gamma delta", "T Cells Regulatory Tregs", 
     "Lymphocytes", "Neutrophils...60", "Eosinophils...61", "Mast Cells", 
     "Dendritic Cells", "Macrophages")
immtable[,ii]=apply(immtable[,ii],2,as.numeric)
data=Final_table
data$Sample=substr(data$Sample,1,12)
data=merge(data,immtable,by='Sample')
ii=c("Intratumor Heterogeneity","Proliferation", "Homologous Recombination Defects","Silent Mutation Rate", "Nonsilent Mutation Rate")
res1=testfn(data,c('scin',ii),'wgii')
res12=testfn(data,c(ii),'scin')

rna_table <-fread("../step8_tme/data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv")
ii = unlist(lapply(rna_table$gene_id, function(x) {
  strsplit(as.character(x), '\\|')[[1]][1]
}))
iidx = which(!duplicated(ii))
ii = ii[iidx]
rna_table = rna_table[iidx, ]
rownames(rna_table) = ii
rna_table = rna_table[which(rownames(rna_table) != "?"), ]
rna_table$gene_id = NULL
rna_table = t(as.data.frame(rna_table))
colnames(rna_table) = ii[-1]
rna_table = data.frame(rna_table)
rna_table$Sample =  substr(rownames(rna_table), 1, 15)
rna_table=rna_table[!duplicated(rna_table$Sample),]
rownames(rna_table)=rna_table$Sample
rna_table=rna_table[,c('Sample',grep('KI67',colnames(rna_table),value = T))]
data=Final_table

data=merge(data,rna_table,by='Sample')
data$MKI67=log2(data$MKI67+1)
ii=c("MKI67")
data$log_wgii=log2(data$wgii+1)
data$log_scin=log2(data$scin+1)
res2=testfn(data,ii,'log_wgii')
res22=testfn(data,ii,'log_scin')



Results=rbind(res0,res1,res2)

dataf=reshape2::dcast(Results, Genome_instability~Cohort, value.var="correlation")
rownames(dataf)=dataf$Genome_instability
dataf=dataf[,-1]
ind=hclust(dist(t(dataf)))$order
ind=colnames(dataf)[ind]
indr=hclust(dist(dataf))$order
indr=rownames(dataf)[indr]
Results$Cohort=factor(Results$Cohort,levels = ind)
Results$Genome_instability=factor(Results$Genome_instability,levels=indr)
Results2 <-Results
Results2$Significant <- ifelse(Results2$pvalue <= 0.05 & Results2$correlation> 0.3 , "Pos", "Not Sig")
Results2$Significant[Results2$pvalue < 0.05 & Results2$correlation< -0.3] <- "Neg"
ii=c("scin","Aneuploidyscore(AS)","Homologous Recombination Defects","Silent Mutation Rate", "Nonsilent Mutation Rate","Intratumor Heterogeneity","Proliferation","MKI67")
Results2$Genome_instability=factor(as.character(Results2$Genome_instability),levels=ii)
pl1=ggplot(Results2, aes(x = Cohort,  y= Genome_instability))
pl1=pl1+geom_point(aes(size = -log(pvalue),fill = correlation),shape=21)+
      scale_fill_gradientn(colours  = seecol(pal = c(rev(pal_seegruen), "white", pal_pinky),n=20))+
  #geom_point(mapping =aes(color=Significant),shape =22,size=5.5) +
  #scale_color_manual(values=c('Pos'="#A89008" ,'Not Sig'="white",'Neg'="#3A90FE"))+
  scale_size(range = c(1, 7)) +  
  theme_cowplot()+ theme(
    axis.text.x = element_text(size=15, angle =-270, vjust=0.5),
    axis.text.y= element_text(size=15),
    legend.text = element_text(size=15),
    legend.title = element_text(size=15),
    #legend.position="none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    panel.background = element_rect(fill="white",colour="white")
  )
pl1=pl1+scale_y_discrete(breaks=ii,labels=c("SCIN","AS",'HRD',"SMR","NSMR", "ITH",'PROLIF','MKI67'))
pl1=pl1+guides(gradientn=guide_legend(title="correlation"),size=guide_legend(title = "-log(pvalue)"))
Results=rbind(res02,res12,res22)

library(reshape2)

dataf=reshape2::dcast(Results, Genome_instability~Cohort, value.var="correlation")
rownames(dataf)=dataf$Genome_instability
dataf=dataf[,-1]
ind=hclust(dist(t(dataf)))$order
ind=colnames(dataf)[ind]
indr=hclust(dist(dataf))$order
indr=rownames(dataf)[indr]
Results$Cohort=factor(Results$Cohort,levels = ind)
Results$Genome_instability=factor(Results$Genome_instability,levels=indr)
Results2 <-Results
Results2$Significant <- ifelse(Results2$pvalue <= 0.05 & Results2$correlation> 0.3 , "Pos", "Not Sig")
Results2$Significant[Results2$pvalue < 0.05 & Results2$correlation< -0.3] <- "Neg"
ii=c("Aneuploidyscore(AS)","Homologous Recombination Defects","Silent Mutation Rate", "Nonsilent Mutation Rate","Intratumor Heterogeneity","Proliferation","MKI67")
Results2$Genome_instability=factor(as.character(Results2$Genome_instability),levels=ii)
pl2 <-ggplot(Results2, aes(x = Cohort,  y = Genome_instability))
pl2=pl2+geom_point(aes(size = -log(pvalue),fill =correlation),shape=21)+
  scale_fill_gradientn(colours  = seecol(pal = c(rev(pal_seegruen), "white", pal_pinky),n=20))+
  #geom_point(mapping =aes(color=Significant),shape =22,size=5.5) +
  #scale_color_manual(values=c('Pos'="#A89008" ,'Not Sig'="white",'Neg'="#3A90FE"))+
  scale_size(range = c(1, 7)) +  
  theme_cowplot()+ theme(
    axis.text.x = element_text(size=15, angle =-270, vjust=0.5),
    axis.text.y= element_text(size=15),
    legend.text = element_text(size=15),
    legend.title = element_text(size=15),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    panel.background = element_rect(fill="white",colour="white")
  )
pl2=pl2+scale_y_discrete(breaks=ii,labels=c("AS",'HRD',"SMR","NSMR", "ITH",'PROLIF','MKI67'))
pl2
#pl1=pl1+theme(plot.margin=margin(r=0.8,t=0.6,unit="cm"))
#pl2=pl2+theme(plot.margin=margin(r=0.8,t=0.6,unit="cm"))
legend <- get_legend(
  pl2 + theme(legend.box.margin = margin(0, 0, 0, 0))
)
#pl2=pl2+theme(legend.position = "none")
prow=plot_grid(pl1,pl2,nrow=2, align = "v", axis = "l",rel_heights  = c(1.2,1),labels = c("WGII score vs GI","SCIN score vs GI"),label_fontface = 'plain')
prow
pdf('plot/pan_can_gi_wgii.pdf',width = 13,height = 4)
print(pl1)
dev.off()
pdf('plot/pan_can_gi_scin.pdf',width = 13,height = 4)
print(pl2)
dev.off()


# 
# 
# 
# data=Final_table
# 
# formula <- y~ x
# Final_table <- merge(Final_table, Prolif_score_2,by='Sample')
# scin_vs_wgii <- ggplot(data, aes(x=wgii, y=scin, colour=Genome_doublings, fill=Genome_doublings))+facet_wrap(.~Cohort,scales = 'free_y') + 
#   xlab("wgii score") + 
#   ylab(paste("scin score")) + 
#   geom_point(size=0.6, alpha=0.6,inherit.aes = T) + 
#   scale_colour_manual("Genome doublings", labels = c("0", "1","2"), values=unname(gcolor[c(1,5,6)]))+
#   scale_fill_manual("Genome doublings", labels = c("0", "1","2"), values=unname(gcolor[c(1,5,6)]))+
#   scale_alpha(range = c(0.1,0.3))+
#   #stat_density2d(aes(alpha=..level..),geom='polygon',bins=20) +
#   geom_smooth(method = "lm", formula = formula, se = F) +
#   stat_cor(aes(label = paste(..r.label.., sep = "~`")),label.y.npc = 'top',inherit.aes = T) +
#   ggtitle('correlation of wgii and scin')+
#   theme(strip.background = element_rect(fill = 'white'),panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), legend.position = 'bottom',axis.text.y=element_text(size=8, colour="black"), axis.text.x=element_text(size=8, colour="black"), axis.title=element_text(size=18, colour="black"), plot.title = element_text(size=18, hjust = 0.5, colour="black"), legend.key=element_rect(fill=NA), legend.text = element_text(size=8, colour="black"), legend.title = element_text(size=16, colour="black"), axis.title.y=element_text(margin=margin(0,10,0,0)), axis.title.x=element_text(margin=margin(10,0,0,0)))
# 
# aneuploidy_vs_wgii <- ggplot(data, aes(x=wgii, y=Aneuploidyscore.AS.))+facet_wrap(.~Cohort,scales = 'free_y') + 
#   xlab("wgii score") + 
#   ylab(paste("aneuploidy score")) + 
#   geom_point(size=0.6, inherit.aes = T, colour="#1396DBFF", fill="#1396DBFF") + 
#   scale_alpha(range = c(0.1,0.3))+
#   stat_density2d(aes(fill = stat(level),alpha=0.5), geom = "polygon") + scale_fill_gradient(low="#1396DBFF", high="white", name="Distribution")+
#   geom_smooth(method = "lm", formula = formula, se = F, colour="#1396DBFF", fill="#1396DBFF") +
#   stat_cor( aes(label = paste(..r.label.., sep = "~`")),label.x.npc = 'left',inherit.aes = T) +
#   ggtitle('correlation of proliferation and wgii')+
#   theme(strip.background = element_rect(fill = 'white'),panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), legend.position = 'bottom',axis.text.y=element_text(size=8, colour="black"), axis.text.x=element_text(size=8, colour="black"), axis.title=element_text(size=18, colour="black"), plot.title = element_text(size=18, hjust = 0.5, colour="black"), legend.key=element_rect(fill=NA), legend.text = element_text(size=8, colour="black"), legend.title = element_text(size=16, colour="black"), axis.title.y=element_text(margin=margin(0,10,0,0)), axis.title.x=element_text(margin=margin(10,0,0,0)))
# 
# 
# aneuploidy_vs_scin <-  ggplot(data, aes(x=scin, y=Aneuploidyscore.AS.))+facet_wrap(.~Cohort,scales = 'free') + 
#   xlab("scin score") + 
#   ylab(paste("aneuploidy score")) + 
#   geom_point(size=0.6, inherit.aes = T, colour="#1396DBFF", fill="#1396DBFF") + 
#   scale_alpha(range = c(0.1,0.3))+
#   stat_density2d(aes(fill = stat(level),alpha=0.5), geom = "polygon") + scale_fill_gradient(low="#1396DBFF", high="white", name="Distribution")+
#   geom_smooth(method = "lm", formula = formula, se = F, colour="#1396DBFF", fill="#1396DBFF") +
#   stat_cor( aes(label = paste(..r.label.., sep = "~`")),label.x.npc = 'left',inherit.aes = T) +
#   ggtitle('correlation of proliferation and wgii')+
#   theme(strip.background = element_rect(fill = 'white'),panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), legend.position = 'bottom',axis.text.y=element_text(size=8, colour="black"), axis.text.x=element_text(size=8, colour="black"), axis.title=element_text(size=18, colour="black"), plot.title = element_text(size=18, hjust = 0.5, colour="black"), legend.key=element_rect(fill=NA), legend.text = element_text(size=8, colour="black"), legend.title = element_text(size=16, colour="black"), axis.title.y=element_text(margin=margin(0,10,0,0)), axis.title.x=element_text(margin=margin(10,0,0,0)))
# 
# 
# proli_vs_wgii <- ggplot(data, aes(x=wgii, y=rates))+facet_wrap(.~Cohort,scales = 'free_y') + 
#   xlab("wgii score") + 
#   ylab(paste("proliferation rates")) + 
#   geom_point(size=0.6, inherit.aes = T, colour="#1396DBFF", fill="#1396DBFF") + 
#   scale_alpha(range = c(0.1,0.3))+
#   stat_density2d(aes(fill = stat(level),alpha=0.5), geom = "polygon") + scale_fill_gradient(low="#1396DBFF", high="white", name="Distribution")+
#   geom_smooth(method = "lm", formula = formula, se = F, colour="#1396DBFF", fill="#1396DBFF") +
#   stat_cor( aes(label = paste(..r.label.., sep = "~`")),label.x.npc = 'left',inherit.aes = T) +
#   ggtitle('correlation of proliferation and wgii')+
#   theme(strip.background = element_rect(fill = 'white'),panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), legend.position = 'bottom',axis.text.y=element_text(size=8, colour="black"), axis.text.x=element_text(size=8, colour="black"), axis.title=element_text(size=18, colour="black"), plot.title = element_text(size=18, hjust = 0.5, colour="black"), legend.key=element_rect(fill=NA), legend.text = element_text(size=8, colour="black"), legend.title = element_text(size=16, colour="black"), axis.title.y=element_text(margin=margin(0,10,0,0)), axis.title.x=element_text(margin=margin(10,0,0,0)))
# 
# proli_vs_scin<- ggplot(data, aes(x=scin, y=rates))+facet_wrap(.~Cohort,scales = 'free') + 
#   xlab("scin score") + 
#   ylab(paste("proliferation rates")) + 
#   geom_point(size=0.6, inherit.aes = T, colour="#1396DBFF", fill="#1396DBFF") + 
#   scale_alpha(range = c(0.1,0.3))+
#   stat_density2d(aes(fill = stat(level),alpha=0.5), geom = "polygon") + scale_fill_gradient(low="#1396DBFF", high="white", name="Distribution")+
#   geom_smooth(method = "lm", formula = formula, se = F, colour="#1396DBFF", fill="#1396DBFF") +
#   stat_cor( aes(label = paste(..r.label.., sep = "~`")),label.x.npc = 'left',inherit.aes = T) +
#   ggtitle('correlation of proliferation and wgii')+
#   theme(strip.background = element_rect(fill = 'white'),panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), legend.position = 'bottom',axis.text.y=element_text(size=8, colour="black"), axis.text.x=element_text(size=8, colour="black"), axis.title=element_text(size=18, colour="black"), plot.title = element_text(size=18, hjust = 0.5, colour="black"), legend.key=element_rect(fill=NA), legend.text = element_text(size=8, colour="black"), legend.title = element_text(size=16, colour="black"), axis.title.y=element_text(margin=margin(0,10,0,0)), axis.title.x=element_text(margin=margin(10,0,0,0)))
# 
# pdf('plot/scin_vs_wgii.pdf', width = 12,height=10)
# print(scin_vs_wgii)
# dev.off()
# pdf('plot/aneuploidy_vs_wgii.pdf', width = 12,height=10)
# print(aneuploidy_vs_wgii)
# dev.off()
# pdf('plot/aneuploidy_vs_scin.pdf', width = 12,height=10)
# print(aneuploidy_vs_scin)
# dev.off()
# pdf('plot/proli_vs_wgii.pdf', width = 12,height=10)
# print(proli_vs_wgii)
# dev.off()
# pdf('plot/proli_vs__scin.pdf', width = 12,height=10)
# print(proli_vs_scin)
# dev.off()
# 
# rna_table <-fread("../step8_tme/data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv")
# ii = unlist(lapply(rna_table$gene_id, function(x) {
#   strsplit(as.character(x), '\\|')[[1]][1]
# }))
# iidx = which(!duplicated(ii))
# ii = ii[iidx]
# rna_table = rna_table[iidx, ]
# rownames(rna_table) = ii
# rna_table = rna_table[which(rownames(rna_table) != "?"), ]
# rna_table$gene_id = NULL
# rna_table = t(as.data.frame(rna_table))
# colnames(rna_table) = ii[-1]
# rna_table = data.frame(rna_table)
# rna_table$Sample =  substr(rownames(rna_table), 1, 15)
# rna_table=rna_table[!duplicated(rna_table$Sample),]
# rownames(rna_table)=rna_table$Sample
# rna_table=rna_table[,c('Sample',grep('KI67',colnames(rna_table),value = T))]
# Final_table=merge(Final_table,rna_table,by='Sample')
# data$log_MKI67=log2(data$MKI67+1)
# data$log_wgii=log2(data$wgii+1)
# data$log_scin=log2(data$scin+1)
# ki67_vs_wgii <- ggplot(data, aes(x=log_wgii, y=log_MKI67))+facet_wrap(.~Cohort,scales = 'free_y') + 
#   xlab("Log2 wgii score") +#scale_y_log10()+scale_x_log10()+
#   ylab(paste("Log2 MKI67 gene expression")) + 
#   geom_point(size=0.6, inherit.aes = T, colour="#1396DBFF", fill="#1396DBFF") + 
#   scale_alpha(range = c(0.1,0.3))+
#   stat_density2d(aes(fill = stat(level),alpha=0.5), geom = "polygon") + scale_fill_gradient(low="#1396DBFF", high="white", name="Distribution")+
#   geom_smooth(method = "lm", formula = formula, se = F, colour="#1396DBFF", fill="#1396DBFF") +
#   stat_cor( aes(label = paste(..r.label.., sep = "~`")),label.x.npc = 'left',inherit.aes = T) +
#   ggtitle('correlation of proliferation and wgii')+
#   theme(strip.background = element_rect(fill = 'white'),panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), legend.position = 'bottom',axis.text.y=element_text(size=8, colour="black"), axis.text.x=element_text(size=8, colour="black"), axis.title=element_text(size=18, colour="black"), plot.title = element_text(size=18, hjust = 0.5, colour="black"), legend.key=element_rect(fill=NA), legend.text = element_text(size=8, colour="black"), legend.title = element_text(size=16, colour="black"), axis.title.y=element_text(margin=margin(0,10,0,0)), axis.title.x=element_text(margin=margin(10,0,0,0)))
# 
# ki67_vs_scin<- ggplot(data, aes(x=log_scin, y=log_MKI67))+facet_wrap(.~Cohort,scales = 'free_y') + 
#   xlab("Log2 scin score") + 
#   ylab(paste("Log2 MKI67 gene expression")) + 
#   geom_point(size=0.6, inherit.aes = T, colour="#1396DBFF", fill="#1396DBFF") + 
#   scale_alpha(range = c(0.1,0.3))+
#   stat_density2d(aes(fill = stat(level),alpha=0.5), geom = "polygon") + scale_fill_gradient(low="#1396DBFF", high="white", name="Distribution")+
#   geom_smooth(method = "lm", formula = formula, se = F, colour="#1396DBFF", fill="#1396DBFF") +
#   stat_cor( aes(label = paste(..r.label.., sep = "~`")),label.x.npc = 'left',inherit.aes = T) +
#   ggtitle('correlation of proliferation and scin')+
#   theme(strip.background = element_rect(fill = 'white'),panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), legend.position = 'bottom',axis.text.y=element_text(size=8, colour="black"), axis.text.x=element_text(size=8, colour="black"), axis.title=element_text(size=18, colour="black"), plot.title = element_text(size=18, hjust = 0.5, colour="black"), legend.key=element_rect(fill=NA), legend.text = element_text(size=8, colour="black"), legend.title = element_text(size=16, colour="black"), axis.title.y=element_text(margin=margin(0,10,0,0)), axis.title.x=element_text(margin=margin(10,0,0,0)))
# pdf('plot/ki67_vs_wgii.pdf', width = 12,height=10)
# print(ki67_vs_wgii)
# dev.off()
# pdf('plot/ki67_vs__scin.pdf', width = 12,height=10)
# print(ki67_vs_scin)
# dev.off()