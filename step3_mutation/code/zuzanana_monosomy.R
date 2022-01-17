rm(list = ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library(data.table)
library(RColorBrewer)
setwd("/data/zhang/pancan_cin/step3_mutation/")
fns = c(
  "DDR gene mutations.csv",
  "DDR gene alterations ONCOPRINT.csv",
  "DDR gene alterations.csv",
  "DDR footprints.csv",
  "DDR epigenetic silencing.csv",
  "DDR deep deletions.csv"
)
pathdir = "../step8_tme/data/TCGA_DDR_Data_Resources_xls"
datalist = list()
for (i in fns) {
  if (i != 'DDR footprints.csv') {
    tmp = fread(file.path(pathdir, i))
    tmp = data.frame(tmp, stringsAsFactors = F)
    coltmp = colnames(tmp)[-c(1:2)]
    rowtmp = as.character(tmp[, 2][-c(1:2)])
    mattmp = tmp[-c(1:2), -c(1:2)]
    rownames(mattmp) = rowtmp
    colnames(mattmp) = coltmp
    datalist[[i]] = mattmp
  }
}


i = "DDR footprints.csv"
tmp = fread(file.path(pathdir, i))
tmp = data.frame(tmp, stringsAsFactors = F)
coltmp = c('Cohort', colnames(tmp)[-c(1:4)])
rowtmp = as.character(tmp[, 2][-c(1:4)])
mattmp = tmp[-c(1:4), -c(1:2, 4)]
rownames(mattmp) = rowtmp
colnames(mattmp) = coltmp
Final_table = mattmp[,c("LOH_n_seg", "LOH_frac_altered", "purity", "ploidy", "genome_doublings","tp53_score")]
Final_table=data.frame(apply(Final_table, 2, as.numeric))
Final_table$Cohort=mattmp$Cohort
rownames(Final_table)=gsub('-','.',rownames(mattmp))
Final_table=Final_table[which(!is.na(Final_table$ploidy)),]
table <- readxl::read_excel("../step2_genome_instability/data/1-s2.0-S1535610818301119-mmc2.xlsx",skip = 1)
table=as.data.frame(table)
rownames(table) <- gsub('-','.',as.character(table$Sample))
table$Sample= gsub('-','.',as.character(table$Sample))
sels=intersect(rownames(table),rownames(Final_table))
data=Final_table[sels,]
data$Sample= gsub('-','.',rownames(data))
Final_table <- merge(data, table, by='Sample')
rownames(Final_table)=Final_table$Sample
Final_table$ploidy_type=cut(Final_table$ploidy, quantile(Final_table$ploidy+runif(length(Final_table$ploidy), -0.006,0.001), c(0,0.06,0.6,1), na.rm=TRUE), right=F,labels=c("monosomy","disomy","polysomy"))
Final_table$ploidy_type[which(is.na(Final_table$ploidy_type))]='polysomy'
sels=intersect(colnames(datalist$`DDR gene mutations.csv`),rownames(Final_table))
Final_table$tp53_mut=as.character(datalist$`DDR gene mutations.csv`['TP53',rownames(Final_table)])
Final_table$tp53_alt=as.character(datalist$`DDR gene alterations.csv`['TP53',rownames(Final_table)])
Final_table$tp53_epi_silence=as.character(datalist$`DDR epigenetic silencing.csv`['TP53',rownames(Final_table)])
Final_table$tp53_deep_del=as.character(datalist$`DDR deep deletions.csv`['TP53',rownames(Final_table)])
Final_table=Final_table[which(Final_table$Genome_doublings==0),]
table(Final_table$ploidy_type)
cli=split(Final_table,as.factor(Final_table$Cohort) )

library(gplots)
library(RColorBrewer)
### this function is used to calcuate the fisher exact test within a matrix
source("code/script/fisher.exact.test.combined.R")
### this function is used to generate the heatmap figures 
source("code/script/heatmap.3.modified.R")
ii = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10",
       "X11", "X12", "X16", "X17", "X18", "X19", "X20")
jj=c("tp53_mut", "tp53_alt", "tp53_epi_silence", "tp53_deep_del")
pdf(file=paste('fisher.single.chromosome.loss.correlation.pdf',sep="."),width=15,height=15)
for(tumor in names(cli)){
  event_del=cli[[tumor]][,ii]
  event_del[event_del==1]=2
  event_del[event_del==-1]=1
  event_del=as.matrix(event_del)
  event_mut=cli[[tumor]][,jj]
  event_mut=apply(event_mut,1:2,function(x) ifelse(x == 1,1,0))
  event=cbind(event_del,event_mut)
  fisher <- fisher.simple.test(event)
  cut.heatmap.default <- 4
  maxp <- max(-log10(unlist(fisher)))
  cut.heatmap <- ifelse(maxp >= cut.heatmap.default,cut.heatmap.default,round(maxp))
  get.mutual.heatmap(fisher[[1]],fisher[[2]],colnames(fisher[[1]]),cut.heatmap,paste("Correlations in single chromosome loss in",tumor),F,F,50)
}
dev.off()


ii = c("X1p", "X1q", "X2p", "X2q", "X3p",
       "X3q", "X4p", "X4q", "X5p", "X5q", "X6p", "X6q", "X7p", "X7q",
       "X8p", "X8q", "X9p", "X9q", "X10p", "X10q", "X11p", "X11q", "X12p",
       "X12q", "X13q", "X14q", "X15q", "X16p", "X16q", "X17p", "X17q",
       "X18p", "X18q", "X19p", "X19q", "X20p", "X20q", "X21q", "X22q")
jj=c("tp53_mut", "tp53_alt", "tp53_epi_silence", "tp53_deep_del")
pdf(file=paste('fisher.arm.chromosome.loss.correlation.pdf',sep="."),width=15,height=15)
for(tumor in names(cli)){
  event_del=cli[[tumor]][,ii]
  event_del[event_del==1]=2
  event_del[event_del==-1]=1
  event_del=as.matrix(event_del)
  event_mut=cli[[tumor]][,jj]
  event_mut=apply(event_mut,1:2,function(x) ifelse(x == 1,1,0))
  event=cbind(event_del,event_mut)
  fisher <- fisher.simple.test(event)
  cut.heatmap.default <- 10
  maxp <- max(-log10(unlist(fisher)))
  cut.heatmap <- ifelse(maxp >= cut.heatmap.default,cut.heatmap.default,round(maxp))
  get.mutual.heatmap(fisher[[1]],fisher[[2]],colnames(fisher[[1]]),cut.heatmap,paste("Correlations in arm chromosome loss in",tumor),F,F,50)
}
dev.off()


Final_table=Final_table[which(!is.na(Final_table$tp53_score)),]
(res=wilcox.test(Final_table$tp53_score[Final_table$ploidy_type %in% "monosomy"], Final_table$tp53_score[Final_table$ploidy_type %in% "disomy"],alternative = 'greater'))
tmp=names(res)
res=lapply(tmp, function(x){res[[x]]})
names(res)=tmp
sink('table/cmp_tp53_score_monosomy_vs_disomy.txt')
print(res)
sink()
(res=wilcox.test(Final_table$tp53_score[Final_table$ploidy_type %in% "monosomy"], Final_table$tp53_score[Final_table$ploidy_type %in% "polysomy"],alternative = 'greater'))
tmp=names(res)
res=lapply(tmp, function(x){res[[x]]})
names(res)=tmp
sink('table/cmp_tp53_score_monosomy_vs_polysomy.txt')
print(res)
sink()
(res=wilcox.test(Final_table$tp53_score[Final_table$ploidy_type %in% "disomy"], Final_table$tp53_score[Final_table$ploidy_type %in% "polysomy"],alternative = 'greater'))
tmp=names(res)
res=lapply(tmp, function(x){res[[x]]})
names(res)=tmp
sink('table/cmp_tp53_score_disomy_vs_polysomy.txt')
print(res)
sink()
reshelp=list()
cnames=unique(Final_table$Cohort)
for(i in cnames){
  res1=NA
  res2=NA
  res3=NA
  tmpmat=subset(Final_table,Final_table$Cohort==i)
  tmp1=tmpmat$tp53_score[tmpmat$ploidy_type %in% "monosomy"]
  tmp2=tmpmat$tp53_score[tmpmat$ploidy_type %in% "disomy"]
  tmp3=tmpmat$tp53_score[tmpmat$ploidy_type %in% "polysomy"]
  if(all(c(length(tmp1),length(tmp2))>2)){
    res1=wilcox.test(tmp1,tmp2)$p.value
  }
  if(all(c(length(tmp1),length(tmp3))>2)){
    res2=wilcox.test(tmp1,tmp3)$p.value
  }
  if(all(c(length(tmp2),length(tmp3))>2)){
    res3=wilcox.test(tmp2,tmp3)$p.value
  }
  reshelp[[i]]=c(res1,res2,res3,i)
}
reshelp=data.frame(do.call(rbind,reshelp),stringsAsFactors = F)
colnames(reshelp)=c('Pvalue_monosomy_vs_disomy','Pvalue_monosomy_vs_polysomy','Pvalue_disomy_vs_polysomy','Cohort')
ploidy_tp53=reshelp
WriteXLS::WriteXLS(c('ploidy_tp53'),'table/tp53_score_vs_ploidy.xlsx',AdjWidth = T)
data_summary <- function(x) {
  m <- median(x)
  ymin <- as.numeric(quantile(x)[2])
  ymax <- as.numeric(quantile(x)[4])
  return(c(y=m,ymin=ymin,ymax=ymax))
}
box1<- ggplot(Final_table, aes(x=factor(ploidy_type), y=tp53_score)) + 
  scale_y_continuous(limits = quantile(Final_table$tp53_score, c(0.01, 0.99))) +
  geom_boxplot(aes(fill=factor(ploidy_type)),width=0.5)+  xlab('sample type')+
  ylab('tp53_score')+
  theme(panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), axis.ticks=element_line(colour="black"),
        axis.text.x=element_text(size=5, colour="black",angle = 60, hjust = 1), axis.text.y=element_text(size=8, colour="black"), axis.title=element_text(size=10, colour="black"),
        axis.title.x = element_text(margin = ggplot2::margin(5,0,0,0)), legend.position = "none") + 
  scale_fill_manual('sample type',values=brewer.pal(n = 5, name = "Greens")) +
  geom_signif(comparisons = list(c("monosomy", "disomy")), annotations="*", y_position =0.85, tip_length = 0.01, size = 0.5, textsize = 3, vjust=0.05) +
  geom_signif(comparisons = list(c("disomy", "polysomy")), annotations="***", y_position =0.9, tip_length = 0.01, size = 0.5, textsize = 3, vjust=0.05) +
  geom_signif(comparisons = list(c("monosomy", "polysomy")), annotations="****", y_position =0.93, tip_length = 0.01, size = 0.5, textsize = 3, vjust=0.05)+
  stat_summary(fun.data=data_summary, color="black", size=0.4)
pdf("plot/tp53_score_vs_ploidy_type_all_pan_cancers.pdf", width =5, height = 4.5)
print(box1)
dev.off()


my_comparisons1 <- list(c("monosomy", "disomy"))
my_comparisons2 <- list( c("disomy", "polysomy"))
my_comparisons3 <- list( c("monosomy", "polysomy"))
cname=c("BLCA","BRCA","CESC","ESCA","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","MESO","SARC","SKCM","STAD","UCEC")
plotdata=subset(Final_table,Final_table$Cohort %in% cname)
p=ggboxplot(plotdata, x = "ploidy_type", y = "tp53_score",#outlier.shape = NA,
          fill = "ploidy_type", palette = c('#99CC99','#CC99CC','#99CCCC'),width = 0.3,bxp.errorbar = TRUE, bxp.errorbar.width = 0.1)+ 
  stat_compare_means(comparisons = my_comparisons1,label.y = 0.1)+ # Add pairwise comparisons p-value
  stat_compare_means(comparisons = my_comparisons2,label.y = 0.2)+
  stat_compare_means(comparisons = my_comparisons3,label.y = 0.85)
p=p+theme(panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), axis.ticks=element_line(colour="black"), axis.text.x=element_text(size=8, colour="black"), axis.text.y=element_text(size=8, colour="black"), axis.title=element_text(size=10, colour="black"), axis.title.x = element_text(margin = ggplot2::margin(5,0,0,0)), legend.key = element_rect(fill="White"), legend.text = element_text(size=14), plot.title = element_text(hjust=0.5)) 
p=p+stat_n_text()
pdf("plot/tp53_score_vs_ploidy_per_cancer.pdf", width =13, height = 10)
facet(p, facet.by = "Cohort",ncol = 5)
dev.off()
WriteXLS::WriteXLS(c('plotdata'),'zuzanne/res_update/table/ssgsea_vs_ploidy_tcga_source_data.xlsx',AdjWidth = T)

write.xlsx(plotdata, file = 'res_update/table/ssgsea_vs_ploidy_tcga_source_data.xlsx',col.names = T,row.names = F)
mat=Final_table[which(Final_table$tp53_alt!='' & Final_table$ploidy_type %in% c('monosomy','polysomy')),]
mat$ploidy_type=as.character(mat$ploidy_type)
cli1=split(mat,as.factor(mat$Cohort) )
cli1=lapply(cli1,function(x){
  tmp=x;
  tmp$ploidy_type=as.character(tmp$ploidy_type);
  tmp$ploidy_type=as.character(tmp$ploidy_type);
  tmp=tmp[which(tmp$ploidy_type %in% c('monosomy','polysomy')),];
  return(tmp)})
rr=lapply(cli1,function(x){return(table(x[,c("ploidy_type","tp53_alt")]))});
iu=names(rr)[which(unlist(lapply(rr,function(x){return(dim(x)[1]==2 & dim(x)[2]==2)})))]


txt=data.frame(pvalue=format(unlist(lapply(iu,function(x){return(fisher.test(rr[[x]])$p.value)})),digits=1),cancer=iu)
