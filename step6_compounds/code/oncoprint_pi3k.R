rm(list=ls())
library("tidyverse") # loads ggplot2, tibble, tidyr, readr, purr, dplyr
library("genefilter")
library("edgeR")
library("ggrepel")
library("gplots")
library("reshape2")
library("data.table")
library("magrittr")
library("GSVA")
library("GSA")
library("limma")
library("edgeR")
library("rio")
library(gdata)
library(gplots)
library(RColorBrewer)
library(impute)
setwd("/data/zhang/pancan_cin/step9_driver_interaction")

colornumber=4
genecol=rev(brewer.pal(colornumber, "RdBu"))[c(1:colornumber)]
scorecol=rev(brewer.pal(colornumber, "PiYG"))[c(1:colornumber)]
setwd("/data/zhang/pancan_cin/step6_compounds/")
load('data/ccle_cin_score.RData')
scores$Genome.doublings=as.factor(scores$Genome.doublings)
data=scores[which(!is.na(scores$tcga_code)),]
data$Cohort=as.character(data$tcga_code)
data=data[which(complete.cases(data$Genome.doublings)),]
Final_table=data[which(data$tcga_code!="UNABLE TO CLASSIFY"),]


pi3k=c("PIK3CA_MUT", "PIK3CA_AMP", "PTEN_MUT","PTEN_DEL")
cnvtable <- fread("data/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct",skip = 2,data.table = F)
rownames(cnvtable)=cnvtable$Name
cnvtable=cnvtable[,-c(1:2)]
cnvtable=subset(cnvtable,rownames(cnvtable) %in% pi3k)
cnvtable=data.frame(t(cnvtable))
colnames(cnvtable)=unlist(lapply(colnames(cnvtable), function(x){tmpx=strsplit(x,'_')[[1]];return(paste0(c(tmpx[2],tmpx[1]),collapse = '.'))}))
cnvtable[,'Sample']=rownames(cnvtable)

rna_table =fread("../step6_compounds/data/CCLE_RNAseq_genes_rpkm_20180929.gct")
rna_table=rna_table[which(!duplicated(rna_table$Description)),]
ii= rna_table$Description
rna_table$Name = NULL
rna_table$Description = NULL
rna_table=t(rna_table)
colnames(rna_table)=ii
rna_table=data.frame(rna_table)
rna_table$Sample = rownames(rna_table)

rna_idx=c('PIK3CA','PTEN')

sels=Reduce(intersect,list(cnvtable$Sample,Final_table$Sample,rna_table$Sample))
cnvtable=subset(cnvtable,Sample %in% sels)
rna_table=subset(rna_table,Sample %in% sels)
Final_table=subset(Final_table,Sample %in% sels)
ii=colnames(cnvtable)[-grep('Sample',colnames(cnvtable))]
z <-merge(Final_table,cnvtable,by='Sample')

rownames(z)=z$Sample
z1=z[,c('wgii','scin')]
z=z[,ii]
z[is.na(z)]=0
z=z[which(rowSums(z)>0),]
z1=z1[rownames(z),]

ii=colnames(z)
col_order=order(colSums(z), decreasing=TRUE)
scoreCol=function(x) {
  score=0
  for(i in 1:length(x)) {
    if(x[i]) {
      score=score + 2^(length(x)-i*1/x[i])
    }
  }
  return(score)
}
scores=apply(z[,col_order], 1, scoreCol)
ii=colnames(cnvtable)[-grep('Sample',colnames(cnvtable))]
z <-merge(Final_table,cnvtable,by='Sample')

rownames(z)=z$Sample
z1=z[,c('wgii','scin')]
z=z[,ii]
z[is.na(z)]=0
z=z[which(rowSums(z)>0),]
z1=z1[rownames(z),]
ii=colnames(z)
col_order=order(colSums(z), decreasing=TRUE)
scoreCol=function(x) {
  score=0
  for(i in 1:length(x)) {
    if(x[i]) {
      score=score + 2^(length(x)-i*1/x[i])
    }
  }
  return(score)
}
scores=apply(z[,col_order], 1, scoreCol)
row_order=order(scores, decreasing=TRUE)
alttype=unique(unlist(lapply(ii,function(x){strsplit(x,'\\.')[[1]][1]})))
altcol=colorRampPalette(yarrr::piratepal(palette="pony",trans=0)[-c(6:7)])(length(alttype))
names(altcol)=alttype
# colMutations =unlist(lapply(ii,function(x){altcol[strsplit(x,'\\.')[[1]][1]]}))
# names(colMutations) <- ii
colMutations=structure(c("#EB5291", "#9DDAF5", "#972C8D", "#972C8D"), .Names = c("AMP.PIK3CA", 
                                                                    "DEL.PTEN", "MUT.PIK3CA", "MUT.PTEN"))

par(bty="n", mgp=c(2,.33,0), mar=rep(3,4), las=1, tcl=-.25, xpd=NA)
plot(NA,NA, xlim=c(0,ncol(z)+2), ylim=c(0,nrow(z)), xaxt="n", yaxt="n", xlab="",ylab="", xaxs="i", yaxs="i")
xtick<-seq(0.5, length(ii)+2+length(c(rna_idx))-0.5, by=1)
axis(side=2, at=xtick, labels=FALSE)
z=z[row_order,col_order]
tmpmat=rna_table[rownames(z),rna_idx]
tmpmat=data.frame(impute::impute.knn(as.matrix(tmpmat))$data)
zmat=sapply(colnames(z), function(i) ifelse(z[,i]>0, colMutations[i], "#FFFFFF"))
zmat=data.frame(zmat)
rownames(zmat)= rownames(z)
colnames(zmat)=colnames(z)
score_mat=z1[row_order,]

score_mat=apply(score_mat, 2, function(y) rank(y) / length(y))
norm_score=score_mat
score_mat=apply(score_mat,2,function(x){scorecol[factor(arules::discretize(x,method='frequency',breaks=colornumber))]})
rownames(score_mat)=rownames(norm_score)
png('plot/onco_pi3k_ccle.png',height =7.3,width=9,units='cm',res=300)
par(bty="n", mgp=c(2,.33,0), mar=c(0.2,7,0.1,0.1), las=1, tcl=-.25, xpd=NA)
#tmpmat=apply(tmpmat, 2, function(y){(y - mean(y)) / sd(y)})
tmpmat=apply(tmpmat, 2, function(y) rank(y) / length(y))
genemat=tmpmat
tmpmat=apply(tmpmat,2,function(x){genecol[factor(arules::discretize(x,method='frequency',breaks=colornumber))]})

imagemat=cbind(score_mat,tmpmat)
imagemat=cbind(imagemat,zmat)
imagemat=t(as.matrix(imagemat))
rownames(imagemat)=c( "WGII",'SCIN',"GEXP.PIK3CA", "GEXP.PTEN", "AMP.PIK3CA", "MUT.PIK3CA", 
                      "MUT.PTEN", "DEL.PTEN")
ii=rownames(imagemat)
textcol=ii[-c(1:4)]
textcol=c(rep('black',4),colMutations[textcol])
textcol=rev(textcol)
plot(NA,NA, xlim=c(0,ncol(imagemat)), ylim=c(0,nrow(imagemat)), xaxt="n", yaxt="n", xlab="",ylab="", xaxs="i", yaxs="i")
rasterImage(imagemat, 0, 0, ncol(imagemat), nrow(imagemat), interpolate=FALSE,)
text(y=xtick,  par("usr")[4], labels=rev(ii), srt=0, pos=2, xpd=T,offset=1.5,cex=1,col=textcol)
dev.off()

pdf('plot/colorbar.pdf',15,19)
par(mar=c(3,3,0,3))
genecol_ba=rev(brewer.pal(colornumber+2, "RdBu"))
scorecol_ba=rev(brewer.pal(colornumber+2, "PiYG"))
plot(NA,NA, xlim=c(2,3), ylim=c(-1,15), xaxt="n", yaxt="n", xlab="",ylab="", xaxs="i", yaxs="i")
n=length(scorecol_ba)
image(y=1:n , x=c(2,2.5), z=matrix(c(1:n), nrow=1), col=scorecol_ba, add=TRUE)
image(y=1:n, x=c(2.5,3), z=matrix(c(1:n), nrow=1), col=genecol_ba, add=TRUE)
jj=altcol
image(y=(n+2):(n+1+length(jj)), x=c(2.5,3), z=matrix(c(1:length(jj)), nrow=1), col=jj, add=TRUE)
axis(side=4, at=seq((n+2),(n+1+length(jj))),  tcl=-.5, label=names(jj), las=1, lwd=.5)
scalelab=unlist(lapply(levels(arules::discretize(norm_score[,2],method='frequency',breaks=colornumber)), function(x){gsub('\\[','',strsplit(x,',')[[1]][1])}))
scalelab=round(as.numeric(scalelab),2)
axis(side=4, at=seq(1,n),  tcl=-.5, label=c(0,scalelab,1), las=1, lwd=.5)
axis(side=2, at=seq(1,n),  tcl=-.5, label=c(0,scalelab,1), las=1, lwd=.5)
text(x=2.5, y=0, "normalised gene expression", pos=4)
text(x=2, y=n+1, "normalised CIN score", pos=4)
text(x=2.5,y=n+5,'alteration',pos=4)
dev.off()
