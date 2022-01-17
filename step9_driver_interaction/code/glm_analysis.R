rm(list = ls())
options(width = 120)
pdf.options(pointsize = 8)
knit_hooks$set(
  smallMar = function(before, options, envir) {
    if (before)
      par(mar = c(3, 3, 1, 1),
          bty = "n",
          mgp = c(2, .5, 0))
  }
)
opts_chunk$set(
  dev = c('my_png', 'pdf'),
  fig.ext = c('png', 'pdf'),
  fig.width = 3,
  fig.height = 3,
  smallMar = TRUE
)
my_png <-  function(file, width, height, pointsize = 12, ...) {
  png(
    file,
    width = 1.5 * width,
    height = 1.5 * height,
    units = "in",
    res = 72 * 1.5,
    pointsize = pointsize,
    ...
  )
}
#opts_knit$set(root.dir = file.path(getwd(),".."))

#' TCGA AML expression data analysis
#' ===========================

#' #### Libraries
library(limma)
library(org.Hs.eg.db)
library(hgu133plus2.db)
library(RColorBrewer)
set1 <- brewer.pal(9, "Set1")
library(cgdsr)
library(CoxHD)
library(mg14)
library(xtable)
library(Hmisc)
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
    mattmp = tmp[-c(1:2),-c(1:2)]
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
mattmp = tmp[-c(1:4),-c(1:2, 4)]
rownames(mattmp) = rowtmp
colnames(mattmp) = coltmp
Final_table = mattmp[,-1]

Final_table = data.frame(apply(Final_table, 2, as.numeric))
Final_table$Cohort = mattmp$Cohort
Final_table$Sample=rownames(mattmp)
Final_table = Final_table[which(!is.na(Final_table$ploidy)), ]
rownames(Final_table)=Final_table$Sample
table <-
  read.delim(
    "../step2_genome_instability/data/TaylorCancerCell_TableS2.txt",
    na.strings = c("", " ", "na", "NA", "n.a.", "#N/A")
  )

sels = intersect(table$Sample,Final_table$Sample)
data = Final_table[sels, ]

Final_table <- merge(data, table, by = 'Sample')
rownames(Final_table) = Final_table$Sample
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
score_table = scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer', 'Primary Solid Tumor'), ]
Final_table <- merge(Final_table, score_table, by = 'Sample')
library(gplots)
library(RColorBrewer)
rna_table <-
  fread("../step8_tme/data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv")
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
sels=intersect(rna_table$Sample,Final_table$Sample)

Final_table <- merge(Final_table, rna_table, by = 'Sample')
Final_table$Cohort = Final_table$Cohort.x
Final_table = Final_table[, which(!(
  colnames(Final_table) %in% c('Cohort.x', 'Cohort.y', 'Type.x', 'Type.y')
))]
onkotablebak=fread('/data/zhang/pancan_cin/step9_driver_interaction/data/onkokb.csv')

onkotable=data.frame(apply(onkotablebak[,-1], 2, function(x){ifelse(x==FALSE,0,1)}),stringsAsFactors = F)
rownames(onkotable)=onkotablebak$SAMPLE_BARCODE
sels=intersect(sels,rownames(onkotable))
sels=intersect(sels,Final_table$Sample[which(Final_table$Co %in% c('COAD'))])
design = cbind(offset=1,onkotable[sels,]) # oncogenic mutations
minF=5 ## Minimal number of alterations
design = design[,colSums(design,na.rm = T)>=minF]
design0 <- design
for(j in 1:ncol(design))
  design[is.na(design[,j]),j] <- mean(design[,j], na.rm=TRUE)

geneExpr=t(rna_table[sels,-grep('Sample',colnames(rna_table))])
load('/data/zhang/pancan_cin/step9_driver_interaction/table/gene_corr.RData')
selg=datalist$scin[which(datalist$scin$Cohort=='COAD'),]
selg=selg[which(abs(selg$Spearman_coef)>=0.35),'Gene']
geneExpr[is.na(geneExpr)]=0
geneExpr=geneExpr[selg,sels]
geneExpr <- t(apply(geneExpr, 1, function(x)(x-min(x))/(max(x)-min(x))))

glm = lmFit(geneExpr, design = design ) 
glm = eBayes(glm)
F.stat <- classifyTestsF(glm[,],fstat.only=TRUE) # remove offset
glm$F <- as.vector(F.stat)
df1 <- attr(F.stat,"df1")
df2 <- attr(F.stat,"df2")
if(df2[1] > 1e6){ # Work around bug in R 2.1
  glm$F.p.value <- pchisq(df1*glm$F,df1,lower.tail=FALSE)
}else
  glm$F.p.value <- pf(glm$F,df1,df2,lower.tail=FALSE)
F.stat <- classifyTestsF(glm[,-1],fstat.only=TRUE) # remove offset
glm$F <- as.vector(F.stat)
df1 <- attr(F.stat,"df1")
df2 <- attr(F.stat,"df2")
if(df2[1] > 1e6){ # Work around bug in R 2.1
  glm$F.p.value <- pchisq(df1*glm$F,df1,lower.tail=FALSE)
}else
  glm$F.p.value <- pf(glm$F,df1,df2,lower.tail=FALSE)

#' #### Random model
#' Compare to a model where all values of the covariates are randomly permuted. If all model assumptions were correct, this wouldn't be needed.
set.seed(42)
rlm <- lmFit(geneExpr[,rownames(design)], apply(design, 2, sample))
rlm <- eBayes(rlm)
F.stat <- classifyTestsF(rlm[,-1],fstat.only=TRUE)
rlm$F <- as.vector(F.stat)
df1 <- attr(F.stat,"df1")
df2 <- attr(F.stat,"df2")
if(df2[1] > 1e6){ # Work around bug in R 2.1
  rlm$F.p.value <- pchisq(df1*rlm$F,df1,lower.tail=FALSE)
}else
  rlm$F.p.value <- pf(rlm$F,df1,df2,lower.tail=FALSE)

#' #### Explained variance by different categories
#' The F-statistic is directly related to the R2.
F.stat <- classifyTestsF(glm[,2:dim(glm)[2]],fstat.only=TRUE) ## All genetics & cytogenetics
df1 <- attr(F.stat,"df1")
df2 <- attr(F.stat,"df2")
F.p.value <- pchisq(df1*F.stat,df1,lower.tail=FALSE)

R.stat <- classifyTestsF(rlm[,2:dim(glm)[2]],fstat.only=TRUE) ## Random

Rall = 1 - 1/(1 + glm$F * (ncol(design)-1)/(nrow(design)-ncol(design)))
Rgenetics = 1 - 1/(1 + F.stat * dim(glm)[2]/(nrow(design)-ncol(design)))
Pgenetics = 1 - 1/(1 + R.stat * dim(glm)[2]/(nrow(design)-ncol(design)))
names(Rgenetics) <- names(Pgenetics) <- names(Rall) <-  rownames(geneExpr)

#' Plot the variance explained by genetics
#+ GE-R2, fig.width=2, fig.height=1.8
pdf('plot/scin_COAD_explained_var.pdf')
par(bty="n", mgp = c(2,.33,0), mar=c(3,2.5,1,1)+.1, las=1, tcl=-.25, xpd=NA)
d <- density(Pgenetics,bw=1e-3)
f <- 1#nrow(gexpr)/512
plot(d$x, d$y * f, col='grey', xlab=expression(paste("Explained variance per gene ", R^2)), main="Differentially expressed SCIN-genes by oncogenic events", lwd=2, type="l", ylab="", xlim=c(0,0.7))
title(ylab="Density", line=1.5)
d <- density(Rgenetics, bw=1e-3)
r <- min(Rgenetics[p.adjust(F.p.value,"BH")<0.05])
x0 <- which(d$x>r)
polygon(d$x[c(x0[1],x0)], c(0,d$y[x0])* f, col=paste(set1[1],"44",sep=""), border=NA)
lines(d$x, d$y* f, col=set1[1], lwd=2)
#points(d$x[x0[1]], d$y[x0[1]]*f, col=set1[1], pch=16)
text(d$x[x0[1]], d$y[x0[1]]*f, pos=4, paste0(sum(Rgenetics > r),"/",length(Rgenetics), " genes q < 0.05"))
arrows(Rgenetics["SMAD4"], par("usr")[4]/7, Rgenetics["SMAD4"], par("usr")[4]/50, length=0.05)
text(Rgenetics["SMAD4"], par("usr")[4]/8, "SMAD4", font=3, pos=3)
dev.off()


testResults <- decideTests(glm, method="hierarchical",adjust.method="BH", p.value=0.05)[,-1]
significantGenes <- sapply(1:ncol(testResults), function(j){
  c <- glm$coefficients[testResults[,j]!=0,j+1]
  table(cut(c, breaks=c(-5,seq(-1.5,1.5,l=7),5)))
})
colnames(significantGenes) <- colnames(testResults)

t <- head(sort(Rgenetics, d=TRUE), 100)
g <- apply(testResults[names(t),], 1, function(x) paste(colnames(testResults)[x!=0], collapse=", "))
print(xtable(cbind(AnnotationDbi::select(org.Hs.eg.db, names(t), c("SYMBOL","GENENAME"), "SYMBOL"),`R^2`=t, `Mutations`=g)), type = 'html' )



glmPrediction <- glm$coefficients %*% t(design)
rlmPrediction <- rlm$coefficients %*% t(design)
colMutations = c(brewer.pal(8,"Set1")[-6], rev(brewer.pal(8,"Dark2")), brewer.pal(7,"Set2"))[c(1:12,16:19,13:15)]

#' #### ABCB7 expression
#+ GE-ABCB7,  fig.width=2, fig.height=1.8
for(w in names(t))
{
    pdf(paste0(c('plot/scin_linear_driver/scin_COAD_explained_variance_',w,'.pdf'),collapse = ''))
    par(bty="n", mgp = c(1.5,.33,0), mar=c(2.5,2.5,1,1)+.1, las=1, tcl=-.25)
      gene <- w
      plot(glmPrediction[w,], geneExpr[w,rownames(design)], ylab=parse(text=paste("Observed  ~ italic(",gene,") ~ expression")), xlab=parse(text=paste("Predicted  ~ italic(",gene,") ~ expression")), pch=16, cex=.8)
      par(xpd=FALSE)
      abline(0,1)
      u <- par("usr")
      par(xpd=NA)
      y <- glm$coefficients[w,-1]+glm$coefficients[w,1]
      u <- par("usr")
      x0 <- rep(u[4]-(u[4]-u[3])/8,ncol(design)-1)
      y0 <- u[4] + 0.05*(u[4]-u[3]) - rank(-y)/length(y) * (u[4]-u[3])/1.2
      d <- density(y)
      lines(d$x, d$y/100+u[3]+0.95, col="black")
      lines(d$x, -d$y/100+u[3]+0.95, col="black")
      points(x=y, y=x0+violinJitter(y, magnitude=0.015)$y,pch=19, col=colMutations)
      text(x=glm$coefficients[w,1], y= u[4], "Model coefficients (logFC)", cex=0.8)
      v <- glm$p.value[w,-1] < 0.01
      if(length(which(v))>0)
      rotatedLabel(y[v], x0[v]+0.1, labels=colnames(design)[-1][v], font=ifelse(grepl("[[:lower:]]", colnames(design)[-1]),1,3)[v], cex=.66, pos=1)
      axis_idxs=round(min(y),2)
      axis_idxe=round(max(y),2)
      axis(at=seq(axis_idxs,axis_idxe,by = 0.1),glm$coefficients[w,1], labels=seq(axis_idxs,axis_idxe,by = 0.1), side=3, cex.axis=.8, line=-1, mgp = c(0.5,.05,0.6), tcl=-.15)
      text(u[2],u[3] + (u[4]-u[3])/2, substitute(paste(R^2==r),list(r=round(Rgenetics[w],2))), pos=2)
    dev.off()
}


