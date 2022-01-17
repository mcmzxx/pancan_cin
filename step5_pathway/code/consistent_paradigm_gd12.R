rm(list=ls())
setwd("/data/zhang/pancan_cin/step5_pathway/")
fls1=list.files('table/paradigm_pathway_gd12',pattern = 'cor_paradigm_wgii_',full.names = T)
fls2=list.files('table/paradigm_pathway_gd12',pattern = 'cor_paradigm_scin_',full.names = T)
results=list()
df=NULL
for(fl in fls1){
  load(fl)
  df=rbind(df,Results)
  tmpname=gsub('.RData','',gsub('cor_paradigm_wgii_','',strsplit(fl,'\\/')[[1]][3]))
  results[[tmpname]]=Results$Pathway[which(abs(Results$Spearman_coef)>=0.3)]
}
consistent=table(unlist(results))
topPathway=names(consistent)[which(consistent>=8)]
plotdf=df[which(df$Pathway %in% topPathway),]
cnames=unique(plotdf$C)
tmpres=lapply(cnames,function(x){
  tmp=subset(df,df$Cohort==x & df$Pathway %in% topPathway);
  res=tmp[,c('Cohort','Pathway','Spearman_coef')];
  rownames(res)=res$Pathway;
  res=res[topPathway,];
  return(res$Spearman_coef)
})
names(tmpres)=cnames
heatdat=data.frame(do.call(cbind,tmpres))
heatdat[is.na(heatdat)]=0
rownames(heatdat)=topPathway
library(pheatmap)
library(RColorBrewer)
pheatmap::pheatmap(heatdat,scale = "none",
                   clustering_method = "ward.D",
                   clustering_distance_cols = "euclidean",
                   clustering_distance_rows = "euclidean",
                   border_color = NA,main ='consistent Pathway of wgii',
                   show_rownames =T,show_colnames = T,treeheight_row = 0,fontsize=3,fontsize_row =3,cellwidth =3,cellheight =3,treeheight_col = 0,
                   file ="/data/zhang/pancan_cin/step5_pathway/plot/consistent_pathway_of_wgii_gd12.pdf",
                   width =10.5, height =28,color = rev(colorRampPalette(brewer.pal(11,"RdBu"))(2000)))


results=list()
df=NULL
for(fl in fls2){
  load(fl)
  Results=Results1
  df=rbind(df,Results)
  tmpname=gsub('.RData','',gsub('cor_paradigm_scin_','',strsplit(fl,'\\/')[[1]][3]))
  results[[tmpname]]=Results$Pathway[which(abs(Results$Spearman_coef)>=0.4)]
}
consistent=table(unlist(results))
topPathway=names(consistent)[which(consistent>=8)]
plotdf=df[which(df$Pathway %in% topPathway),]
cnames=unique(plotdf$C)
tmpres=lapply(cnames,function(x){
  tmp=subset(df,df$Cohort==x & df$Pathway %in% topPathway);
  res=tmp[,c('Cohort','Pathway','Spearman_coef')];
  rownames(res)=res$Pathway;
  res=res[topPathway,];
  return(res$Spearman_coef)
})
names(tmpres)=cnames
heatdat=data.frame(do.call(cbind,tmpres))
heatdat[is.na(heatdat)]=0
rownames(heatdat)=topPathway
library(pheatmap)
library(RColorBrewer)
pheatmap::pheatmap(heatdat,scale = "none",
                   clustering_method = "ward.D",
                   clustering_distance_cols = "euclidean",
                   clustering_distance_rows = "euclidean",
                   border_color = NA,main ='consistent Pathway of scin',
                   show_rownames =T,show_colnames = T,treeheight_row = 0,fontsize=3,fontsize_row =3,cellwidth =3,cellheight =3,treeheight_col = 0,
                   file ="/data/zhang/pancan_cin/step5_pathway/plot/consistent_pathway_of_scin_gd12.pdf",
                   width =10.5, height =18,color = rev(colorRampPalette(brewer.pal(11,"RdBu"))(2000)))