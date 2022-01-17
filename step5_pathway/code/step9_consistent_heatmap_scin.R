rm(list=ls())
df=read.table('/data/zhang/pancan_cin/step5_pathway/table/paradigm_pathway_corr_results_wgii.txt',header = T,sep='\t',stringsAsFactors = F)
colnames(df)=c('cancer','pathway','cor_spearman','p_spearman','cor_pearson','p_pearson')
cnames=unique(df$cancer)
results=lapply(cnames, function(x){
  tmp=subset(df,df$cancer==x)
  return(tmp$pathway[which(abs(tmp$cor_spearman)>=0.3)])
})
consistent=table(unlist(results))
toppathway=names(consistent)[which(consistent>=8)]
plotdf=df[which(df$pathway %in% toppathway),]
cnames=unique(plotdf$cancer)
tmpres=lapply(cnames,function(x){
  tmp=subset(df,df$cancer==x & df$pathway %in% toppathway);
  res=tmp[,c('cancer','pathway','cor_spearman')];
  rownames(res)=res$pathway;
  res=res[toppathway,];
  return(res$cor_spearman)
})
names(tmpres)=cnames
heatdat=data.frame(do.call(cbind,tmpres))
heatdat[is.na(heatdat)]=0
rownames(heatdat)=toppathway
library(pheatmap)
library(RColorBrewer)
pheatmap::pheatmap(heatdat,scale = "none",
                   clustering_method = "ward.D",
                   clustering_distance_cols = "euclidean",
                   clustering_distance_rows = "euclidean",
                   border_color = NA,main ='consistent pathway of wgii',
                   show_rownames =T,show_colnames = T,treeheight_row = 0,fontsize=3,fontsize_row =3,cellwidth =3,cellheight =3,treeheight_col = 0,
                   file ="/data/zhang/pancan_cin/step5_pathway/plot/consistent_pathway_of_wgii.pdf",
                   width =10.5, height =28,color = rev(colorRampPalette(brewer.pal(11,"RdBu"))(2000)))