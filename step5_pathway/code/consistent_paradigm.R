rm(list=ls())
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(ggplot2)
setwd("/data/zhang/pancan_cin/step5_pathway/")
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer','Primary Solid Tumor'),]
maxmin_norm=function(x){(x-min(x))/(max(x)-min(x)) }
Final_table$Sample <- substr(Final_table$Sample, 1, 12)
rownames(Final_table)=Final_table$Sample
Final_table$Genome_doublings <- as.factor(Final_table$`Genome doublings`)
imm=fread('data/merge_merged_reals.txt')
ii=imm$Gene
rownames(imm)=imm$Gene
imm$Gene=NULL
imm=t(as.data.frame(imm))
colnames(imm)=ii
imm=data.frame(imm)
imm$Sample=rownames(imm)
# Final_table=merge(Final_table,til,by='Sample')
Final_table=merge(Final_table,imm,by='Sample')
data=Final_table[which(complete.cases(Final_table$Genome_doublings)),]
data=data[with(data,order(Cohort,scin)),]
data=data[which(data$Cohort!='SKCM'),]
counts=table(data$Cohort)
xindex=integer(sum(counts)) # preallocate to required length with zero values
ci <- 1 # cumulative count for vector indexing
for (i in counts) {
  tmp=seq(i)
  xindex[ci:(ci+i-1)]=tmp 
  ci <- ci+i
}
data$NO <- xindex
data1=split(data,data$Cohort)
data$NO_ncin=unlist(lapply(data1,function(x){sample(x$NO)}))
dataf=ddply(data, .(Cohort), summarise,pscin=median(scin),pwgii=median(wgii),px=median(NO))
dataf$NO=dataf$NO
dataf$Cohort=as.character(dataf$Cohort)

ii=dataf$Cohort[order(dataf$pwgii)]
dataf$Cohort=factor(dataf$Cohort,levels = ii)

fls1=list.files('table/paradigm_pathway',pattern = 'cor_paradigm_wgii_',full.names = T)
fls2=list.files('table/paradigm_pathway',pattern = 'cor_paradigm_scin_',full.names = T)
results=list()
df=NULL
for(fl in fls1){
  load(fl)
  df=rbind(df,Results)
  tmpname=gsub('.RData','',gsub('cor_paradigm_wgii_','',strsplit(fl,'\\/')[[1]][3]))
  results[[tmpname]]=Results$Pathway[which(Results$Spearman_coef>=0.3)]
}
consistent=table(unlist(results))
topPathway=names(consistent)[which(consistent>=7)]
topPathway=topPathway[-grep('_',topPathway)][1:15]

#consistent[topPathway]
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

#colind=colnames(heatdat)[hclust(dist(t(heatdat)))$order]
colind=dataf[order(dataf$pwgii),'Cohort']
rowind=rownames(heatdat)[hclust(dist(heatdat))$order]
heatdat$Pathway=rownames(heatdat)
df2=gather(heatdat, Tissue, Corr_wgii, ACC:UVM, factor_key=TRUE)
df2$Pathway=factor(df2$Pathway,levels=rowind)
df2$Tissue=factor(df2$Tissue,levels=colind)
df3=data[,c(topPathway,"Cohort")]
rownames(df3)=rownames(data)
df3=ddply(df3, .(Cohort), numcolwise(median))
rownames(df3)=df3$Cohort
df3=apply(df3[,-1],2,maxmin_norm)
df3=data.frame(t(df3))
df3$Pathway=rownames(df3)
df3=gather(df3, Tissue, Paradigm, ACC:UVM,factor_key=TRUE)
df3$Pathway=factor(df3$Pathway,levels=rowind)
df3$Tissue=factor(df3$Tissue,levels=colind)
pwgii = ggplot() + geom_point(
  aes(
    x = px,
    y = pwgii
  ),data = dataf,
  colour = 'red'
)
pwgii=pwgii+ylab('WGII score')
pwgii =pwgii + facet_grid(. ~ Cohort, scales = 'free', space = 'fix')+ theme_cowplot() + theme(
  axis.text.x = element_text(size =1, angle = -270,colour = 'white'),
  axis.title.x = element_blank(),axis.title.y = element_text(size =10),
  axis.ticks.x = element_blank(),panel.background = element_rect(fill="white",colour="white")
)
pwgii= pwgii + theme(strip.background = element_rect(fill = 'white'),strip.text.x = element_blank(),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pwgii=pwgii+scale_y_continuous(labels = scales::label_number(accuracy = 0.1))
pwgii
heatplot=ggplot(df2, aes(y=Pathway,x=Tissue)) + 
  geom_tile(aes(fill = Corr_wgii), colour = "white") + 
  geom_text(aes(label = round(Corr_wgii,2)), size=2) +
  scale_fill_distiller(palette = "RdBu",type ='div',direction =-1 ) +
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"),
        panel.grid = element_blank(),axis.title.x =element_blank(),axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))

heatplot=ggplot(df3, aes(y=Pathway,x=Tissue)) + 
  geom_tile(aes(fill = Paradigm), colour = "white") + 
  geom_text(aes(label = round(Paradigm,2)), size=2) +
  scale_fill_distiller(palette = "RdBu",type ='div',direction =-1 ) +
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"),
        panel.grid = element_blank(),axis.text.y = element_text(color = '#006747FF'),axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1),axis.title = element_blank())
heatplot=heatplot+theme(legend.position = "none")
mainplot_wgii=plot_grid(pwgii,heatplot,ncol = 1,rel_heights =c(1,4),align = 'v',axis = "l")
print(mainplot_wgii)
pdf("plot/paradigm_order_by_wgii.pdf",width =6,height = 5)
print(mainplot_wgii)
dev.off()


results=list()
df=NULL
for(fl in fls2){
  load(fl)
  Results=Results1
  df=rbind(df,Results)
  tmpname=gsub('.RData','',gsub('cor_paradigm_scin_','',strsplit(fl,'\\/')[[1]][3]))
  results[[tmpname]]=Results$Pathway[which(Results$Spearman_coef>=0.3)]
}
consistent=table(unlist(results))
topPathway=names(consistent)[which(consistent>=7)]
topPathway=topPathway[-grep('_',topPathway)][1:15]
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

#colind=colnames(heatdat)[hclust(dist(t(heatdat)))$order]
colind=dataf[order(dataf$pscin),'Cohort']
rowind=rownames(heatdat)[hclust(dist(heatdat))$order]
heatdat$Pathway=rownames(heatdat)
df2=gather(heatdat, Tissue, Corr_scin, ACC:UVM, factor_key=TRUE)
df2$Pathway=factor(df2$Pathway,levels=rowind)
df2$Tissue=factor(df2$Tissue,levels=colind)
ii=dataf$Cohort[order(dataf$pscin)]
dataf$Cohort=factor(dataf$Cohort,levels = ii)
pscin = ggplot() + geom_point(
  aes(
    x = px,
    y = pscin
  ),data = dataf,
  colour = 'red'
)
pscin=pscin+ylab('SCIN score')
pscin =pscin + facet_grid(. ~ Cohort, scales = 'free', space = 'fix')+ theme_cowplot() + theme(
  axis.text.x = element_text(size =1, angle = -270,colour = 'white'),
  axis.title.x = element_blank(),axis.title.y = element_text(size =10),
  axis.ticks.x = element_blank(),panel.background = element_rect(fill="white",colour="white")
)
pscin= pscin + theme(strip.background = element_rect(fill = 'white'),strip.text.x = element_blank(),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pscin
heatplot=ggplot(df2, aes(y=Pathway,x=Tissue)) + 
  geom_tile(aes(fill = Corr_scin), colour = "white") + 
  geom_text(aes(label = round(Corr_scin,2)), size=2) +
  scale_fill_distiller(palette = "RdBu",type ='div',direction =-1 ) +
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"),
        panel.grid = element_blank(),axis.title.x =element_blank(),axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))

df3=data[,c(topPathway,"Cohort")]
rownames(df3)=rownames(data)
df3=ddply(df3, .(Cohort), numcolwise(median))
rownames(df3)=df3$Cohort
df3=apply(df3[,-1],2,maxmin_norm)
df3=data.frame(t(df3))
df3$Pathway=rownames(df3)
df3=gather(df3, Tissue, Paradigm, ACC:UVM,factor_key=TRUE)
df3$Pathway=factor(df3$Pathway,levels=rowind)
df3$Tissue=factor(df3$Tissue,levels=colind)

heatplot=ggplot(df3, aes(y=Pathway,x=Tissue)) + 
  geom_tile(aes(fill = Paradigm), colour = "white") + 
  geom_text(aes(label = round(Paradigm,2)), size=2) +
  scale_fill_distiller(palette = "RdBu",type ='div',direction =-1 ) +
  theme(legend.position = "bottom", legend.key.width = unit(2, "cm"),
        panel.grid = element_blank(),axis.text.y = element_text(color = '#006747FF'),axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1),axis.title = element_blank())
guides(color=guide_legend(title="PARADIGM pathway activity"))
legend=get_legend(heatplot+ theme(legend.key = element_rect(colour = NA, fill = NA),legend.box.margin = margin(-1, -1, -1, 0)))
heatplot=heatplot+theme(legend.position = "none")
mainplot_scin=plot_grid(pscin,heatplot,ncol = 1,rel_heights =c(1,4),align = 'v',axis = "l")
print(mainplot_scin)

pdf("plot/paradigm_order_by_scin.pdf",width =6,height = 5)
print(mainplot_scin)
dev.off()


prow=plot_grid(mainplot_wgii,mainplot_scin, align = "v", axis = "l",ncol=1,rel_heights= c(1,1), labels = c("A","C"),label_x = -0.005)
pdf('plot/paradigm_order_by_cin.pdf',width=9.3,height = 10.7)
plot_grid(prow,legend, align = "v",axis='l',ncol=1,rel_heights =c(10,0.7))
dev.off()