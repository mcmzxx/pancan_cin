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
gcolor=yarrr::piratepal("appletv") 
setwd("/data/zhang/pancan_cin/step6_compounds/")
load('data/ccle_cin_score.RData')
scores$Genome_doublings=as.factor(scores$Genome.doublings)
data=scores[which(!is.na(scores$tcga_code)),]
data$Cohort=as.character(data$tcga_code)
data=data[which(complete.cases(data$Genome_doublings)),]
data=data[with(data,order(Cohort,scin)),]
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
dataf=ddply(data, .(Cohort), summarise,py=median(scin),px=median(NO))
dataf$Cohort=as.character(dataf$Cohort)
ii=dataf$Cohort[order(dataf$py)]
dataf$Cohort=factor(dataf$Cohort,levels = ii)
data$Cohort=factor(data$Cohort,levels = ii)
pscin = ggplot(data, aes(x = NO, y = scin))+geom_boxplot( coef = 0,outlier.shape=NA,color='red',lwd=0.1) + geom_point(size = 0.01)
pscin = pscin + geom_point(
  aes(
    x = px,
    y = py
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
bar_scin=ggplot(data, aes(x=scin,fill=Genome_doublings) ) +
  geom_histogram() + scale_fill_manual(values = unname(gcolor[c(1, 5, 6)])) + theme_cowplot()+
  theme(axis.text.x = element_text(size = 10, angle = -270,vjust=0.5),legend.text = element_text(size =10),
        axis.title.x = element_blank(),axis.title.y= element_blank(),
        panel.background = element_rect(fill="white",colour="white"),
        axis.text.y = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_flip()
bar_wgii = ggplot(data, aes(x = wgii, fill = Genome_doublings)) +
  geom_histogram() + scale_fill_manual(name='WGD',values = unname(gcolor[c(1, 5, 6)])) + theme_cowplot()+
  theme(axis.text.x = element_text(size = 10, angle = -270,vjust=0.5),legend.text = element_text(size = 10),legend.title = element_text(size=10),axis.title.x = element_blank(),axis.title.y= element_blank(),
        panel.background = element_rect(fill="white",colour="white"),
        axis.text.y = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_flip()
pwgii = ggplot(data, aes(x = Cohort, y = wgii))
pwgii=pwgii + geom_quasirandom(aes(color = Genome_doublings),method = "smiley",varwidth = F,size = 0.01)+scale_color_manual(name='WGD',values=unname(gcolor[c(1,5,6)]))
pwgii=pwgii+ylab('WGII score')
pwgii = pwgii+ facet_grid(. ~ Cohort, scales = 'free', space = 'fix') + theme_cowplot()+theme(
  axis.text.x = element_text(size = 10, angle = -270,vjust=0.5),axis.title.y = element_text(size =10),
  axis.title.x = element_blank(),
  axis.ticks.x = element_blank(),
  panel.background = element_rect(fill="white",colour="white")
)
pwgii = pwgii + theme(strip.background = element_rect(fill = 'white'),strip.text.x =element_blank(),legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())

legend <- get_legend(
  bar_wgii + theme(legend.box.margin = margin(0, -3, 0, 0))
)
wgiiplot=plot_grid(pwgii,bar_wgii+theme(legend.position = "none"), align = "h", axis = "b",ncol=2,rel_widths = c(3, .2))
scinplot=plot_grid(pscin,bar_scin+theme(legend.position = "none"), align = "h", axis = "b",ncol=2,rel_widths = c(3, .2))
scinplot=scinplot+theme(plot.margin=margin(b=-0.5,unit="cm"))
prow=plot_grid(scinplot, wgiiplot, ncol=1, align = "v", axis = "l",rel_heights = c(2, 3))
pdf('plot/pan_can_cin_ccle.pdf',width = 13,height = 5)
plot_grid(prow,legend,rel_widths = c(3, .15))
dev.off()