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
library(arules)
gcolor=yarrr::piratepal("appletv") 
setwd("/data/zhang/pancan_cin/step1_calculate_cin/")
load('data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer','Primary Solid Tumor'),]
Final_table$Sample <- substr(Final_table$Sample, 1, 15)
rownames(Final_table)=Final_table$Sample
Final_table$Genome_doublings <- as.factor(Final_table$`Genome doublings`)
data=Final_table[which(complete.cases(Final_table$Genome_doublings)),]
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
pscin =pscin+
  facet_grid(. ~ Cohort, scales = 'free', space = 'fix',switch = "x")+ theme_cowplot() + theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),axis.title.y = element_text(size =15),
    axis.ticks.x = element_blank(),panel.background = element_rect(fill="white",colour="white")
  )
pscin= pscin + theme(strip.background = element_rect(fill = 'white'),strip.placement = "outside",
                     strip.text.x =element_text(size =15, angle = -270,vjust=0.1),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank())
bar_scin=ggplot(data, aes(x=scin,fill=Genome_doublings) ) +
  geom_histogram() + scale_fill_manual(values = unname(gcolor[c(1, 5, 6)])) + theme_cowplot()+
  theme(axis.text.x = element_text(size =15, angle =-270,vjust=0.5),legend.text = element_text(size =12),
        legend.title = element_text(size=12),legend.key.size = unit(0.2, "lines"),
        axis.title.x = element_text(size =15),axis.title.y= element_blank(),
        panel.background = element_rect(fill="white",colour="white"),
        axis.text.y = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_flip()+guides(fill=guide_legend(title="WGD"))
scinplot=plot_grid(pscin,bar_scin+theme(legend.position =c(0.1,0.8)), align = "h", axis = "b",ncol=2,rel_widths = c(3, .25))
print(scinplot)
pdf('plot/pan_can_scin.pdf',width = 13,height = 3)
print(scinplot)
dev.off()