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

plist=split(data,data$Cohort)
plist=lapply(plist,function(x){
  x$wgii_type=discretize(x$wgii,method = 'cluster',categories = 2)
  x$wgii_type=factor(x$wgii_type)
  levels(x$wgii_type) = c('low','high')
  #x$fraction_hi=length(which(x$wgii_type=='high'))/dim(x)[1]
  x$fraction_hi=length(which(x$Genome_doublings!=0))/dim(x)[1]
  #x$fraction_hi=median(x$wgii)
  x$fraction_hi=format(round(x$fraction_hi, digits=3), nsmall = 3) 
  return(x)
})
data_w=data.frame(do.call(rbind,plist))
#pro_hi=unlist(lapply(plist, function(x){length(which(x$wgii_type=='high'))/dim(x)[1]}))
#pro_hi=unlist(lapply(plist, function(x){length(which(x$Genome_doublings!=0))/dim(x)[1]}))
pro_hi=unlist(lapply(plist, function(x){median(x$wgii)}))
kk=names(pro_hi)[order(pro_hi)]
data_w$Cohort_w=factor(as.character(data_w$Cohort),levels = kk)

bar_wgii = ggplot(data_w, aes(x = wgii, fill = Genome_doublings)) +
  geom_histogram() + scale_fill_manual(name='WGD',values = unname(gcolor[c(1, 5, 6)])) + theme_cowplot()+
  theme(axis.text.x = element_text(size =15, angle = -270,vjust=0.5),legend.text = element_text(size =12),
        legend.title = element_text(size=12),legend.key.size = unit(0.2, "lines"),
        axis.title.x = element_text(size =15),axis.title.y= element_blank(),
        panel.background = element_rect(fill="white",colour="white"),
        axis.text.y = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_flip()+guides(fill=guide_legend(title="WGD"))
bar_wgii
pwgii = ggplot(data_w, aes(x =fraction_hi , y = wgii))
pwgii=pwgii + geom_quasirandom(aes(color = Genome_doublings),method = "smiley",varwidth = F,size = 0.01)+scale_color_manual(name='WGD',values=unname(gcolor[c(1,5,6)]))
pwgii=pwgii+ylab('WGII score')

pwgii = pwgii+ facet_grid(. ~ Cohort_w, scales = 'free',space = 'fix',switch = "x") + theme_cowplot()+theme(
  axis.text.x = element_text(size =15, angle = 90,vjust=0.5),axis.title.y = element_text(size =15),
  axis.title.x = element_blank(),
  axis.ticks.x = element_blank(),
  legend.position = "none",
  panel.background = element_rect(fill="white",colour="white")
)


pwgii = pwgii +theme(strip.background = element_rect(fill = 'white'),strip.placement = "outside",
                       strip.text.x =element_text(size =15, angle =90,vjust=0.5),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pwgii
wgiiplot=plot_grid(pwgii,bar_wgii+theme(legend.position =c(0.4,0.65)), align = "h", axis = "b",ncol=2,rel_widths = c(3, .25))

print(wgiiplot)
pdf('plot/pan_can_wgii.pdf',width = 13,height = 4.5)
print(wgiiplot)
dev.off()