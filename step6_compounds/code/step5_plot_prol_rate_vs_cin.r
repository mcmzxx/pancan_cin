rm(list=ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
setwd("/data/zhang/pancan_cin/step6_compounds/")
#Load Table S2 from data from Taylor et al, Cancer Cell 2018 (https://www.sciencedirect.com/science/article/pii/S1535610818301119)
load('data/ccle_cin_score.RData')
Final_table=scores
Final_table$rates=1/Final_table$Doubling.Time.Calculated.hrs

res1=cor.test(Final_table$wgii, Final_table$rates, method = "spearman")
tmp=names(res1)
res1=lapply(tmp, function(x){res1[[x]]})
names(res1)=tmp
res2=cor.test(Final_table$ncin, Final_table$rates, method = "spearman")
tmp=names(res2)
res2=lapply(tmp, function(x){res2[[x]]})
names(res2)=tmp
res3=cor.test(Final_table$scin, Final_table$rates, method = "spearman")
tmp=names(res3)
res3=lapply(tmp, function(x){res3[[x]]})
names(res3)=tmp
sink('table/wgii_vs_proliferation_cor.txt')
print(res1)
sink()

sink('table/ncin_vs_proliferation_cor.txt')
print(res2)
sink()

sink('table/scin_vs_proliferation_cor.txt')
print(res3)
sink()

scatter1 <- ggplot(Final_table, aes(x=wgii, y=rates)) + #scale_x_log10() +
  #scale_y_log10() +
  xlab("wgii score") + 
  ylab(paste("proliferation rates [1/h]")) + 
  ggtitle(paste0("")) +
  geom_point(colour="grey70", size=1, alpha=1) +
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(low="white",high="grey20") +
  scale_alpha(range = c(0.1,0.3)) +
  guides(alpha="none", fill="none") +
  theme(panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), axis.text.y=element_text(size=13, colour="black"), axis.text.x=element_text(size=13, colour="black"), axis.title=element_text(size=16, colour="black"), plot.title = element_text(size=17, hjust = 0.5, colour="black"), legend.key=element_rect(fill=NA), legend.text = element_text(size=13, colour="black"), legend.title = element_text(size=16, colour="black"), axis.title.x=element_text(margin=margin(10,0,0,0)))


scatter2 <- ggplot(Final_table, aes(x=ncin, y=rates))+ #scale_x_log10() +
  #scale_y_log10() +
  xlab("ncin score") + 
  ylab(paste("proliferation rates [1/h]")) + 
  ggtitle(paste0("")) +
  geom_point(colour="grey70", size=1, alpha=1) +
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(low="white",high="grey20") +
  scale_alpha(range = c(0.1,0.3)) +
  guides(alpha="none", fill="none") +
  theme(panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), axis.text.y=element_text(size=13, colour="black"), axis.text.x=element_text(size=13, colour="black"), axis.title=element_text(size=16, colour="black"), plot.title = element_text(size=17, hjust = 0.5, colour="black"), legend.key=element_rect(fill=NA), legend.text = element_text(size=13, colour="black"), legend.title = element_text(size=16, colour="black"), axis.title.x=element_text(margin=margin(10,0,0,0)))

scatter3 <- ggplot(Final_table, aes(x=scin, y=rates))+ #scale_x_log10() +
  #scale_y_log10() +
  xlab("scin score") + 
  ylab(paste("proliferation rates [1/h]")) + 
  ggtitle(paste0("")) +
  geom_point(colour="grey70", size=1, alpha=1) +
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(low="white",high="grey20") +
  scale_alpha(range = c(0.1,0.3)) +
  guides(alpha="none", fill="none") +
  theme(panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), axis.text.y=element_text(size=13, colour="black"), axis.text.x=element_text(size=13, colour="black"), axis.title=element_text(size=16, colour="black"), plot.title = element_text(size=17, hjust = 0.5, colour="black"), legend.key=element_rect(fill=NA), legend.text = element_text(size=13, colour="black"), legend.title = element_text(size=16, colour="black"), axis.title.x=element_text(margin=margin(10,0,0,0)))

pdf("plot/cin_score_vs_proliferation.pdf", height  = 4,width=12)
plot_grid(scatter1,scatter2,scatter3, labels = c("A", "B","C"),ncol=3)
dev.off()

