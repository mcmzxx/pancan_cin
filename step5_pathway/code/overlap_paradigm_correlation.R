rm(list=ls())
library(ggplot2)
library(ggsignif)
library(matrixStats)
library(ggrepel)
library(cowplot)
library(data.table)
setwd("/data/zhang/pancan_cin/step5_pathway/")
table <- read.delim("../step2_genome_instability/data/TaylorCancerCell_TableS2.txt", na.strings = c("", " ", "na", "NA", "n.a.", "#N/A"))
rownames(table) <- as.character(table$Sample)
load('../step1_calculate_cin/data/pan_can_cin_score.RData')
Final_table=scores[scores$sample_type_detail %in% c('Primary Blood Derived Cancer','Primary Solid Tumor'),]
rownames(Final_table)=Final_table$Sample
sels=intersect(rownames(table),rownames(Final_table))
Final_table=Final_table[sels,]
Final_table$Sample=substr(Final_table$Sample,1,12)
rownames(Final_table)=Final_table$Sample
# til=gdata::read.xls('data/TILMap_TableS1.xlsx')
# til$Sample=as.character(til$ParticipantBarcode)
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
colnames(Final_table)[41:length(colnames(Final_table))]=ii
library(GSA)
gset=GSA.read.gmt('paradigm_pathway/paradigm_pathway/pid_120912_genesets.gmt')
names(gset$genesets)=gset$geneset.names
tcga_pathway=Final_table[,1:40]
for(i in gset$geneset.names){
  intersect_gene=intersect(gset$genesets[[i]],colnames(Final_table))
  if(length(intersect_gene)>1){
    tcga_pathway[,i]=apply(Final_table[,intersect_gene], 1, sum)
  }else if(length(intersect_gene)==1){
    tcga_pathway[,i]=Final_table[,intersect_gene]
  }
}
tmpc=unique(tcga_pathway$Cohort)
ii=colnames(tcga_pathway)[41:length(colnames(tcga_pathway))]
reslist=list()
for(j in tmpc){
  Results <- data.frame()
  tmpmat=subset(tcga_pathway,tcga_pathway$Cohort==j)
  for(i in ii){
    if(length(unique(tmpmat[,i]))>5){
      cell_type=i
      cor_s <- cor.test(tmpmat$wgii, tmpmat[,i], method = "spearman")
      Results=rbind(Results, data.frame(j,cell_type, as.numeric(cor_s$estimate), as.numeric(cor_s$p.value)))
     
    }}
  Results$FDR_s <- p.adjust(Results$as.numeric.cor_s.p.value., method = "fdr")
  names(Results) <- c("Cohort","Cell_type","Spearman_coef", "Spearman_pvalue","Spearman_FDR")
  reslist[[j]]=Results
}
Results=data.frame(do.call(rbind,reslist))
names(Results) <- c("Cohort","Cell_type","Spearman_coef", "Spearman_pvalue","Spearman_FDR")
write.table(Results, "table/merge_paradigm_pathway_corr_results_wgii.txt", quote=F, sep="\t", row.names = F)
Results2 <-Results[order(Results$Spearman_pvalue),]
Results2$text <- "no"
Results2$text[Results2$Spearman_FDR < 0.05 & Results2$Spearman_coef>0.4] <- "yes"
Results2$text <- as.factor(Results2$text)
Results2$Significant <- ifelse(Results2$Spearman_FDR < 0.05 & Results2$Spearman_coef > 0, "Pos", "Not Sig")
Results2$Significant[Results2$Spearman_FDR < 0.05 & Results2$Spearman_coef < 0] <- "Neg"
plot_drugs<- ggplot(Results2, aes(x = Spearman_coef, y = -log10(Spearman_pvalue))) +facet_wrap(~Cohort)+
  geom_point(aes(color = Significant),size=1) +
  xlab("cor_spearman_wgii_merged_paradigm_pathway_score") +
  ylab("-log10(p-value)") +
  scale_color_manual(values = c(tropiccolor[1], tropiccolor[1], tropiccolor[3])) +
  theme_bw(base_size = 12) + theme(axis.text.x=element_text(size=13, colour="black"), axis.text.y=element_text(size=13, colour="black"), axis.title=element_text(size=16, colour="black"), legend.key = element_rect(fill="White"), legend.text = element_text(size=13), legend.title = element_text(size=14), legend.position = "bottom") +
  geom_text_repel(
    data = subset(Results2, text == "yes"),
    aes(label = Cell_type),
    size = 1,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) +
  guides(color=guide_legend(title="Significant (FDR < 0.05)"))
pdf("plot/volcano_wgii_vs_merged_paradigm_pathway.pdf", width=28,height=10)
print(plot_drugs)
dev.off()
save(tcga_pathway,file = 'data/tcga_pathway.RData')
write.table(tcga_pathway, "table/tcga_pathway.txt", quote=F, sep="\t", row.names = F)



m <-
  matrix(
    nrow = length(gset$genesets),
    ncol = length(gset$genesets),
    data = NA
  )
for (i in 1:length(gset$genesets)) {
  for (j in 1:length(gset$genesets)) {
    intersection <- intersect(gset$genesets[[i]], gset$genesets[[j]])
    m[i, j] <- length(intersection)
  }
}
rownames(m) <- names(gset$genesets)
colnames(m) <- names(gset$genesets)


m <- na.omit(melt(m))

colnames(m) <- c("x", "y", "value")


g <- ggplot(m, aes(x, y))
g <- g + theme_bw() + xlab("") + ylab("")
g <- g + ggtitle('overlap gene in paradigm') + theme(plot.title = element_text(hjust = 0.5))
g <- g + geom_tile(aes(fill = value), color = 'white')
nz.data <- m[m$value > 0, ]
g <- g + geom_text(data = nz.data, aes(x = x, y = y, label = value))
g <- g + scale_fill_gradient(
    low = 'white',
    high = 'darkblue',
    space = 'Lab',
    guide = guide_colorbar(title = "# Genes")
  )
g <- g + theme(
  axis.text.x = element_text(angle = 90),
  axis.ticks = element_blank(),
  axis.line = element_blank(),
  panel.border = element_blank(),
  panel.grid.major = element_line(color = '#eeeeee')
)

pdf("plot/overlap_gene_paradigm_pathway.pdf", width=28,height=28)
print(g)
dev.off()

pdf("plot/overlap_gene_paradigm_pathway.pdf", width=28,height=28)
library(pheatmap)
library(RColorBrewer)
pheatmap::pheatmap(m,scale = "row",
                   # clustering_method = "ward.D",
                   # clustering_distance_cols = "euclidean",
                   # clustering_distance_rows = "euclidean",
                   cluster_cols = F,
                   cluster_rows = F,
                   border_color = NA,main ='overlap gene in paradigm',
                   show_rownames =T,show_colnames = T,treeheight_row = 0,fontsize=3,fontsize_row =3,treeheight_col = 0,
                   file ="/data/zhang/pancan_cin/step5_pathway/plot/overlap_gene_paradigm_pathway.pdf",
                   width =28, height =28,color = rev(colorRampPalette(brewer.pal(11,"RdBu"))(2000)))
