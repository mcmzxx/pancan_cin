rm(list=ls())
setwd("/data/zhang/pancan_cin/step1_calculate_cin/")
TCGA_cohorts_table <- read.delim("../TCGA_cohorts.txt")
TCGA_cohorts <- TCGA_cohorts_table$Cohort
system("mkdir Firebrowse_download")
setwd("Firebrowse_download")
Table_final_results <- data.frame()
for(i in 1:length(TCGA_cohorts)) {
  cohort_TCGA <- TCGA_cohorts[i]
  if(length(grep(paste("gdac.broadinstitute.org_", cohort_TCGA,".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0", sep=""), list.files())) ==0 ){
    system(paste("wget http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/", cohort_TCGA,"/20160128/gdac.broadinstitute.org_", cohort_TCGA,".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0.tar.gz", sep=""))
    system(paste("tar -zxvf gdac.broadinstitute.org_", cohort_TCGA,".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0.tar.gz", sep=""))
    system(paste("rm gdac.broadinstitute.org_", cohort_TCGA,".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0.tar.gz", sep=""))
  }
  if(length(grep(paste("wget http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/", cohort_TCGA, "/20160128/gdac.broadinstitute.org_", cohort_TCGA,".Merge_Clinical.Level_1.2016012800.0.0.tar.gz", sep=""), list.files())) ==0 ){
    system(paste("wget http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/", cohort_TCGA, "/20160128/gdac.broadinstitute.org_", cohort_TCGA,".Merge_Clinical.Level_1.2016012800.0.0.tar.gz", sep=""))
    system(paste("tar -zxvf gdac.broadinstitute.org_", cohort_TCGA,".Merge_Clinical.Level_1.2016012800.0.0.tar.gz", sep=""))
    system(paste("rm gdac.broadinstitute.org_", cohort_TCGA,".Merge_Clinical.Level_1.2016012800.0.0.tar.gz", sep=""))
    
  }
  Expression_table <- read.delim(paste("gdac.broadinstitute.org_", cohort_TCGA,".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0/", cohort_TCGA,".rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt", sep=""), na.strings = c("null", "", " ", "NA", "na"))
  Expression_table_2 <- Expression_table[-c(1:30),c(1,grep("raw_count", as.character(unlist(Expression_table[1, ]))))]
  Expression_table_2$Hybridization.REF = sapply(strsplit(as.character(Expression_table_2$Hybridization.REF),"[[:punct:]]"), `[`, 1)
  
  # mean across repeated genes
  Expression_table_3 <- as.data.frame(Expression_table_2)
  Expression_table_3[,2:ncol(Expression_table_3)] <- lapply(Expression_table_3[,2:ncol(Expression_table_3)], function(x) as.numeric(as.character(x)))
  
  Expression_table_3 <- aggregate(Expression_table_3[,2:ncol(Expression_table_3)],by=list(name=Expression_table_3$Hybridization.REF),data=Expression_table_3,FUN=mean)
  
  rownames(Expression_table_3)=Expression_table_3[,1]
  Expression_table_3=Expression_table_3[,-1]
  Expression_table_3=as.data.frame(t(Expression_table_3))
  
  Final_Expression_table=Expression_table_3
  Final_Expression_table$Sample=substr(rownames(Final_Expression_table),1,15)
  rownames(Final_Expression_table)=substr(rownames(Final_Expression_table),1,15)
  
  Final_Expression_table$sample_type <- NA
  Final_Expression_table$sample_type[grep(".01$", rownames(Final_Expression_table))] <- "Tumor"
  Final_Expression_table$sample_type[grep(".02$", rownames(Final_Expression_table))] <- "Tumor"
  Final_Expression_table$sample_type[grep(".03$", rownames(Final_Expression_table))] <- "Tumor"
  Final_Expression_table$sample_type[grep(".04$", rownames(Final_Expression_table))] <- "Tumor"
  Final_Expression_table$sample_type[grep(".05$", rownames(Final_Expression_table))] <- "Tumor"
  Final_Expression_table$sample_type[grep(".06$", rownames(Final_Expression_table))] <- "Tumor"
  Final_Expression_table$sample_type[grep(".07$", rownames(Final_Expression_table))] <- "Tumor"
  Final_Expression_table$sample_type[grep(".08$", rownames(Final_Expression_table))] <- "Tumor"
  Final_Expression_table$sample_type[grep(".09$", rownames(Final_Expression_table))] <- "Tumor"
  Final_Expression_table$sample_type[grep(".10$", rownames(Final_Expression_table))] <- "Normal"
  Final_Expression_table$sample_type[grep(".11$", rownames(Final_Expression_table))] <- "Normal"
  Final_Expression_table$sample_type[grep(".12$", rownames(Final_Expression_table))] <- "Normal"
  Final_Expression_table$sample_type[grep(".13$", rownames(Final_Expression_table))] <- "Normal"
  Final_Expression_table$sample_type[grep(".14$", rownames(Final_Expression_table))] <- "Normal"
  
  Final_Expression_table$sample_type_detail <- NA
  Final_Expression_table$sample_type_detail[grep(".01$", rownames(Final_Expression_table))] <- "Primary Solid Tumor"
  Final_Expression_table$sample_type_detail[grep(".02$", rownames(Final_Expression_table))] <- "Recurrent Solid Tumor"
  Final_Expression_table$sample_type_detail[grep(".03$", rownames(Final_Expression_table))] <- "Primary Blood Derived Cancer"
  Final_Expression_table$sample_type_detail[grep(".04$", rownames(Final_Expression_table))] <- "Recurrent Blood Derived Cancer"
  Final_Expression_table$sample_type_detail[grep(".05$", rownames(Final_Expression_table))] <- "Additional - New Primary"
  Final_Expression_table$sample_type_detail[grep(".06$", rownames(Final_Expression_table))] <- "Metastatic"
  Final_Expression_table$sample_type_detail[grep(".07$", rownames(Final_Expression_table))] <- "Additional Metastatic"
  Final_Expression_table$sample_type_detail[grep(".08$", rownames(Final_Expression_table))] <- "Human Tumor Original Cells"
  Final_Expression_table$sample_type_detail[grep(".09$", rownames(Final_Expression_table))] <- "Primary Blood Derived Cancer - Bone Marrow"
  Final_Expression_table$sample_type_detail[grep(".10$", rownames(Final_Expression_table))] <- "Blood Derived Normal"
  Final_Expression_table$sample_type_detail[grep(".11$", rownames(Final_Expression_table))] <- "Solid Tissue Normal"
  Final_Expression_table$sample_type_detail[grep(".12$", rownames(Final_Expression_table))] <- "Buccal Cell Normal"
  Final_Expression_table$sample_type_detail[grep(".13$", rownames(Final_Expression_table))] <- "EBV Immortalized Normal"
  Final_Expression_table$sample_type_detail[grep(".14$", rownames(Final_Expression_table))] <- "Bone Marrow Normal"
  
  Final_Expression_table$Cohort <- cohort_TCGA
  
  Table_final_results=rbind(Table_final_results, Final_Expression_table)
  
  print(cohort_TCGA)
  
} # TCGA cohort

saveRDS(Table_final_results, "../TCGA_expression_table_reads.rds")
setwd("../")

list_cohorts <- c("BLCA", "BRCA", "CESC",
                  "COADREAD", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", 
                  "LAML", "LGG", "LIHC", "LUAD", "LUSC","OV", "PAAD", "PRAD", "SKCM", "STAD", "TGCT", "THCA","UCEC")

library(limma)
Table_final_results <- readRDS("TCGA_expression_table_reads.rds")
expression_matrix <- t(as.matrix(Table_final_results[,1:20281]))
pdf("plot/Voom.pdf", width = 10)
voom.analysis=voom(expression_matrix, plot=TRUE, normalize.method = "quantile")
dev.off()
Voom_table <- merge(as.data.frame(t(voom.analysis$E)), Table_final_results[,20282:20285], by=0)
rownames(Voom_table) <- Voom_table[,1]
Voom_table = Voom_table[,-1]
# Proliferation rates from paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5186797/pdf/fphys-07-00644.pdf
# The table used is in https://github.com/cdiener/proliferation/blob/master/results/pred_rates.csv

Prolif_score <- read.csv("TCGA_rol_rates.csv")
Prolif_score$sample_id <- gsub("-",".",Prolif_score$patient_barcode)
Prolif_score <- Prolif_score[Prolif_score$tumor %in% "TRUE",]

# summarize by mean
Prolif_score_2 <- aggregate(Prolif_score[,3],by=list(sample_id=Prolif_score$sample_id),data=Prolif_score,FUN=mean)
names(Prolif_score_2)[2] <- "rates"

Final_expression_table <- read.delim("CA20_expression_table_TCGA_CA20allCohorts.txt")
Final_expression_table_tum <- Final_expression_table[Final_expression_table$sample_type_detail %in% "Primary Solid Tumor",]
Final_expression_table_tum$sample_id <- substr(Final_expression_table_tum$Sample,1,12)

Final_table <- merge(Final_expression_table_tum, Prolif_score_2, by.x=27, by.y=1)

cor.test(Final_table$CA20, Final_table$rates, method = "spearman")

scatter <- ggplot(Final_table, aes(x=CA20, y=rates)) + 
  xlab("CA20 score") + 
  ylab(paste("Predicted proliferation rates [1/h]")) + 
  ggtitle(paste0("")) +
  geom_point(colour="grey70", size=1, alpha=1) +
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(low="white",high="grey20") +
  scale_alpha(range = c(0.1,0.3)) +
  guides(alpha="none", fill="none") +
  scale_x_continuous(breaks=seq(-100,100,10)) +
  scale_y_continuous(breaks = seq(-100,100,0.01)) +
  theme(panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), axis.text.y=element_text(size=13, colour="black"), axis.text.x=element_text(size=13, colour="black"), axis.title=element_text(size=16, colour="black"), plot.title = element_text(size=17, hjust = 0.5, colour="black"), legend.key=element_rect(fill=NA), legend.text = element_text(size=13, colour="black"), legend.title = element_text(size=16, colour="black"), axis.title.x=element_text(margin=margin(10,0,0,0)))

pdf("FigS1a.pdf", height  = 4)

print(scatter)

dev.off()


library(dplyr)
library(plyr)
library(reshape)

Final_expression_table <- read.delim("CA20_expression_table_TCGA_CA20allCohorts.txt")
Final_expression_table_tum <- Final_expression_table[Final_expression_table$sample_type %in% "Tumor",]

test_table <- ddply(Final_expression_table_tum, .(Cohort), mutate, Median_mut_n = homen(CA20))
test_table = sort_df(test_table, "Median_mut_n")
test_table$Cohort = factor(test_table$Cohort, levels = unique(test_table$Cohort))

library(ggplot2)
library(ggsignif)
library(scater)

data_summary <- function(x) {
  m <- homen(x)
  ymin <- as.numeric(quantile(x)[2])
  ymax <- as.numeric(quantile(x)[4])
  return(c(y=m,ymin=ymin,ymax=ymax))
}

plot <- ggplot(test_table, aes(x=Cohort, y=CA20)) + 
  geom_dotplot(aes(colour=Cohort, fill=Cohort), binaxis='y', stackdir='center', binwidth = .5, stackratio = .7, dotsize=1) +
  scale_y_continuous(breaks = seq(-100,100,10)) +
  labs(x="", y = "CA20 score")+
  theme(panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), axis.text.x=element_text(size=14, angle=45, hjust=1, colour="black"), axis.text.y=element_text(size=12, colour="black"), axis.title.y=element_text(size=16, colour="black", margin=margin(0,10,0,0)), legend.key = element_rect(fill="White"), legend.text = element_text(size=14), plot.title = element_text(size=17, hjust=0.5, colour="black")) +
  stat_summary(fun.data=data_summary, color="black") +
  guides(fill=F, colour=F)

pdf("Fig1b.pdf", width = 12, height = 6)

print(plot)

dev.off()


Final_expression_table$Sample <- substr(Final_expression_table$Sample, 1, 12)
Final_expression_table$Sample_ID <- rownames(Final_expression_table)

list_cohorts <- c("BLCA", "BRCA", "COADREAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC")

Expression_table_CA20 <- Final_expression_table[Final_expression_table$Cohort %in% list_cohorts,]

test_table <- ddply(Expression_table_CA20[Expression_table_CA20$sample_type %in% "Tumor",], .(Cohort), mutate, Median_mut_n = homen(CA20))
test_table = sort_df(test_table, "Median_mut_n")
test_table$Cohort = factor(test_table$Cohort, levels = unique(test_table$Cohort))

Expression_table_CA20$Cohort = factor(Expression_table_CA20$Cohort, levels = unique(test_table$Cohort))

pdf("Fig1c.pdf", width = 12, height = 6)

# geom spli violin
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

plot_sample_type_violin <- ggplot(Expression_table_CA20, aes(x=Cohort, y=CA20, fill=factor(sample_type))) + 
  geom_split_violin(width=1.13) +
  scale_y_continuous(breaks = seq(-100,100,10)) +
  labs(x="", y = "CA20 score")+
  theme(panel.background = element_rect(fill="white",colour="white"), axis.line=element_line(colour="black"), axis.text.x=element_text(size=14, angle=45, hjust=1, colour="black"), axis.text.y=element_text(size=12, colour="black"), axis.title.y=element_text(size=16, colour="black", margin=margin(0,10,0,0)), legend.key = element_rect(fill="White"), legend.text = element_text(size=14), plot.title = element_text(size=17, hjust=0.5, colour="black")) +
  guides(fill=F, colour=F) +
  scale_fill_manual(values = c("Normal"=tropiccolor[1], "Tumor"=tropiccolor[3]))

print(plot_sample_type_violin)

dev.off()


