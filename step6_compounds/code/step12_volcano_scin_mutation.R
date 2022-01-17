rm(list=ls())
library(ggplot2)
library(ggsignif)
library(plyr)
library(ggrepel)
library(cowplot)
library(data.table)
setwd("/data/zhang/pancan_cin/step6_compounds/")
#Load Table S2 from data from Taylor et al, Cancer Cell 2018 (https://www.sciencedirect.com/science/article/pii/S1535610818301119)
load('data/ccle_cin_score.RData')
Final_table=scores
Final_table$rates=1/Final_table$Doubling.Time.Calculated.hrs
Final_table$Cohort=Final_table$tcga_code
Final_table=Final_table[which(!is.na(Final_table$Cohort)),]
Final_table=Final_table[which(Final_table$tcga_code!="UNABLE TO CLASSIFY"),]
Mut_table <- fread("data/CCLE_DepMap_18q3_maf_20180718.txt")
Mut_table$Sample=substr(Mut_table$Tumor_Sample_Barcode,1,15)
Results <- read.delim("table/scin_allMutations_linear_model_results.txt")
Results <- Results[order(Results$Pvalue),]
Results <- Results[!Results$Gene %in% ".",]
Results$text <- "no"
Results$text[c(1:10,1235,1424)] <- "yes"
Results$text <- as.factor(Results$text)
Results$Significant <- ifelse(Results$FDR < 0.05 & Results$MUTminusWT_diff > 0, "Pos", "Not Sig")
Results$Significant[Results$FDR < 0.05 & Results$MUTminusWT_diff < 0] <- "Neg"

Results$Significant2 <- ifelse(Results$FDR < 0.05 & Results$Coef > 0, "Pos", "Not Sig")
Results$Significant2[Results$FDR < 0.05 & Results$Coef < 0] <- "Neg"

pdf("plot/lm_scin_mutation.pdf")
plot <- ggplot(Results, aes(x = Coef, y = -log10(Pvalue))) +
  geom_point(aes(color = Significant2)) +
  xlab("Coefficient (linear model)") +
  ylab("-log10(p-value)") +
  scale_color_manual(values = c(tropiccolor[1], tropiccolor[1], tropiccolor[3])) +
  theme_bw(base_size = 12) + theme(axis.text.x=element_text(size=13, colour="black"), axis.text.y=element_text(size=13, colour="black"), axis.title=element_text(size=16, colour="black"), legend.key = element_rect(fill="White"), legend.text = element_text(size=13), legend.title = element_text(size=14), legend.position = "bottom") +
  geom_text_repel(
    data = subset(Results, text == "yes"),
    aes(label = Gene),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  ) +
  guides(color=guide_legend(title="Significant (FDR < 0.05)"))
print(plot)
dev.off()
Final_expression_table_tum=Final_table
Final_expression_table_tum$TP53 <- "wt"
Final_expression_table_tum$TP53[Final_expression_table_tum$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% "TP53"]] <- "mut"
Final_expression_table_tum$TP53 <- factor(Final_expression_table_tum$TP53, levels = c("wt", "mut"))

Results_TP53 <- data.frame()

for(i in 1:length(unique(Final_expression_table_tum$Cohort))){
  
  cohort <- paste(unique(Final_expression_table_tum$Cohort)[i])
  
  if(length(which(Final_expression_table_tum$TP53 %in% "mut" & Final_expression_table_tum$Cohort %in% cohort)) >=20 & length(which(Final_expression_table_tum$TP53 %in% "mut" & Final_expression_table_tum$Cohort %in% cohort)) >=20){
    
    Final_expression_table_tum_cohort <- Final_expression_table_tum[Final_expression_table_tum$Cohort %in% cohort,]
    
    fit <- lm(scin ~ TP53, Final_expression_table_tum_cohort)
    
    pvalue = summary(fit)$coefficients[2,4]  
    coef = summary(fit)$coefficients[2,1]
    
    wt <- mean(Final_expression_table_tum_cohort$scin[Final_expression_table_tum_cohort$TP53 %in% "wt"])
    mut <- mean(Final_expression_table_tum_cohort$scin[Final_expression_table_tum_cohort$TP53 %in% "mut"])
    
    WT_samples <- length(which(Final_expression_table_tum$TP53 %in% "wt" & Final_expression_table_tum$Cohort %in% cohort))
    MUT_samples <- length(which(Final_expression_table_tum$TP53 %in% "mut" & Final_expression_table_tum$Cohort %in% cohort))
    
    Results_TP53=rbind(Results_TP53, data.frame(cohort, coef, pvalue, wt, mut, mut-wt, WT_samples, MUT_samples))
    
  }
}

Results_TP53$FDR <- p.adjust(Results_TP53$pvalue, method = "fdr")
Results_TP53 <- Results_TP53[,c(1:3,9,4:8)]
names(Results_TP53) <- c("Cohort", "Coef", "Pvalue", "FDR", "WTmean", "MUTmean", "MUTminusWT_diff", "WT_samples", "MUT_samples")


### plot for all cohorts

Results_TP53 <- Results_TP53[order(Results_TP53$Coef, decreasing = T),]

order_names=as.character(Results_TP53$Cohort)
Results_TP53$Cohort <- factor(Results_TP53$Cohort, levels=order_names)

Results_TP53$Significant <- ifelse(Results_TP53$FDR < 0.05 & Results_TP53$Coef > 0, "Pos", "Not Sig")
Results_TP53$Significant[Results_TP53$FDR < 0.05 & Results_TP53$Coef < 0] <- "Neg"
Results_TP53$Significant=factor(Results_TP53$Significant,levels=c('Not Sig','Pos'))
Plot <- ggplot(Results_TP53, aes(Results_TP53$Cohort, Results_TP53$Coef)) +
  geom_bar(width=0.9, stat = "identity", aes(fill=factor(Significant)), colour=NA) + 
  labs(x = "", y = "Coefficient (linear model)") +
  theme_bw(base_size = 12) + theme(axis.text.y=element_text(size=8, colour="black"), axis.text.x=element_text(size=4, angle=45, hjust=1, colour="black"), axis.title=element_text(size=9, colour="black"), legend.key = element_rect(fill="White"), legend.text = element_text(size=13), legend.title = element_text(size=14), legend.position = "bottom") +
  guides(fill=FALSE) +
  scale_fill_manual(values = c(tropiccolor[1], tropiccolor[3])) 

name= paste("plot/scin_tp53_mut.pdf")
pdf(name, height = 2, width = 4)
print(Plot)
dev.off()


# Driver mutations from Catalog of Validated Oncogenic Mutations (https://www.cancergenomeinterpreter.org/mutations)
Driver <- read.delim("../step3_mutation/data/catalog_of_validated_oncogenic_mutations.tsv")

## which ones are driver?
Mut_table$ID <- paste0("chr", Mut_table$Chromosome, ":g.", Mut_table$Start_position, Mut_table$Reference_Allele, ">", Mut_table$Tumor_Seq_Allele1)
Mut_table$driver <- Mut_table$ID %in% Driver$gdna

Results <- data.frame()

for(gene in as.character(unique(Mut_table$Hugo_Symbol[Mut_table$driver %in% "TRUE"]))){
  # at least 10 samples per group
  if(length(which(Final_expression_table_tum$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene & Mut_table$driver %in% "TRUE"])) >=10 & length(which(!Final_expression_table_tum$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene & Mut_table$driver %in% "TRUE"])) >=10){
    
    Final_expression_table_tum$Mutation_selected <- "WT"
    Final_expression_table_tum$Mutation_selected[Final_expression_table_tum$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene]] <- "passenger"
    Final_expression_table_tum$Mutation_selected[Final_expression_table_tum$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene & Mut_table$driver %in% "TRUE"]] <- "driver"
    
    #only driver vs wt
    Final_expression_table_tum <- Final_expression_table_tum[Final_expression_table_tum$Mutation_selected %in% c("WT", "driver"),]
    Final_expression_table_tum$Mutation_selected <- droplevels(factor(Final_expression_table_tum$Mutation_selected))
    Final_expression_table_tum$Mutation_selected <- factor(Final_expression_table_tum$Mutation_selected, levels=c("WT", "driver"))
    
    fit <- lm(scin ~ Mutation_selected + Cohort, Final_expression_table_tum)
    
    pvalue = summary(fit)$coefficients[2,4]  
    coef = summary(fit)$coefficients[2,1]
    
    wt <- mean(Final_expression_table_tum$scin[Final_expression_table_tum$Mutation_selected %in% "WT"])
    mut <- mean(Final_expression_table_tum$scin[Final_expression_table_tum$Mutation_selected %in% "driver"])
    
    WT_samples <- length(which(!Final_expression_table_tum$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene]))
    MUT_samples <- length(which(Final_expression_table_tum$Sample %in% Mut_table$Sample[Mut_table$Hugo_Symbol %in% gene]))
    
    Results=rbind(Results, data.frame(gene, coef, pvalue, wt, mut, mut-wt, WT_samples, MUT_samples))
    
  }
}

Results$FDR <- p.adjust(Results$pvalue, method = "fdr")
Results <- Results[order(Results$pvalue),c(1:3,9,4:8)]
names(Results) <- c("Gene", "Coef", "Pvalue", "FDR", "WTmean", "MUT_driver_mean", "MUTminusWT_diff", "WT_samples", "MUT_driver_samples")
write.table(Results, "table/scin_DriverMutations_linear_model_results.txt", quote=F, sep="\t", row.names = F)


## plot

Results <- read.delim("table/scin_DriverMutations_linear_model_results.txt")
Results <- Results[order(Results$Pvalue),]
Results <- Results[!Results$Gene %in% ".",]
Results$text <- "no"
Results$text[1:3] <- "yes"
Results$text <- as.factor(Results$text)
Results$Significant <- ifelse(Results$FDR < 0.05 & Results$MUTminusWT_diff > 0, "Pos", "Not Sig")
Results$Significant[Results$FDR < 0.05 & Results$MUTminusWT_diff < 0] <- "Neg"

Results$Significant2 <- ifelse(Results$FDR < 0.05 & Results$Coef > 0, "Pos", "Not Sig")
Results$Significant2[Results$FDR < 0.05 & Results$Coef < 0] <- "Neg"

pdf("plot/scin_driver_mut.pdf", width = 4.5, height = 5)

plot <- ggplot(Results, aes(x = Coef, y = -log10(Pvalue))) +
  geom_point(aes(color = Significant2), size=2) +
  xlab("Coefficient (linear model)") +
  ylab("-log10(p-value)") +
  scale_color_manual(values = c(tropiccolor[1], tropiccolor[1], "red")) +
  theme_bw(base_size = 11) + theme(axis.text=element_text(size=10, colour="black"), axis.title=element_text(size=13, colour="black"), legend.key = element_rect(fill="White"), legend.text = element_text(size=11), legend.title = element_text(size=12), legend.position = "bottom") +
  geom_text_repel(
    data = subset(Results, text == "yes"),
    aes(label = Gene),
    size = 4,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.2, "lines")
  ) +
  guides(color=guide_legend(title="Significant (FDR < 0.05)"))

print(plot)

dev.off()


