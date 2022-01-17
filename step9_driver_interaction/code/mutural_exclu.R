rm(list=ls())
library(ComplexHeatmap)
library(data.table)
onkotablebak=fread('/data/zhang/pancan_cin/step9_driver_interaction/data/onkokb.csv')
onkotable=data.frame(apply(onkotablebak[,-1], 2, function(x){ifelse(x==FALSE,0,1)}),stringsAsFactors = F)
onkotable[is.na(onkotable)]=0
rownames(onkotable)=onkotablebak$SAMPLE_BARCODE
cg=colnames(onkotable)
alca=unlist(lapply(cg,function(x){strsplit(x,'\\.')[[1]][1]}))
alcg=unlist(lapply(cg,function(x){strsplit(x,'\\.')[[1]][2]}))

tmpmat=onkotable
tmpmat[,]=""
for(i in cg){
tmpa=strsplit(i,'\\.')[[1]][1]
tmpg=strsplit(i,'\\.')[[1]][2]
tmpmat[,i]=ifelse(onkotable[,i]==1,tmpa,'')
}
resmat=list()
for(i in unique(alcg)){
  idx=which(alcg==i)
  if(length(idx)>1)
  {
    tmpres=tmpmat[,idx]
    resmat[[i]]=apply(tmpres, 1, function(x){paste0(as.character(x)[which(x!='')],collapse = ';')})
  }
  else{
    resmat[[i]]=tmpres
  }
}
resdata=do.call(cbind,resmat)
finaldata=apply(resdata, c(1,2),function(x){tmpstr=x;
if(nchar(tmpstr)>1)
  {return( paste0(tmpstr,';'))}
else{
  return(tmpstr)
}})
mat = read.table(system.file("extdata", package = "ComplexHeatmap", 
                             "tcga_lung_adenocarcinoma_provisional_ras_raf_mek_jnk_signalling.txt"), 
                 header = TRUE, stringsAsFactors = FALSE, sep = "\t")
mat[is.na(mat)] = ""
rownames(mat) = mat[, 1]
mat = mat[, -1]
mat=  mat[, -ncol(mat)]
mat = t(as.matrix(mat))
mat[1:3, 1:3]
col = c("FUSION" = "blue", "AMP" = "red", "MUT" = "#008000","EPISIL"='cyan')
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
 FUSION = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["FUSION"], col = NA))
  },
  # bug red
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["AMP"], col = NA))
  },
  # small green
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["MUT"], col = NA))
  },
 DEL = function(x, y, w, h) {
   grid.rect(x, y, w-unit(0.3, "mm"), h*0.23, 
             gp = gpar(fill = col["DEL"], col = NA))
 },
 
 EPISIL = function(x, y, w, h) {
   grid.rect(x, y, w-unit(0.5, "mm"), h*0.63, 
             gp = gpar(fill = col["EPISIL"], col = NA))
 }
)
column_title = "OncoPrint for TCGA Lung Adenocarcinoma, genes in Ras Raf MEK JNK signalling"
heatmap_legend_param = list(title = "Alternations", at = c("HOMDEL", "AMP", "MUT"),labels = c("Deep deletion", "Amplification", "Mutation"))

oncoPrint(mat,
          alter_fun = alter_fun, col = col, 
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                                             foo1 = 1:172,
                                             bar1 = anno_points(1:172)),
          left_annotation = rowAnnotation(foo2 = 1:26),
          right_annotation = rowAnnotation(bar2 = anno_barplot(1:26)),
          column_title = column_title, heatmap_legend_param = heatmap_legend_param)