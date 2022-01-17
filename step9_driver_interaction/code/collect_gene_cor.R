rm(list = ls())
datalist=list()
ii=list.files('/data/zhang/pancan_cin/step5_pathway/table',full.names = T)
ii=grep('gene._corr_results_ncin.txt',ii,value = T)
datalist[['ncin']]=do.call(rbind,lapply(ii, function(x){read.table(x,header = T,  sep = '\t',stringsAsFactors = F)}))
ii=list.files('/data/zhang/pancan_cin/step5_pathway/table',full.names = T)
ii=grep('gene._corr_results_scin.txt',ii,value = T)
datalist[['scin']]=do.call(rbind,lapply(ii, function(x){read.table(x,header = T,  sep = '\t',stringsAsFactors = F)}))
ii=list.files('/data/zhang/pancan_cin/step5_pathway/table',full.names = T)
ii=grep('gene._corr_results_wgii.txt',ii,value = T)
datalist[['wgii']]=do.call(rbind,lapply(ii, function(x){read.table(x,header = T,  sep = '\t',stringsAsFactors = F)}))
save(datalist,file = '/data/zhang/pancan_cin/step9_driver_interaction/table/gene_corr.RData')


