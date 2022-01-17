
### This code was developed by Hailei Zhang and Jaegil Kim from Broad GDAC group 
#### this function is used to read the feature table and obtain the copy number and mutation info and then do preprocess the matrix, 
#####calculate the fisher exact test and then draw the covarance figure

### featuretable: the input file of *.samplefeatures.txt which was created from Aggregate_AnalysisFeatures
### tumor: cancer type
### cut.arm: the threshod of q value for the copy number arm level event 
### cut.focal: the threshod of q value for the copy number focal event 
### fold.arm.gain: the threshold of copy number gain (log2 ratio) for arm level event
### fold.arm.loss: the threshold of copy number loss (log2 ratio) for arm level event
### fold.focal.gain: the threshold of copy number gain (log2 ratio) for focal event
### fold.focal.loss: the threshold of copy number loss (log2 ratio) for focal event
### maxevent: the nimimal sample number of events


createcorfigure <- function(feature.file,tumor,cut.arm=0.1,cut.focal=0.1,fold.arm.gain=0.4,fold.arm.loss=-0.4,maxevent=5,fold.focal.gain=0.4,fold.focal.loss=-0.4){


##################################################################################
##################################### figure 1 muation cooccurence ################
##################################################################################

source("./script/plot_mutation_cooccurrence.R")
event.mutation <- createmutationfigure(feature.file="./CESC-TP.samplefeatures.txt",tumor="CESC-TP")

###################################################################################################
################################# figure 2 coocurrece of copy number arm level events ##############
###################################################################################################
source("./script/plot_cn_arm_cooccurrence.R")
event.arm <- creatCNarmfigure(feature.file="./CESC-TP.samplefeatures.txt",tumor="CESC-TP",cut.arm=0.1,fold.arm.gain=0.4,fold.arm.loss=-0.4,maxevent=5)

###################################################################################################
################################# figure 3 coocurrece of copy number focal level events ##############
###################################################################################################
source("./script/plot_cn_focal_cooccurrence.R")
event.focal <- creatCNfocalfigure(feature.file="./CESC-TP.samplefeatures.txt",tumor="CESC-TP",cut.focal=0.1,fold.focal.gain=0.4,fold.focal.loss=-0.4,maxevent=5)
  
###################################################################################################
################################# figure 4 coocurrece of all genomic events ##############
###################################################################################################

source("./script/plot_all_event.R")
creatalleventsfigure(event.mutation,event.arm,event.focal,tumor)
}


createcorfigure(feature.file="./CESC-TP.samplefeatures.txt",tumor="CESC-TP",cut.arm=0.1,cut.focal=0.1,fold.arm.gain=0.4,fold.arm.loss=-0.4,maxevent=5,fold.focal.gain=0.4,fold.focal.loss=-0.4)

