setwd("~/Research projects/Sierra dendro/sierra-dendro") # Derek on Latim-GIS-S
source("R/summarize_cluster.R")

#####################################################################################
## The summarize.cluster function does four things:
## 1) Displays errors in .pos files (CooRecorder files) in the specified cluster
##      These files should be fixed at some early stage of crossdating
## 2) Produces 2 plots of correlations among cores:
##      2a) A plot of all cores in the cluster
##      2b) A plot of all the cores included in the reference chronology
## 3) Updates the file Dendro/Ancillary data/cluster_summary_XX.csv
##      This file is intended to help prioritize cores to inspect when crossdating
## 4) Displays the timespan of the reference chronology
##      Cores cannot be cross-dated to years later than the one displayed
####################################################################################

list.plots()

summarize.cluster("SS81B",type="plot") # 

##SS81B RR33