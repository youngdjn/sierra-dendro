setwd("~/UC Davis/Research Projects/Sierra dendro/sierra-dendro") # Derek on Derek's computer
setwd("~/Research projects/Sierra dendro/sierra-dendro") # Derek on Latim-GIS-S

tree.ring.folder = "data/dendro/coorecorder-measurements" # on repository
#tree.ring.folder = "S:/FacultyData/LATIMER/LATIMERShared/DYoung Dendro/CooRecorder measurements" # on networked computer

# load the summarize.cluster() function and several other convenience functions
source("scripts/dendro/dendro-functions/summarize_cluster.R")

#####################################################################################
## The summarize.cluster() function does 5 things:
## 1) Displays errors in .pos files (CooRecorder files) in the specified cluster
##      These files should be fixed at some early stage of crossdating
## 2) Produces 2 plots of correlations among cores:
##      2a) A plot of all cores in the cluster
##      2b) A plot of all the cores included in the reference chronology
## 3) Updates the summary CSV file in "data/dendro/crossdating-summaries"
##      This file is intended to help prioritize cores to inspect when crossdating
## 4) Updates the reference chronology in "data/dendro/reference-chronologies"
## 5) Displays the timespan of the reference chronology
##      Cores cannot be cross-dated to years later than the one displayed
####################################################################################


## Example of running for a specific plot:
list.plots()
summarize.cluster(cluster = "RS58",type = "plot", name="RS58")

## Example of running for a specific group of trees:
focal.trees = c("1402","4020","1403","4021","4038","1471","1406","4002","1435","4042","1448","1474","1094","2119","2125","2123","2137","4047","1453","1465")
focal.trees = "1213"
summarize.cluster(cluster = focal.trees, type = "tree", name = "TreeGroup1") # The "name" will be used for the filename of the output reference chronology

## Example of running for a cluster:
summarize.cluster(cluster = "SH",type = "cluster",clean.ref.chron=FALSE) # clean.ref.chron (which is true by default) truncates the core when corssdating records indicate: a ring removal or addition did not make sense in the image; crossdating was limited by reference chronology length 

## Example of running for ALL CORES that we have plot data and tree ring data for
summarize.cluster(cluster = "ALL",type = "cluster") 
