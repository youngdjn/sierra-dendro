setwd("~/UC Davis/Research Projects/Sierra dendro/sierra-dendro")

library(tidyverse)


#####
##### Find cores with crossdating problems #####
#####

# no pith comment or pith comment unrecognized
# "pith hit" or "pith offset" but no measurement
# no pos files matching
# cores truncated at uncs


source("scripts/dendro/dendro-functions/read_pos_extended_bai.R")

a <- read.csv("data/dendro/crossdating-records/Crossdating progress records (Responses).csv",header=TRUE,stringsAsFactors=FALSE)
cores <- unique(a$Core.number)

chron <- open.chron(cores)

write.csv(pith.problems,"data/dendro/crossdating-records/pith_problems.csv",row.names=FALSE)
write.csv(pos.problems,"data/dendro/crossdating-records/pos_problems_2.csv",row.names=FALSE)
write.csv(unc.problems,"data/dendro/crossdating-records/unc_problems.csv",row.names=FALSE)
write.csv(dated.problems,"data/dendro/crossdating-records/dated_problems.csv",row.names=FALSE)




#####
##### Find pos files without crossdating records #####
#####

## Get all pos files

directory = "data/dendro/coorecorder-measurements"
pos_files = list.files(directory) %>% toupper()

# remove any "t" and ".pos"

meas_cores = str_replace(pos_files,"T.POS","")
meas_cores = str_replace(meas_cores,".POS","")

cores_without_cd = setdiff(meas_cores,cores %>% toupper())
## good, all of our measured cores have been crossdated



####
#### Are there any plot trees without cores? ####
####


# plot trees
plot_data = read.csv("data/plot-and-tree/processed/trees_loc_simple.csv",header=TRUE,stringsAsFactors=FALSE)
plot_trees = plot_data$tree.id %>% toupper()
secondary_plot_trees = plot_data$former.id %>% toupper()
all_plot_trees = union(secondary_plot_trees,plot_trees) %>% unique()


cores_no_suffix = str_replace(meas_cores,"A","")
cores_no_suffix = str_replace(cores_no_suffix,"B","")
cores_no_suffix = str_replace(cores_no_suffix,"Z","")

plot_trees_w_no_cores = setdiff(plot_trees,cores_no_suffix) # 88 field-mesaured trees missing cores
secondary_plot_trees_w_no_cores = setdiff(secondary_plot_trees,cores_no_suffix) # no secondary trees missing cores
# there were 21 secondary trees
# so there were 67 field-measured trees missing cores


### Are there measured cores without plot data?

cores_with_no_plot = setdiff(cores_no_suffix,all_plot_trees)







