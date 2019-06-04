#### This script: ####
# - loads compiled tree growth and climate data to explore what is there 
# just starting here to get a handle on the data and its structure . . . 

#### Load and merge data ####
years <- read.csv("./data/compiled-for-analysis/years.csv")
plots <- read.csv("./data/compiled-for-analysis/plots.csv")
trees <- read.csv("./data/compiled-for-analysis/trees.csv")

# check structure
dim(years); head(years)
dim(plots); head(plots)
dim(trees); head(trees)

# Merge summary stats for individual tree size and growth into tree-level data for overall descriptive plots of data

# Using tree.id, merge tree-level info (species, plot.id, dbh, elev, cluster, voronoi.area) into yearly increment data set

# Using plot.id, merge plot-level info (ppt.norm, tmean.norm, x.plot, y.plot, rad.tot,rad.03, rad.06, voronoi.area.mean) into yearly increment data set. 


