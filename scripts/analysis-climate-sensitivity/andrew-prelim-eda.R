#### This script: ####
# - loads compiled tree growth and climate data to explore what is there 
# just starting here to get a handle on the data and its structure . . . 

library(tidyverse)

#### Load and merge data ####
years <- read.csv("./data/compiled-for-analysis/years.csv")
plots <- read.csv("./data/compiled-for-analysis/plots.csv")
trees <- read.csv("./data/compiled-for-analysis/trees.csv")

# check structure
dim(years); head(years)
dim(plots); head(plots)
dim(trees); head(trees)

# Merge summary stats for individual tree size and growth into tree-level data for overall descriptive plots of data
years$growth <- years$raw_width # Choose growth metric to look at

# Make a data matrix containing growth for trees on rows and years on columns
n.trees <- nrow(trees)
year.ids <- sort(unique(years$year))
n.years <- length(year.ids)
growth_tree_by_year <- matrix(NA, n.trees, n.years)
for (i in 1:n.trees) {
  z <- filter(years, tree.id==trees$tree.id[i])
  for (j in 1:n.years) {
    zz <- z$growth[which(z$year==year.ids[j])]
    if (!is.null(zz)) growth_tree_by_year[i,j] <- zz else growth_tree_by_year[i,j] <- NA
  }
}

# check autocorrelation in tree growth time series
trees_ts <- apply(growth_tree_by_year, 1, acf, na.action=na.exclude)
trees_AR1 <- unlist(lapply(trees_ts, f <- function(x) {return(x$acf[2])}))
hist(trees_AR1)
trees$AR1 <- trees_AR1
ggplot(trees) + geom_histogram(aes(AR1)) + facet_wrap(~species)



# Using tree.id, merge tree-level info (species, plot.id, dbh, elev, cluster, voronoi.area) into yearly increment data set

# Using plot.id, merge plot-level info (ppt.norm, tmean.norm, x.plot, y.plot, rad.tot,rad.03, rad.06, voronoi.area.mean) into yearly increment data set. 
