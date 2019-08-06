#### This script: ####
# - loads compiled tree growth and climate data to explore what is there 
# just starting here to get a handle on the data and its structure . . . 

library(tidyverse)
library(ggExtra)

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
tree.ids <- sort(unique(trees$tree.id))
n.years <- length(year.ids)
year.ids <- sort(unique(years$year))
growth_tree_by_year <- matrix(NA, n.trees, n.years)
for (i in 1:n.trees) {
  z <- filter(years, tree.id==tree.ids[i])
  for (j in 1:n.years) {
    zz <- z$growth[which(z$year==year.ids[j])]
    if (!is.null(zz)) growth_tree_by_year[i,j] <- zz else growth_tree_by_year[i,j] <- NA
  }
}
# Turn into data frame with year labels on columns and tree labels on rows
growth_tree_by_year <- as.data.frame(growth_tree_by_year)
colnames(growth_tree_by_year) <- as.character(year.ids)
rownames(growth_tree_by_year) <- tree.ids

# Calculate lag1 autocorrelation in tree growth time series
trees_ts <- apply(growth_tree_by_year, 1, acf, na.action=na.exclude, plot=FALSE)
trees_AR1 <- unlist(lapply(trees_ts, f <- function(x) {return(x$acf[2])}))
hist(trees_AR1)
trees$AR1 <- trees_AR1
# Display histogram of autocorrelation values by species
ggplot(trees) + geom_histogram(aes(AR1)) + facet_wrap(~species) + theme_bw() + removeGrid()
# And by cluster 
ggplot(trees) + geom_histogram(aes(AR1)) + facet_wrap(~cluster) + theme_bw() + removeGrid()
# Scatterplot vs elevation
plot(AR1~elev, trees) # almost all the variation is within site
summary(lm(AR1~species+cluster, trees)) # stronger autocorrelation in the low-elevation clusters

# Using tree.id, merge tree-level info (species, plot.id, dbh, elev, cluster, voronoi.area) into yearly increment data set
d_long <- merge(years, trees[,c("tree.id", "species", "plot.id", "x", "y", "elev", "cluster", "voronoi.area", "radius", "radius.external", "AR1")], by="tree.id", )

# Using plot.id, merge plot-level info (ppt.norm, tmean.norm, x.plot, y.plot, rad.tot,rad.03, rad.06, voronoi.area.mean) into yearly increment data set. 
d_long <- merge(d_long, plots[, c("plot.id", "ppt.norm", "tmean.norm", "x.plot", "y.plot", "rad.tot", "rad.03", "rad.06", "voronoi.area.mean", "n.trees")], by="plot.id")
head(d_long)


