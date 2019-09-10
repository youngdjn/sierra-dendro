#### This script: ####
# - loads compiled tree growth and climate data to explore what is there 
# just starting here to get a handle on the data and its structure . . . 

library(tidyverse)
library(ggExtra)
library(broom)

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

##### Run time series model for each tree separately ####

# Goal here is to partition variation for each tree into long-term trend, autocorrelation, and sensitivity to precipitation and temperature. 

# QUESTIONS: 

# What subset of data to use -- should we look only at trees with a record going back a certain number of years? How far back should we look? For consistency, should we look no further back than a few decades?  Keep all species? 
# Proposal: for now, keep all trees, do analyses on the past 50 years (1965-2014) but check for conflation of patterns with age of tree and length of record. 

# How many trees have ring widths available for each year?
table(d_long$year[!is.na(d_long$raw_width)])
plot(table(d_long$year[!is.na(d_long$raw_width)]), type="l", ylab="Number of trees")

# What's the distribution of length of records across trees? 
hist(table(d_long$tree.id[!is.na(d_long$raw_width)]))

# Do we know the age of trees? It seems like we need to include age or size in the model? 
# Show raw_width time series
ggplot(data=d_long, aes(y=raw_width, x=year)) + geom_line(color="darkgray", aes(group=tree.id)) + geom_smooth(se=FALSE) + theme_classic()
by_tree <- group_by(d_long, tree.id)
z <- do(by_tree, tidy(lm(raw_width~year, .)))
trend_coefs <- filter(z, term=="year")$estimate
hist(trend_coefs) # leans negative, but not by much
# try looking at resids
tree_augmented <- do(by_tree, augment(lm(raw_width~year, .)))
ggplot(tree_augmented, aes(year, .resid)) + geom_line(aes(group = tree.id), alpha = 1 / 3) +  geom_smooth(se = FALSE)

# Note >100 trees seem to be missing density info (voronoi.area)
sum(is.na(trees$voronoi.area))

## Subset data for analysis. 
d <- filter(d_long, year > 1964) # keep last 50 years
tree_record_length_table <- table(d_long$tree.id[!is.na(d_long$raw_width)])
hist(tree_record_length_table, main="Histogram of tree-ring record lengths")
min_record_length <- 30 # drop trees with shorter record than this number of years
d <- filter(d_long, tree.id %in% names(tree_record_length_table)[tree_record_length_table >= min_record_length])
dim(d)
length(unique(d$tree.id))

# Keep only the necessary columns
d <- select(d, cluster, plot.id, tree.id, species, year, rwi, raw_width, voronoi.area, radius.external, ppt.norm, tmean.norm, rad.tot, starts_with("ppt.z"), starts_with("tmean.z"))
head(d)




