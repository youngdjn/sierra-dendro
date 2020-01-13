# This script: 
# 1) combines Derek's processed data into a long-form data frame for analysis;
# 2) filters the data based on species, year, and number of data points per tree; 

library(dplyr)
library(ggplot2)

#### Load and merge data ####
years <- read.csv("./data/compiled-for-analysis/years.csv")
plots <- read.csv("./data/compiled-for-analysis/plots.csv")
trees <- read.csv("./data/compiled-for-analysis/trees.csv")

# Using tree.id, merge tree-level info (species, plot.id, dbh, elev, cluster, voronoi.area) into yearly increment data set
d_long <- left_join(years, trees[,c("tree.id", "species", "plot.id", "x", "y", "elev", "cluster", "voronoi.area", "radius", "radius.external", "AR1")], by="tree.id", )

# Using plot.id, merge plot-level info (ppt.norm, tmean.norm, x.plot, y.plot, rad.tot,rad.03, rad.06, voronoi.area.mean) into yearly increment data set. 
d_long <- merge(d_long, plots[, c("plot.id", "ppt.norm", "tmean.norm", "x.plot", "y.plot", "rad.tot", "rad.03", "rad.06", "voronoi.area.mean", "n.trees")], by="plot.id")
head(d_long)

## Subset data for analysis. 
# Check length of each tree record
tree_record_length_table <- table(d_long$tree.id[!is.na(d_long$raw_width)])
hist(tree_record_length_table, main="Histogram of tree-ring record lengths")
min_record_length <- 30 # drop trees with shorter record than this number of years
# Assemble data
d <- filter(d_long, tree.id %in% names(tree_record_length_table)[tree_record_length_table >= min_record_length])
dim(d)
length(unique(d$tree.id))
d <- filter(d, year > 1961) # keep last 53 years (50 plus 3 lags)
d <- filter(d, species == "PSME") # keep just one species
#d <- filter(d, species %in% c("ABCO", "PIPO", "PSME")) # keep the 3 common species
dim(d)
length(unique(d$tree.id))

# Keep only the necessary columns
d <- dplyr::select(d, cluster, plot.id, tree.id, species, year, rwi, raw_width, voronoi.area, radius.external, ppt.norm, tmean.norm, rad.tot, starts_with("ppt.z"), starts_with("tmean.z"))
head(d)


# Center and scale unstandardized columns 
d <- mutate(d, ppt.norm.std = scale(ppt.norm), tmean.norm.std = scale(tmean.norm), year.std = scale(year), rad.tot.std = scale(rad.tot), voronoi.area.std = scale(voronoi.area))

# Create a new variable that represents "within-plot long-term growth rate z-score": how much faster or slower does each tree grow than the average for its plot? 
d$plot_mean_rw <- with(d, ave(raw_width, plot.id, FUN=function(x) mean(x, na.rm=T)))
d$tree_mean_rw <- with(d, ave(raw_width, tree.id, FUN=function(x) mean(x, na.rm=T)))
d$plot_sd_rw <- with(d, ave(tree_mean_rw, plot.id, FUN=function(x) sd(x, na.rm=T)))
d <- mutate(d, tree.growth.z = (tree_mean_rw - plot_mean_rw)/plot_sd_rw)