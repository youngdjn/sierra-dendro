#### This script: ####
# - loads compiled tree growth and climate data to explore what is there 
# just starting here to get a handle on the data and its structure . . . 
 

library(tidyverse)
library(ggExtra)
library(broom)
library(lme4)
library(fields)

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
year.ids <- sort(unique(years$year))
n.years <- length(year.ids)
growth_tree_by_year <- matrix(NA, n.trees, n.years)
for (i in 1:n.trees) {
  z <- filter(years, tree.id==tree.ids[i])
  for (j in 1:n.years) {
    zz <- z$growth[which(z$year==year.ids[j])]
    if (!is.null(zz)) growth_tree_by_year[i,j] <- zz else growth_tree_by_year[i,j] <- NA
  }
}

# Visualize whole data matrix
image.plot(growth_tree_by_year)



# Write these data in form readable by dplr as a rwl object: 
dat <- as.data.frame(t(growth_tree_by_year))
colnames(dat) <- tree.ids
rownames(dat) <- year.ids
write.csv(dat, "./working-data/tree_rings_out.csv")



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
summary(lm(AR1~species+cluster, trees)) # slightly stronger autocorrelation in the low-elevation clusters

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

# For now, look only at the common 3 species: 
d_long <- filter(d_long, species %in% c("ABCO", "PIPO", "PSME"))

# How many trees have ring widths available for each year?
table(d_long$year[!is.na(d_long$raw_width)])
#plot(table(d_long$year[!is.na(d_long$raw_width)]), type="l", ylab="Number of trees")
# By cluster: 
plot(table(d_long$year[!is.na(d_long$raw_width) & d_long$cluster.x == "SL"]), type="l", col="orange", ylab="Number of trees")
lines(table(d_long$year[!is.na(d_long$raw_width) & d_long$cluster.x == "NH"]), type="l", col="blue")
lines(table(d_long$year[!is.na(d_long$raw_width) & d_long$cluster.x == "NL"]), type="l", col="cyan")
lines(table(d_long$year[!is.na(d_long$raw_width) & d_long$cluster.x == "SH"]), type="l", col="red")

# By species
plot(table(d_long$year[!is.na(d_long$raw_width) & d_long$species == "PSME"]), type="l", col="darkgreen", ylab="Number of trees")
lines(table(d_long$year[!is.na(d_long$raw_width) & d_long$species == "PIPO"]), type="l", col="orange")
lines(table(d_long$year[!is.na(d_long$raw_width) & d_long$species == "ABCO"]), type="l", col="blue")
#lines(c(1975,1975), c(0, 500))
legend("topleft", c("PSME", "PIPO", "ABCO"), col=c("darkgreen", "orange", "blue"), lwd=rep(3, 3), cex=0.75)


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

# Quick check with lmer 
m1 <- lmer(raw_width ~ ppt.z*tmean.z + ppt.z1*tmean.z1 + cluster + (1|tree.id) , data=d)
summary(m1)

# all main variables with raw ring widths
m_raw <- lmer(raw_width ~ ppt.norm.std + tmean.norm.std  + year.std + ppt.z + ppt.z1 + ppt.z2 + ppt.z3 + tmean.z + tmean.z1 + tmean.z2 + tmean.z3 + ppt.z*tmean.z + (1|tree.id)  + (1+ppt.z|cluster), data=d)
summary(m_raw)
coef(m_raw)$cluster$ppt.z

# all main variables with ring width index
m_rwi <- lm(rwi ~ ppt.norm.std + tmean.norm.std + year.std + ppt.z*cluster + ppt.z1 + ppt.z2 + tmean.z + tmean.z1 , data=d)
summary(m_rwi)

# check influence, if any, of local variation in tree growth rate on climate sensitivity
m_rwi <- lm(rwi ~ ppt.norm.std + tmean.norm.std + ppt.z*tree.growth.z + tmean.z*tree.growth.z + plot.id, data=d)
summary(m_rwi, digits=2)

m_rwi <- lm(rwi ~ ppt.norm.std + tmean.norm.std + ppt.z*tree.growth.z + tmean.z*tree.growth.z + plot.id, data=d, subset = d$ppt.z <= -1) # sensitivity to ppt higher in drought 
summary(m_rwi, digits=2)

m_rwi <- lm(rwi ~ ppt.norm.std + tmean.norm.std + ppt.z*tree.growth.z + tmean.z*tree.growth.z + plot.id, data=d, subset = d$ppt.z >- 1) # ppt sensitivity lower in wet years; tree growth rate more positively associated with ppt sensitivity in wet years



#######################
#### Next steps: 
#### 1) Figure out best (stiff spline?) detrending method 
#### 2) Fit regressions to individual trees to explore how coeffieient variation is partitioned among plots and regions, and how well this variation can be explained by normal climate or other plot-level variables. 
#### 3) Explore nonlinear response shapes for tree growth relation to precipitation
#### 4) Refit individual tree regressions for the subset of dry years only, and the subset of wet/normal years only. Do climatic relationships within trees, and among plots, change a lot in wet vs dry years? 


# Checking nonlinearity of growth response 

# At the cluster level -- this combines responses across trees and sites within a cluster. So nonlinearity may or may not be present within individual trees.
library(gam)
m_gam <- gam(rwi ~ ppt.norm.std + tmean.norm.std + s(ppt.z) + s(tmean.z), data=d[d$cluster=="NH",])
plot(m_gam)
# Across the board there seems a steeper sensitivity to ppt in drought -- in the ppt.z values more than 1 sd below average. In the north high sites ONLY, there's a clear flattening or even reversal of sensitivity at >1 sd above mean ppt. Potentially this is a snowpack effect? 
# Temperature isn't that important, and the splines are all over the place. Hard to interpret. 

# At the plot level within the "NH" cluster -- this combines responses across trees and sites within a cluster. So nonlinearity may or may not be present within individual trees.
