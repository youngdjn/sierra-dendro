# This script: 
# 1) combines Derek's processed data into a long-form data frame for analysis;
# 2) filters the data based on species, year, and number of data points per tree; 

library(dplyr)
library(ggplot2)
library(brms)
library(rstan)

#### Load and merge data ####
years <- read.csv("./data/compiled-for-analysis/years.csv")
plots <- read.csv("./data/compiled-for-analysis/plots.csv")
trees <- read.csv("./data/compiled-for-analysis/trees.csv")

# Keep only relevant rows of tree data
trees_select <- trees[ , c("tree.id", "species", "plot.id", "x", "y", "elev", "cluster", "voronoi.area", "radius", "radius.external")]

# Merge
d_long <- left_join(trees_select, years, by="tree.id") %>% 
    left_join(plots, by="plot.id")

#### Subset data for analysis ####

# Filter by number of data points per tree 
tree_record_length_table <- table(d_long$tree.id[!is.na(d_long$raw_width)]) # Check length of each tree record
hist(tree_record_length_table, main="Histogram of tree-ring record lengths")
min_record_length <- 30 # drop trees with shorter record than this number of years
d <- filter(d_long, tree.id %in% names(tree_record_length_table)[tree_record_length_table >= min_record_length]) # Drop trees with less than minimum length

# Filter by year -- drop all data before a given date 
# For now, keep all years back to 1956
#length(unique(d$tree.id))
#d <- filter(d, year > 1955) # keep last 53 years (50 plus 3 lags)

# Filter by species -- keep just the big 3
d <- filter(d, species %in% c("ABCO", "PIPO", "PSME")) # keep the 3 common species
dim(d)
length(unique(d$tree.id))

# Save file 
#write.csv(d, "./working-data/tree-ring-data-long.csv")



#### Set up data for analysis #### 

d <- read.csv("./working-data/tree-ring-data-long.csv")

# Keep only the necessary columns
d <- dplyr::select(d, cluster.x, cluster.y, plot.id, tree.id, species, year, rwi, raw_width, ba, bai, voronoi.area, radius.external, ppt.norm, tmean.norm, rad.tot, starts_with("ppt.z"), starts_with("tmean.z"))
head(d)

#### NOTE: Get TMAX as well -- Derek will extract

#### NOTE: Get spline values as well 

#### NOTE: ALSO get lag1 rwi values

# Create lag1 variable for RWI
#rwi_wide <- select(d, c(tree.id, rwi, year)) %>%
#            reshape(idvar = "tree.id", v.names="rwi", timevar="year", direction="wide")
#rwi_mat <- as.matrix(select(rwi_wide, -tree.id))
#rwi_wide_lag1 <- as.data.frame(t(apply(rwi_mat, 1, lag)))
#rwi_wide_lag1$tree.id <- rwi_wide$tree.id
#rwi_long_lag1 <- reshape(rwi_wide_lag1, direction="long", idvar="tree.id", timevar="year", varying=list(1:ncol(rwi_mat)), sep=".") # Broken -- doesn't preserve year

# Center and scale unstandardized columns 
d <- mutate(d, ppt.norm.std = scale(ppt.norm), tmean.norm.std = scale(tmean.norm), year.std = scale(year), rad.tot.std = scale(rad.tot), voronoi.area.std = scale(voronoi.area))


#### Specify model ####

# Simplest model at tree level has one one weather variable (ppt.z) and its short-term lag (ppt.z1), one autoregressive term (lag(growth, 1)). 
# Simplest model at plot level -- explaining variation in mean growth, and in sensitivity of growth to precipitation, among plots, using long-term normal precipitation. 

mod.form <- brmsformula(rwi ~ ppt.z*ppt.norm.std + ppt.z1 + lag(rwi, 1) + (1 + lag(rwi, 1)| plot.id))

brms.data <- dplyr::filter(d, species=="PSME")

prior <- c(set_prior("normal(0, 2)", class = "Intercept", coef = ""), set_prior("normal(0, 2)", class = "b", coef = ""), set_prior("cauchy(0, 10)", class = "sd", coef = ""))

m0 <- brm(mod.form, data=brms.data, prior=prior, family=gaussian(), chains=3, cores=3, iter=100)

summary(m0)
