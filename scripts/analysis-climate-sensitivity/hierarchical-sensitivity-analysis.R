# This script: 

# 1) Based on results from hierarchical-sensitivity-analysis-exploration.R, sets up hierarchical analysis of tree-ring growth sensitivity, using brms 

#### NOTES: 
#### - Get TMAX as well -- Derek will extract
#### - Get spline values as well to use as covariate?  

library(car)
library(dplyr)
library(reshape2)
library(readr)
library(ggplot2)
library(gridExtra)
library(brms)
library(rstan)
library(bayesplot)
library(arm)
library(purrr)
library(tidybayes)
library(brmstools)


# Set stan options 
options(mc.cores = parallel::detectCores() - 2)
rstan_options(auto_write = TRUE)

#### Set up data for analysis #### 

d <- read.csv("./working-data/tree-ring-data-long.csv")

# Keep only the columns needed for analysis
d <- dplyr::select(d, cluster.x, cluster.y, plot.id, tree.id, species, year, starts_with("rwi"), raw_width, ba, bai, voronoi.area, radius.external, ppt.norm, tmean.norm, rad.tot, starts_with("ppt"), starts_with("tmean"))
head(d)

# Center and scale unstandardized columns 
# Using a custom function because using scale() sometimes produces bad side-effects with stan. 
stdize <- function(x) {return((x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE))}
d <- mutate(d, ppt.std = stdize(ppt), ppt1.std = stdize(ppt1), ppt2.std = stdize(ppt2), ppt3.std = stdize(ppt3), tmean.std=stdize(tmean), tmean1.std=stdize(tmean1), tmean2.std=stdize(tmean2), tmean3.std=stdize(tmean3), ppt.norm.std = stdize(ppt.norm), tmean.norm.std = stdize(tmean.norm), year.std = stdize(year), rad.tot.std = stdize(rad.tot), voronoi.area.std = stdize(voronoi.area))

# Separate data into species 
d_psme <- dplyr::filter(d, species=="PSME")
d_pipo <- dplyr::filter(d, species=="PIPO")
d_abco <- dplyr::filter(d, species=="ABCO")


# Run linear regression model for precipitation sensitivity only, one species at a time

brms_data <- d_psme

# Optionally subset the data for testing
sub <- 0.3
brms_data <- brms_data[sample(1:round(nrow(brms_data)*sub)),]

brms_formula <- as.formula("raw_width ~ rwi1 + rwi2 + ppt.std * ppt.norm.std + ppt1.std + tmean.std + rad.tot.std + voronoi.area.std + voronoi.area.std:ppt.std + (1 + rwi1 + rwi2 + ppt.std + ppt.std:voronoi.area.std+ ppt1.std + tmean.std|cluster.y/plot.id)")

set_prior("normal(0, 1)", class="b")
#set_prior("lkj(2)", class = "cor")
#set_prior("cauchy(0,2)", class = "sd")

m1_psme <- brm(brms_formula, data=brms_data, iter=2000)
forest(m1_psme, grouping = "cluster.y", par="ppt.std:voronoi.area.std")

summary(m1_psme)
conditional_effects(m1_psme)
loo(m1_psme) # with ppt and lagged rwi as random slopes, looic= -1966.4

# rad.tot.std:ppt.std has little effect and little variation among clusters
#

get_variables(m1_psme)


m1_psme %>% 
  spread_draws(r_cluster.y[,"ppt.std"]) %>%
  mutate(cluster.y_mean = b_ppt.std + r_cluster.y[,"ppt.std"]) %>%
  ggplot(aes(y = cluster.y, x = cluster.y_mean)) + 
  stat_halfeyeh() # Change point Rank: Plumas, Sierra, Tahoe, Yose


# Try spline fit 
brms_formula <- as.formula("rwi ~ rwi1 + rwi2 + s(ppt.std) + ppt.norm.std + ppt1.std + tmean.std + rad.tot.std + voronoi.area.std + (1|cluster.y/plot.id)")
set_prior("normal(0, 1)", class="b")
m2_psme <- brm(brms_formula, data=brms_data, iter=2000)
conditional_effects(m2_psme)  
loo(m1_psme, m2_psme) # nonhierarchical spline fit is quite nonlinear, but nonetheless, the resulting model minus random slopes is worse.

# Next steps: 
# Test for among-cluster and among-plot variation in other variables
# Test for interactions between ppt.std and various environmental variables (rad.tot, voronoi, ppt.norm, tmean.norm, tmean, rwi1)

# To compare species changes in sensitivity, maybe plot sensitivity in each cluster for each species? Or predict sensitivity under wet vs dry conditions for each species, and compare? 

# Try predicting how trees from Plumas would respond to the Sierra environment, and vice versa. Check how different the prediction is using nonlinear vs linear fits. 


