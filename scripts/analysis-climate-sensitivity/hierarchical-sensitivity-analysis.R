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

lin_formula <- as.formula("rwi ~ rwi1 + rwi2 + rwi3 + ppt.std + ppt1.std + tmean.std + tmean1.std + rad.tot.std + voronoi.area.std + (1 + ppt.std|cluster.y/plot.id)")

set_prior("normal(0, 1)", class="b")
set_prior("lkj(2)", class = "cor")
set_prior("cauchy(0,2)", class = "sd")

brm(lin_formula, d_psme)


