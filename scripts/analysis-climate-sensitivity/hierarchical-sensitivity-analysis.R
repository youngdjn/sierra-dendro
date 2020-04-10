# This script: 
# 1) combines Derek's processed data into a long-form data frame for analysis;
# 2) filters the data based on species, year, and number of data points per tree; 

library(dplyr)
library(readr)
library(ggplot2)
library(brms)
library(rstan)

# Set stan options 
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

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

# Keep only the columns needed for analysis
d <- dplyr::select(d, cluster.x, cluster.y, plot.id, tree.id, species, year, rwi, raw_width, ba, bai, voronoi.area, radius.external, ppt.norm, tmean.norm, rad.tot, ppt, tmean, starts_with("ppt.z"), starts_with("tmean.z"))
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




#### Test run a mixed model on the data in stan 

model.path <- ("./scripts/analysis-climate-sensitivity/lmm_random_slopes.stan")
d_psme <- dplyr::filter(d, species=="PSME")
d_psme <- d_psme[!is.na(d_psme$rwi) & !is.na(d_psme$ppt.z),]
d_psme$plot_index <- match(d_psme$plot.id, unique(d_psme$plot.id))
d_psme_plot <- summarise(group_by(d_psme, plot.id), raw_width_mean = mean(raw_width, na.rm=T), ppt.norm = mean(ppt.norm, na.rm=T))
  
d_psme$plot_index <- match(d_psme$plot.id, unique(d_psme$plot.id))
stan.data <- list(N=nrow(d_psme), N_groups = max(d_psme$plot_index), y = d_psme$rwi, x = d_psme$ppt.z, group_index = d_psme$plot_index)

m <- stan(file=model.path, data=stan.data)
print(m)
check_hmc_diagnostics(m)
stan_dens(m, pars = c("sigma_y", "sigma_a", "sigma_b")

#### Specify random slope model in model ####

# Simplest model at tree level has one one weather variable (ppt.z) and its short-term lag (ppt.z1), one autoregressive term (lag(growth, 1)). 
# Simplest model at plot level -- allowing sensitivity of growth to precipitation to vary among plots (random slopes). 

mod.form <- brmsformula(rwi ~ ppt.z + ppt.norm.std + ppt.z1 + lag(rwi, 1) + (1 + ppt.z | plot.id))

brms.data <- dplyr::filter(d, species=="PSME")

prior <- c(set_prior("normal(0, 2)", class = "Intercept", coef = ""), set_prior("normal(0, 2)", class = "b", coef = ""), set_prior("cauchy(0, 2)", class = "sd", coef = ""))

m0 <- brm(mod.form, data=brms.data, prior=prior, family=gaussian(), chains=3, cores=3, iter=1000)

# Try fitting simple changepoint model in brms
# This is based on code by Paul Buerkner at https://discourse.mc-stan.org/t/piecewise-linear-mixed-models-with-a-random-change-point/5306 

bform <- bf(
  rwi ~ b0 + b1 * (ppt.z - alpha) * step(alpha - ppt.z) + 
    b2 * (ppt.z - alpha) * step(ppt.z - alpha) + b3 * ppt.z1, 
  b0 + b1 + b2 + b3 + alpha ~ 1 + (1|plot.id/tree.id),
  # to keep omega within the age range of 600 and 1400
  #nlf(omega ~ inv_logit(alpha) * 2 - 1),
  nl = TRUE
)

df <- filter(d_psme, cluster.y == "Plumas")

bprior <- prior(normal(0, 3), nlpar = "b0") +
  prior(normal(0, 3), nlpar = "b1") +
  prior(normal(0, 3), nlpar = "b2") +
  prior(normal(0, 3), nlpar = "b3") +
  #prior(normal(0, 3), nlpar = "b4") +
  #prior(normal(0, 3), nlpar = "b5") +
  prior(normal(0, 0.2), nlpar = "alpha")

make_stancode(bform, data = df, prior = bprior)

fit <- brm(bform, data = df, prior = bprior, chains = 3, iter=1000)
summary(fit)

# you need the github version of brms for this to run
marginal_effects(fit)

ranef(fit)

