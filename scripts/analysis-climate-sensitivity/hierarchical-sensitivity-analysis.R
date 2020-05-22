# This script: 
# 1) combines Derek's processed data into a long-form data frame for analysis;
# 2) filters the data based on species, year, and number of data points per tree; 

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

#### Load and merge data ####
years <- read.csv("./data/compiled-for-analysis/years.csv")
plots <- read.csv("./data/compiled-for-analysis/plots.csv")
trees <- read.csv("./data/compiled-for-analysis/trees.csv")

# Keep only relevant rows of tree data
trees_select <- trees[ , c("tree.id", "species", "plot.id", "x", "y", "elev", "cluster", "voronoi.area", "radius", "radius.external")]

# Create lagged values of rwi
years <- years %>% 
  split(.$tree.id) %>%
  map_dfr(f <- function(x) {return(mutate(x, rwi1=lag(rwi,1), rwi2=lag(rwi,2), rwi3=lag(rwi,3)))})

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
d <- dplyr::select(d, cluster.x, cluster.y, plot.id, tree.id, species, year, starts_with("rwi"), raw_width, ba, bai, voronoi.area, radius.external, ppt.norm, tmean.norm, rad.tot, starts_with("ppt"), starts_with("tmean"))
head(d)

#### NOTE: Get TMAX as well -- Derek will extract

#### NOTE: Get spline values as well 

# Center and scale unstandardized columns 
stdize <- function(x) {return((x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE))}
d <- mutate(d, ppt.std = stdize(ppt), ppt1.std = stdize(ppt1), ppt2.std = stdize(ppt2), ppt3.std = stdize(ppt3), tmean.std=stdize(tmean), tmean1.std=stdize(tmean1), tmean2.std=stdize(tmean2), tmean3.std=stdize(tmean3), ppt.norm.std = stdize(ppt.norm), tmean.norm.std = stdize(tmean.norm), year.std = stdize(year), rad.tot.std = stdize(rad.tot), voronoi.area.std = stdize(voronoi.area))




#### Test run a mixed model on the data in stan 

# SINGLE_LEVEL VERSION

# Set up data for stan
d_test <- dplyr::filter(d, species=="PSME")
d_test <- d_test[!is.na(d_test$rwi) & !is.na(d_test$ppt.z) & !is.na(d_test$ppt.std),]
# optionally subsample for testing speed
d_test <- filter(d_test, cluster.y == "Plumas") # cluster.y == "Yose" |

#stan.data <- list(N=nrow(d_test),  y = d_test$rwi, x = d_test$ppt.z)
# test with standardized absolute precip data
stan.data <- list(N=nrow(d_test),  y = d_test$rwi, x = as.vector(d_test$ppt.std))

# Run model 
model.path <- ("./scripts/analysis-climate-sensitivity/simple_piecewise_linear_regression.stan")
m <- stan(file=model.path, data=stan.data, iter=2000, chains=3, cores = 6, verbose=T)

# Check results
print(m)
check_hmc_diagnostics(m)
stan_dens(m, pars = c("a","b1", "b2","sigma_y", "cp"))
plot(m, pars=c("a","b1", "b2","sigma_y", "cp"))

# Check results for cp one cluster at a time using ppt.std as x 
# Sierra: -0.27
# Yose: 0.22
# Tahoe: -0.42
# Plumas: -0.24
# OK I don't understand this. It seems inconsistent with the plots of the raw data. Maybe we need to look down into the plot level, and that's where most of the variation is. 


# Single level version for all plots, one plot at a time
d_test <- dplyr::filter(d, species=="PSME")
d_test <- d_test[!is.na(d_test$rwi) & !is.na(d_test$ppt.std),]
d_test$plot_index <- match(d_test$plot.id, unique(d_test$plot.id))
plot_info <- summarise(group_by(d_test, plot_index), plot.id = first(plot.id), cluster.y = first(cluster.y), ppt.norm = first(ppt.norm), rad.tot = first(rad.tot))

stan.data <- list(N=nrow(d_test), N_groups = max(d_test$plot_index), y = d_test$rwi, x = as.vector(d_test$ppt.z), group_index = d_test$plot_index)

model.path <- ("./scripts/analysis-climate-sensitivity/simple_piecewise_linear_regression_by_plot.stan")
m <- stan(file=model.path, data=stan.data, iter=1000, chains=3)

#save(m, file = "./working-data/piecewise_plot_no_pool.RData")
#save(m, file = "./working-data/piecewise_plot_no_pool_ppt-z.RData")
#save(m, file = "./working-data/piecewise_plot_no_pool_fixed_a.RData")

# Check results
summary(m, pars="cp")$summary

a <- as.data.frame(m, pars="a")
cp <- as.data.frame(m, pars="cp")
b1 <- as.data.frame(m, pars="b1")
b2 <- as.data.frame(m, pars="b2")
plot_info <- cbind(plot_info, cp = apply(cp, 2, mean), a = apply(a, 2, mean), b1 = apply(b1, 2, mean), b2 = apply(b2, 2, mean))

cp_plot <- ggplot(plot_info, aes(x = ppt.norm, y = cp, color = cluster.y)) + geom_point() + theme_classic()
a_plot <- ggplot(plot_info, aes(x = ppt.norm, y = a, color = cluster.y)) + geom_point() + theme_classic()
b1_plot <- ggplot(plot_info, aes(x = ppt.norm, y = b1, color = cluster.y)) + geom_point() + theme_classic()
b2_plot <- ggplot(plot_info, aes(x = ppt.norm, y = b2, color = cluster.y)) + geom_point() + theme_classic()
grid.arrange(a_plot, cp_plot, b1_plot, b2_plot, nrow=2)


cor(plot_info[,c("a", "b1", "b2", "cp", "ppt.norm", "rad.tot")])

names(a) <- as.character(unique(d_test$plot.id))
mcmc_intervals(a)
names(cp) <- as.character(unique(d_test$plot.id))
mcmc_intervals(cp)

plot(cp~a, plot_info)
summary(lm(b1~scale(ppt.norm), plot_info))
summary(lm(b2~scale(ppt.norm), plot_info))
summary(lm(cp~scale(ppt.norm), plot_info))
summary(lm(a~scale(ppt.norm), plot_info))
# within plumas, changepoint location is negatively correlated with ppt.norm and with rad.tot





# HIERARCHICAL VERSION 

# Set up data for stan
d_test <- dplyr::filter(d, species=="PSME")
d_test <- d_test[!is.na(d_test$rwi) & !is.na(d_test$ppt.z),]

# Subsample data to speed up testing 
#d_test <- filter(d_test, cluster.y == "Sierra" | cluster.y == "Yose")
d_test$plot_index <- match(d_test$plot.id, unique(d_test$plot.id))
d_test$cluster_index <- match(d_test$cluster.y, unique(d_test$cluster.y))
d_plot_test <- summarise(group_by(d_test, plot.id), cluster.y = first(cluster.y), raw_width_mean = mean(raw_width, na.rm=T), ppt.norm = mean(ppt.norm, na.rm=T), ppt.norm.std = mean(ppt.norm.std, na.rm=T), ppt.mean.within = mean(ppt, na.rm=T), ppt.sd.within = sd(ppt, na.rm=T))
d_plot_test
d_cluster_test <- summarise(group_by(d_test, cluster.y), cluster.x = first(cluster.x), raw_width_mean = mean(raw_width, na.rm=T), ppt.norm = mean(ppt.norm, na.rm=T), ppt.norm.std = mean(ppt.norm.std, na.rm=T), ppt_mean = mean(ppt, na.rm=T))
# Add in the average within-cluster sd of ppt 
d_within_test <- summarise(group_by(d_plot_test, cluster.y), ppt.mean.within = mean(ppt.mean.within, na.rm=T), ppt.sd.within = mean(ppt.sd.within, na.rm=T))
# Hack -- reorder the summary to match the cluster index -- not sure how to do this better 
d_within_test <- d_within_test[c(3, 2, 4, 1),]

# Clusters as groups 
#stan.data <- list(N=nrow(d_test), N_groups = max(d_test$cluster_index), y = d_test$rwi, x = as.vector(d_test$ppt.z), group_index = d_test$cluster_index, x_group_mean = d_within_test$ppt.mean.within, x_group_sd = d_within_test$ppt.sd.within)

# Plots as groups 


# Run model 

# using ppt.z
#model.path <- ("./scripts/analysis-climate-sensitivity/hierarchical_piecewise_linear_regression.stan")
#m <- stan(file=model.path, data=stan.data, iter=2000, chains=3)
stan.data <- list(N=nrow(d_test), N_groups = max(d_test$plot_index), y = d_test$rwi, x = as.vector(d_test$ppt.std), group_index = d_test$plot_index, x_group_mean = d_plot_test$ppt.mean.within, x_sd = d_plot_test$ppt.sd.within)

# using ppt
stan.data <- list(N=nrow(d_test), N_groups = max(d_test$plot_index), y = d_test$rwi, x = as.vector(d_test$ppt), group_index = d_test$plot_index) #x_mean = mean(d_test$ppt, na.rm=T), x_sd = sd(d_test$ppt, na.rm=T)

model.path <- ("./scripts/analysis-climate-sensitivity/hierarchical_piecewise_absolute_ppt.stan")
m <- stan(file=model.path, data=stan.data, iter=2000, chains=3)



# Check results
msum <- summary(m)
head(msum$summary, 20)
check_hmc_diagnostics(m)
stan_dens(m, pars = c("mu_a", "mu_b1", "mu_b2", "mu_cp", "sigma_y", "sigma_b1", "sigma_b2", "sigma_cp"))
plot(m, pars=c("b1", "b2", "cp"))


# 2-level model with pptnorm predicting group-level coefficients

# Set up data for stan
d_test <- dplyr::filter(d, species=="PSME")
d_test <- d_test[!is.na(d_test$rwi) & !is.na(d_test$ppt.z),]

# Subsample data to speed up testing 
#d_test <- filter(d_test, cluster.y == "Sierra" | cluster.y == "Yose")
d_test$plot_index <- match(d_test$plot.id, unique(d_test$plot.id))
d_test$cluster_index <- match(d_test$cluster.y, unique(d_test$cluster.y))
d_plot_test <- summarise(group_by(d_test, plot.id), cluster.y = first(cluster.y), raw_width_mean = mean(raw_width, na.rm=T), ppt.norm = mean(ppt.norm, na.rm=T), ppt.norm.std = mean(ppt.norm.std, na.rm=T), ppt.mean.within = mean(ppt, na.rm=T), ppt.sd.within = sd(ppt, na.rm=T))
d_plot_test

# Run model 

# using ppt.z
#model.path <- ("./scripts/analysis-climate-sensitivity/hierarchical_piecewise_linear_regression.stan")
#m <- stan(file=model.path, data=stan.data, iter=2000, chains=3)
stan.data <- list(N=nrow(d_test), N_groups = max(d_test$plot_index), y = d_test$rwi, x = as.vector(d_test$ppt.std), group_index = d_test$plot_index, x_group_mean = d_plot_test$ppt.mean.within, x_group_sd = d_plot_test$ppt.sd.within, z = d_plot_test$ppt.norm.std)

model.path <- ("./scripts/analysis-climate-sensitivity/piecewise_linear_regression_2-level.stan")
m <- stan(file=model.path, data=stan.data, iter=2000, chains=3)



# Check results
msum <- summary(m)
head(msum$summary, 20)
check_hmc_diagnostics(m)
stan_dens(m, pars = c("mu_a", "mu_b1", "mu_b2", "mu_cp", "sigma_y", "sigma_b1", "sigma_b2", "sigma_cp"))
plot(m, pars=c("b1", "b2", "cp", "b_group_cp"))


# Function to plot distribution of parameters while relabeling groups
parplot <- function(model.name, pars=NULL, labels=NULL) { 
  require(bayesplot)
  if (is.null(pars)) df <- as.data.frame(m) else df <- as.data.frame(m, pars=pars)
  if (!is.null(labels)) names(df) = labels
  mcmc_intervals(df)
}

parplot(m, pars="cp_abs", labels = as.character(d_plot_test$plot.id))
parplot(m, pars="cp", labels = as.character(d_plot_test$plot.id))
parplot(m, pars="a", labels = as.character(d_plot_test$plot.id))
parplot(m, pars="b1", labels = as.character(d_plot_test$plot.id))
parplot(m, pars="b2", labels = as.character(d_plot_test$plot.id))

cp_abs <- as.data.frame(m, pars="cp_abs")
names(cp_abs) <- as.character(d_plot_test$plot.id)
mcmc_intervals(cp_abs)

plot(m, pars = "cp_abs") + scale_y_discrete(labels=as.character(d_plot_test$plot.id))
plot(m, pars = "b1")
plot(m, pars = "b2")
stan_trace(m, pars=c("a", "mu_b1", "mu_b2", "mu_cp"))

# Note order of clusters: Tahoe  Sierra Yose   Plumas
# Using ppt.z, tahoe sierra and yose in expected order, but plumas has steepest low-precip sensitivity. Should 'unstandardize" these changepoints for interpretation 


# Next steps: 
# - add lag1 rwi as explanatory variable


          
#### Specify random slope model ####



# Simplest model at tree level has one one weather variable (ppt.z) and its short-term lag (ppt.z1), one autoregressive term (lag(growth, 1)). 
# Simplest model at plot level -- allowing sensitivity of growth to precipitation to vary among plots (random slopes). 

mod.form <- brmsformula(rwi ~ ppt.z + ppt.norm.std + ppt.z1 + lag(rwi, 1) + (1 + ppt.z | plot.id))

brms_data <- dplyr::filter(d, species=="PSME")

prior <- c(set_prior("normal(0, 2)", class = "Intercept", coef = ""), set_prior("normal(0, 2)", class = "b", coef = ""), set_prior("cauchy(0, 2)", class = "sd", coef = ""))

m0 <- brm(mod.form, data=brms.data, prior=prior, family=gaussian(), chains=3, cores=3, iter=1000)

# Try fitting simple changepoint model in brms
# This is based on code by Paul Buerkner at https://discourse.mc-stan.org/t/piecewise-linear-mixed-models-with-a-random-change-point/5306 

bform <- bf(
  rwi ~ b0 + b1 * (ppt.z - alpha) * step(alpha - ppt.z) + 
    b2 * (ppt.z-alpha) * step(ppt.z - alpha) + b3 * ppt.z1, 
  #b0 + b1 + b2 + b3 + alpha ~ 1 + (1|plot.id/tree.id), #(1|cluster.y/plot.id)
  #rwi ~ (b0 + b1 * ppt.z )*step(alpha - ppt.z) + 
  #(b2 +  b3 * ppt.z) * step(ppt.z - alpha),
  b0 + b1 + b2 + b3 + alpha ~ 1 + (1|plot.id),
  nl = TRUE
)

df <- filter(brms_data, cluster.y == "Sierra")

bprior <- prior(normal(0, 3), nlpar = "b0") +
  prior(normal(0, 3), nlpar = "b1") +
  prior(normal(0, 3), nlpar = "b2") +
  prior(normal(0, 3), nlpar = "b3") +
  #prior(normal(0, 3), nlpar = "b4") +
  #prior(normal(0, 3), nlpar = "b5") +
  prior(normal(0, 0.5), nlpar = "alpha")

make_stancode(bform, data = df, prior = bprior)

fit <- brm(bform, data = df, prior = bprior, chains = 3, iter=500)
summary(fit)

# you need the github version of brms for this to run
marginal_effects(fit)

ranef(fit)

## Check single plots for ppt vs ppt.z 
bform0 <- bf(
  rwi ~ b0 + b1 * ppt.std, 
  b0  + b1 ~ (1|tree.id),
  nl = TRUE
)

bform1 <- bf(
  rwi ~ b0 + b1 * (ppt.std - alpha) * step(alpha - ppt.std) + 
    b2 * (ppt.std-alpha) * step(ppt.std - alpha), 
  b0  + b1 + b2 + alpha ~ (1|tree.id),
  nl = TRUE
)

bprior0 <- prior(normal(0, 3), nlpar="b0") +
  prior(normal(0, 3), nlpar = "b1")

bprior1 <- prior(normal(0, 3), nlpar="b0") +
  prior(normal(0, 3), nlpar = "b1") +
  prior(normal(0, 3), nlpar = "b2") +
  prior(normal(0, 0.5), nlpar = "alpha")

#Plumas plots SS81  SS81B SSC1  SSC5B SS82  SSA76 TR15  TULLB TR01 
#brms_data <- dplyr::filter(d, species=="PSME" & plot.id=="TR01")

# Sierra plots RR20  RR33  RR7   RT02C RR71  RT02D RT01  RT02B RR32 
brms_data <- dplyr::filter(d, species=="PSME" & plot.id=="RT01")

#make_stancode(bform1, data = brms_data, prior = bprior)
m0 <- brm(bform0, data=brms_data, prior=bprior0, family=gaussian(), chains=3, cores=3, iter=1000)
m1 <- brm(bform1, data=brms_data, prior=bprior1, family=gaussian(), chains=3, cores=3, iter=1000)

loo(m0, m1, reloo=TRUE)


marginal_effects(m0) 

ggplot(brms_data, aes(x=ppt, y=rwi, group=tree.id)) + geom_point() + legend_none()

# Basically, most Plumas plots don't see the accellerating sensitivity at lower preciptation. They are shallow and linear. Not sure what to make of that! TR01 seems to be an exception with a decline around 1000 mm ppt 

# Idea: rather than focusing on changepoints and slopes, what about focusing on the total amount of lost productivity due to enhanced sensitivity at low precipitation? This could be something like: take shallow slope for normal/high precipitation conditions and extrapolate low-precip growth rate (to the left of the changepoint if any). Then subtract the fitted growth rate according to the steeper "left-hand" slope, and sum. We could calculate total area of that triangle (down to some standard minimum like 400mm). AND we could calculate the actual loss, which combines that area with the observed rainfall distribution. You could even calculate the relative proportion variation in "lost growth" attributable to tree sensitivity, and to rainfall distribution. This would be a tree-specific property but could be summarized for plots and clusters, or correlated with environmental variables. 


##### I think we need to go through and see how many plots have substantially nonlinear response.

# to do this, we can just loop through all the plots and save results, then set up another function to extract the bits we need from the model objects. 
## Check single plots for ppt vs ppt.z 
bform0 <- bf(
  rwi ~ b0 + b1 * ppt.std, 
  b0  + b1 ~ (1|tree.id),
  nl = TRUE
)

bform1 <- bf(
  rwi ~ b0 + b1 * (ppt.std - alpha) * step(alpha - ppt.std) + 
    b2 * (ppt.std-alpha) * step(ppt.std - alpha), 
  b0  + b1 + b2 + alpha ~ (1|tree.id),
  nl = TRUE
)

bprior0 <- prior(normal(0, 3), nlpar="b0") +
  prior(normal(0, 3), nlpar = "b1")

bprior1 <- prior(normal(0, 3), nlpar="b0") +
  prior(normal(0, 3), nlpar = "b1") +
  prior(normal(0, 3), nlpar = "b2") +
  prior(normal(0, 0.5), nlpar = "alpha")


brms_data <- filter(d, species=="PSME" & cluster.y == "Sierra")
brms_data$plot.id <- droplevels(brms_data$plot.id)
plot_m1 <- brms_data %>% 
  split(.$plot.id) %>%
  map(~ brm(bform1, data=., prior=bprior1, family=gaussian(), chains=3, cores=3, iter=1000))

lapply(plot_m1, FUN=f<-function(m){return(summary(plot_m1[[1]])$fixed[4])})

# Compare models for all plots
plot_loo <- list() 
for (i in 1:length(plot_m0)) plot_loo <- c(plot_loo, loo(plot_m0[[i]], plot_m1[[i]]))
plot_loo
marginal_effects(m1_sierra)
# most plots better with linear model 

# All plots together for Sierra
bform0 <- bf(
  rwi ~ b0 + b1 * ppt.std, 
  b0  + b1 ~ (1|plot.id),
  nl = TRUE
)
bform1 <- bf(
  rwi ~ b0 + b1 * (ppt.std - alpha) * step(alpha - ppt.std) + 
    b2 * (ppt.std-alpha) * step(ppt.std - alpha), 
  b0  + b1 + b2 + alpha ~ (1|plot.id),
  nl = TRUE
)
brms_data <- filter(d, species=="PSME")# & cluster.y == "Yose")
#brms_data <- filter(brms_data, tree.id == "1405")
m0_Yose <- brm(bform1, data=brms_data, prior=bprior1, family=gaussian(), chains=3, cores=3, iter=1000)

# Test fitting smoothed response for each cluster
# Q does brms pool across clusters at all withthis syntax?
test1 <- brm(rwi~s(ppt.std, by=cluster.y), data=brms_data, chains=3, cores=3, iter=1000)
marginal_smooths(test1)






marginal_effects(m1_tahoe)
m0_sierra <- brm(bform0, data=brms_data, prior=bprior0, family=gaussian(), chains=3, cores=3, iter=1000)
loo(m0_sierra, m1_sierra)
# For the whole cluster together, nonlinear is better (dloo = -10.7)

# Graphically compare results for 4 clusters

alphas <- data.frame(sierra = as.data.frame(m1_sierra, pars="alpha")$b_alpha_Intercept, yose = as.data.frame(m1_yose, pars="alpha")$b_alpha_Intercept, tahoe = as.data.frame(m1_tahoe, pars="alpha")$b_alpha_Intercept, plumas = as.data.frame(m1_plumas, pars="alpha")$b_alpha_Intercept)

b1s <- data.frame(sierra = as.data.frame(m1_sierra, pars="b1")$b_b1_Intercept, yose = as.data.frame(m1_yose, pars="b1")$b_b1_Intercept, tahoe = as.data.frame(m1_tahoe, pars="b1")$b_b1_Intercept, plumas = as.data.frame(m1_plumas, pars="b1")$b_b1_Intercept)

b2s <- data.frame(sierra = as.data.frame(m1_sierra, pars="b2")$b_b2_Intercept, yose = as.data.frame(m1_yose, pars="b2")$b_b2_Intercept, tahoe = as.data.frame(m1_tahoe, pars="b2")$b_b2_Intercept, plumas = as.data.frame(m1_plumas, pars="b2")$b_b2_Intercept)
mcmc_intervals(alphas)
mcmc_intervals(b1s)
mcmc_intervals(b2s)

loo(m1_sierra, m0_sierra)
loo(m1_yose, m0_yose)
loo(m1_tahoe, m0_tahoe)
loo(m1_plumas, m0_plumas)
# Much better fit with nonlinear model in Yose, Tahoe, and Plumas. Less advantage for Sierra (borderline equallly good to do linear fit). 

# Now let's compare what it looks like if we fit all clusters together
# (could be an issue with shrinkage among clusters?)
brms_data <- filter(d, species=="PSME")
bform0 <- bf(
  rwi ~ b0 + b1 * ppt.std, 
  b0  + b1 ~ (1|cluster.y/plot.id),
  nl = TRUE
)
bform1 <- bf(
  rwi ~ b0 + b1 * (ppt.std - alpha) * step(alpha - ppt.std) + 
    b2 * (ppt.std-alpha) * step(ppt.std - alpha), 
  b0  + b1 + b2 + alpha ~ (1|cluster.y/plot.id),
  nl = TRUE
)
bprior0 <- prior(normal(0, 3), nlpar="b0") +
  prior(normal(0, 3), nlpar = "b1")
bprior1 <- prior(normal(0, 3), nlpar="b0") +
  prior(normal(0, 3), nlpar = "b1") +
  prior(normal(0, 3), nlpar = "b2") +
  prior(normal(0, 0.5), nlpar = "alpha")
m1_all <- brm(bform1, data=brms_data, prior=bprior1, family=gaussian(), chains=3, cores=3, iter=2000)
m0_all <- brm(bform0, data=brms_data, prior=bprior0, family=gaussian(), chains=3, cores=3, iter=2000)
marginal_effects(m1_all)

# Model comparisons: ~ ppt.std + (1|cluster.y/plot.id) - all Pareto k estimates ok
loo(m0_all, m1_all)
#elpd_diff se_diff
#m1_all    0.0       0.0 
#m0_all -170.8      18.4 

# Model comparison with covariates
#brms_data <- filter(d, species=="PSME")
# subset plots randomly from the clusters with more (9 plots per cluster)
plots <- c(as.character(unique(d$plot.id[d$cluster.y=="Sierra"])), as.character(sample(unique(d$plot.id[d$cluster.y=="Yose"]), 9)), as.character(sample(unique(d$plot.id[d$cluster.y=="Tahoe"]), 9)), as.character(unique(d$plot.id[d$cluster.y=="Plumas"])))
brms_data <- filter(d, species=="PSME" & plot.id %in% plots & year > 1976)

bform0 <- bf(
  rwi ~ b0 + b1*ppt.std + b2*ppt1.std + b3*tmean.std + b4 * rad.tot.std + b5*rwi1 + b6*rwi2, 
  b0 + b1 + b2 ~ (1|clister.y),
  b3 + b4 + b5 + b6   ~ (1|cluster.y/plot.id),
  nl = TRUE
)
bform1 <- bf(
  rwi ~ b0 + b1 * (ppt.std - alpha) * step(alpha - ppt.std) + 
    b2 * (ppt.std-alpha) * step(ppt.std - alpha) + b3*ppt1.std + b4*tmean.std + b5*rad.tot.std + b6*rwi1 + b7*rwi2,
   alpha + b1 + b2 ~ (1|cluster.y/plot.id),
   b0 + b3 + b4 + b5 + b6 + b7 ~ (1|cluster.y),
  nl = TRUE
)
bprior0 <- prior(normal(0.5, 0.5), nlpar="b0") +
  prior(normal(0, 1), nlpar = "b1") +   
  prior(normal(0, 1, nlpar = "b2") +
  prior(normal(0, 1), nlpar = "b3") +
  prior(normal(0, 1), nlpar = "b4") +
  prior(normal(0, 1), nlpar = "b5") +
  prior(normal(0, 1), nlpar = "b6")
  
bprior1 <- prior(normal(0.5, 0.5), nlpar="b0") +
  prior(normal(0, 1), nlpar = "b1") +
  prior(normal(0, 1), nlpar = "b2") +
  prior(normal(0, 1), nlpar = "b3") +
  prior(normal(0, 1), nlpar = "b4") +
  prior(normal(0, 1), nlpar = "b5") +
  prior(normal(0, 1), nlpar = "b6") +
  prior(normal(0, 1), nlpar = "b7") +
  prior(normal(-1, 0.5), nlpar = "alpha")

# check autocorrelation of predictors
cor(brms_data[,c("ppt.std", "ppt1.std","tmean.std","rad.tot.std", "rwi1", "rwi2")], use="pairwise.complete")
vif(lm(rwi~ppt.std + I(ppt^2) +  ppt1.std + tmean.std + rad.tot.std + rwi1+rwi2, data=brms_data)) # low vif for linear terms, high if we add ppt^2

make_stancode(bform1, data=brms_data, prior=bprior1, family=gaussian())

m0_all <- brm(bform0, data=brms_data, prior=bprior0, family=gaussian(), chains=3, cores=3, iter=1000)
system.time(m1_all <- brm(bform1, data=brms_data, prior=bprior1, family=gaussian(), chains=3, cores=3, iter=2000))
# 1000 iterationsof 3 parallel chains took ~10 hours
save(m1_all, m0_all, file="./working-data/current_regression_models.RData")

summary(m1_all)
summary(m0_all)
loo(m0_all, m1_all)

marginal_effects(m1_all)
pp_check(m1_all)
pp_check(m0_all)
stanplot(m1_all)
bayes_R2(m1_all)

get_variables(m1_all)

# plot parameters by cluster
m1_all %>% 
  spread_draws(b_alpha_Intercept, r_cluster.y__alpha[cluster.y,]) %>%
  mutate(cluster.y_mean = b_alpha_Intercept + r_cluster.y__alpha) %>%
  ggplot(aes(y = cluster.y, x = cluster.y_mean)) +
  stat_halfeyeh() # Change point Rank: Plumas, Sierra, Tahoe, Yose
m1_all %>%
  spread_draws(b_b1_Intercept, r_cluster.y__b1[cluster.y,]) %>%
  mutate(cluster.y_mean = b_b1_Intercept + r_cluster.y__b1) %>%
  ggplot(aes(y = cluster.y, x = cluster.y_mean)) +
  stat_halfeyeh() # steeper slopes all pretty similar
m1_all %>%
  spread_draws(b_b2_Intercept, r_cluster.y__b2[cluster.y,]) %>%
  mutate(cluster.y_mean = b_b2_Intercept + r_cluster.y__b2) %>%
  ggplot(aes(y = cluster.y, x = cluster.y_mean)) +
  stat_halfeyeh() # To the right of inflection point, slopes shallowest in Plumas, steepest in Yose and Sierra. (what does this mean? note it corresponds to those having the steepest slopes in a linear fit)
m1_all %>%
  spread_draws(b_b0_Intercept, r_cluster.y__b0[cluster.y,]) %>%
  mutate(cluster.y_mean = b_b0_Intercept + r_cluster.y__b0) %>%
  ggplot(aes(y = cluster.y, x = cluster.y_mean)) +
  stat_halfeyeh() 

ranef(m1_all)
    
    
# Results: really no variation at the tree level (sds for random slopes and intercepts very near zero), so removing that level 
# Convergence is poor when both intercept and b1 slope are allowed to vary by plot. So I'll allow it to vary only by cluster. 
# Of the other explanatory variables, only b6 (lag1 rwi) is important, but it also shows not too much variation by cluster or plot. So for now, allow those to vary onlyby cluster. 
 
 

# check individual tree response shapes
ggplot(d[d$cluster.y == "Sierra",], aes(ppt.std, rwi))  + theme_bw() + geom_point() + geom_smooth(method=loess, se=FALSE, fullrange=FALSE) + facet_wrap(~tree.id) + theme(legend.position="none")
ggplot(d[d$cluster.y == "Plumas",], aes(ppt.std, rwi))  + theme_bw() + geom_point() + geom_smooth(method=loess, se=FALSE, fullrange=FALSE) + facet_wrap(~tree.id) + theme(legend.position="none")
ggplot(d[d$cluster.y == "Yose"& d$cluster.x == "SL",], aes(ppt.std, rwi))  + theme_bw() + geom_point() + geom_smooth(method=loess, se=FALSE, fullrange=FALSE) + facet_wrap(~plot.id) + theme(legend.position="none")
ggplot(d[d$cluster.y == "Tahoe"& d$cluster.x == "NL",], aes(ppt.std, rwi))  + theme_bw() + geom_point() + geom_smooth(method=loess, se=FALSE, fullrange=FALSE) + facet_wrap(~plot.id) + theme(legend.position="none")
# cluster-wide response shapes
ggplot(d[d$cluster.y == "Sierra",], aes(ppt.std, rwi))  + theme_bw() + geom_point() + geom_smooth(method=loess, se=FALSE, fullrange=FALSE) + theme(legend.position="none")
ggplot(d[d$cluster.y == "Plumas",], aes(ppt.std, rwi))  + theme_bw() + geom_point() + geom_smooth(method=loess, se=FALSE, fullrange=FALSE) + theme(legend.position="none")
ggplot(d[d$cluster.y == "Yose" & d$cluster.x == "SL",], aes(ppt.std, rwi))  + theme_bw() + geom_point() + geom_smooth(method=loess, se=FALSE, fullrange=FALSE) + theme(legend.position="none")
ggplot(d[d$cluster.y == "Tahoe" & d$cluster.x == "NL",], aes(ppt.std, rwi))  + theme_bw() + geom_point() + geom_smooth(method=loess, se=FALSE, fullrange=FALSE) + theme(legend.position="none")

# Other ideas 

# Is there still autocorrelation in the residuals? 

# Are there species differences in coefficients? Species differences by cluster? 

# Strength of autocorrelation by tree -- any systematic differences across precip gradients? Could look at both raw autocor in raw ring width, and at rwi, and at residuals of model? 

# Some quantification of "additional lost growth" due to nonlinearity of response to low ppt. This is a combination of sensitivity threshold, low-precip slope, and the distribution of observed precip. Could map areas with high sensitivity, and also areas with high realized loss (by time period?)

# Can variation in threshold (alpha) be explained by latitude, or temperature variation at the same level of precip? 
