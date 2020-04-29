# This script: 
# 1) combines Derek's processed data into a long-form data frame for analysis;
# 2) filters the data based on species, year, and number of data points per tree; 

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
d <- mutate(d, ppt.std = scale(ppt), ppt.norm.std = scale(ppt.norm), tmean.norm.std = scale(tmean.norm), year.std = scale(year), rad.tot.std = scale(rad.tot), voronoi.area.std = scale(voronoi.area))




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
m <- stan(file=model.path, data=stan.data, iter=2000, chains=3)

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
stan_trace(m, pars=c("a", "mu_b1", "mu_b2", "mu_cp")

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

ggplot(brms_data, aes(x=ppt, y=rwi, colour=tree.id)) + geom_point() + legend_none()

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
plot_m0 <- brms_data %>% 
  split(.$plot.id) %>%
  map(~ brm(bform0, data=., prior=bprior0, family=gaussian(), chains=3, cores=3, iter=1000))


plots <- unique(brms_data$plot.id)

fit_compare_piecewise <- function(brms_data, iter=1000, chains=3, cores=3){ 
  require(brms); 
  bform0 <- bf(rwi ~ b0 + b1 * ppt.std, b0  + b1 ~ (1|tree.id), nl = TRUE)
  bform1 <- bf(rwi ~ b0 + b1 * (ppt.std - alpha) * step(alpha - ppt.std) + b2 * (ppt.std-alpha) * step(ppt.std - alpha), b0  + b1 + b2 + alpha ~ (1|tree.id), nl = TRUE)
  bprior0 <- prior(normal(0, 3), nlpar="b0") + prior(normal(0, 3), nlpar = "b1")
  bprior1 <- prior(normal(0, 3), nlpar="b0") + prior(normal(0, 3), nlpar = "b1") + prior(normal(0, 3), nlpar = "b2") + prior(normal(0, 0.5), nlpar = "alpha")
  m0 <- brm(bform0, data=brms_data, prior=bprior0, family=gaussian(), chains=chains, cores=cores, iter=iter)
  m1 <- brm(bform1, data=brms_data, prior=bprior1, family=gaussian(), chains=chains, cores=cores, iter=iter)
  loo_out <- loo(m0, m1, cores=cores)
  return(list(m0, m1, loo_out))
}

m0_model_list <-list()
m1_model_list <- list()
loo_result_list <- list()
for (i in plots[1:2]) {
  brms_sub <- dplyr::filter(brms_data, plot.id == i)
  out <- fit_compare_piecewise(brms_sub)
  m0_model_list <- c(m0_model_list, out[[1]])
  m1_model_list <- c(m1_model_list, out[[2]])
  loo_result_list <- c(loo_result_list, out[[3]])
}

