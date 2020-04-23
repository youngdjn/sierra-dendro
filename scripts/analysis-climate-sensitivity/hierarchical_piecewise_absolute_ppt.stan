
data {
  int<lower=0> N;   // number of data items
  int<lower=0> N_groups; // number of groups
  vector[N] x;   // individual-level predictor vector
  vector[N] y;      // outcome vector
  //vector[N_groups] x_group;   // group-level predictor vector 
  int<lower=1,upper=N_groups> group_index[N]; // group index
}
parameters {
  real a[N_groups];	// needs group-specific intercept because ppt not centered at 1 for each group any more 
  real mu_a;           // overall intercept #1
  real b1[N_groups];       // random slopes for groups #1 
  real mu_b1;         // overall slope #1 
  real b2[N_groups];       // random slopes #2 for groups
  real mu_b2;         // overall slope #2
  real cp[N_groups];	// random changepoints for groups
  real mu_cp;           // overall changepoint
  //real b_group; // group-level coefficient
  real<lower=0> sigma_y;  // error scale -- individual level
  real<lower=0> sigma_a;  // error scale -- group level intercept #1
  real<lower=0> sigma_b1;  // error scale -- group level slope #1
  real<lower=0> sigma_b2;  // error scale -- group level slope #2
  real<lower=0> sigma_cp;  // error scale -- group level changepoint
}
transformed parameters { 
  vector[N] conditional_mean; 
  vector[N_groups] cp_abs;
  vector[N] x_std;
  real x_mean;
  real x_sd;
  
  // standardize x
  x_mean = mean(x);
  x_sd = sd(x);
  for (i in 1:N) {
    x_std[i] = (x[i]-x_mean) / x_sd;
  }
    
  // conditional_mean depends on whether x is <> cp
  for (i in 1:N) {
    if (x_std[i] < cp[group_index[i]]) {
      conditional_mean[i] = a[group_index[i]] + b1[group_index[i]] * (x_std[i] - cp[group_index[i]]);
    } else {
      conditional_mean[i] = a[group_index[i]] + b2[group_index[i]] * (x_std[i] - cp[group_index[i]]);
    }
  }
  
  // back-calculate changepoint values for each group based on its mean and sd of ppt 
  for (j in 1:N_groups) {
    cp_abs[j] = cp[j] * x_sd + x_mean;
  }
}
model {
  // group-level priors
  for (n in 1:N_groups) {
  	a[n] ~ normal(mu_a, sigma_a);
    b1[n] ~ normal(mu_b1, sigma_b1); // + b_group_b1 * x_group[n]
    b2[n] ~ normal(mu_b2, sigma_b2); // + b_group_b2 * x_group[n]
    cp[n] ~ normal(mu_cp, sigma_cp); // + b_group_cp * x_group[n]
  }
  // overall priors
  mu_a ~ normal(0, 2); 
  mu_b1 ~ normal(0, 2);
  mu_b2 ~ normal(0, 2);
  mu_cp ~ normal(0, 2);
  //b_group ~ normal(0, 2);
  sigma_y ~ normal(0,2) T[0,]; 
  sigma_a ~ normal(0,2) T[0,]; 
  sigma_b1 ~ normal(0,2) T[0,];
  sigma_b2 ~ normal(0,2) T[0,];
  sigma_cp ~ normal(0,2) T[0,];
  
  // Likelihood
  for (i in 1:N)
    y[i] ~ normal(conditional_mean[i], sigma_y);  // likelihood
}

