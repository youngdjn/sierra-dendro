
data {
  int<lower=0> N;   // number of data items
  int<lower=0> N_groups; // number of groups
  vector[N] x;   // individual-level predictor vector
  vector[N] y;      // outcome vector
  vector[N_groups] x_group;   // group-level predictor vector 
  int<lower=1,upper=N_groups> group_index[N]; // group index
}
parameters {
  real a1[N_groups];	// random intercepts #1 for groups
  real mu_a1;           // overall intercept #1
  real b1[N_groups];       // random slopes for groups #1 
  real mu_b1;         // overall slope #1 
  real a2[N_groups];	// random intercepts #2 for groups
  real mu_a2;           // overall intercept #2
  real b2[N_groups];       // random slopes #2 for groups
  real mu_b2;         // overall slope #2
  real cp[N_groups];	// random changepoints for groups
  real mu_cp;           // overall changepoint
  //real b_group; // group-level coefficient
  real<lower=0> sigma_y;  // error scale -- individual level
  real<lower=0> sigma_a1;  // error scale -- group level intercept #1
  real<lower=0> sigma_b1;  // error scale -- group level slope #1
  real<lower=0> sigma_a2;  // error scale -- group level intercept #2
  real<lower=0> sigma_b2;  // error scale -- group level slope #2
  real<lower=0> sigma_cp;  // error scale -- group level changepoint
}
transformed parameters { 
  vector[N] conditional_mean; 
  real slope_difference;      // the difference between slope_after and slope_before

  slope_difference = mu_b2 - mu_b1;  

  // conditional_mean depends on whether x is <> cp
  for (i in 1:N) {
    if (x[i] < cp[group_index[i]]) {
      conditional_mean[i] = a1[group_index[i]] + b1[group_index[i]] * (x[i] - cp[group_index[i]]);
    } else {
      conditional_mean[i] = a2[group_index[i]] + b2[group_index[i]] * (x[i] - cp[group_index[i]]);
    }
  }
}
model {
  mu_a1 ~ normal(0, 2); // priors
  mu_b1 ~ normal(0, 2);
  mu_a2 ~ normal(0, 2); // priors
  mu_b2 ~ normal(0, 2);
  mu_cp ~ normal(0, 0.5);
  //b_group ~ normal(0, 2);
  sigma_y ~ normal(0,2) T[0,]; 
  sigma_a1 ~ normal(0,2) T[0,]; 
  sigma_b1 ~ normal(0,2) T[0,];
  sigma_a2 ~ normal(0,2) T[0,]; 
  sigma_b2 ~ normal(0,2) T[0,];
  sigma_cp ~ normal(0,2) T[0,];
  for (n in 1:N_groups) 
  	a1[n] ~ normal(mu_a1, sigma_a1);
  for (n in 1:N_groups)
    b2[n] ~ normal(mu_b1, sigma_b1); // + b_group * x_group[n]
  for (n in 1:N_groups) 
    a2[n] ~ normal(mu_a2, sigma_a2);
  for (n in 1:N_groups) 
    b2[n] ~ normal(mu_b2, sigma_b2); // + b_group * x_group[n]
  for (n in 1:N_groups) 
    cp[n] ~ normal(mu_cp, sigma_cp); // + b_group * x_group[n]
  }
  for (i in 1:n)
    y[i] ~ normal(conditional_mean[i], sigma_y);  // likelihood
}

