
data {
  int<lower=0> N;   // number of data items
  int<lower=0> N_groups; // number of groups
  vector[N] x;   // individual-level predictor vector
  vector[N] y;      // outcome vector
  int<lower=1,upper=N_groups> group_index[N]; // group index
}

parameters {
  real a[N_groups];	// needs group-specific intercept because ppt not centered at 1 for each group any more 
  real b1[N_groups];       // slopes #1 for groups
  real b2[N_groups];       // slopes #2 for groups
  real cp[N_groups];	// changepoints for groups
  real<lower=0> sigma_y;  // error scale -- individual level
}

transformed parameters { 
  vector[N] conditional_mean; 

  // conditional_mean depends on whether x is <> cp
  for (i in 1:N) {
    if (x[i] < cp[group_index[i]]) {
      conditional_mean[i] = a[group_index[i]] + b1[group_index[i]] * (x[i] - cp[group_index[i]]);
    } else {
      conditional_mean[i] = a[group_index[i]] + b2[group_index[i]] * (x[i] - cp[group_index[i]]);
    }
  }
}

model {
  // group-level priors -- here no pooling 
  for (n in 1:N_groups) {
  	a[n] ~ normal(0, 2);
    b1[n] ~ normal(0, 2); 
    b2[n] ~ normal(0, 2);
    cp[n] ~ normal(0, 2);
  }
  sigma_y ~ normal(0,2) T[0,]; 
  
  // Likelihood
  for (i in 1:N)
    y[i] ~ normal(conditional_mean[i], sigma_y);  // likelihood
}

