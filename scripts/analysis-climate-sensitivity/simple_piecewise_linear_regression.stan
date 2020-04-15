
data {
  int<lower=0> N;   // number of data items
  vector[N] x;   // individual-level predictor vector
  vector[N] y;      // outcome vector
}
parameters {
  real a1;	// random intercepts #1 for groups
  real b1;       // random slopes for groups #1 
  real a2;	// random intercepts #2 for groups
  real b2;       // random slopes #2 for groups
  real cp;	// random changepoints for groups
  //real b_group; // group-level coefficient
  real<lower=0> sigma_y;  // error scale -- individual level
}
transformed parameters { 
  vector[N] conditional_mean; 
  real slope_difference;      // the difference between slope_after and slope_before

  slope_difference = b2 - b1;  

  // conditional_mean depends on whether x is <> cp
  for (i in 1:N) {
    if (x[i] < cp) {
      conditional_mean[i] = a1 + b1 * (x[i] - cp);
    } else {
      conditional_mean[i] = a2 + b2 * (x[i] - cp);
    }
  }
}
model {
  // priors
  sigma_y ~ normal(0,2) T[0,]; 
  a1 ~ normal(0, 2);
  a2 ~ normal(0, 2);  
  b1 ~ normal(0, 2);  
  b2 ~ normal(0, 2);  
  cp ~ normal(0, 0.5);
  // likelihood
  for (i in 1:N)
    y[i] ~ normal(conditional_mean[i], sigma_y);  // likelihood
}

