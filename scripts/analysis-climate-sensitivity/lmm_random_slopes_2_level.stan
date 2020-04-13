
data {
  int<lower=0> N;   // number of data items
  int<lower=0> N_groups; // number of groups
  vector[N] x;   // individual-level predictor vector
  vector[N] y;      // outcome vector
  vector[N_groups] x_group;   // group-level predictor vector 
  int<lower=1,upper=N_groups> group_index[N]; // group index
}
parameters {
  real a[N_groups];	// random intercepts
  real mu_a;           // overall intercept
  real b[N_groups];       // coefficients for predictors
  real mu_b;
  real b_group; // group-level coefficient
  real<lower=0> sigma_y;  // error scale -- individual level
  real<lower=0> sigma_a;  // error scale -- group level
  real<lower=0> sigma_b;  // error scale -- group level slopes
}
model {
  mu_a ~ normal(0, 2); // priors
  mu_b ~ normal(0, 2);
  b_group ~ normal(0, 2);
  sigma_y ~ normal(0,2) T[0,]; 
  sigma_a ~ normal(0,2) T[0,]; 
  sigma_b ~ normal(0,2) T[0,]; // gamma(2, 0.5); 
  for (n in 1:N_groups)
	  a[n] ~ normal(mu_a, sigma_a);
  for (n in 1:N_groups)
	  b[n] ~ normal(mu_b + b_group * x_group[n], sigma_b);
  for (n in 1:N)
  	y[n] ~ normal(x[n] * b[group_index[n]] + a[group_index[n]], sigma_y);  // likelihood
}

