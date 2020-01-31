
data {
  int<lower=0> N;   // number of data items
  int<lower=0> N_groups; // number of groups
  vector[N] x;   // predictor vector
  vector[N] y;      // outcome vector
  int<lower=1,upper=N_groups> group_index[N]; // group index
}
parameters {
  real a[N_groups];	// random intercepts
  real mu_a;           // overall intercept
  real b[N_groups];       // coefficients for predictors
  real mu_b;
  real<lower=0> sigma_y;  // error scale -- individual level
  real<lower=0> sigma_a;  // error scale -- group level
  real<lower=0> sigma_b;  // error scale -- group level intercepts
}
model {
  mu_a ~ normal(0, 2); // priors
  mu_b ~ normal(0, 2);
  sigma_y ~ normal(0,2) T[0,]; 
  sigma_a ~ normal(0,2) T[0,]; 
  sigma_b ~ normal(0,2) T[0,]; // gamma(2, 0.5); 
  for (n in 1:N_groups)
	a[n] ~ normal(mu_a, sigma_a);
  for (n in 1:N_groups)
	b[n] ~ normal(mu_b, sigma_b);
  for (n in 1:N)
  	y[n] ~ normal(x[n] * b[group_index[n]] + a[group_index[n]], sigma_y);  // likelihood
}

