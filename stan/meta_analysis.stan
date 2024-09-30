functions {
#include gaussian_process.stan
#include convolve.stan
#include observed_in_window.stan
}

data {
  int ut;
  int t;
  int n_indicators;
  matrix[t, n_indicators] mean_r;
  matrix[t, n_indicators] sd_r;
  real lengthscale_alpha; // alpha for gp lengthscale prior
  real lengthscale_beta;  // beta for gp lengthscale prior
  int <lower = 1> M; // approximate gp dimensions
  real L; // approximate gp boundary
  real inc_zero;
  real r_prior_sd_mean;
  real r_prior_sd_sd;
 
}

transformed data {
  // set up approximate gaussian process
  matrix[t, M] PHI = setup_gp(M, L, t);
}

parameters {
  real<lower = 0> rho; // length scale
  real<lower = 0> alpha; // scale
  vector[M] eta; // eta
  vector<lower =0>[t] sigma_r;
  matrix[t, n_indicators] ri_z;
   vector[t] rz;
  
   real<lower=0> r_prior_sd;
  
}

transformed parameters {
  vector[t] gp;
  vector[t] r;
  matrix[t, n_indicators] ri;
  // update gaussian process
 
  gp = update_gp(PHI, M, L, alpha, rho, eta, 0);
  r = gp + r_prior_sd * rz;
  
  // relative probability of infection
 // r=inc_zero*gp;
  for(i in 1:t){
    for(n in 1:n_indicators){
        ri[i,n] = r[i] + sigma_r[i] * ri_z[i,n];
    }
  }
  // calculate detectable cases
 
}

model {
  // gaussian process priors
  rho ~ inv_gamma(lengthscale_alpha, lengthscale_beta);
  alpha ~ std_normal() T[0,];
  eta ~ std_normal();
  sigma_r ~ normal(0, 0.5);
  //r ~ normal(0,10);// r_prior_sd);
  r_prior_sd ~ normal(r_prior_sd_mean, r_prior_sd_sd);
  rz ~std_normal();
  for(i in 1:t){
    for(n in 1:n_indicators){
     
        mean_r[i, n] ~ normal(ri[i,n], sd_r[i,n]);
        ri_z[i, n] ~normal(0,1);

    }
  }
 
}
