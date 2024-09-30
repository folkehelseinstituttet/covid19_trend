functions {
#include gaussian_process.stan
#include convolve.stan
#include observed_in_window.stan
}

data {
  int ut;
  int t;
  int n_obs;
  int<lower=0> A[n_obs];
  real lengthscale_alpha; // alpha for gp lengthscale prior
  real lengthscale_beta;  // beta for gp lengthscale prior
  int <lower = 1> M; // approximate gp dimensions
  real L; // approximate gp boundary
  real inc_zero;
}

transformed data {
  // set up approximate gaussian process
  matrix[t, M] PHI = setup_gp(M, L, t);
}

parameters {
  real<lower = 0> rho; // length scale
  real<lower = 0> alpha; // scale
  vector[M] eta; // eta
  
}

transformed parameters {
  vector[t] gp;
  vector[t] infections;
  real<lower=0> I_W[n_obs];
  // update gaussian process
  gp = update_gp(PHI, M, L, alpha, rho, eta, 0);
  // relative probability of infection
  infections =inc_zero*inv_logit(gp);
  // calculate detectable cases
  for (i in 1:n_obs){
    I_W[i] <- sum(infections[(i*7-6):(i*7)]);
  }

}

model {
  // gaussian process priors
  rho ~ inv_gamma(lengthscale_alpha, lengthscale_beta);
  alpha ~ std_normal() T[0,];
  eta ~ std_normal();
  A ~ poisson(I_W)	;
}
