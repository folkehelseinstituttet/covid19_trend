functions {
#include gaussian_process.stan
#include convolve.stan
#include observed_in_window.stan
}

data {
  int ut;
  int t;
  int obs;
  vector[obs] prev;
  vector[obs] prev_sd2;
  int prev_stime[obs];
  int prev_etime[obs];
  int pbt;
  vector[pbt] prob_detect_mean;
  vector[pbt] prob_detect_sd;
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
  real<lower = 0> sigma;
  vector<lower = 0, upper = 1>[pbt] prob_detect;
}

transformed parameters {
  vector[t] gp;
  vector[t] infections;
  vector[t] dcases;
  vector[obs] odcases;
  vector[obs] combined_sigma;
  // update gaussian process
  gp = update_gp(PHI, M, L, alpha, rho, eta, 0);
  // relative probability of infection
  infections = inv_logit(inc_zero + gp);
  // calculate detectable cases
  dcases = convolve(infections, prob_detect);
  // calculate observed detectable cases
  odcases = observed_in_window(dcases, prev_stime, prev_etime, ut, obs);
  //combined standard error
  combined_sigma = sqrt(square(sigma) + prev_sd2);
}

model {
  // gaussian process priors
  rho ~ inv_gamma(lengthscale_alpha, lengthscale_beta);
  alpha ~ std_normal() T[0,];
  eta ~ std_normal();

  // prevalence observation model
  for (i in 1:pbt) {
    prob_detect[i] ~ normal(prob_detect_mean[i], prob_detect_sd[i]) T[0, 1];
  }
  sigma ~ normal(0.005, 0.0025) T[0,];
  prev ~ normal(odcases, combined_sigma);
}

generated quantities {
  // vector[t - ut] R;
  // vector[t - 1] r;
  // vector<lower = 0>[t]  cumulative_infections;
  // real est_prev[obs];
  // // cumulative incidence
  // cumulative_infections = cumulative_sum(infections);

  // // sample estimated prevalence
  // est_prev = normal_rng(odcases, combined_sigma);
  // // sample generation time
  // real gtm_sample = normal_rng(gtm[1], gtm[2]);
  // real gtsd_sample = normal_rng(gtsd[1], gtsd[2]);
  // // calculate Rt using infections and generation time
  // R = calculate_Rt(infections, ut, gtm_sample, gtsd_sample, gtmax, 1);
  // // calculate growth
  // r = calculate_growth(infections, 1);
}
