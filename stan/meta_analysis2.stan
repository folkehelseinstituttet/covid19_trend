functions {
#include gaussian_process.stan
#include convolve.stan
#include observed_in_window.stan
}

data {
  
  int t;
  int n_indicators;
  matrix[t, n_indicators] mean_r;
  matrix[t, n_indicators] sd_r;
 real r_prior_sd;
 
}



parameters {
  vector<lower =0>[t] sigma_r;
  vector[t] r;
  matrix[t, n_indicators] ri_z;
 // real <lower=0> r_prior_sd;
  
}

transformed parameters {
    matrix[t, n_indicators] ri;
  //for(i in 1:t){
  //  for(n in 1:n_indicators){
  //      ri[i,n] = r[i] + sigma_r[i] * ri_z[i,n];
  //  }
  //}
  ri = rep_matrix(r, n_indicators) + rep_matrix(sigma_r, n_indicators) .* ri_z;
  // calculate detectable cases
 
}

model {
  // gaussian process priors
 // r_prior_sd ~ normal(0, 1);
  sigma_r ~ normal(0, 0.5);
  r ~ normal(0, r_prior_sd);
 // for(i in 1:t){
    for(n in 1:n_indicators){
        mean_r[, n] ~ normal(ri[,n], sd_r[,n]);
        ri_z[, n] ~normal(0,1);

    }
  //}
  

  
}
