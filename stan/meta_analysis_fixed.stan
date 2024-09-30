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

  vector[t] r;
  
}

model {
  // gaussian process priors


  r ~ normal(0, r_prior_sd);
  for(i in 1:t){
    for(n in 1:n_indicators){
        mean_r[i, n] ~ normal(r[i], sd_r[i,n]);
    }
  }
  

  
}
