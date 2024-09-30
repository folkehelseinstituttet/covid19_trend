
data {
  
  int t;
  int n_indicators;
  matrix[t, n_indicators] mean_r;
  matrix[t, n_indicators] sd_r;
  real r_prior_sd;
 
}



parameters {
   vector[t] r;
  vector<lower =0>[t] sigma_r;
 // vector[t] rz;
  matrix[t, n_indicators] ri_z;
  real alpha_rw;
  real beta_rw;
  real<lower=0> sigma_rw;
  
}

transformed parameters {
   
  //  r[1] = rz[1]*r_prior_sd;
   // for(i in 2:t){
   //     r[i] = alpha_rw + beta_rw * r[i-1] + sigma_rw * rz[i];
   // }
    
   
    matrix[t, n_indicators] ri;
  //for(i in 1:t){
   // for(n in 1:n_indicators){
   //     ri[i,n] = r[i] + sigma_r[i] * ri_z[i,n];
    //}

   ri = rep_matrix(r, n_indicators) + rep_matrix(sigma_r, n_indicators) .* ri_z;

 // }
  // calculate detectable cases
 
}

model {
  // gaussian process priors

  sigma_r ~ normal(0, 0.5);
  r[1] ~ normal(0, r_prior_sd);
  sigma_rw ~ normal(0.1, 0.05);
 // rz ~ normal(0, 1);
  r[2:t] ~ normal(alpha_rw + beta_rw * r[1:(t - 1)], 0.1);

  
  
  

 // for(i in 1:t){
    for(n in 1:n_indicators){
        mean_r[,n ] ~ normal(ri[,n], sd_r[,n]);
        ri_z[,n] ~normal(0,1);
  //      ri_z[i, n] ~normal(0,1);

    }
 // }
  

  
}
