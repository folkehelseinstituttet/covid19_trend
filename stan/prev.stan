
data {
  int<lower=1> N;
  int fork[N];
  int tot[N];
  int pos[N];
  int tot_tested[N];
}

transformed data{

}


parameters{
  real<lower=0, upper=1> p_fork[N];
  real<lower=0, upper=1> p_pos[N];
}

transformed parameters{
}

model {
      fork ~ binomial(tot, p_fork);
      pos ~ binomial(tot_tested, p_pos);	
}

generated quantities {
  real <lower=0, upper=1> p_covid[N];
  for(i in 1:N){
    p_covid[i] = p_fork[i]*p_pos[i];
  }	       

}
