dt <- user(0.1)
initial(time) <- 0
update(time) <- (step + 1) * dt
steps_per_day <- 1/dt

update(I) <- I + nI - nR
update(P) <- P + nI - nP
update(logbeta) <- if(beta_in[1]!=0) log(beta_in[step]) else ( if(time%%beta_time==0) logbeta + rnorm(0, sigma_log_beta)else logbeta)
update(inc)<-if(time %% steps_per_day==0) nI else inc+nI
initial(I) <- I_ini
initial(P) <- P_ini
initial(logbeta) <- beta_ini
initial(inc)<-0
beta <- exp(logbeta)

nI <- rbinom(N, 1 - exp(-beta*I/N*dt))
nR <- rbinom(I, 1 - exp(-1/Ti*dt))
nP <- rbinom(P, 1- exp(-1/Tp*dt))

N <- user(5.4e6)
Ti <- user(5)
Tp <- user(6)
I_ini <- user(10)
P_ini <- user(10)
sigma_log_beta <- user(0.01)
beta_ini <- user(0.1)
N_steps <- user()
beta_in[] <- user(0)
dim(beta_in) <- N_steps
beta_time<-user(1)
