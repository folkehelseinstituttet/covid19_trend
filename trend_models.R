library(ggplot2)
library(data.table)
library(dplyr)

#odin_model <- odin.dust::odin_dust("prevmod.R")
prev_mod <- rstan::stan_model("stan/prev.stan", auto_write=TRUE)
inc2prev <- rstan::stan_model("stan/inc2prev_simple.stan", auto_write=TRUE)
disagg_mod <- rstan::stan_model("stan/disagg.stan", auto_write=TRUE)




inc_from_prev <- function(prev, prev_sd, interval="day", start_dates, end_dates, prob_detectable=NULL, n_parts=400, n_threds=2, gp_m=0.3){

  start_times <- as.numeric(difftime(start_dates, start_dates[1], units=interval))+1
  end_times <- as.numeric(difftime(end_dates, start_dates[1], units=interval))+1

  stan_data <- list(
    ut=0,
    prev=prev,
    prev_sd2=prev_sd^2,
    t=max(end_times),
    obs=length(prev),
    prev_stime=start_times,
    prev_etime=end_times,
    pbt=nrow(prob_detectable),
    prob_detect_mean=rev(prob_detectable$mean),
    prob_detect_sd=rev(prob_detectable$sd),
    lengthscale_alpha=0.5,
    lengthscale_beta=0.5,
    M=round(gp_m*max(end_times)),
    L=2,
    inc_zero=0.02
  )

  res <- rstan::sampling(inc2prev, 
                data=stan_data,
                iter=1000,
                chains=2,
                include=TRUE,
                  cores = 4)
  

  infs <- rstan::extract(res, "infections")$infections

  return(list(estimated_incidence=infs, 
  dates=seq(start_dates[1], end_dates[length(end_dates)], by=interval)))
}


plot_inc <- function(inc){

  df <- as.data.frame(t(inc$estimated_incidence))
  df$date <- inc$dates
  df <- melt(df, id.vars="date")
  ggplot(df) + geom_line(aes(x=date, y=value, group=variable)) + ylab("Incidence")
}

estimate_r_t_direct <- function(cases, len, weekends, N=1000, offsets=FALSE){
  if(weekends){
    r_t <- estimate_r_t(cases$I, (len+1):(nrow(cases)-len), len, cases$date,lubridate::wday(cases$date) %in% c(7,1), N=N, offsets = offsets)
  }else{
    r_t <- estimate_r_t(cases$I, (len+1):(nrow(cases)-len), len, cases$date,N=N, offsets=offsets)
  }
}


inner_r_model <- function(I, len, weekends,N=100){
      if(length(weekends)>1){
        m <- lm(log(I)~t + weekend,data=data.frame(I=I,t=1:length(I), weekend=weekends))
      }else{
            m <- lm(log(I)~t,data=data.frame(I=I,t=1:length(I))) #MASS::glm.nb(I~t ,data=data.frame(I=I,t=1:len))
      }
      growth <- coef(m)[2]
      growth_sd <- summary(m)$coefficients[2,2]
      return(rnorm(N, mean=growth, sd=growth_sd))
}

inner_r_model_nb <- function(I, len, weekends,N=500, offset=FALSE){
      if(length(weekends)>1){
         if(any(offset)){
         m <-MASS::glm.nb(I~t + weekend + offset(log(tot_N)) ,data=data.frame(I=I,t=1:length(I),weekend=weekends, tot_N=offset))
          }else{
        m <- MASS::glm.nb(I~t + weekend,data=data.frame(I=I,t=1:length(I), weekend=weekends))
          }
      }else{
          if(any(offset)){
         m <-MASS::glm.nb(I~t + offset(log(tot_N)) ,data=data.frame(I=I,t=1:length(I), tot_N=offset))
      }else{
        m <-MASS::glm.nb(I~t ,data=data.frame(I=I,t=1:length(I)))
      }
      }
      growth <- coef(m)[2]
      growth_sd <- summary(m)$coefficients[2,2]
      return(rnorm(N, mean=growth, sd=growth_sd))
}

estimate_r_t <- function(hist, times, len, dates,  weekends=0, cores=1, N=100, type="nb", offsets=FALSE, rand=FALSE){
  # Estimate r_t centred on each date with "len" days before and after
  if(length(dim(hist))==0 | is.null(hist)) hist <- matrix(hist, ncol=length(hist), nrow=1)
  
  

  if(rand){
    I <- hist
    new_Is <- matrix(0, ncol=dim(I)[2], nrow=rand)
    for(i in 1:rand){
      a <- cbind(sample(1:dim(I)[1], dim(I)[2], replace=T), 1:dim(I)[2])
      new_Is[i, ] <- I[a]
    }
    hist <- new_Is
      }
  n_sims <- dim(hist)[1]


  r_t <- list()
  for(t in times){
    if(type=="nb"){
      rs <- unlist(parallel::mclapply(1:n_sims, function(i){inner_r_model_nb(hist[i,(t-len):(t+len)], len, if(length(weekends)>1) weekends[(t-len):(t+len)] else 0, N=N, offset=if(length(offsets)> 1) offsets[(t-len):(t+len)] else FALSE)}, mc.cores=cores))
    }else if(type=="lm"){
      rs <- unlist(parallel::mclapply(1:n_sims, function(i){inner_r_model(hist[i,(t-len):(t+len)], len, if(length(weekends)>1) weekends[(t-len):(t+len)] else 0, N=N)}, mc.cores=cores))
    }
    r_t[[t]] <- data.table(t=t,
                           r=median(rs),
                           r_l=quantile(rs, 0.2),
                           r_u=quantile(rs,0.8),
                           r_l95=quantile(rs, 0.025),
                           r_u95=quantile(rs,0.975),
                           date=dates[t],
                           I=median(hist[,t]),
                           Imin=quantile(hist[,t], 0.2),
                           Imax=quantile(hist[,t], 0.8)
                           )
  }
  rbindlist(r_t)
}


est_prev_double <- function(N_tot, N_sub_1, tot_tested, N_pos){
  stan_data <- list(
    N = length(N_tot),
    fork = N_sub_1,
    tot = N_tot,
    pos = N_pos,
    tot_tested = tot_tested
  )
  
  fit = rstan::sampling(prev_mod,
               data=stan_data,
               iter=2500,
               chains=2,
               include=TRUE,
               cores = 4)

  p_cov <- rstan::extract(fit, "p_covid")$p_covid
  
  return(list(prev=colMeans(p_cov), prev_sd=matrixStats::colSds(p_cov)))
}

interp_week <- function(P, dates){
  t <- dates - dates[1] 
  spline <- smooth.spline(t, P)
  predict(spline, (t[1]+1):t[length(t)])$y
}


plot_inc <- function(inc){

  dfs  <- list()
  for(i in 1:dim(inc$estimated_incidence)[1]){
    dfs[[i]] <- data.frame(i=inc$estimated_incidence[i,], 
                            date=inc$dates, sim=i)
  }
  comb <- rbindlist(dfs)

  ggplot(comb) + geom_line(aes(x=date, y=i, group=sim))

}


plot_histories <- function(hists, var="I", labels=NULL, norm=1){

  if(is.null(labels)) labels<- 1:length(hists)
  j <- 1
  all_dfs <- list()
  for(hist in hists){
    dfs<-list()
    sims <- dim(hist)[2]
    for(sim in 1:sims){
      dfs[[sim]] <-data.frame(t=hist[1, sim,],
                              I=hist[2, sim,],
                              P=hist[3, sim,],
                              i=hist[5,sim,],
                              beta=exp(hist[4, sim,]),
                              R=exp(hist[4, sim,])*5,
                              sim=paste0(labels[j], sim),
                              lab=labels[j])
    }
    
    all_dfs[[length(all_dfs) + 1]] <- rbindlist(dfs)
    j <- j+1
  }
  df <- rbindlist(all_dfs)
  df[, factor_group:=paste(sim)]
  
  ggplot(df) + geom_line(aes(x=t, y=get(var)/norm, group=sim, color=factor(lab))) + ylab(var)
  
}




diasagg <- function(A, interval="day", start_dates, end_dates, n_parts=400, n_threds=2, gp_m=1, inc_zero=100, L=5){

  start_times <- as.numeric(difftime(start_dates, start_dates[1], units=interval))+1
  end_times <- as.numeric(difftime(end_dates, start_dates[1], units=interval))+1

  stan_data <- list(
    ut=0,
    A=A,
    t=max(end_times),
    n_obs=length(A),
    prev_stime=start_times,
    prev_etime=end_times,
    lengthscale_alpha=0.5,
    lengthscale_beta=0.5,
    M=round(gp_m*max(end_times)),
    L=L,
    inc_zero=inc_zero
  )

  res <- rstan::sampling(disagg_mod, 
                data=stan_data,
                iter=1000,
                chains=2,
                include=TRUE,
                  cores = 4)
  

  infs <- rstan::extract(res, "infections")$infections

  return(list(estimated_incidence=infs, 
  dates=seq(start_dates[1], end_dates[length(end_dates)], by=interval)))
}




meta_analysis<- function(rs, model=meta_mod_random_rw, stan_control=list(max_treedepth=12, adapt_delta=0.9), iter=1500, chains=4, cores=4, r_prior_sd=0.3){
  rs <- cbind(rs, get_mean_sd(rs))

  groups <- unique(rs$grp)
  n_inds <- length(groups)
  dates <- sort(unique(rs$date))
  t <- length(dates)

  means <- matrix(0, nrow=t, ncol=n_inds)
  sd <- matrix(1e6, nrow=t, ncol=n_inds)
  for(i in 1:t){
    for(g in 1:n_inds){
      group <- groups[g]
      x_date <- dates[i]
      if(group %in% rs[date==x_date, grp]){
         means[i,g] <- rs[grp==group & date==x_date, mean]
        sd[i,g] <- rs[grp==group & date==x_date, sd]
      }
    }
  }


  stan_data <- list(
    mean_r=means,
    sd_r=sd,
    t=dim(means)[1],
    n_indicators=n_inds,
    r_prior_sd=r_prior_sd

  )

  res <- rstan::sampling(model,
                data=stan_data,
                iter=iter,
                chains=chains,
                include=TRUE,
                control=stan_control,
                  cores = cores)


  r <- rstan::extract(res, "r")$r
  r_est <-   data.table(date=dates, r=apply(r, 2, mean), r_l95=apply(r, 2, function(x) quantile(x, 0.025)), r_u95=apply(r, 2, function(x) quantile(x, 0.975)))



  

  
 if("sigma_r[1]" %in% names(res)){
    I2 <- calc_H2_I2(res, means,sd)
 } else{
     I2 <- NULL
 }
  return(list(r_est=r_est, I2=I2, r_samples=r, full_results=res))
}





rest <- function(){

  r_fixed <- rstan::extract(res_fixed, "r")$r
  r_random <- rstan::extract(res_random, "r")$r
  r_random_gp <- rstan::extract(res_random_gp, "r")$r
  r_random_gp2 <- rstan::extract(res_random_gp2, "r")$r
  


  r_out <-rbind(data.table(date=dates, r=apply(r_fixed, 2, mean), r_l95=apply(r_fixed, 2, function(x) quantile(x, 0.025)), r_u95=apply(r_fixed, 2, function(x) quantile(x, 0.975))) %>% mutate(model="fixed"),
                data.table(date=dates, r=apply(r_random, 2, mean), r_l95=apply(r_random, 2, function(x) quantile(x, 0.025)), r_u95=apply(r_random, 2, function(x) quantile(x, 0.975))) %>% mutate(model="random"),
           #     data.table(date=dates, r=apply(r_random_gp, 2, mean), r_l95=apply(r_random_gp, 2, function(x) quantile(x, 0.025)), r_u95=apply(r_random_gp, 2, function(x) quantile(x, 0.975))) %>% mutate(model="random_gp"),
          #      data.table(date=dates, r=apply(r_random_gp2, 2, mean), r_l95=apply(r_random_gp2, 2, function(x) quantile(x, 0.025)), r_u95=apply(r_random_gp2, 2, function(x) quantile(x, 0.975))) %>% mutate(model="random_gp2"),
              
                )
                
  r_sigma <- rstan::extract(res_random_rw, "sigma_r")$sigma_r
  
  w_k <- sd^2/( sd^2 + colMeans(r_sigma)^2)
  sd[sd>1e5] <- NaN
  rms <- rowMeans(sd, na.rm = TRUE)

  plot(colMeans(r_sigma)/rms)


  ggplot(rs) + geom_ribbon(aes(x=date, ymin=r_l95, ymax=r_u95, fill=grp), alpha=0.6)
  ggsave("prob.png")
  ggplot(r_out, aes(x=date, y=r, ymin=r_l95, ymax=r_u95, fill=model, color=model)) + geom_line() + geom_ribbon(alpha=0.7) 
  ggsave("sol1.png")
  ggplot(r_out, aes(x=date, y=r, ymin=r_l95, ymax=r_u95)) + geom_line() + geom_ribbon(alpha=0.5) + geom_ribbon(aes(ymin=r_l95, ymax=r_u95, fill=grp), alpha=0.6, data=rs) + facet_wrap(~model) + ylab("r")
  ggsave("sol2.png")  




  return(list(estimated_incidence=infs))
}

calc_H2_I2 <- function(res, means,sd, probs=c(0.025, 0.5, 0.975)){
  sd[sd>1e5] <- NaN 
  r_sigma <- rstan::extract(res, "sigma_r")$sigma_r
  r <- rstan::extract(res, "r")$r

  n  <- rowSums(!is.na(sd))
  geom_avg_var <-  (matrixStats::rowProds(sd^2, na.rm=T))^(1/n)
  geom_avg_var  <- matrix(rep(geom_avg_var, each=nrow(r_sigma)), ncol=dim(r_sigma)[2], nrow=dim(r_sigma)[1])
  alt_I2 <- r_sigma^2/(r_sigma^2 + geom_avg_var)

#  m

  # means_k <- means
  # means_k[means_k==0] <- NA
  # theta_k  <- t(matrix(rep(means_k, nrow(r_sigma)), ncol=nrow(r_sigma), nrow=length(means)))

  # w <-1 / sd**2
  # Qs <- matrix(0, nrow=dim(means_k)[1], ncol=dim(r_sigma)[1])
 
  # Is <- matrix(0, nrow=dim(means_k)[1], ncol=dim(r_sigma)[1])
  # for(i in 1:dim(means_k)[1]){
  #   theta_k <- means_k[i,]#t(matrix(rep(means_k[i,], nrow(r_sigma)), ncol=nrow(r_sigma), nrow=length(means_k[i,])))
  #   for(j in 1:dim(r)[1]){
  #     Q_i <- sum(1/(1/w[i,] + r_sigma[j,i]^2)*(theta_k-  r[j,i])^2, na.rm = TRUE)
  #     Q_i[Q_i < 0 ] <- 0
  #     K <- sum(!is.na(theta_k))
  #     I_i <- (Q_i - (K - 1))/Q_i
  #     I_i[I_i<0] <- 0
  #     Qs[i,j] <- Q_i
  #     Is[i,j] <- I_i
  #   }
  # }

  
  return(matrixStats::colQuantiles(alt_I2, probs=probs))
}



get_mean_sd <- function(x1){
  return(data.frame(mean=x1$r, sd=(x1$r - x1$r_l95)/1.96))

}

sample_cor <- function(cor, n, N=1){
  z <- atanh(cor)
  se <- 1/sqrt(n-3)

  transformed_sims <- rnorm(N, mean=z, sd=se)
  return(tanh(transformed_sims))
}


get_cor <- function(x1, x2, N=5000){
  
  date_overlap <- as.Date(intersect(x1$date, x2$date), origin=as.Date("1970-01-01"))
  if(length(date_overlap)< 2){
     return(data.frame(cor=NA, cor_sd=NA, 
      l50=NA, u50=NA,
      l80=NA, u80=NA,
       l95=NA, u95=NA))
  }
  mean_sd_1 <- get_mean_sd(x1 %>% filter(date %in% date_overlap))
  mean_sd_2 <- get_mean_sd(x2 %>% filter(date %in% date_overlap))
 
  cors <- c()
  for(i in 1:N){
    x1_i <- rnorm(nrow(mean_sd_1), mean=mean_sd_1$mean, sd=mean_sd_1$sd)
    x2_i <- rnorm(nrow(mean_sd_2), mean=mean_sd_2$mean, sd=mean_sd_2$sd)
    cors <- c(cors, sample_cor(cor(x1_i, x2_i), length(x1_i), N=1))
  } 
 
  return(data.frame(cor=mean(cors), cor_sd=sd(cors), 
  l50=quantile(cors, c(0.25)), u50=quantile(cors, c(0.75)),
  l80=quantile(cors, c(0.01)), u80=quantile(cors, c(0.9)),
  l95=quantile(cors, c(0.0025)), u95=quantile(cors, c(0.975))))
}



est_delays_max_cor <- function(data){
  delays <- list()
    for(lab in unique(data$label)){
    res <- data.frame()
     for(delay in -15:20){ 
       hosps <- data %>% filter(label=="hosp")
      cases <- data %>% filter(label==lab)
      cases$date <- cases$date + delay
      date_overlap <- as.Date(lubridate::intersect(hosps$date, cases$date), origin="1970-01-01")
       res <- rbind(res, data.frame(delay=delay, label=lab, cor=cor(hosps %>% filter(date %in% date_overlap) %>% select(r), cases %>% filter(date %in% date_overlap) %>% select(r), use="pairwise.complete.obs")))
    }
     delays[[lab]] <- list(delay= res %>% arrange(desc(r)) %>% slice(1) %>% pull(delay), max_cor=res %>% arrange(desc(r)) %>% slice(1) %>% pull(r))
  }
  shifted <- copy(data)
  setDT(shifted)
  for(key in names(delays)){
    shifted[label==key, date:=date+delays[[key]]$delay] 
  }
  return(list(delays=delays, shifted=shifted))

}

cor_matrix <- function(res){
  cor_mat <- list()
  for(lab in unique(res$label)){
    for(lab2 in unique(res$label)){
      if(lab==lab2) {
      cor_mat[[length(cor_mat) + 1]] <- data.frame(cor=1, cor_sd=0, l50=1, u50=1,
      u80=1, l80=1, u95=1, l95=1) %>% mutate(label1=lab, label2=lab2)

      }else{
      x1 <- res %>% filter(label==lab)
      x2 <- res %>% filter(label==lab2)
      coro <- get_cor(x1, x2)
      cor_mat[[length(cor_mat) + 1]] <- coro %>% mutate(label1=lab, label2=lab2)
      }
    }
  }
  return(rbindlist(cor_mat, use.names=TRUE))
}



get_gt <- function(combined, r_sim){
  # Data sources listed in supplementary materials
  gts <- list(
    "wuhan" = list(m=4.95, cv=1.2),
    "alpha" = list(m=4.35, cv=0.73),
    "delta" = list(m = 3.65, cv=0.7),
    "omicron"= list(m=2.99, cv=0.65)
  )

  # Taken form   #variants <- fread("https://raw.githubusercontent.com/EU-ECDC/Respiratory_viruses_weekly_data/main/data/variants.csv")
  transitions <- list(
    list(from="wuhan", to="alpha", date=as.Date("2021-02-15"), width=10),
    list(from="alpha", to="delta", date=as.Date("2021-06-28"), width=10),
    list(from="delta", to="omicron", date=as.Date("2021-12-20"), width=10)
    
  )
  i <- 1
  transition <- transitions[[i]]
  scales <- c()
  rates <- c()
  for(d in as.list(unique(combined$date))){
    kern <- as.numeric((d - transition$date)/transition$width)
    frac  <- 1/(1 + exp(-kern))
    mean <- gts[[transition$from]]$m*(1-frac) + gts[[transition$to]]$m*frac
    cv <- gts[[transition$from]]$cv*(1-frac) + gts[[transition$to]]$cv*frac
    tmp <- epitrix::gamma_mucv2shapescale(mean,cv)# sd/mean)
    scales <- c(scales, tmp$scale)
    rates <- c(rates, 1/tmp$shape)
    if(frac > 0.999 & i < length(transitions)){
      i <- i + 1
      transition <- transitions[[i]]
    }
  }
  combined$scale <- rep(scales, each=dim(r_sim)[1])
  combined$rate <- rep(rates, each=dim(r_sim)[1])
  return(combined)

}
