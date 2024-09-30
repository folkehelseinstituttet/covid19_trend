
source("trend_models.R")

metas <- list()
stan_controls <- list(list(max_treedepth=12, adapt_delta=0.9), list(max_treedepth=12, adapt_delta=0.95), list(max_treedepth=15, adapt_delta=0.98))
meta_mod_random_rw <- rstan::stan_model("stan/meta_analysis_rw.stan", auto_write=TRUE)
models <- list(meta_mod_random_rw)
for(r_len in c(8, 10, 12, 15)){
    shifted <- readRDS(glue::glue("trend_estimates_shifted_{r_len}.rds"))
    short_shifted <- shifted #%>% filter(date < as.Date("2020-12-01"))
    metas <- lapply(1:length(models), function(i) meta_analysis(short_shifted, model=models[[i]], stan_control=stan_controls[[i]], r_prior_sd=0.1, iter=2000))
#metas <- lapply(1:3, function(i) meta_analysis(short_shifted, model=models[[i]], stan_control=stan_controls[[i]], r_prior_sd=0.1, iter=2000))
    saveRDS(metas, glue::glue("meta_analysis_{r_len}.rds"))
}