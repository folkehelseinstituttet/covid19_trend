source("data/prepare_data.R")
source("trend_models.R")

est_growth_rates <- function(data, r_len){

       r_msis <- estimate_r_t_direct(data$msis, len=r_len, weekend=TRUE, N=1000)
       r_perct_positive <- estimate_r_t_direct(data$perc_pos, len=r_len, weekends=TRUE, offsets=data$perc_pos$tot_N)
       r_hosp_main_cause <- estimate_r_t_direct(data$hosp_main_cause, len=r_len, weekend=TRUE, N=1000)# SHUS
       r_selv <- estimate_r_t_direct(data$self_tests, r_len,weekend=TRUE)
       r_sympto <- estimate_r_t(data$sympto_inc$estimated_incidence, (r_len +1):(length(data$sympto_inc$dates)-r_len) ,r_len, data$sympto_inc$dates, cores=1, type="lm", rand=100)
       r_with_covid19 <- estimate_r_t(data$with_covid_inc$estimated_incidence, (r_len +1):(length(data$with_covid_inc$dates)-r_len) ,r_len, data$with_covid_inc$dates, cores=1, type="lm", rand=100)
       r_norsyss <- estimate_r_t(data$inc_norsyss$estimated_incidence, (r_len +1):(length(data$inc_norsyss$dates)-r_len) ,r_len, data$inc_norsyss$dates, cores=1, type="lm", rand=100)
       r_waste <- estimate_r_t(data$inc_wastewater$estimated_incidence, (r_len+1):(length(data$inc_wastewater$dates)-r_len) , r_len,data$inc_wastewater$dates , type="lm", rand=100)
       r_deaths <- estimate_r_t(data$inc_deaths$estimated_incidence, (r_len +1):(length(data$inc_deaths$dates)-r_len) ,r_len, data$inc_deaths$dates, cores=1, type="lm", rand=100)%>% filter(date > as.Date("2020-03-17"))
       #Moba data comes in two separate waves
       moba_phase1 <- data$moba[1:86]
       moba_phase2 <- data$moba[113:nrow(data$moba)]
       r_moba_phase1 <- estimate_r_t_direct(moba_phase1, len=r_len, weekends=FALSE, offsets=moba_phase1$tot_N)
       r_moba_phase2 <- estimate_r_t_direct(moba_phase2, len=r_len, weekends=FALSE, offsets=moba_phase2$tot_N)
       r_moba <- rbind(r_moba_phase1, r_moba_phase2)

       to_plot <- rbindlist(list(r_msis %>% mutate(label="cases"),
                            r_perct_positive %>% mutate(label="perc_pos"),
                            r_hosp_main_cause %>% mutate(label="hosp"),
                            r_sympto %>% mutate(label="sympto"),
                            r_with_covid19 %>% mutate(label="with_covid"),
                            r_waste %>% mutate(label="wastewater"), 
                            r_moba %>% mutate(label="moba", date=as.Date(date)), 
                            r_norsyss %>% mutate(label="norsyss"),
                            r_selv %>% mutate(label="rdt"), 
                            r_deaths %>% mutate(label="deaths")), fill=TRUE) #%>% filter(date > as.Date("2022-12-01"))
       return(to_plot %>% mutate(r_len=r_len))
}

data <- get_data()


#Prevalence to incidence
sympto_prev <- est_prev_double(data$sympto$N,
                        data$sympto$symptoms,
                        data$sympto$tested,
                        data$sympto$positive)
data$sympto_prev <- sympto_prev
data$sympto_inc <- inc_from_prev(sympto_prev$prev,  sympto_prev$prev_sd, start_dates = data$sympto$date -7, end_dates=data$sympto$date, prob_detectable=data$prob_detectable_pcr)
data$with_covid_inc <- inc_from_prev(data$with_covid$prev, data$with_covid$prev_sd, start_dates=data$with_covid$date -1, end_dates=data$with_covid$date, prob_detectable=data$prob_detectable_pcr)
data$inc_norsyss <- inc_from_prev(data$norsyss$frac, data$norsyss$se, start_dates=data$norsyss$date -6, end_dates=data$norsyss$date, prob_detectable=data$prob_detectable_pcr)
data$inc_wastewater <- inc_from_prev(data$wastewater$prev, data$wastewater$prev_sd, start_dates=data$wastewater$date -1, end_dates=data$wastewater$date, prob_detectable=data$prob_detectable_waste)
data$inc_deaths <- diasagg(data$deaths$I, start_dates=data$deaths$date , end_dates=data$deaths$date+6,gp_m=0.3)

saveRDS(data, "data_with_incidence.rds")


for(r_len in c(8, 10, 12, 15)){
       gr <- est_growth_rates(data, r_len)
       saveRDS(gr, glue::glue("trend_estimates_{r_len}.rds"))
       del_shift <- est_delays_max_cor(gr)
       delays <- del_shift$delays
       shifted <- del_shift$shifted
       saveRDS(delays, glue::glue("delays_{r_len}.rds"))
       shifted[, year:=year(date)]
       shifted[, day:=lubridate::yday(date) + as.Date("2024-01-01")]
       shifted[, grp:=label]
       shifted[label=="moba" & date < as.Date("2022-04-15"), grp:="moba1"]
       shifted[label=="moba" & date > as.Date("2022-04-15"), grp:="moba2"]

       saveRDS(shifted, glue::glue("trend_estimates_shifted_{r_len}.rds"))
}





