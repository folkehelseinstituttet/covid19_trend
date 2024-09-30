library(data.table)
library(dplyr)

max_date <- as.Date("2024-01-01")


convolve <- function(x1, x2){
  max_n <- length(x1) + length(x2) -1
  y <- numeric(max_n)
  x1 <- c(x1, rep(0, max_n - length(x1)))
  x2 <- c(x2, rep(0, max_n - length(x2)))
  for(i in 1:max_n){
    for(j in 1:i){
      y[i] <- y[i] + x1[j] * x2[i-j+1]
    }
   
  }
  y
}



get_data <- function(max_date=as.Date("2024-01-01")){
  # Lab tests
  msis <- fread("data/msis_data.csv")[, date:=as.Date(date)]
  msis <- msis[date<max_date, .(I=sum(n)), by=.(date)]
  test_events <- readRDS("data/testhendelser.rds")

  perc_pos <- msis %>% left_join(test_events %>% filter(date <max_date)) %>% mutate(p=I/n, tot_N=n)





  hospital_data <- fread("data/covid_data.csv")
  hospital_data[, date:=as.Date(reference_date)]
  hospital_data[, n_hospital_main_cause:=N]
  nat_hospital <- fread("data/hosp_data.csv")
  nat_hospital <- nat_hospital %>% tidyr::complete(date=seq(min(date), max(date), by="1 day")) %>% tidyr::replace_na(list(I=0))
  
# Main cause hosps and hosp prev
sympto <- fread("data/symptometer.csv")
sympto[, WK:=paste0(substr(as.character(UKE),1,4), "-W",substr(as.character(UKE),5,6), "-1")]
sympto[, date:=ISOweek::ISOweek2date(WK)]
sympto <- sympto[, .(N=antallBesvarelser, symptoms=antallForkj, tested=antallMedForkjTestetPos + antallMedForkjTestetNeg,
              positive=antallMedForkjTestetPos,date)]

with_covid <- fread("data/with_covid.csv")
with_covid <- with_covid %>% tidyr::replace_na(list(main=0, not_main=0))%>% mutate( not_main=COVID.19 - not_main) %>%
  select(
    date,
    I=not_main,
   tot_N=N) %>% filter(date < as.Date("2023-09-26")) %>% mutate(I=ifelse(I<0, 0, I), date=as.Date(date)) 
  with_covid <- with_covid%>% mutate(prev=I/tot_N,
       prev_sd=sqrt(prev*(1-prev))/sqrt(tot_N))%>%tidyr::replace_na(list(prev_sd=0)) %>%filter(!is.na(prev) & date >= as.Date("2020-03-03"))
with_covid[prev_sd<=0, prev_sd:=0.001]



  # Syndromic surveillance
  norsyss <- readRDS("data/202401-ro-2024-03-07.RDS")
  combined_norsyss <- norsyss[, .(I=sum(consultations_icpc2group_n), tot_N=mean(consultations_all_n)), by=.(isoyearweek, date)][date >=as.Date("2020-03-08")]
  combined_norsyss[, frac:=I/tot_N]
  combined_norsyss[, se:=sqrt(frac*(1-frac))/sqrt(tot_N)]


  #Wastewater
  waste <- fread("data/wastewater.csv")

  waste[, WK:=paste0(gsub("-","-W",isoweek),"-4")]
  waste[, date:=ISOweek::ISOweek2date(WK)]
  waste[, cov_fc_ratio:=get("SARS-CoV-2 RNA")]
  waste[, day:= as.numeric(date - min(date))]
  waste <- waste[, .(date, amount=cov_fc_ratio, day)]
  waste$prev <- waste$amount / (10*max(waste$amount))
  waste$prev_sd <- waste$prev*0.2

  #MoBA cohort
  moba <- fread("data/MoBa_insidens.csv")
  moba <- moba %>% select(date=Dato, I=Pos, tot_N=Freq)


  # Self-reported RATs
  rat <- fread("data/selvtest.csv")
    rat$date <- as.Date(rat$date)
  

  #Mortality
  deaths <- readRDS("data/daar.rds") %>% mutate(WK=paste0(gsub("-","-W",isoyearweek),"-3"),
                                                date=ISOweek::ISOweek2date(WK), 
                                                I=n)

  #Detection probability for PCR tests by day since infection. Data from : https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-021-01982-x
  prob_detectable <- fread("data/prob_detectable.csv")

  prob_detectable <- melt(
      copy(prob_detectable),
      value.name = "p", id.vars = "sample"
  )
  prob_detectable[, time := as.numeric(as.character(variable))]
  prob_detectable <- prob_detectable[, .(
    median = median(p),
    mean = mean(p),
    sd = sd(p)
  ),
  by = time
  ]
  

  # Time between infection and Viral load in gastrointestinal tract. Distribution values from: https://www.sciencedirect.com/science/article/pii/S0043135421004504
  shape_rate_onset <- epitrix::gamma_mucv2shapescale(5.3, 3.2/5.3)
  dens_sympt_onset <- distcrete::distcrete("gamma", 1, shape=shape_rate_onset$shape, scale=shape_rate_onset$scale)$d(0:50)
  shape_rate_vl <- epitrix::gamma_mucv2shapescale(6.7, 6.9/6.7)
  dens_vl <- distcrete::distcrete("gamma", 1, shape=shape_rate_vl$shape, scale=shape_rate_vl$scale)$d(0:50)
  combined <- convolve(dens_sympt_onset, dens_vl)

  waste_sld <- data.frame(time=1:length(combined), mean=combined, sd=combined*0.1)


  #sld <- CVXR::value(CVXR::conv(dens_sympt_onset, dens_vl))



  data <- list(msis=msis,
            hosp_main_cause=nat_hospital,
            sympto=sympto,
            with_covid=with_covid,
            wastewater=waste,
            norsyss=combined_norsyss,
            moba=moba,
            deaths=deaths,
            perc_pos=perc_pos,
            self_tests=rat,
            prob_detectable_pcr=prob_detectable,
            prob_detectable_waste=waste_sld)

  return(data)
}








