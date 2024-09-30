
library(ggplot2)
library(data.table)
library(dplyr)
source("trend_models.R")

common_theme <- theme_minimal() #+ theme(text = element_text(size=20))

# Labels
t_cases <- "Cases"
t_perc_pos <- "Proportion Positive"
t_hosp <- "Hospital Admissions"
t_sympto <- "Proportion Positive - Survey"
t_with_covid <- "Hospital Prevalence"
t_wastewater <- "Wastewater"
t_moba <- "Proportion positive - Cohort"
t_norsyss <- "GP consultations"
t_rdt <- "Positive RAT tests"
t_deaths <- "Deaths"
order <- c(t_wastewater, t_sympto, t_moba, t_norsyss, t_rdt, t_with_covid, t_cases, t_perc_pos, t_hosp, t_deaths)

ylabels <- list(
  "wastewater"="Normalised Concentration",
  "sympto"="Proportion positive",
  "moba"="Proportion testing positive",
  "norsyss"="Proportion of consultations",
  "self_tests"="Positive antigen tests",
  "with_covid"="Prevalence",
  "msis"="Cases",
  "perc_pos"="Percent positive",
  "hosp_main_cause"="Weekly admissions",
  "deaths"="Weekly Deaths"
)

data <- readRDS("data_with_incidence.rds")

 # Figure 1 - All data sources
frac_indicators <- c("sympto", "with_covid", "norsyss", "moba", "perc_pos")
frac_data <- list(
  sympto=data.frame(date=data$sympto$date, p=data$sympto_prev$prev, pmin=data$sympto_prev$prev-1.96*data$sympto_prev$prev_sd, pmax=data$sympto_prev$prev+1.96*data$sympto_prev$prev_sd),
  norsyss=data.frame(date=data$norsyss$date, p=data$norsyss$frac, pmin=data$norsyss$frac - 1.96*data$norsyss$se, pmax=data$norsyss$frac + 1.96*data$norsyss$se)
)
plot_data <- data
plot_data$wastewater <- plot_data$wastewater %>% mutate(I=amount)
daily <- c("msis", "hosp_main_cause", "with_covid", "moba", "self_tests", "perc_pos")
plots <- list()
for(key in names(ylabels)){
  if(key %in% frac_indicators){
      if(key %in% daily){
        plot_data[[key]]$week <- ISOweek::ISOweek(plot_data[[key]]$date)
        pl <- plot_data[[key]] %>% group_by(week) %>% summarise(I=sum(I), tot_N=sum(tot_N), date= min(date) + 3)
        pl <- pl %>% mutate(p=I/tot_N, se=sqrt(p*(1-p))/sqrt(tot_N), pmin=p - 1.96*se, pmax=p+1.96*se)
       
      }else{
        pl <- frac_data[[key]]
      }
    q <- ggplot(pl, aes(x=date, y=p, ymin=pmin, ymax=pmax)) + geom_line() + geom_ribbon(alpha=0.3) + ylab(ylabels[[key]]) + xlab("Date") + scale_y_continuous(labels = scales::percent, trans="log10") 
    }else{
      if(key %in% daily){
        plot_data[[key]]$week <- ISOweek::ISOweek(plot_data[[key]]$date)
        pl <- plot_data[[key]] %>% group_by(week) %>% summarise(I=sum(I),date= min(date) + 3)
        }else{
          pl <- plot_data[[key]]
        }
   q <- ggplot(pl) + geom_line(aes(x=date, y=I)) + ylab(ylabels[[key]]) + xlab("Date") + scale_y_continuous(trans="log10")
  }
  plots[[length(plots) + 1]] <- q + xlim(c(as.Date("2020-03-01"), as.Date("2024-01-01"))) + common_theme + theme(text = element_text(size=14))
}
plot_grid <- purrr::partial(cowplot::plot_grid, nrow=2, ncol=5, labels=c(t_wastewater,t_sympto,t_moba,t_norsyss,   t_rdt,  t_with_covid, t_cases,t_perc_pos, t_hosp, t_deaths), hjust=0, vjust=1.5, scale=0.95, label_size=16)
do.call(plot_grid, plots)
ggsave("paper/figures/all_indicators.png", width=18, height=10)
ggsave("paper/figures/all_indicators.pdf", width=18, height=10)






plot_rs <- function(r_len){
  unshifted_rs <- readRDS(glue::glue("trend_estimates_{r_len}.rds"))
  shifted <- readRDS(glue::glue("trend_estimates_shifted_{r_len}.rds"))
  delays <- readRDS(glue::glue("delays_{r_len}.rds"))
  metas <- readRDS(glue::glue("meta_analysis_{r_len}.rds"))
  unshifted_rs <- unshifted_rs %>% mutate(nice_label = recode(label, cases=t_cases, perc_pos=t_perc_pos,"hosp"=t_hosp, "sympto"=t_sympto, with_covid=t_with_covid, wastewater=t_wastewater, moba=t_moba, "norsyss"=t_norsyss, rdt=t_rdt, deaths=t_deaths))
  shifted <- shifted %>% mutate(nice_label = recode(label, cases=t_cases, perc_pos=t_perc_pos,"hosp"=t_hosp, "sympto"=t_sympto, with_covid=t_with_covid, wastewater=t_wastewater, moba=t_moba, "norsyss"=t_norsyss, rdt=t_rdt, deaths=t_deaths))
  unshifted_rs$nice_label <- factor(unshifted_rs$nice_label, levels=order)
  shifted$nice_label <- factor(shifted$nice_label, levels=order)
  ggplot(shifted[date < as.Date("2024-01-01")], aes(x=day, y=r, ymin=r_l, ymax=r_u, group=grp, color=nice_label, fill=nice_label)) + geom_line() + geom_ribbon(alpha=0.6)  + 
          scale_y_continuous("Growth Rate", labels = scales::percent) + common_theme + facet_wrap(.~year) + theme(text = element_text(size=30), legend.position="bottom")  + scale_x_date(date_labels="%b") + labs(fill="Indicator", color="Indicator") + xlab("Date")
  ggsave(glue::glue("paper/figures/fig_growth_rates_{r_len}.png"), width=25, height=14)
  ggsave(glue::glue("paper/figures/fig_growth_rates_{r_len}.pdf"), width=25, height=14)



  # Shifting by yearly time-shifts and preparing data for plot
  yearly_shifted <- list()
  delay_year <- list()
  cors <- list()
  for(x_year in 2020:2023){
    del_shift <- est_delays_max_cor(unshifted_rs %>% filter(year(date)==x_year))
    delays_i <- del_shift$delays
    shifted_i <- del_shift$shifted
    yearly_shifted[[length(yearly_shifted) +1]] <- shifted_i
    delay_year[[length(delay_year) + 1]] <- data.frame(year=x_year, delay=as.numeric(lapply(delays_i, function(x) x$delay)), col=names(delays_i))
    df <- data.frame()
    for(col in unique(shifted_i$label)){
      if(col!="hosp"){
        cor <- get_cor(shifted_i %>% filter(label=="hosp"), shifted_i %>% filter(label==col))
        df <- rbind(df, data.frame(year=x_year, col=col, cor= cor$cor, cor_l=cor$l95, cor_u=cor$u95))
      }
    }
    cors[[length(cors)+1]] <- df
  }
  delay_all <- rbindlist(delay_year)
  cors <- rbindlist(cors)
  overall_delay <- data.frame(col=names(delays), delay=as.numeric(lapply(delays, function(x) x$delay)), year="All years")
  comb_shifted <- rbindlist(yearly_shifted)
  new <- comb_shifted %>% group_by(label) %>% tidyr::complete(date=seq(min(date), max(date), by="1 day")) %>% tidyr::fill(r, r_l, r_u, r_l95, r_u95)
  new_shifted <- new %>% group_by(label, date) %>% summarise(across(c(r, r_l, r_u, r_l95, r_u95), mean))
  df <- data.frame()
    for(col in unique(shifted$label)){
      if(col!="hosp"){
        cor <- get_cor(shifted %>% filter(label=="hosp"), shifted %>% filter(label==col))
        df <- rbind(df, data.frame(year=x_year, col=col, cor= cor$cor, cor_l=cor$l95, cor_u=cor$u95))
      }
    }

  df$year <- "All years"

  delay_all <- rbind(delay_all, overall_delay, fill=TRUE)
  delay_all[, size:=4]
  delay_all[year=="All years", size:=5]

  cors <- rbind(cors, df)
  cors$year <-  factor(cors$year, levels=c(2020, 2021, 2022, 2023,  "All years"))
  delay_all$year <- factor(delay_all$year, levels=c(2020, 2021, 2022, 2023,  "All years"))

  cors[, size:=4]
  cors[year=="All years", size:=5]

  delay_all <- delay_all %>% mutate(nice_label = recode(col, cases=t_cases, perc_pos=t_perc_pos, "hosp"=t_hosp, "sympto"=t_sympto, with_covid=t_with_covid, wastewater=t_wastewater, moba=t_moba, "norsyss"=t_norsyss, rdt=t_rdt, deaths=t_deaths))
  cors <- cors %>% mutate(nice_label = recode(col, cases=t_cases, perc_pos=t_perc_pos, hosp=t_hosp, "sympto"=t_sympto, with_covid=t_with_covid, wastewater=t_wastewater, moba=t_moba, "norsyss"=t_norsyss, rdt=t_rdt, deaths=t_deaths))
  df <- df %>% mutate(nice_label = recode(col, cases=t_cases, perc_pos=t_perc_pos,"hosp"=t_hosp, "sympto"=t_sympto, with_covid=t_with_covid, wastewater=t_wastewater, moba=t_moba, "norsyss"=t_norsyss, rdt=t_rdt, deaths=t_deaths))

  delay_all$nice_label <- factor(delay_all$nice_label, levels=order)
  df$nice_label <- factor(df$nice_label, levels=order)
  cors$nice_label <- factor(cors$nice_label, levels=order)
  colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "black")
  cols <- delay_all %>% count(col) %>% filter(n==2) %>% pull(col)

  delay_all_pl <- delay_all %>% filter(!(col %in% cols & year!="All years"))
  cors_pl <- cors %>% filter(!(col %in% cols & year!="All years") & !(col=="sympto" & year==2020))

  q1 <- ggplot(delay_all_pl %>% filter(col!="hosp"), aes(x=nice_label)) + geom_point(aes(y=delay, group=factor(year), color=factor(year), size=size), position = position_dodge(width = 0.60)) +  scale_size_identity() + common_theme + ylab("Delay(Days)") + scale_color_manual("Year", values=colors) + xlab("") +  theme(text = element_text(size=20)) + scale_x_discrete(labels = function(x) stringr::str_wrap(x, width=10))
  #+  geom_point(aes(x=col, y=delays, color="All years"), data=overall_delay, shape=95, size=20)# + geom_linerange(aes(ymin=cor_l, ymax=cor_u), position = position_dodge(width = 0.60)) 
  q2 <- ggplot(cors_pl, aes(x=nice_label)) + geom_point(aes(y=cor, group=factor(year), color=factor(year)), size=cors_pl$size, position = position_dodge(width = 0.60)) + geom_linerange(aes(ymin=cor_l, ymax=cor_u, group=factor(year), color=factor(year)), size=1.2, position = position_dodge(width = 0.60))  +  common_theme + ylab("Correlation") + scale_color_manual("Year", values=colors) + xlab("Indicator") +  theme(text = element_text(size=20)) + scale_x_discrete(labels = function(x) stringr::str_wrap(x, width=10)) #+  geom_point(aes(x=col, y=delays, color="Global"), data=overall_delay, shape=95, size=20)# + geom_linerange(aes(ymin=cor_l, ymax=cor_u), position = position_dodge(width = 0.60)) 

  main <- cowplot::plot_grid(q1 + theme(legend.position="none"), q2  + theme(legend.position="bottom"),  ncol=1, labels=c("A)", "B)"), scale=0.98)
  ggsave(glue::glue("paper/figures/fig_delay_cor_{r_len}.png"), width=18, height=12)
  ggsave(glue::glue("paper/figures/fig_delay_cor_{r_len}.pdf"), width=18, height=12)



  # Fig 4 - Correlation Matrix
  a <- cor_matrix(new_shifted)

  a <- a %>% mutate(label1_n = recode(label1, cases=t_cases, perc_pos=t_perc_pos,"hosp"=t_hosp, "sympto"=t_sympto, with_covid=t_with_covid, wastewater=t_wastewater, moba=t_moba, "norsyss"=t_norsyss, rdt=t_rdt, deaths=t_deaths),
                    label2_n = recode(label2, cases=t_cases, perc_pos=t_perc_pos,"hosp"=t_hosp, "sympto"=t_sympto, with_covid=t_with_covid, wastewater=t_wastewater, moba=t_moba, "norsyss"=t_norsyss, rdt=t_rdt, deaths=t_deaths))

  a$label1_n <- factor(a$label1_n, levels=order)
  a$label2_n <- factor(a$label2_n, levels=rev(order))

  a$cor_cat <- cut(a$cor, breaks=c(0, 0.2, 0.4, 0.7, 0.9, 1), labels=c( "0-0.2", "0.2-0.4", "0.4-0.7", "0.7-0.9", "0.9-1.0"))

  ggplot(a) +  geom_raster(aes(x=label1_n, y=label2_n, fill=cor_cat))+geom_text(aes(x=label1_n, y=label2_n, label=ifelse(l95!=u95, glue::glue("{round(cor,2)} ({round(l95, 2)} - {round(u95, 2)})"), 1)), size=6) + labs(fill="Correlation") +scale_x_discrete(limits=rev) +theme_minimal()+ theme(text = element_text(size=24), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("") + ylab("")  + scale_x_discrete(labels = function(x) stringr::str_wrap(x, width=10)) + scale_fill_brewer(palette="Blues", labels = function(breaks) {breaks[is.na(breaks)] <- "No overlap"; breaks})
  ggsave(glue::glue("paper/figures/fig_corr_{r_len}.png"), width=28, height=14, bg="white")
  ggsave(glue::glue("paper/figures/fig_corr_{r_len}.pdf"), width=28, height=14, bg="white")



  selected_index <- 1

  r_sim <- metas[[selected_index]]$r_samples

  hetrogeneity <- as.data.frame(metas[[selected_index]]$I2) %>% mutate(date=metas[[selected_index]]$r_est$date)
  colnames(hetrogeneity) <- c("lower", "median", "upper", "date")
  combined <- as.data.frame(t(r_sim)) %>% mutate(date=metas[[selected_index]]$r_est$date) %>% tidyr::pivot_longer(- date, values_to="r", names_to="sim") 
  combined <- get_gt(combined, r_sim)
  combined <- combined %>% mutate(R= (1 + r/rate)^(scale))
  combined <- combined %>% mutate(R=ifelse(is.na(R), 0, R))

  combined_sum <- combined %>% group_by(date) %>% summarise(
    gmean=mean(r),mean_l=quantile(r, 0.025), mean_u=quantile(r, 0.975),
    gR=mean(R),R_l=quantile(R, 0.025), R_u=quantile(R, 0.975)
  )
  combined <- combined %>% mutate(I=cumprod(exp(r)))

  plot_f <- rbind(combined_sum %>% select(date, mean=gmean, lower=mean_l, higher=mean_u) %>% mutate(type="Combined", name="Combined"), shifted %>% filter(label=="hosp")%>% mutate(type=t_hosp, name="test") %>% select(type, date, mean=r, lower=r_l, higher=r_u, name))
  plot_f <- rbind(plot_f, shifted %>% filter(label=="cases")%>% mutate(type="Cases", name="test") %>% select(type, date, mean=r, lower=r_l, higher=r_u, name))
  plot_f$type <- factor(plot_f$type, c("Combined", t_hosp ,t_cases, t_deaths, "Fraction"))
  q3 <- ggplot(combined_sum %>% filter(date > as.Date("2020-03-10") & date < as.Date("2024-01-01")), aes(x=date, y=gR, ymin=R_l, ymax=R_u)) + geom_line() + geom_ribbon(alpha=0.4) + labs(fill="", color="") + common_theme + scale_color_brewer(palette="Dark2") +  scale_fill_brewer(palette="Dark2") + xlab("") + ylab("Reproduction Number") + theme(legend.position="None") #+ scale_y_continuous(sec.axis = sec_axis(~ (1 + kappa*.*mean)^(1/kappa), name = "Reproduction number")) + xlab("")
  q1 <- ggplot(plot_f %>% filter(date > as.Date("2020-03-10") & date < as.Date("2024-01-01")), aes(x=date, y=mean, ymin=lower, ymax=higher, color=type, fill=type)) + geom_line() + geom_ribbon(alpha=0.8) + labs(fill="", color="") + common_theme + scale_color_brewer(palette="Dark2") +  scale_fill_brewer(palette="Dark2") + xlab("") + ylab("Growth Rate") + theme(legend.position="bottom") #+ scale_y_continuous(sec.axis = sec_axis(~ (1 + kappa*.*mean)^(1/kappa), name = "Reproduction number")) + xlab("")
  q2 <- ggplot(hetrogeneity, aes(x=date, y=median, ymin=lower, ymax=upper)) + geom_line() + geom_ribbon(alpha=0.4) + common_theme + ylab("Heterogeneity (I^2)") + xlab("")
  I <- combined %>% group_by(sim) %>% mutate(I=1*cumprod(exp(r))) %>% mutate(I=I/max(I))#data.frame(I=10*cumprod(exp(combined$mean)), t=combined_sum$date) %>% mutate(cum=cumsum(I))
  q4 <- ggplot(I %>% filter(date> as.Date("2020-03-10") & date < as.Date("2024-01-01"))) + geom_line(aes(x=date, y=I, group=sim)) + scale_y_continuous(trans="log10") + common_theme + ylab("Relative Incidence") + xlab("Date")
  top_row <- cowplot::plot_grid(NULL, q1 +  theme(text = element_text(size=25),legend.spacing.x=grid::unit(20, "pt"), legend.key.spacing=grid::unit(20, "pt")), rel_widths=c(0.015, 1))
  middle_row1 <- cowplot::plot_grid(NULL, q2 +  theme(text = element_text(size=25)), rel_widths=c(0.015, 1))
  middle_row2 <- cowplot::plot_grid(NULL, q3 +  theme(text = element_text(size=25)), rel_widths=c(0.015, 1))
  bottom_row <- cowplot::plot_grid(q4 +  theme(text = element_text(size=25)),NULL, rel_widths=c(1, 0.045))
  cowplot::plot_grid(top_row,middle_row1,middle_row2, bottom_row, nrow=4, labels=c("A)", "B)", "C)", "D)"))
  
  ggsave(glue::glue("paper/figures/combined_r_{r_len}.png"), height=20, width=12, bg="white")
  ggsave(glue::glue("paper/figures/combined_r_{r_len}.pdf"), height=20, width=12, bg="white")

}
parallel::mclapply(c(8,10,12,15), plot_rs, mc.cores=4)


Rs <- data.frame(dates=as.Date(c("2020-02-17", "2020-03-14", "2020-03-15", "2020-04-19", "2020-04-20", "2020-05-10", "2020-05-11", "2020-06-30", "2020-07-01", "2020-07-31", "2020-08-01", "2020-08-31")),
                  R=c(3.24, 3.24, 0.49, 0.49, 0.66, 0.66, 0.65, 0.65, 0.95, 0.95, 1.09, 1.07),
                  R_min=c(2.47, 2.47, 0.41, 0.41, 0.29, 0.29, 0.12, 0.12, 0.14, 0.14, 0.73, 0.73),
                  R_max=c(3.98, 3.98, 0.58, 0.58, 1.1, 1.1, 1.09, 1.09, 1.63, 1.63, 1.36, 1.36))


q3 <- ggplot(combined_sum %>% filter( date < as.Date("2020-09-01")))  + geom_line(aes(x=dates, y=R), data=Rs) + geom_ribbon(aes(x=dates, ymin=R_min, ymax=R_max), data=Rs, alpha=0.3)  +geom_line( aes(x=date, y=gR, ymin=R_l, ymax=R_u)) + geom_ribbon( aes(x=date, y=gR, ymin=R_l, ymax=R_u),alpha=0.4) + labs(fill="", color="") + common_theme + scale_color_brewer(palette="Dark2") +  scale_fill_brewer(palette="Dark2") + xlab("") + ylab("Reproduction Number") + theme(legend.position="None") #+ scale_y_continuous(sec.axis = sec_axis(~ (1 + kappa*.*mean)^(1/kappa), name = "Reproduction number")) + xlab("")
ggsave("test.png")


# SUPP

r_sups <- list()
for(r_l_sens in c(8, 10, 12, 15)){
       r_sups[[length(r_sups) + 1]] <- estimate_r_t_direct(data$hosp_main_cause, len=r_l_sens, weekend=TRUE, N=1000) %>% mutate(r_len=r_l_sens, label=t_hosp)
       r_sups[[length(r_sups) + 1]] <- estimate_r_t_direct(data$msis, len=r_l_sens, weekend=TRUE, N=1000) %>% mutate(r_len=r_l_sens, label=t_cases)
}

r_sups <- rbindlist(r_sups)

ggplot(r_sups %>% filter(r_len!=3 & label==t_cases), aes(x=date, y=r, ymin=r_l95, ymax=r_u95)) + geom_line(aes(color=factor(r_len*2 +1))) + geom_ribbon(aes(fill=factor(r_len*2+1)), alpha=0.7) + facet_wrap(.~label) +common_theme + scale_color_brewer(palette="Dark2") +  scale_fill_brewer(palette="Dark2") + xlab("Date") + ylab("Growth Rate") + theme(legend.position="bottom")  + labs(fill="Window Length", color="Window Length")   +  theme(text = element_text(size=25))
ggsave("paper/figures/S_sensitivity_r_len.png", width=14, height=8, bg="white")
ggsave("paper/figures/S_sensitivity_r_len.pdf", width=14, height=8, bg="white", device=cairo_ps)
