########################
# Plot descriptive statistics
########################

#source("analyses/run_process_country_data.R")
rogan_gladen<-function(obs_prev,sens,spec) (obs_prev + spec -1)/(spec+sens-1)

studies<-c("ESP","DNK","NLD", "CHE","IRN")
i<-studies[1]
x<-readRDS(paste0("data/derived/",i,"/",i,"_agebands.RDS"))

#####################################################
## Align serology and deaths by age where possible. Document assumptions where no perfect alignment
#####################################################
### Spain. Perfect alignment of deaths and seroprevalence age groups.
i<-"ESP"
x<-readRDS(paste0("data/derived/",i,"/",i,"_agebands.RDS"))
curr_sero<-x$seroprev_group
curr_sero<-curr_sero %>%
  dplyr::mutate(ageband_deaths=cut(curr_sero$age_high, breaks=c(0,x$deaths_group$age_high))) %>%
  dplyr::group_by(ageband_deaths) %>%
  dplyr::summarise(ObsDaymin=mean(ObsDaymin),
                   ObsDaymax=mean(ObsDaymax),
                   n_tested=sum(n_tested),
                   n_positive=sum(n_positive)) %>%
  dplyr::mutate(seroprevalence=n_positive/n_tested) %>%
  dplyr::ungroup()

curr_sero$age_mid<-0.5*(x$deaths_group$age_low + x$deaths_group$age_high)
curr_sero$age_mid[nrow(curr_sero)]<-95
curr_sero$inf_pop_crude<-curr_sero$seroprevalence * x$popN*x$prop_pop$pop_prop
curr_sero$ifr_age_crude<-x$deaths_group$deaths_at_sero/curr_sero$inf_pop_crude
curr_sero$seroprev_adj_ss<-rogan_gladen(curr_sero$seroprevalence, x$sero_sens,x$sero_spec)
curr_sero$inf_pop_adj_ss<-curr_sero$seroprev_adj_ss * x$popN*x$prop_pop$pop_prop
curr_sero$ifr_age_adj_ss<-x$deaths_group$deaths_at_sero/curr_sero$inf_pop_adj_ss
esp_res<-curr_sero

### Netherlands. Not perfectly aligned.
# Assumptions.
# 1) 31-40 seroprevalence will be equivalent to the 30-39 age group, same for other 10 year age bands 30-70.
# 2) 18-30 seroprevalence = 20-29 seroprevalence
# 3) 60-72 seroprevalence = 60-69 seroprevalence
# 4) all other seroprevalence = national average.
i<-"NLD"
x<-readRDS(paste0("data/derived/",i,"/",i,"_agebands.RDS"))
curr_sero<-x$deaths_group
curr_sero$age_mid<-0.5*(curr_sero$age_low + curr_sero$age_high)
curr_sero$age_mid[which(curr_sero$ageband=='94-999')]<-98  ## TODO check what average is in this age group.
## enter national average by default.
curr_sero$seroprevalence<-x$seroprev$seroprev
#######
# enter for specific age groups
i<-which(x$seroprev_group$ageband=='18-30')
curr_sero$seroprevalence[which(curr_sero$ageband=='19-24' | curr_sero$ageband=='24-29')]<-
      x$seroprev_group$seroprevalence[i]
i<-which(x$seroprev_group$ageband=='30-40')
curr_sero$seroprevalence[which(curr_sero$ageband=='29-34' | curr_sero$ageband=='34-39')]<-
  x$seroprev_group$seroprevalence[i]
i<-which(x$seroprev_group$ageband=='40-50')
curr_sero$seroprevalence[which(curr_sero$ageband=='39-44' | curr_sero$ageband=='44-49')]<-
  x$seroprev_group$seroprevalence[i]
i<-which(x$seroprev_group$ageband=='50-60')
curr_sero$seroprevalence[which(curr_sero$ageband=='49-54' | curr_sero$ageband=='54-59')]<-
  x$seroprev_group$seroprevalence[i]
i<-which(x$seroprev_group$ageband=='60-72')
curr_sero$seroprevalence[which(curr_sero$ageband=='59-64' | curr_sero$ageband=='64-69')]<-
  x$seroprev_group$seroprevalence[i]

x$prop_pop<-x$prop_pop[which(!is.na(x$prop_pop$ageband)),]

curr_sero$inf_pop_crude<-curr_sero$seroprevalence * x$popN*x$prop_pop$pop_prop
curr_sero$ifr_age_crude<-x$deaths_group$deaths_at_sero/curr_sero$inf_pop_crude
curr_sero$seroprev_adj_ss<-rogan_gladen(curr_sero$seroprevalence, x$sero_sens,x$sero_spec)
curr_sero$inf_pop_adj_ss<-curr_sero$seroprev_adj_ss * x$popN*x$prop_pop$pop_prop
curr_sero$ifr_age_adj_ss<-x$deaths_group$deaths_at_sero/curr_sero$inf_pop_adj_ss
nld_res<-curr_sero


### Denmark
i<-"DNK"
x<-readRDS(paste0("data/derived/",i,"/",i,"_agebands.RDS"))
curr_sero<-x$deaths_group
curr_sero$age_mid<-0.5*(curr_sero$age_low + curr_sero$age_high)
curr_sero$age_mid[which(curr_sero$ageband=='89-999')]<-95
## enter national average by default.
curr_sero$seroprevalence<-x$seroprev$seroprev
curr_sero$seroprevalence[which(curr_sero$ageband=='59-69')]<-
  x$seroprev_group$seroprevalence[which(x$seroprev_group$ageband=='59-69')]

curr_sero$inf_pop_crude<-curr_sero$seroprevalence * x$popN*x$prop_pop$pop_prop
curr_sero$ifr_age_crude<-x$deaths_group$deaths_at_sero/curr_sero$inf_pop_crude
curr_sero$seroprev_adj_ss<-rogan_gladen(curr_sero$seroprevalence, x$sero_sens,x$sero_spec)
curr_sero$inf_pop_adj_ss<-curr_sero$seroprev_adj_ss * x$popN*x$prop_pop$pop_prop
curr_sero$ifr_age_adj_ss<-x$deaths_group$deaths_at_sero/curr_sero$inf_pop_adj_ss
dnk_res<-curr_sero


plot(esp_res$age_mid,100*esp_res$ifr_age_crude,xlab="age group (years)",ylab="Crude IFR",xlim=c(0,100),
     ylim=c(0,17))
points(esp_res$age_mid,100*esp_res$ifr_age_adj_ss,col="blue")
plot(nld_res$age_mid,100*nld_res$ifr_age_crude,xlab="age group years",ylab="Crude IFR",xlim=c(0,100),
     ylim=c(0,17))
points(nld_res$age_mid,100*nld_res$ifr_age_adj_ss,col="blue")
plot(dnk_res$age_mid,100*dnk_res$ifr_age_crude,xlab="age group years",ylab="Crude IFR",xlim=c(0,100),
     ylim=c(0,17))
points(dnk_res$age_mid,100*dnk_res$ifr_age_adj_ss,col="blue")
