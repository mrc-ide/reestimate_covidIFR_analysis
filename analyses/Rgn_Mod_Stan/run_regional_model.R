### FIT REGIONAL AND AGE SEROPREVALENCE TO RE-ESTIMATE SPECIFICITY - USING RSTAN.
library(rstan)
library(shinystan)
library(tidyr)
library(plotrix)
library(reshape)

################################
# COMPILE RSTAN MODEL
################################
model_reg_age_full<-stan_model("analyses/Rgn_Mod_Stan/s_IFR_estimate_region_age_full_CC5.stan")
model_reg_age_full_sepSero<-stan_model("analyses/Rgn_Mod_Stan/s_IFR_estimate_region_age_full_sep_seroCC5.stan")

###################
# PROCESS DATA
###################

dat0<-readRDS("results/descriptive_results/descriptive_results_datamap.RDS")

dat_age <- dat0 %>%
  dplyr::filter(breakdown == "ageband")
dat_reg <- dat0 %>%
  dplyr::filter(breakdown == "region")



######
# SPAIN
######
####### Extract data & filter to most recent survey
curr_study_id<-"ESP1-2"
curr_dat_age <- dat_age$plotdat[[which(dat_age$study_id==curr_study_id)]] %>%
  dplyr::filter(seromidpt == obsday & obsdaymax==max(obsdaymax))
seromidpt<-curr_dat_age$seromidpt[1]
N_deaths_age<-round(curr_dat_age$cumdeaths)

####### Filter out duplicated seroprevalence data
inds<-!duplicated(curr_dat_age$n_positive,curr_dat_age$n_tested)
x_sero_age<-curr_dat_age$n_positive[inds]
N_sero_age<-curr_dat_age$n_tested[inds]

agebrks_d<-c(seq(0,90,10),999)

# seroassay validation data
x_sens_validat<-dat_age$data[[which(dat_age$study_id==curr_study_id)]]$sero_sens$npos
N_sens_validat<-dat_age$data[[which(dat_age$study_id==curr_study_id)]]$sero_sens$ntest
N_spec_validat<-dat_age$data[[which(dat_age$study_id==curr_study_id)]]$sero_spec$ntest
x_spec_validat<-dat_age$data[[which(dat_age$study_id==curr_study_id)]]$sero_spec$nneg


# pop by age
pop_age<-curr_dat_age$popn

### regional data - seroprevalence and deaths
i<-which(dat_reg$study_id==curr_study_id)
curr_dat_reg<-dat_reg$data[[i]]
sero_reg<- curr_dat_reg$seroprev_group %>%
  dplyr::filter(ObsDaymax==max(ObsDaymax)) %>%
  dplyr::mutate(n_positive=round(n_positive)) %>%
  dplyr::arrange(region)
x_sero_reg<-round(sero_reg$n_positive)
N_sero_reg<-sero_reg$n_tested
# death data
# total deaths
deaths_at_sero <- curr_dat_reg$deaths_TSMCMC %>%
  dplyr::mutate(cumdeaths=cumsum(deaths)) %>%
  dplyr::filter(ObsDay == seromidpt)
deaths_reg<-curr_dat_reg$deaths_propMCMC %>%
  dplyr::arrange(region)
N_deaths_reg<-round(deaths_at_sero$cumdeaths * deaths_reg$death_prop)


# pop by age and region
pop_reg_age<-dplyr::select(curr_dat_reg$prop_pop, -popN) %>%
  dplyr::arrange(region)
pop_reg <- curr_dat_reg$prop_pop %>%
  dplyr::arrange(region) %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(popN = sum(popN)) %>%
  dplyr::pull(popN)

pop_reg_age<-spread(pop_reg_age, key = ageband, value = pop_prop)
pop_reg_age<-as.matrix(pop_reg_age[2:ncol(pop_reg_age)])
colnames(pop_reg_age) <- NULL

## region plot to check all is well.
plot(x_sero_reg/N_sero_reg, 100000*N_deaths_reg/pop_reg)

########################
# FIT RSTAN MODEL
########################

nIter<-20000

options(mc.cores = 2) # parallel::detectCores())

t1<-Sys.time()  #### RUN TAKES APPROX 19 MINS ON 2 CORES, less on 4.
fit_reg_age_full <- sampling(model_reg_age_full,list(nr=length(pop_reg),
                                                     na=length(pop_age),
                                                     x_seror=x_sero_reg,
                                                     N_seror=N_sero_reg,
                                                     x_seroa = x_sero_age,
                                                     N_seroa = N_sero_age,
                                                     N_deathsr=N_deaths_reg,
                                                     N_deathsa=N_deaths_age,
                                                     tot_obsd = sum(N_deaths_age),
                                                     popr=pop_reg,
                                                     prop_pop_reg = pop_reg/sum(pop_reg),
                                                     popa=pop_age,
                                                     prop_pop_age = pop_age/sum(pop_age),
                                                     pop_reg_age = round(pop_reg_age*sum(pop_reg)),
                                                     prop_pop_reg_age = pop_reg_age,
                                                     x_sens_validat=x_sens_validat,
                                                     N_sens_validat=N_sens_validat,
                                                     x_spec_validat=x_spec_validat,
                                                     N_spec_validat=N_spec_validat),
                             iter=nIter, chains=4,control = list(adapt_delta = 0.98, max_treedepth = 15))
t2<-Sys.time()
t2-t1
#print(fit_reg_age_full)
## file too big for github so write elsewhere.
if(write2file) saveRDS(fit_reg_age_full, "C:/Users/Lucy/Dropbox (SPH Imperial College)/IFR update/rgn_mod_results/fit_spain_reg_age_full.rds")

params<-rstan::extract(fit_reg_age_full)
#plot(density(params$specificity))
if(write2file) write.csv(params$specificity, file="analyses/Rgn_Mod_Stan/results/spain_spec_reg_age.csv",row.names = F,col.names = NULL)
if(write2file) write.csv(params$sensitivity, file="analyses/Rgn_Mod_Stan/results/spain_sens_reg_age.csv",row.names = F,col.names = NULL)



#################
# SWEDEN
#################
curr_study_id<-"SWE1"
####### Extract data & filter to most recent survey
curr_dat_age <- dat_age$plotdat[[which(dat_age$study_id==curr_study_id)]] %>%
  dplyr::filter(seromidpt == obsday & obsdaymax==max(obsdaymax))
seromidpt<-curr_dat_age$seromidpt[1]
N_deaths_age<-round(curr_dat_age$cumdeaths)

####### Filter out duplicated seroprevalence data
inds<-!duplicated(curr_dat_age$seroprev)
seroprev<-curr_dat_age$seroprev[inds]
N_sero_age<-500  # inferred from CI
x_sero_age<-round(seroprev*N_sero_age)

agebrks_d<-c(0,seq(9,89,10),999)

# seroassay validation data
x_sens_validat<-dat_age$data[[which(dat_age$study_id==curr_study_id)]]$sero_sens$npos
N_sens_validat<-dat_age$data[[which(dat_age$study_id==curr_study_id)]]$sero_sens$ntest
N_spec_validat<-dat_age$data[[which(dat_age$study_id==curr_study_id)]]$sero_spec$ntest
x_spec_validat<-dat_age$data[[which(dat_age$study_id==curr_study_id)]]$sero_spec$nneg


# pop by age
pop_age<-curr_dat_age$popn

### regional data - seroprevalence and deaths
curr_dat_reg<-read.csv("data/raw/SWE1_regions.csv") %>%
  dplyr::arrange(region)
N_sero_reg<-curr_dat_reg$n_inferred
x_sero_reg<-round(N_sero_reg*curr_dat_reg$seroprev)
# death data
deaths_reg<-curr_dat_reg$deaths_propMCMC %>%
  dplyr::arrange(region)
N_deaths_reg<-round(deaths_at_sero$cumdeaths * $death_prop)

# pop by age and region
# demography (non-US Census data and BRA on its own)
pop_reg_age <- readr::read_tsv("data/raw/non_usa_non_bra_population.tsv") %>%
  dplyr::select(-c("reference")) %>%
  dplyr::filter(study_id=="SWE1")
%>%
  dplyr::select(curr_dat_reg$prop_pop, -popN)

pop_reg <- curr_dat_reg$prop_pop %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(popN = sum(popN)) %>%
  dplyr::pull(popN)

pop_reg_age<-spread(pop_reg_age, key = ageband, value = pop_prop)
pop_reg_age<-as.matrix(pop_reg_age[2:ncol(pop_reg_age)])
colnames(pop_reg_age) <- NULL

