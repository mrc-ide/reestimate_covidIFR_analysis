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
model_reg_only_full<-stan_model("analyses/Rgn_Mod_Stan/s_IFR_estimate_region_age_full_no_age_seroCC5.stan")


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
# FIT RSTAN MODEL - SPAIN
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

agebrks_d<-c(0,seq(9,79,10),999)

# seroassay validation data
x_sens_validat<-dat_age$data[[which(dat_age$study_id==curr_study_id)]]$sero_sens$npos
N_sens_validat<-dat_age$data[[which(dat_age$study_id==curr_study_id)]]$sero_sens$ntest
N_spec_validat<-dat_age$data[[which(dat_age$study_id==curr_study_id)]]$sero_spec$ntest
x_spec_validat<-dat_age$data[[which(dat_age$study_id==curr_study_id)]]$sero_spec$nneg


### regional data - seroprevalence and deaths
curr_dat_reg<-dat_reg$data[[which(dat_reg$study_id==curr_study_id)]]
sero_reg<- curr_dat_reg$seroprev_group %>%
  dplyr::filter(ObsDaymax==max(ObsDaymax)) %>%
  dplyr::mutate(n_positive=round(n_positive)) %>%
  dplyr::arrange(region)
N_sero_reg<-sero_reg$n_tested
x_sero_reg<-sero_reg$n_positive
# death data
# total deaths
deaths_at_sero <- curr_dat_reg$deaths_TSMCMC %>%
  dplyr::mutate(cumdeaths=cumsum(deaths)) %>%
  dplyr::filter(ObsDay == seromidpt)
### regional deaths data.
# read in full set of regions in sweden to work out proportion of deaths in those 3 regions where we have data:
# cumulative deaths
SWEdeathsdf <- readr::read_tsv("data/raw/cumulative_deaths.tsv") %>%
  dplyr::filter(study_id=="SWE1" & region!="all") %>%
  dplyr::mutate(date_start_survey = lubridate::ymd(date_start_survey), # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format
                date_end_survey = lubridate::ymd(date_end_survey),
                prop_deaths=n_deaths/sum(n_deaths)) %>%
  dplyr::filter(region %in% curr_dat_reg$deaths_propMCMC$region) %>%
  dplyr::arrange(region)

N_deaths_reg<-round(deaths_at_sero$cumdeaths * SWEdeathsdf$prop_deaths)

## only include total deaths from the regions included
deaths_age<-curr_dat_age %>%
  dplyr::ungroup() %>%
  dplyr::mutate(age_mid=ifelse(ageband=='79-89' | ageband=='89-999',90,age_mid)) %>%
  dplyr::select(age_mid,cumdeaths) %>%
  dplyr::group_by(age_mid) %>%
  dplyr::summarise(cumdeaths=sum(cumdeaths))
N_deaths_age<-round(sum(N_deaths_reg) * deaths_age$cumdeaths/sum(deaths_age$cumdeaths))

#### check total deaths are the same by age and region
sum(N_deaths_reg)
sum(N_deaths_age)

# pop by age and region
pop_reg_age<-dplyr::select(curr_dat_reg$prop_pop, -pop_prop) %>%
  dplyr::arrange(region)
### careful with ordering:
pop_reg_age <- pop_reg_age %>%
  dplyr::mutate(age_low = as.numeric(stringr::str_extract(ageband, "[0-9]+?(?=-)")),
                age_low2=ifelse(ageband=='79-89' | ageband=='89-999',90,age_low)) %>%
  dplyr::group_by(region,age_low2) %>%
  dplyr::summarise(popN=sum(popN)) %>%
  dplyr::arrange(region, age_low2)

pop_reg_age<-spread(pop_reg_age, key = age_low2, value = popN)
pop_reg_age<-as.matrix(pop_reg_age[2:ncol(pop_reg_age)])
colnames(pop_reg_age) <- NULL

pop_reg <- rowSums(pop_reg_age)
pop_age <- colSums(pop_reg_age)
## check total population the same by region and by age
if(sum(pop_reg)!=sum(pop_age)) print("Error, different population totals by age versus by region")



## region plot to check all is well.
plot(x_sero_reg/N_sero_reg, 100000*N_deaths_reg/pop_reg)
## age plot to check all is well.
plot(1:length(pop_age),N_deaths_age/pop_age)

########################
# FIT RSTAN MODEL - SWEDEN
########################
nIter<-20000

options(mc.cores = 4) # parallel::detectCores())

t1<-Sys.time()  #### RUN TAKES APPROX 19 MINS ON 2 CORES, less on 4.
fit_reg_age_full <- sampling(model_reg_only_full,list(nr=length(pop_reg),
                                                     na=length(pop_age),
                                                     x_seror=x_sero_reg,
                                                     N_seror=N_sero_reg,
                                                     N_deathsr=N_deaths_reg,
                                                     N_deathsa=N_deaths_age,
                                                     tot_obsd = sum(N_deaths_reg),
                                                     popr=pop_reg,
                                                     prop_pop_reg = pop_reg/sum(pop_reg),
                                                     popa=pop_age,
                                                     prop_pop_age = pop_age/sum(pop_age),
                                                     pop_reg_age = pop_reg_age,
                                                     prop_pop_reg_age = pop_reg_age/sum(pop_reg_age),
                                                     x_sens_validat=x_sens_validat,
                                                     N_sens_validat=N_sens_validat,
                                                     x_spec_validat=x_spec_validat,
                                                     N_spec_validat=N_spec_validat),
                             iter=nIter, chains=4,control = list(adapt_delta = 0.99, max_treedepth = 15))
t2<-Sys.time()
t2-t1
#print(fit_reg_age_full)
## file too big for github so write elsewhere.
if(write2file) saveRDS(fit_reg_age_full, "C:/Users/Lucy/Dropbox (SPH Imperial College)/IFR update/rgn_mod_results/fit_sweden_reg_age_full.rds")

params<-rstan::extract(fit_reg_age_full)
#plot(density(params$specificity))
if(write2file) write.csv(params$specificity, file="analyses/Rgn_Mod_Stan/results/sweden_spec_reg_age.csv",row.names = F,col.names = NULL)
if(write2file) write.csv(params$sensitivity, file="analyses/Rgn_Mod_Stan/results/sweden_sens_reg_age.csv",row.names = F,col.names = NULL)


#################
# New York State
#################
###################
# PROCESS DATA
###################
####### Extract data
# seroprevalence data
curr_study_id<-"NYS1"
####### Extract data & filter to most recent survey
curr_dat_age <- dat_age$plotdat[[which(dat_age$study_id==curr_study_id)]] %>%
  dplyr::filter(seromidpt == obsday & obsdaymax==max(obsdaymax))
seromidpt<-curr_dat_age$seromidpt[1]

agebrks_d<-c(0,seq(9,79,10),999)

# seroassay validation data
x_sens_validat<-dat_age$data[[which(dat_age$study_id==curr_study_id)]]$sero_sens$npos
N_sens_validat<-dat_age$data[[which(dat_age$study_id==curr_study_id)]]$sero_sens$ntest
N_spec_validat<-dat_age$data[[which(dat_age$study_id==curr_study_id)]]$sero_spec$ntest
x_spec_validat<-dat_age$data[[which(dat_age$study_id==curr_study_id)]]$sero_spec$nneg


### regional data - seroprevalence and deaths
curr_dat_reg<-dat_reg$data[[which(dat_reg$study_id==curr_study_id)]]
sero_reg<- curr_dat_reg$seroprev_group %>%
  dplyr::filter(ObsDaymax==max(ObsDaymax)) %>%
  dplyr::mutate(n_positive=round(n_positive)) %>%
  dplyr::arrange(region)
N_sero_reg<-sero_reg$n_tested
x_sero_reg<-sero_reg$n_positive
# death data
# total deaths
deaths_at_sero <- curr_dat_reg$deaths_TSMCMC %>%
  dplyr::mutate(cumdeaths=cumsum(deaths)) %>%
  dplyr::filter(ObsDay == seromidpt)

# regional level
curr_sero<-curr_dat_reg$data[[1]]$seroprevMCMC
x_sero_reg<- $n_positive
N_sero_reg<-dat_reg$n_tested

agebrks_d<-c(0,seq(9,79,10),999)  ## join 90+ with 80-90 as population data does not have 90+ category, only 85+
sero<-dplyr::filter(dat_age,!duplicated(n_positive,n_tested))
x_sero_age<-sero$n_positive
N_sero_age<-sero$n_tested

# seroassay validation data
x_sens_validat<-204
N_sens_validat<-234
N_spec_validat<-  92
x_spec_validat<-92

sens<-x_sens_validat/N_sens_validat
spec<-x_spec_validat/N_spec_validat


# national level
N_sero_nat<-sum(N_sero_age)
obs_prev<-(0.125*sens) + (1-0.125)*(1-spec)  ## calculate backwards using sens and spec.
x_sero_nat<-round(N_sero_nat*obs_prev) #cannot derive from regional/age data as is weighted.


# death data
N_deaths_reg<-round(dat_reg$cumdeaths)

dat_age$ageband<-c(1:9,9)
aged<-dat_age %>%
  dplyr::group_by(ageband) %>%
  dplyr::summarise(cumdeaths=sum(cumdeaths),
                   popn=sum(popn))
N_deaths_age<-round(aged$cumdeaths)

# pop by age and region
pop_reg0<-dat_reg$popn
pop_age0<-aged$popn

## groupings of seroprevalence in the study:
wr<-c("New York_Westchester|New York_Rockland")
long_island<-c("New York_Nassau|New York_Suffolk")
nyc<-"New York_New-York|New York_Kings|New York_Bronx|New York_Queens|New York_Richmond"
upstate<-c("Albany|Allegany|Broome|Cattaraugus|Cayuga|Chautauqua|Chemung|Chenango|Clinton|Columbia|Cortland|Delaware|Dutchess|Erie|Essex|Franklin|Fulton|Genessee|Greene|Hamilton|Herkimer|Jefferson|Lewis|Livingston|Madison|Monroe|Montgomery|Niagara|Oneida|Onondaga|Ontario|Orange|Orleans|Oswego|Otsego|Putnam|Rensselaer|St. Lawrence|Saratoga|Schenectady|Schoharie|Schuyler|Seneca|Steuben|Sullivan|Tioga|Tompkins|Ulster|Warren|Washington|Wayne|Wyoming|Yates")

#dat_pop<-read.csv("C:/Users/Lucy/Documents/GitHub/reestimate_covidIFR_analysis/data/raw/USA_County_Demographic_Data.csv")

populationdf <- readr::read_csv("C:/Users/Lucy/Documents/GitHub/reestimate_covidIFR_analysis/data/raw/USA_County_Demographic_Data.csv") %>%
  tidyr::gather(., key = "strata", value = "population", 3:ncol(.)) %>%
  dplyr::filter(stringr::str_detect(strata, "Both_", negate = TRUE)) %>%
  dplyr::filter(stringr::str_detect(strata, "_Total", negate = TRUE)) %>%
  dplyr::mutate(
    country = "USA",
    Countysp = gsub(" County", "", County),
    Countysp = gsub(" ", "-", Countysp),
    region = paste0(State, "_", Countysp),
    ageband = stringr::str_split_fixed(strata, "[A-Za-z]_", n = 2)[,2],
    ageband = ifelse(stringr::str_detect(ageband, "\\+"),
                     paste0(stringr::str_extract_all(ageband, "[0-9]+", simplify = TRUE), "-", 999),
                     ageband),
    age_low = as.numeric( stringr::str_split_fixed(ageband, "-[0-9]+", n = 2)[,1] ),
    age_high = as.numeric( stringr::str_split_fixed(ageband, "[0-9]-", n = 2)[,2] ),
    gender = stringr::str_extract_all(strata, "[A-Za-z]+", simplify = TRUE)[,1],
    age_breakdown = 1,
    for_regional_analysis = 1,
    gender_breakdown = 1
  ) %>%
  dplyr::select(c("country", "age_low", "age_high", "region", "gender", "population", "age_breakdown", "for_regional_analysis", "gender_breakdown")) %>%
  dplyr::left_join(., readr::read_csv("C:/Users/Lucy/Documents/GitHub/reestimate_covidIFR_analysis/data/raw/usa_study_id_county_key.csv"), by = "region")


NYSpopdf <- populationdf %>%
  dplyr::filter(grepl("New York", region))
wrp <- NYSpopdf %>%
  dplyr::filter(grepl(wr,region)) %>%
  dplyr::mutate(region="Westchester and Rockland")
long_islandp<-NYSpopdf %>%
  dplyr::filter(grepl(long_island,region)) %>%
  dplyr::mutate(region="Long Island")
nycp <- NYSpopdf %>%
  dplyr::filter(grepl(nyc,region)) %>%
  dplyr::mutate(region="New York City")
upstatep <- NYSpopdf %>%
  dplyr::filter(grepl(upstate,region)) %>%
  dplyr::mutate(region="Upstate New York")

NYpopdf <- rbind(wrp,long_islandp,nycp,upstatep)
NYpopdf<- NYpopdf %>%
  dplyr::mutate(age2=cut(age_high,agebrks_d)) %>%
  dplyr::group_by(region, age2) %>%
  dplyr::summarise(population = sum(population)) %>%
  dplyr::ungroup()

pop_reg_age<-melt(NYpopdf)
levels(pop_reg_age$age2)<-1:10
summ_pop<-spread(pop_reg_age, key = age2, value = value)

## checked regions in same order in the rgn RDS

pop_reg_age<-as.matrix(summ_pop[3:ncol(summ_pop)])
colnames(pop_reg_age) <- NULL
## convert to proportion of population in each age group per region
prop_pop_reg_age_regSum1<-pop_reg_age / rowSums(pop_reg_age)
prop_pop_reg_age_ageSum1<-pop_reg_age
for(i in 1:ncol(pop_reg_age)) prop_pop_reg_age_ageSum1[,i]<-pop_reg_age[,i] / sum(pop_reg_age[,i])

pop_age<-colSums(pop_reg_age)
pop_reg<-rowSums(pop_reg_age)

d_i<-c(rep(1,6),rep(2,3))
na_s<-length(x_sero_age)
pop_age_sero<-rep(0,na_s)
for(i in 1:length(N_deaths_age)) pop_age_sero[d_i[i]]<-pop_age_sero[d_i[i]]+pop_age[i]
nr<-length(N_deaths_reg)
pop_reg_age_sero<-matrix(rep(0,na_s*nr),nrow=nr,ncol=na_s)
for(i in 1:length(N_deaths_age)) {
  for(r in 1:nr) {
    pop_reg_age_sero[r,d_i[i]]<-pop_reg_age_sero[r,d_i[i]] + pop_reg_age[r,i]
  }
}



## quick region plot to check all is well.
plot(x_sero_reg/N_sero_reg, 100000*N_deaths_reg/pop_reg0)


########################
# FIT RSTAN MODELS
########################

nIter<-20000


### REGION & AGE, ASSUME RR of sero constant by age and region (i.e. no interaction)
fit_reg_age_full <- sampling(model_reg_age_full,list(nr=length(pop_reg),
                                                     na=length(pop_age),
                                                     na_s=length(x_sero_age),
                                                     d_i=d_i,
                                                     x_seror=x_sero_reg,
                                                     N_seror=N_sero_reg,
                                                     x_seroa = x_sero_age,
                                                     N_seroa = N_sero_age,
                                                     N_deathsr=N_deaths_reg,
                                                     N_deathsa=N_deaths_age,
                                                     tot_obsd = sum(N_deaths_reg),
                                                     popr=pop_reg,
                                                     prop_pop_reg = pop_reg/sum(pop_reg),
                                                     popa=pop_age,
                                                     popas= pop_age_sero, #// population in age groups for the serosurvey.
                                                     prop_pop_age = pop_age/sum(pop_age),
                                                     pop_reg_age = pop_reg_age,
                                                     prop_pop_reg_age = pop_reg_age/sum(pop_reg_age),
                                                     prop_pop_reg_age_sero=pop_reg_age_sero/sum(pop_reg_age_sero),
                                                     x_sens_validat=x_sens_validat,
                                                     N_sens_validat=N_sens_validat,
                                                     x_spec_validat=x_spec_validat,
                                                     N_spec_validat=N_spec_validat),
                             iter=nIter, chains=4,control = list(adapt_delta = 0.98, max_treedepth = 13))  # each chain has 200 iterations. Default is to use half of iterations for burn in.

print(fit_reg_age_full)


#################
# ITALY
#################
curr_study_id<-"ITA1"
