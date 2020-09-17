library(drjacoby)
library(tidyverse)
setwd("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/analyses/Rgn_Age_Mod_Compromise/temp_ESP")
#setwd("~/Documents/GitHub/reestimate_covidIFR_analysis/analyses/Rgn_Age_Mod_Compromise/temp_ESP")

#............................................................
# likelihoods
#...........................................................
cpp_loglike <- readLines("ESP_natcubicspline_loglike.cpp")
# remove comments which cause issues when coercing to string in this format
commlines <- grep("//", cpp_loglike)
cpp_loglike <- cpp_loglike[! 1:length(cpp_loglike) %in% commlines]
cpp_loglike <- gsub("\\\n", "", cpp_loglike)
cpp_loglike <- capture.output(cat(cpp_loglike))


cpp_logprior <- readLines("ESP_logprior.cpp")
# remove comments which cause issues when coercing to string in this format
commlines <- grep("//", cpp_logprior)
cpp_logprior <- cpp_logprior[! 1:length(cpp_logprior) %in% commlines]
cpp_logprior <- gsub("\\\n", "", cpp_logprior)
cpp_logprior <- capture.output(cat(cpp_logprior))



#............................................................
# data read in
#...........................................................
ESPrgn <- readRDS("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/data/derived/ESP1-2/ESP1-2_agebands.RDS")
#ESPage <- readRDS("~/Documents/GitHub/reestimate_covidIFR_analysis/data/derived/ESP1-2/ESP1-2_agebands.RDS")
ESPrgn <- readRDS("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/data/derived/ESP1-2/ESP1-2_regions.RDS")
#ESPrgn <- readRDS("~/Documents/GitHub/reestimate_covidIFR_analysis/data/derived/ESP1-2/ESP1-2_regions.RDS")

# store region names
rgnnames <- unique(ESPrgn$prop_pop$region)

 #......................
# wrangle data
#......................
# major liftover
obs_deaths <- ESPrgn$deaths_TSMCMC$deaths
age_prop_strata_obs_deaths <- ESPage$deaths_propMCMC$death_prop
rgn_prop_strata_obs_deaths <- ESPrgn$deaths_propMCMC$death_prop

datinput <- list(obs_deaths = obs_deaths,
                 age_prop_strata_obs_deaths = age_prop_strata_obs_deaths,
                 rgn_prop_strata_obs_deaths = rgn_prop_strata_obs_deaths,
                 age_obs_serologypos = round(ESPage$seroprevMCMC$n_positive),
                 age_obs_serologyn = round(ESPage$seroprevMCMC$n_tested),
                 rgn_obs_serologypos = round(ESPrgn$seroprevMCMC$n_positive),
                 rgn_obs_serologyn = round(ESPrgn$seroprevMCMC$n_tested))

#..................
# inputs
#..................
# get demog_pa
demog_pa <- ESPrgn$prop_pop %>%
  dplyr::group_by(region) %>%
  dplyr::mutate(totpop = sum(popN),
                demog_pa = popN/totpop) %>%
  dplyr::pull(demog_pa)
# demog age
demog_age <- ESPage$prop_pop %>%
  dplyr::group_by(ageband) %>%
  dplyr::summarise(popN = sum(popN)) %>%
  dplyr::pull(popN)

# demog rgn
demog_rgn <- ESPrgn$prop_pop %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(popN = sum(popN)) %>%
  dplyr::pull(popN)



misc_list = list(rcensor_day = .Machine$integer.max,
                 days_obsd = max(ESPrgn$deaths_TSMCMC$ObsDay),
                 n_sero_obs = 2,
                 max_seroday_obsd = 153,
                 sero_survey_start = c(118, 139),
                 sero_survey_end = c(132, 153),
                 n_knots = 5,
                 demog_pa = demog_pa,
                 Ademog = demog_age,
                 Rdemog = demog_rgn,
                 agestratlen = 10,
                 rgnstratlen = 17)

# param df
ifr_paramsdf <- tibble::tibble(name = paste0("ma", 1:10),
                               min  = rep(0, 10),
                               init = rep(0.1, 10),
                               max = rep(0.4, 10))

Ane_paramsdf <- tibble::tibble(name = paste0("Ane", 1:10),
                               min  = rep(0.5, 10),
                               init = rep(1, 10),
                               max = rep(1.5, 10))

Rne_paramsdf <- tibble::tibble(name = paste0("Rne", 1:17),
                               min  = rep(0.1, 17),
                               init = rep(1, 17),
                               max = rep(10, 17))

knot_paramsdf <- tibble::tibble(name = paste0("x", 1:4),
                                min  = c(rep(0, 3), 216),
                                init = c(rep(0.5, 3), 220),
                                max =  c(rep(1, 3), 230))

infxn_paramsdf <- tibble::tibble(name = paste0("y", 1:5),
                                 min  = c(0, 0, 0, 0, 0),
                                 init = c(0.05, 0.05, 9, 0.05, 0.05),
                                 max =  c(1, 1, 17.66, 1, 1))

tod_paramsdf <- tibble::tibble(name = c("mod", "sod", "sero_con_rate"),
                               min  = c(17,      0,     16),
                               init = c(19.66,   0.99,   18.3),
                               max =  c(22,      1,     21))

sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.83,    0.8),
                                init =  c(0.85,    0.99),
                                max =   c(1.00,    1.00))

df_params <- rbind.data.frame(ifr_paramsdf, Ane_paramsdf, Rne_paramsdf,
                              infxn_paramsdf, knot_paramsdf,
                              tod_paramsdf, sens_spec_tbl)
#......................
# mod
#......................
n_chains <- 5
n_cores <- parallel::detectCores()

if (n_cores < n_chains) {
  mkcores <- n_cores - 1
} else {
  mkcores <- n_chains
}

cl <- parallel::makeCluster(mkcores)
mcmcout <- drjacoby::run_mcmc(data = datinput,
                              df_params = df_params,
                              misc = misc_list,
                              loglike = cpp_loglike,
                              logprior = cpp_logprior,
                              burnin = 1e4,
                              samples = 1e4,
                              chains = 5,
                              rungs = 25,
                              GTI_pow = 3,
                              cluster = cl)
parallel::stopCluster(cl)
gc()

#............................................................
# out
#...........................................................
# save out pieces that we will need for functions later

out <- list(
  rgnnames = rgnnames,
  mcmcout = mcmcout,
  datinput = datinput,
  misc_list = misc_list
)

saveRDS(out, "mcmcout_newrgnl_mod.RDS")


