####################################################################################
## Purpose: run spain on LL
##
## Notes:
####################################################################################
setwd("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis")
library(rslurm)
library(COVIDCurve)
library(tidyverse)
source("R/covidcurve_helper_functions.R")

#............................................................
# read in AGE DATA and experiment
#...........................................................
age <- readRDS("data/derived/ESP/ESP_agebands.RDS")


#......................
# make IFR models
#......................
wrap_make_IFR_model <- function(x) {
  # make dfs
  ifr_paramsdf <- make_ma_reparamdf(num_mas = 10)
  knot_paramsdf <- make_splinex_reparamdf(max_xvec = list("name" = "x4", min = 125, init = 131, max = 137, dsc1 = 125, dsc2 = 137),
                                          num_xs = 4)

  infxn_paramsdf <- make_spliney_reparamdf(max_yvec = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                           num_ys = 5)

  sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1"),
                                  min =   c(0.83,     0.50,   10,         117),
                                  init =  c(0.85,     0.99,   10,         125),
                                  max =   c(0.87,     1.00,   10,         131),
                                  dsc1 =  c(850,      10,     5,         117),
                                  dsc2 =  c(150,      3,      15,        131))
  noise_paramsdf <- make_noiseeff_reparamdf(num_Nes = 10, min = 0, init = 5, max = 10)

  # bring together
  df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, knot_paramsdf, sens_spec_tbl, noise_paramsdf)
  #......................
  # format data
  #......................
  dictkey <- tibble::tibble(ageband = age$seroprev_group$ageband, Strata = paste0("ma", 1:10))
  # deaths
  age$deaths <- age$deaths %>%
    dplyr::left_join(., dictkey) %>%
    dplyr::select(c("ObsDay", "Strata", "Deaths"))
  # seroprev
  age$obs_serology <- dplyr::left_join(age$seroprev_group, dictkey) %>%
    dplyr::mutate(SeroDay = "sero_day1") %>%
    dplyr::rename(SeroPrev = seroprevalence) %>%
    dplyr::select(c("SeroDay", "Strata", "SeroPrev"))

  inputdata <- list(obs_deaths = age$deaths,
                    obs_serology = age$obs_serology)

  demog <- age$prop_pop %>%
    dplyr::left_join(., dictkey) %>%
    dplyr::select(c("Strata", "popN")) %>%
    dplyr::mutate(popN = round(popN))


  # make mod
  mod1 <- make_IFRmodel_agg$new()
  mod1$set_MeanOnset(18.8)
  mod1$set_CoefVarOnset(0.45)
  mod1$set_level("Time-Series")
  mod1$set_IFRparams(paste0("ma", 1:10))
  mod1$set_maxMa("ma10")
  mod1$set_Knotparams(paste0("x", 1:4))
  mod1$set_relKnot("x4")
  mod1$set_Infxnparams(paste0("y", 1:5))
  mod1$set_relInfxn("y5")
  mod1$set_Serotestparams(c("sens", "spec", "sero_rate"))
  mod1$set_Serodayparams(c("sero_day1"))
  mod1$set_Noiseparams(paste0("Ne", 1:10))
  mod1$set_data(inputdata)
  mod1$set_demog(demog)
  mod1$set_paramdf(df_params)
  mod1$set_rho(rep(1, 10))
  mod1$set_rcensor_day(.Machine$integer.max)
  # out
  mod1
}

#......................
# make maps
#......................
map <- tibble::as_tibble(expand.grid(rungs = c(10, 25, 50),
                                     GTI_pow = c(2, 2.5, 3, 3.5, 4.0, 4.5, 5.0, 5.5, 6),
                                     burnin = 1e3,
                                     samples = 1e3))


map$modelobj <- purrr::map(map$rungs, wrap_make_IFR_model)

#......................
# wrapper for run
#......................
run_wrapper <- function(modelobj, rungs, GTI_pow, burnin, samples) {
  fit <- COVIDCurve::run_IFRmodel_agg(IFRmodel = modelobj,
                                      reparamIFR = TRUE,
                                      reparamInfxn = TRUE,
                                      reparamKnots = TRUE,
                                      chains = 10,
                                      burnin = burnin,
                                      samples = samples,
                                      rungs = rungs,
                                      GTI_pow = GTI_pow,
                                      silent = FALSE)
  mc_accept_mean <- mean(fit$mcmcout$diagnostics$mc_accept$value)
  mc_accept_min <- min(fit$mcmcout$diagnostics$mc_accept$value)
  out <- list(fit = fit,
              mc_accept_mean = mc_accept_mean,
              mc_accept_min = mc_accept_min)
  return(out)
}

#......................
# send out on slurm
#......................
ntry <- 30
sjob <- rslurm::slurm_apply(f = run_wrapper,
                            params = map,
                            jobname = 'newRandomNetuneMCMC_ESPage',
                            nodes = ntry,
                            cpus_per_node = 1,
                            submit = T,
                            slurm_options = list(mem = "24g",
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "2-00:00:00"))



#............................................................
# read in REGION DATA and experiment
#...........................................................
rgn <- readRDS("data/derived/ESP/ESP_regions.RDS")


#......................
# make IFR models
#......................
wrap_make_IFR_model <- function(x) {
  # make dfs
  ifr_paramsdf <- make_ma_reparamdf(num_mas = 17)
  knot_paramsdf <- make_splinex_reparamdf(max_xvec = list("name" = "x4", min = 125, init = 131, max = 137, dsc1 = 125, dsc2 = 137),
                                          num_xs = 4)

  infxn_paramsdf <- make_spliney_reparamdf(max_yvec = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                           num_ys = 5)

  sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1"),
                                  min =   c(0.83,     0.50,   10,         117),
                                  init =  c(0.85,     0.99,   10,         125),
                                  max =   c(0.87,     1.00,   10,         131),
                                  dsc1 =  c(850,      10,     5,         117),
                                  dsc2 =  c(150,      3,      15,        131))
  noise_paramsdf <- make_noiseeff_reparamdf(num_Nes = 17, min = 0, init = 5, max = 10)

  # bring together
  df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, knot_paramsdf, sens_spec_tbl, noise_paramsdf)
  #......................
  # format data
  #......................
  dictkey <- tibble::tibble(region = rgn$seroprev_group$region, Strata = paste0("ma", 1:17))
  # deaths
  rgn$deaths <- rgn$deaths %>%
    dplyr::left_join(., dictkey) %>%
    dplyr::select(c("ObsDay", "Strata", "Deaths"))
  # seroprev
  rgn$obs_serology <- dplyr::left_join(rgn$seroprev_group, dictkey) %>%
    dplyr::mutate(SeroDay = "sero_day1") %>%
    dplyr::rename(SeroPrev = seroprevalence) %>%
    dplyr::select(c("SeroDay", "Strata", "SeroPrev"))

  inputdata <- list(obs_deaths = rgn$deaths,
                    obs_serology = rgn$obs_serology)

  demog <- rgn$prop_pop %>%
    dplyr::left_join(., dictkey) %>%
    dplyr::select(c("Strata", "popN")) %>%
    dplyr::mutate(popN = round(popN))



  # make mod
  mod1 <- make_IFRmodel_agg$new()
  mod1$set_MeanOnset(18.8)
  mod1$set_CoefVarOnset(0.45)
  mod1$set_level("Time-Series")
  mod1$set_IFRparams(paste0("ma", 1:17))
  mod1$set_maxMa("ma14") # madrid
  mod1$set_Knotparams(paste0("x", 1:4))
  mod1$set_relKnot("x4")
  mod1$set_Infxnparams(paste0("y", 1:5))
  mod1$set_relInfxn("y5")
  mod1$set_Serotestparams(c("sens", "spec", "sero_rate"))
  mod1$set_Serodayparams(c("sero_day1"))
  mod1$set_Noiseparams(paste0("Ne", 1:17))
  mod1$set_data(inputdata)
  mod1$set_demog(demog)
  mod1$set_paramdf(df_params)
  mod1$set_rho(rep(1, 17))
  mod1$set_rcensor_day(.Machine$integer.max)
  # out
  mod1
}

#......................
# make maps
#......................
map <- tibble::as_tibble(expand.grid(rungs = c(10, 25, 50),
                                     GTI_pow = c(2, 2.5, 3, 3.5, 4.0, 4.5, 5.0, 5.5, 6),
                                     burnin = 1e4,
                                     samples = 1e4))


map$modelobj <- purrr::map(map$rungs, wrap_make_IFR_model)

#......................
# wrapper for run
#......................
run_wrapper <- function(modelobj, rungs, GTI_pow, burnin, samples) {
  fit <- COVIDCurve::run_IFRmodel_agg(IFRmodel = modelobj,
                                      reparamIFR = TRUE,
                                      reparamInfxn = TRUE,
                                      reparamKnots = TRUE,
                                      chains = 10,
                                      burnin = burnin,
                                      samples = samples,
                                      rungs = rungs,
                                      GTI_pow = GTI_pow,
                                      silent = FALSE)
  mc_accept_mean <- mean(fit$mcmcout$diagnostics$mc_accept$value)
  mc_accept_min <- min(fit$mcmcout$diagnostics$mc_accept$value)
  out <- list(fit = fit,
              mc_accept_mean = mc_accept_mean,
              mc_accept_min = mc_accept_min)
  return(out)
}

#......................
# send out on slurm
#......................
ntry <- 30
sjob <- rslurm::slurm_apply(f = run_wrapper,
                            params = map,
                            jobname = 'newRandomNetuneMCMC_ESPRGN',
                            nodes = ntry,
                            cpus_per_node = 1,
                            submit = T,
                            slurm_options = list(mem = "24g",
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "2-00:00:00"))
