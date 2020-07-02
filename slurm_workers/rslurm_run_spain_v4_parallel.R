####################################################################################
## Purpose: This is a script that will wrap a DRAKE make process to
##          run multiple ESP data iterations to fine to the MCoupling
##          params
##
## Notes:
####################################################################################
setwd("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis")
library(drake)
library(parallel)
library(COVIDCurve)
library(tidyverse)
source("R/covidcurve_helper_functions.R")

#............................................................
# Make Paramset and write to disk for input into MCMC
#...........................................................
make_IFR_model_spain <- function(num_mas, maxMa, groupvar, dat) {
  # make dfs
  ifr_paramsdf <- make_ma_reparamdf(num_mas = num_mas)
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
  noise_paramsdf <- make_noiseeff_reparamdf(num_Nes = num_mas, min = 0, init = 5, max = 10)

  # bring together
  df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, knot_paramsdf, sens_spec_tbl, noise_paramsdf)
  #......................
  # format data
  #......................
  dictkey <- tibble::tibble(groupvar = unlist(dat$seroprev_group[, groupvar]), "Strata" = paste0("ma", 1:num_mas))
  colnames(dictkey) <- c(paste(groupvar), "Strata")
  # deaths
  dat$deaths <- dplyr::left_join(dat$deaths, dictkey) %>%
    dplyr::select(c("ObsDay", "Strata", "Deaths"))
  # seroprev
  dat$obs_serology <- dplyr::left_join(dat$seroprev_group, dictkey) %>%
    dplyr::mutate(SeroDay = "sero_day1") %>%
    dplyr::rename(SeroPrev = seroprevalence) %>%
    dplyr::select(c("SeroDay", "Strata", "SeroPrev"))

  inputdata <- list(obs_deaths = dat$deaths,
                    obs_serology = dat$obs_serology)

  demog <- dat$prop_pop %>%
    dplyr::left_join(., dictkey) %>%
    dplyr::select(c("Strata", "popN")) %>%
    dplyr::mutate(popN = round(popN))


  # make mod
  mod1 <- make_IFRmodel_agg$new()
  mod1$set_MeanOnset(18.8)
  mod1$set_CoefVarOnset(0.45)
  mod1$set_level("Time-Series")
  mod1$set_IFRparams(paste0("ma", 1:num_mas))
  mod1$set_maxMa(maxMa)
  mod1$set_Knotparams(paste0("x", 1:4))
  mod1$set_relKnot("x4")
  mod1$set_Infxnparams(paste0("y", 1:5))
  mod1$set_relInfxn("y5")
  mod1$set_Serotestparams(c("sens", "spec", "sero_rate"))
  mod1$set_Serodayparams(c("sero_day1"))
  mod1$set_Noiseparams(paste0("Ne", 1:num_mas))
  mod1$set_data(inputdata)
  mod1$set_demog(demog)
  mod1$set_paramdf(df_params)
  mod1$set_rho(rep(1, num_mas))
  mod1$set_rcensor_day(.Machine$integer.max)
  # out
  mod1
}

# read raw data and make param maps
rawage <- readRDS("data/derived/ESP/ESP_agebands.RDS")
rawrgn <- readRDS("data/derived/ESP/ESP_regions.RDS")
agemap <- tibble::as_tibble(expand.grid(rungs = c(10, 25, 50),
                                        GTI_pow = c(2, 2.5, 3, 3.5, 4.0, 4.5, 5.0, 5.5, 6),
                                        burnin = 1e2,
                                        samples = 1e2)) %>%
  dplyr::mutate(lvl = "age",
                num_mas = 10,
                maxMa = "ma10",
                groupvar = "ageband",
                dat = list(rawage))

rgnmap <- tibble::as_tibble(expand.grid(rungs = c(10, 25, 50),
                                        GTI_pow = c(2, 2.5, 3, 3.5, 4.0, 4.5, 5.0, 5.5, 6),
                                        burnin = 1e2,
                                        samples = 1e2)) %>%
  dplyr::mutate(lvl = "region",
                num_mas = 17,
                maxMa = "ma14",
                groupvar = "region",
                dat = list(rawrgn))

param_map <- dplyr::bind_rows(agemap, rgnmap)
param_map$modelobj <- purrr::pmap(param_map[, c("num_mas", "maxMa", "groupvar", "dat")], make_IFR_model_spain)



#............................................................
# MCMC Object
#...........................................................
run_MCMC <- function(modelobj, rungs, GTI_pow, burnin, samples) {
  start <- Sys.time()
  n_chains <- 5
  n_cores <- parallel::detectCores()

  if (n_cores < n_chains) {
    mkcores <- n_cores/2
  } else {
    mkcores <- n_chains/2
  }
  #......................
  # make cluster object to parallelize chains
  #......................
  cl <- parallel::makeCluster(mkcores)
  fit <- COVIDCurve::run_IFRmodel_agg(IFRmodel =modelobj,
                                      reparamIFR = TRUE,
                                      reparamInfxn = TRUE,
                                      reparamKnots = TRUE,
                                      chains = n_chains,
                                      burnin = burnin,
                                      samples = samples,
                                      rungs = rungs,
                                      GTI_pow = GTI_pow,
                                      cluster = cl)
  mc_accept_mean <- mean(fit$mcmcout$diagnostics$mc_accept$value)
  mc_accept_min <- min(fit$mcmcout$diagnostics$mc_accept$value)
  time_elapse <- Sys.time() - start
  # out
  out <- list(fit = fit,
              mc_accept_mean = mc_accept_mean,
              mc_accept_min = mc_accept_min,
              time_elapse = time_elapse)

  return(out)
}

#......................
# send out on slurm
#......................
clstparam_map <- param_map %>%
  dplyr::select(c("modelobj", "rungs", "GTI_pow", "burnin", "samples"))

ntry <- 60
sjob <- rslurm::slurm_apply(f = run_MCMC,
                            params = clstparam_map,
                            jobname = 'rslurm_parallel_test',
                            nodes = ntry,
                            cpus_per_node = 1,
                            preschedule_cores = FALSE,
                            submit = TRUE,
                            slurm_options = list(mem = "24g",
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 mail-type = "END",
                                                 mail-user = "nbrazeau@med.unc.edu",
                                                 time = "36:00:00"))
