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

  sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day"),
                                  min =   c(0.83,     0.97,   10,         117),
                                  init =  c(0.85,     0.99,   10,         125),
                                  max =   c(0.87,     1.00,   10,         131),
                                  dsc1 =  c(8500,     9900,    5,         117),
                                  dsc2 =  c(1500,     100,     15,        131))

  # bring together
  df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, knot_paramsdf, sens_spec_tbl)
  # format data
  inputdata <- list(obs_deaths = age$deaths,
                    obs_serologyrate = age$seroprev$seroprev)
  # make pa
  pa <- age$prop_pop$pop_prop * age$seroprev_group$seroprev
  pa <- pa/sum(pa)
  # make mod
  mod1 <- COVIDCurve::make_IFRmodel_agg$new()
  mod1$set_MeanOnset(18.8)
  mod1$set_CoefVarOnset(0.45)
  mod1$set_level("Time-Series")
  mod1$set_data(inputdata)
  mod1$set_IFRparams(paste0("ma", 1:10))
  mod1$set_maxMa("ma10")
  mod1$set_Knotparams(paste0("x", 1:4))
  mod1$set_relKnot("x4")
  mod1$set_Infxnparams(paste0("y", 1:5))
  mod1$set_relInfxn("y3")
  mod1$set_Seroparams(c("sens", "spec", "sero_rate", "sero_day"))
  mod1$set_popN(age$popN)
  mod1$set_paramdf(df_params)
  mod1$set_pa(pa)
  mod1$set_rcensor_day(.Machine$integer.max)
  # out
  mod1
}

#......................
# make maps
#......................
map <- tibble::as_tibble(expand.grid(rungs = c(10, 20, 25, 30, 40, 50),
                         GTI_pow = c(2, 2.5, 3, 3.5, 4.0, 4.5, 5.0, 5.5, 6),
                         burnin = c(1e3, 1e4, 1e5),
                         samples = c(1e3, 1e4, 1e5)))


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
                                      GTI_pow = GTI_pow)
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
ntry <- 360
sjob <- rslurm::slurm_apply(f = run_wrapper,
                            params = map,
                            jobname = 'tuneMCMC_ESPage',
                            nodes = ntry,
                            cpus_per_node = 1,
                            submit = T,
                            slurm_options = list(mem = "24g",
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "2-00:00:00"))


