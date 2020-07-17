####################################################################################
## Purpose: This is a script that will wrap a DRAKE make process to
##          run multiple ESP data iterations to fine tune the MCoupling
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
tod_paramsdf <- tibble::tibble(name = c("mod", "sod"),
                               min  = c(10,    0.01),
                               init = c(14,    0.7),
                               max =  c(20,    1.00),
                               dsc1 = c(2.7,   -0.23),
                               dsc2 = c(0.05,   0.05))
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1"),
                                min =   c(0.83,    0.50,     0,           117),
                                init =  c(0.85,    0.99,     0.1,          125),
                                max =   c(0.87,    1.00,     1,          131),
                                dsc1 =  c(886.5,   520.5,     70,        117),
                                dsc2 =  c(114.5,   10.5,      30,        131))
#......................
# regions
#......................
rawrgn <- readRDS("data/derived/ESP/ESP_regions.RDS")
ESP_rgn_mod <- make_IFR_model_fit(num_mas = 17, maxMa = "ma14",
                                  groupvar = "region",  dat = rawrgn,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 125, init = 131, max = 137, dsc1 = 125, dsc2 = 137),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")

#......................
# agebands
#......................
rawage <- readRDS("data/derived/ESP/ESP_agebands.RDS")
ESP_age_mod <- make_IFR_model_fit(num_mas = 10, maxMa = "ma10",
                                  groupvar = "ageband",  dat = rawage,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 125, init = 131, max = 137, dsc1 = 125, dsc2 = 137),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")

#......................
# make param maps
#......................
param_map <- tibble::as_tibble(expand.grid(rungs = c(10, 25, 50),
                                           GTI_pow = c(2, 2.5, 3, 3.5, 4.0, 4.5, 5.0, 5.5, 6),
                                           burnin = 1e4,
                                           samples = 1e4))
# age bands
age_mod_map <- tibble::tibble(name = c("ESP_agebands"),
                              modelobj = list(ESP_age_mod))
age_mod_map <- dplyr::bind_cols(age_mod_map, param_map)

# regions
rgn_mod_map <- tibble::tibble(name = c("ESP_rgns"),
                              modelobj = list(ESP_rgn_mod))
rgn_mod_map <- dplyr::bind_cols(rgn_mod_map, param_map)

# bring together
mod_param_map <- dplyr::bind_rows(age_mod_map, rgn_mod_map)


# select what we need for fits and make outpaths
dir.create("data/param_map/ESP_MixPower_Fits/", recursive = T)
lapply(split(mod_param_map, 1:nrow(mod_param_map)), function(x){
  saveRDS(x, paste0("data/param_map/ESP_MixPower_Fits/",
                    x$name, "_GTI", x$GTI_pow, "_rung", x$rungs, "_burn", x$burnin, "_smpl", x$samples, ".RDS"))
})


#............................................................
# MCMC Object
#...........................................................
run_MCMC <- function(path) {
  mod <- readRDS(path)
  #......................
  # make cluster object to parallelize chains
  #......................
  start <- Sys.time()
  n_chains <- 10
  n_cores <- parallel::detectCores()

  if (n_cores < n_chains) {
    mkcores <- n_cores - 1
  } else {
    mkcores <- n_chains
  }

  cl <- parallel::makeCluster(mkcores)
  fit <- COVIDCurve::run_IFRmodel_agg(IFRmodel = mod$modelobj[[1]],
                                      reparamIFR = TRUE,
                                      reparamInfxn = TRUE,
                                      reparamKnots = TRUE,
                                      chains = n_chains,
                                      burnin = mod$burnin,
                                      samples = mod$samples,
                                      rungs = mod$rungs,
                                      GTI_pow = mod$GTI_pow,
                                      cluster = cl)
  parallel::stopCluster(cl)
  gc()

  mc_accept_mean <- mean(fit$mcmcout$diagnostics$mc_accept$value)
  mc_accept_min <- min(fit$mcmcout$diagnostics$mc_accept$value)
  time_elapse <- Sys.time() - start
  # out
  dir.create("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/ESP_MixPower_Fits/", recursive = TRUE)
  outpath = paste0("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/ESP_MixPower_Fits/",
                   mod$name, "_GTI", mod$GTI_pow, "_rung", mod$rungs, "_burn", mod$burnin, "_smpl", mod$samples, ".RDS")
  out <- list(fit = fit,
              mc_accept_mean = mc_accept_mean,
              mc_accept_min = mc_accept_min,
              time_elapse = time_elapse)

  saveRDS(out, file = outpath)

  return(0)
}


#............................................................
# Make Drake Plan
#...........................................................
# due to R6 classes being stored in environment https://github.com/ropensci/drake/issues/961
# Drake can't find <environment> in memory (obviously).
# Need to either wrap out of figure out how to nest better

# read files in after sleeping to account for file lag
Sys.sleep(60)
file_param_map <- list.files(path = "data/param_map/ESP_MixPower_Fits/",
                             pattern = "*.RDS",
                             full.names = TRUE)
file_param_map <- tibble::tibble(path = file_param_map)


#............................................................
# Make Drake Plan
#...........................................................
plan <- drake::drake_plan(
  fits = target(
    run_MCMC(path),
    transform = map(
      .data = !!file_param_map
    )
  )
)


#......................
# call drake to send out to slurm
#......................
options(clustermq.scheduler = "slurm",
        clustermq.template = "drake_workers/slurm_clustermq_LL.tmpl")
make(plan, parallelism = "clustermq", jobs = nrow(file_param_map),
     log_make = "ESP_PowerFits_drake.log", verbose = 2,
     log_progress = FALSE,
     log_build_times = FALSE,
     recoverable = FALSE,
     history = FALSE,
     session_info = FALSE,
     lock_envir = FALSE, # unlock environment so parallel::clusterApplyLB in drjacoby can work
     lock_cache = FALSE)


