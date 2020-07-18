####################################################################################
## Purpose: Use Drake to fit IFR Models
##
## Notes:
####################################################################################
setwd("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis")
library(drake)
library(parallel)
library(COVIDCurve)
library(tidyverse)
source("R/covidcurve_helper_functions.R")


#...................................................................................
# Make Paramset and write to disk for input into MCMC
#.................................................................................
# global
tod_paramsdf <- tibble::tibble(name = c("mod", "sod"),
                               min  = c(10,    0.01),
                               init = c(18,    0.45),
                               max =  c(25,    1.00),
                               dsc1 = c(2.9,   -0.78),
                               dsc2 = c(0.05,   0.05))

#............................................................
# ESP
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1"),
                                min =   c(0.83,    0.50,     5,           117),
                                init =  c(0.85,    0.99,     10,          125),
                                max =   c(0.87,    1.00,     15,          131),
                                dsc1 =  c(886.5,   52.5,     2.15,        117),
                                dsc2 =  c(114.5,   1.5,      0.05,        131))
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


#............................................................
# Denmark
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1"),
                                min =   c(0.50,    0.50,     5,           96),
                                init =  c(0.85,    0.99,     10,          100),
                                max =   c(1.00,    1.00,     15,          107),
                                dsc1 =  c(128.5,   647.5,     2.15,        96),
                                dsc2 =  c(27.5,    4.5,      0.05,        107))
#......................
# regions
#......................
rawrgn <- readRDS("data/derived/DNK/DNK_regions.RDS")
DNK_rgn_mod <- make_IFR_model_fit(num_mas = 3, maxMa = "ma1",
                                  groupvar = "region",  dat = rawrgn,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 125, init = 131, max = 138, dsc1 = 125, dsc2 = 138),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")
#......................
# agebands
#......................
rawage <- readRDS("data/derived/DNK/DNK_agebands.RDS")
DNK_age_mod <- make_IFR_model_fit(num_mas = 5, maxMa = "ma5",
                                  groupvar = "ageband",  dat = rawage,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 125, init = 131, max = 138, dsc1 = 125, dsc2 = 138),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")

#............................................................
# Netherlands
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1"),
                                min =   c(0.50,    0.50,     5,           91),
                                init =  c(0.85,    0.99,     10,          100),
                                max =   c(1.00,    1.00,     15,          105),
                                dsc1 =  c(128.5,   647.5,     2.15,        91),
                                dsc2 =  c(27.5,    4.5,      0.05,        105))
#......................
# regions
#......................
rawrgn <- readRDS("data/derived/NLD/NLD_regions.RDS")
NLD_rgn_mod <- make_IFR_model_fit(num_mas = 25, maxMa = "ma1",
                                  groupvar = "region",  dat = rawrgn,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 125, init = 131, max = 138, dsc1 = 125, dsc2 = 138),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")

#......................
# agebands
#......................
rawage <- readRDS("data/derived/NLD/NLD_agebands.RDS")
NLD_age_mod <- make_IFR_model_fit(num_mas = 10, maxMa = "ma10",
                                  groupvar = "ageband",  dat = rawage,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 125, init = 131, max = 138, dsc1 = 125, dsc2 = 138),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")

#............................................................
# Brazil
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1"),
                                min =   c(0.50,    0.50,     5,           64),
                                init =  c(0.85,    0.99,     10,          68),
                                max =   c(1.00,    1.00,     15,          71),
                                dsc1 =  c(850.5,   99.5,     2.15,        64),
                                dsc2 =  c(150.5,    1.5,     0.05,        71))
#......................
# regions
#......................
rawrgn <- readRDS("data/derived/BRA/BRA_regions.RDS")
BRA_rgn_mod <- make_IFR_model_fit(num_mas = 5, maxMa = "ma5",
                                  groupvar = "region",  dat = rawrgn,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 78, init = 85, max = 88, dsc1 = 78, dsc2 = 88),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")
#......................
# agebands
#......................
rawage <- readRDS("data/derived/BRA/BRA_agebands.RDS")
BRA_age_mod <- make_IFR_model_fit(num_mas = 10, maxMa = "ma10",
                                  groupvar = "ageband",  dat = rawage,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 78, init = 85, max = 88, dsc1 = 78, dsc2 = 88),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")

#............................................................
# Come Together
#...........................................................
fit_map <- tibble::tibble(
  name = c("ESP_age", "ESP_rgn", "DNK_age", "DNK_rgn", "NLD_age", "NLD_rgn", "BRA_age", "BRA_rgn"),
  modelobj = list(ESP_age_mod, ESP_rgn_mod, DNK_age_mod, DNK_rgn_mod, NLD_age_mod, NLD_rgn_mod, BRA_age_mod, BRA_rgn_mod),
  rungs = 50,
  GTI_pow = c(2.5, 4.5, 3, 3, 3, 3, 3, 3),
  burnin = 1e4,
  samples = 1e4
)



# select what we need for fits and make outpaths
dir.create("data/param_map/ModFits/", recursive = T)
lapply(split(fit_map, 1:nrow(fit_map)), function(x){
  saveRDS(x, paste0("data/param_map/ModFits/",
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
    mkcores <- n_chains/2
  }

  cl <- parallel::makeCluster(mkcores)
  fit <- COVIDCurve::run_IFRmodel_agg(IFRmodel = mod$modelobj[[1]],
                                      reparamIFR = TRUE,
                                      reparamInfxn = TRUE,
                                      reparamKnots = TRUE,
                                      reparamSpec = TRUE,
                                      chains = n_chains,
                                      burnin = mod$burnin,
                                      samples = mod$samples,
                                      rungs = mod$rungs,
                                      GTI_pow = mod$GTI_pow,
                                      cluster = cl)

  parallel::stopCluster(cl)
  gc()

  # out
  dir.create("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/ModFits/", recursive = TRUE)
  outpath = paste0("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/ModFits/",
                   mod$name, "_GTI", mod$GTI_pow, "_rung", mod$rungs, "_burn", mod$burnin, "_smpl", mod$samples, ".RDS")
  saveRDS(fit, file = outpath)

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
file_param_map <- list.files(path = "data/param_map/ModFits/",
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
     log_progress = TRUE,
     log_build_times = FALSE,
     recoverable = FALSE,
     history = FALSE,
     session_info = FALSE,
     lock_envir = FALSE, # unlock environment so parallel::clusterApplyLB in drjacoby can work
     lock_cache = FALSE)



cat("************** Drake Finished **************************")




