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
# onset to deaths
tod_paramsdf <- tibble::tibble(name = c("mod", "sod"),
                               min  = c(10,     0.01),
                               init = c(14,     0.7),
                               max =  c(20,     1.00),
                               dsc1 = c(2.657,  -0.236),
                               dsc2 = c(0.01,   0.01))


#............................................................
# ESP
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1"),
                                min =   c(0.83,    0.50,     0,          117),
                                init =  c(0.85,    0.99,     0.7,        125),
                                max =   c(0.87,    1.00,     1,          131),
                                dsc1 =  c(123.5,   156.5,    700,        117),
                                dsc2 =  c(30.5,    0.5,      300,        131))
# https://www.thelancet.com/cms/10.1016/S0140-6736(20)31483-5/attachment/25c80941-a8c5-470e-a6a8-fde7397b9547/mmc1.pdf
# based on supp table 3
#......................
# regions
#......................
rawrgn <- readRDS("data/derived/ESP/ESP_regions.RDS")
ESP_rgn_mod <- make_IFR_model_fit(num_mas = 17, maxMa = "ma14",
                                  groupvar = "region",  dat = rawrgn,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 123, init = 131, max = 137, dsc1 = 123, dsc2 = 137),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")

#......................
# agebands
#......................
rawage <- readRDS("data/derived/ESP/ESP_agebands.RDS")
ESP_age_mod <- make_IFR_model_fit(num_mas = 10, maxMa = "ma10",
                                  groupvar = "ageband",  dat = rawage,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 123, init = 131, max = 137, dsc1 = 123, dsc2 = 137),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")


#............................................................
# Denmark
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1"),
                                min =   c(0.50,    0.50,     0,           96),
                                init =  c(0.85,    0.99,     0.7,          100),
                                max =   c(1.00,    1.00,     1,          107),
                                dsc1 =  c(128.5,   647.5,    700,        96),
                                dsc2 =  c(27.5,    4.5,      300,        107))
#......................
# regions
#......................
rawrgn <- readRDS("data/derived/DNK/DNK_regions.RDS")
DNK_rgn_mod <- make_IFR_model_fit(num_mas = 3, maxMa = "ma1",
                                  groupvar = "region",  dat = rawrgn,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 124, init = 131, max = 138, dsc1 = 124, dsc2 = 138),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")
#......................
# agebands
#......................
rawage <- readRDS("data/derived/DNK/DNK_agebands.RDS")
DNK_age_mod <- make_IFR_model_fit(num_mas = 5, maxMa = "ma5",
                                  groupvar = "ageband",  dat = rawage,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 124, init = 131, max = 138, dsc1 = 124, dsc2 = 138),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")

#............................................................
# Netherlands
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1"),
                                min =   c(0.50,    0.50,     0,          91),
                                init =  c(0.85,    0.99,     0.7,        100),
                                max =   c(1.00,    1.00,     1,          105),
                                dsc1 =  c(171.5,   281.5,    700,        91),
                                dsc2 =  c(3.5,     1.5,      300,        105))
#......................
# regions
#......................
rawrgn <- readRDS("data/derived/NLD/NLD_regions.RDS")
NLD_rgn_mod <- make_IFR_model_fit(num_mas = 25, maxMa = "ma1",
                                  groupvar = "region",  dat = rawrgn,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 124, init = 131, max = 138, dsc1 = 124, dsc2 = 138),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")

#......................
# agebands
#......................
rawage <- readRDS("data/derived/NLD/NLD_agebands.RDS")
NLD_age_mod <- make_IFR_model_fit(num_mas = 10, maxMa = "ma10",
                                  groupvar = "ageband",  dat = rawage,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 124, init = 131, max = 138, dsc1 = 124, dsc2 = 138),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")

#............................................................
# Brazil
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1"),
                                min =   c(0.50,    0.50,     0,          64),
                                init =  c(0.85,    0.99,     0.7,        68),
                                max =   c(1.00,    1.00,     1,          71),
                                dsc1 =  c(850.5,   990.5,    700,        64),
                                dsc2 =  c(150.5,    10.5,    300,        71))
#......................
# regions
#......................
rawrgn <- readRDS("data/derived/BRA/BRA_regions.RDS")
BRA_rgn_mod <- make_IFR_model_fit(num_mas = 5, maxMa = "ma5",
                                  groupvar = "region",  dat = rawrgn,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 74, init = 81, max = 88, dsc1 = 74, dsc2 = 88),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")
#......................
# agebands
#......................
rawage <- readRDS("data/derived/BRA/BRA_agebands.RDS")
BRA_age_mod <- make_IFR_model_fit(num_mas = 10, maxMa = "ma10",
                                  groupvar = "ageband",  dat = rawage,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 74, init = 81, max = 88, dsc1 = 74, dsc2 = 88),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")


#............................................................
# New York City
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1"),
                                min =   c(0.50,    0.50,     0,          109),
                                init =  c(0.85,    0.99,     0.7,        114),
                                max =   c(1.00,    1.00,     1,          118),
                                dsc1 =  c(204.5,   990.5,    700,        109),
                                dsc2 =  c(234.5,   10.5,     300,        118))
#......................
# regions
#......................
rawage <- readRDS("data/derived/USA/NYC_NY_1_agebands.RDS")
NYC_age_mod <- make_IFR_model_fit(num_mas = 5, maxMa = "ma5",
                                  groupvar = "ageband",  dat = rawage,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 174, init = 181, max = 188, dsc1 = 174, dsc2 = 188),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")


#............................................................
# Great Britian
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1"),
                                min =   c(0.50,    0.50,     0,          109),
                                init =  c(0.99,    0.99,     0.7,        114),
                                max =   c(1.00,    1.00,     1,          118),
                                dsc1 =  c(990.5,   990.5,    700,        109),
                                dsc2 =  c(10.5,   10.5,     300,        118))
#......................
# regions
#......................
rawrgn <- readRDS("data/derived/UK/GBR_regions.RDS")
GBR_rgn_mod <- make_IFR_model_fit(num_mas = 7, maxMa = "ma1",
                                  groupvar = "region",  dat = rawrgn,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 174, init = 181, max = 188, dsc1 = 174, dsc2 = 188),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")

#......................
# agebands
#......................
rawage <- readRDS("data/derived/UK/GBR_agebands.rds")
GBR_age_mod <- make_IFR_model_fit(num_mas = 10, maxMa = "ma10",
                                  groupvar = "ageband",  dat = rawage,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 174, init = 181, max = 188, dsc1 = 174, dsc2 = 188),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")



#............................................................
# Come Together
#...........................................................
fit_map <- tibble::tibble(
  name = c("ESP_age", "ESP_rgn", "DNK_age", "DNK_rgn",
           "NLD_age", "NLD_rgn", "BRA_age", "BRA_rgn",
           "NYC_age", "GBR_age", "GBR_rgn"),
  modelobj = list(ESP_age_mod, ESP_rgn_mod, DNK_age_mod, DNK_rgn_mod,
                  NLD_age_mod, NLD_rgn_mod, BRA_age_mod, BRA_rgn_mod,
                  NYC_age_mod, GBR_age_mod, GBR_rgn_mod),
  rungs = 50,
  GTI_pow = c(2.5, 4.5, 3, 3, 3, 3, 3, 3, 3, 3, 3),
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
    mkcores <- n_chains
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




