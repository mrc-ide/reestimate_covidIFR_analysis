#........................................................................................
## Purpose: Use Drake to fit IFR Models
##
## Notes:
#........................................................................................
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
#---- BRA1 #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1"),
                                min =   c(0.50,    0.50,     0,          136),
                                init =  c(0.85,    0.99,     0.9,        140),
                                max =   c(1.00,    1.00,     1,          143),
                                dsc1 =  c(850.5,   990.5,    900,        136),
                                dsc2 =  c(150.5,    10.5,    100,        143))
#......................
# regions
#......................
rawrgn <- readRDS("data/derived/BRA/BRA_regions.RDS")
BRA_rgn_mod <- make_IFR_model_fit(num_mas = 5, maxMa = "ma5",
                                  groupvar = "region",  dat = rawrgn,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 174, init = 181, max = 188, dsc1 = 174, dsc2 = 188),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")
#......................
# agebands
#......................
rawage <- readRDS("data/derived/BRA/BRA_agebands.RDS")
BRA_age_mod <- make_IFR_model_fit(num_mas = 10, maxMa = "ma10",
                                  groupvar = "ageband",  dat = rawage,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 174, init = 181, max = 188, dsc1 = 174, dsc2 = 188),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")


#............................................................
#---- CHE1 #-----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1", "sero_day2", "sero_day3", "sero_day4", "sero_day5"),
                                min =   c(0.50,    0.50,     0,         97,          104,         111,         118,         124),
                                init =  c(0.85,    0.99,     0.9,       100,         108,         115,         121,         127),
                                max =   c(1.00,    1.00,     1,         104,         111,         118,         124,         131),
                                dsc1 =  c(156.5,   176.5,    900,       97,          104,         111,         118,         124),
                                dsc2 =  c(25.5,    0.5,      100,       104,         111,         118,         124,         131))
#......................
# regions
#......................
rawrgn <- readRDS("data/derived/CHE/CHE_region.RDS")
CHE_rgn_mod <- make_IFR_model_fit(num_mas = 1, maxMa = "ma1",
                                  groupvar = "region",  dat = rawrgn,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 134, init = 140, max = 146, dsc1 = 134, dsc2 = 146),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = c("sero_day1", "sero_day2", "sero_day3", "sero_day4", "sero_day5"))
#......................
# agebands
#......................
rawage <- readRDS("data/derived/CHE/CHE_agebands.RDS")
CHE_age_mod <- make_IFR_model_fit(num_mas = 9, maxMa = "ma9",
                                  groupvar = "ageband",  dat = rawage,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 134, init = 140, max = 146, dsc1 = 134, dsc2 = 146),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = c("sero_day1", "sero_day2", "sero_day3", "sero_day4", "sero_day5"))


#............................................................
#---- DNK1 #-----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1"),
                                min =   c(0.50,    0.50,     0,          97),
                                init =  c(0.85,    0.99,     0.9,        110),
                                max =   c(1.00,    1.00,     1,          124),
                                dsc1 =  c(128.5,   647.5,    900,        97),
                                dsc2 =  c(27.5,    4.5,      100,        124))
#......................
# regions
#......................
rawrgn <- readRDS("data/derived/DNK/DNK_regions.RDS")
DNK_rgn_mod <- make_IFR_model_fit(num_mas = 3, maxMa = "ma1",
                                  groupvar = "region",  dat = rawrgn,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 124, init = 131, max = 139, dsc1 = 124, dsc2 = 139),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")
#......................
# agebands
#......................
rawage <- readRDS("data/derived/DNK/DNK_agebands.RDS")
DNK_age_mod <- make_IFR_model_fit(num_mas = 5, maxMa = "ma5",
                                  groupvar = "ageband",  dat = rawage,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 124, init = 131, max = 139, dsc1 = 124, dsc2 = 139),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")


#............................................................
#---- ESP1-2 #-----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1", "sero_day2"),
                                min =   c(0.83,    0.50,     0,          118,        139),
                                init =  c(0.85,    0.99,     0.7,        125,        146),
                                max =   c(1.00,    1.00,     1,          132,        153),
                                dsc1 =  c(123.5,   156.5,    700,        118,        139),
                                dsc2 =  c(30.5,    0.5,      300,        132,        153))
# https://www.thelancet.com/cms/10.1016/S0140-6736(20)31483-5/attachment/25c80941-a8c5-470e-a6a8-fde7397b9547/mmc1.pdf
# based on supp table 3
#......................
# regions
#......................
# NB in ecdc data, DNK was one day "ahead" of ESP
rawrgn <- readRDS("data/derived/ESP/ESP_regions.RDS")
ESP_rgn_mod <- make_IFR_model_fit(num_mas = 17, maxMa = "ma14",
                                  groupvar = "region",  dat = rawrgn,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 193, init = 200, max = 205, dsc1 = 193, dsc2 = 205),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = c("sero_day1", "sero_day2"))

#......................
# agebands
#......................
rawage <- readRDS("data/derived/ESP/ESP_agebands.RDS")
ESP_age_mod <- make_IFR_model_fit(num_mas = 10, maxMa = "ma10",
                                  groupvar = "ageband",  dat = rawage,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 193, init = 200, max = 205, dsc1 = 193, dsc2 = 205),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = c("sero_day1", "sero_day2"))


#...........................................................
#---- GBR2 #-----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1"),
                                min =   c(0.50,    0.50,     0,          94),
                                init =  c(0.99,    0.99,     0.7,        107),
                                max =   c(1.00,    1.00,     1,          120),
                                dsc1 =  c(990.5,   990.5,    700,        94),
                                dsc2 =  c(10.5,   10.5,     300,         120))
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
rawage <- readRDS("data/derived/UK/GBR_agebands.RDS")
GBR_age_mod <- make_IFR_model_fit(num_mas = 10, maxMa = "ma10",
                                  groupvar = "ageband",  dat = rawage,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 174, init = 181, max = 188, dsc1 = 174, dsc2 = 188),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                  serodayparams = "sero_day1")



#............................................................
#---- NLD1 #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1"),
                                min =   c(0.50,    0.50,     0,          92),
                                init =  c(0.85,    0.99,     0.9,        100),
                                max =   c(1.00,    1.00,     1,          106),
                                dsc1 =  c(171.5,   281.5,    900,        92),
                                dsc2 =  c(3.5,     1.5,      100,        106))
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



# #............................................................
# # New York City
# #...........................................................
# sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1", "sero_day2"),
#                                 min =   c(0.50,    0.50,     0,          82,          109),
#                                 init =  c(0.85,    0.99,     0.7,        87,          114),
#                                 max =   c(1.00,    1.00,     1,          91,          118),
#                                 dsc1 =  c(204.5,   990.5,    700,        82,          109),
#                                 dsc2 =  c(234.5,   10.5,     300,        91,          118))
# #......................
# # agebands
# #......................
# rawage <- readRDS("data/derived/USA/NYC_NY_1_cdc1_agebands.RDS")
# NYC_age_mod <- make_IFR_model_fit(num_mas = 5, maxMa = "ma5",
#                                   groupvar = "ageband",  dat = rawage,
#                                   num_xs = 4, max_xveclist = list("name" = "x4", min = 174, init = 181, max = 188, dsc1 = 174, dsc2 = 188),
#                                   num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
#                                   sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
#                                   serodayparams = c("sero_day1", "sero_day2"))



#............................................................
#---- Come Together #----
#...........................................................
fit_map <- tibble::tibble(
  name = c("BRA_age", "BRA_rgn",
           "CHE1_age",
           "DNK_age", "DNK_rgn",
           "ESP_age", "ESP_rgn",
           "GBR_age", "GBR_rgn",
           "NLD_age", "NLD_rgn"),
  modelobj = list(BRA_age_mod, BRA_rgn_mod,
                  CHE_age_mod,
                  DNK_age_mod, DNK_rgn_mod,
                  ESP_age_mod, ESP_rgn_mod,
                  GBR_age_mod, GBR_rgn_mod,
                  NLD_age_mod, NLD_rgn_mod),
  rungs = 50,
  GTI_pow = 3,
  burnin = 1e4,
  samples = 1e4
)



# select what we need for fits and make outpaths
dir.create("data/param_map/Modfits/", recursive = T)
lapply(split(fit_map, 1:nrow(fit_map)), function(x){
  saveRDS(x, paste0("data/param_map/Modfits/",
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

  # out
  dir.create("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/Modfits/", recursive = TRUE)
  outpath = paste0("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/Modfits/",
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
file_param_map <- list.files(path = "data/param_map/Modfits/",
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
     log_make = "Modfits_drake.log", verbose = 2,
     log_progress = TRUE,
     log_build_times = FALSE,
     recoverable = FALSE,
     history = FALSE,
     session_info = FALSE,
     lock_envir = FALSE, # unlock environment so parallel::clusterApplyLB in drjacoby can work
     lock_cache = FALSE)



cat("************** Drake Finished **************************")




