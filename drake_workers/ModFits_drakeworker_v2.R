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
                               min  = c(10,     0),
                               init = c(14,     0.7),
                               max =  c(30,     1),
                               dsc1 = c(2.66,   50),
                               dsc2 = c(0.25,   50))

#............................................................
#---- BRA1 #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate"),
                                min =   c(0.50,    0.50,     10),
                                init =  c(0.85,    0.99,     15),
                                max =   c(1.00,    1.00,     30),
                                dsc1 =  c(850.5,   990.5,    2.8),
                                dsc2 =  c(150.5,    10.5,    0.1))
#......................
# regions
#......................
rawrgn <- readRDS("data/derived/BRA/BRA_regions.RDS")
BRA_rgn_mod <- make_IFR_model_fit(num_mas = 5, maxMa = "ma2", # Northern Region
                                  groupvar = "region",  dat = rawrgn,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 174, init = 181, max = 188, dsc1 = 174, dsc2 = 188),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 19.17, dsc1 = 0, dsc2 = 19.17),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)
#......................
# agebands
#......................
rawage <- readRDS("data/derived/BRA/BRA_agebands.RDS")
BRA_age_mod <- make_IFR_model_fit(num_mas = 10, maxMa = "ma10",
                                  groupvar = "ageband",  dat = rawage,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 174, init = 181, max = 188, dsc1 = 174, dsc2 = 188),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 19.17, dsc1 = 0, dsc2 = 19.17),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)


#............................................................
#---- CHE1 #-----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate"),
                                min =   c(0.50,    0.50,     10),
                                init =  c(0.85,    0.99,     15),
                                max =   c(1.00,    1.00,     30),
                                dsc1 =  c(156.5,   176.5,    2.8),
                                dsc2 =  c(25.5,    0.5,      0.1))
#......................
# regions
#......................
rawrgn <- readRDS("data/derived/CHE/CHE_region.RDS")
CHE_rgn_mod <- make_IFR_model_fit(num_mas = 1, maxMa = "ma1", # Geneva (only 1 region)
                                  groupvar = "region",  dat = rawrgn,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 134, init = 140, max = 146, dsc1 = 134, dsc2 = 146),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 13.12, dsc1 = 0, dsc2 = 13.12),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)
#......................
# agebands
#......................
rawage <- readRDS("data/derived/CHE/CHE_agebands.RDS")
CHE_age_mod <- make_IFR_model_fit(num_mas = 9, maxMa = "ma9",
                                  groupvar = "ageband",  dat = rawage,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 134, init = 140, max = 146, dsc1 = 134, dsc2 = 146),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 13.12, dsc1 = 0, dsc2 = 13.12),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)


#............................................................
#---- DNK1 #-----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate"),
                                min =   c(0.50,    0.50,     10),
                                init =  c(0.85,    0.99,     15),
                                max =   c(1.00,    1.00,     30),
                                dsc1 =  c(128.5,   647.5,    2.8),
                                dsc2 =  c(27.5,    4.5,      0.1))
#......................
# regions
#......................
rawrgn <- readRDS("data/derived/DNK/DNK_regions.RDS")
DNK_rgn_mod <- make_IFR_model_fit(num_mas = 3, maxMa = "ma1", # Capital Region
                                  groupvar = "region",  dat = rawrgn,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 192, init = 199, max = 206, dsc1 = 192, dsc2 = 206),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 15.57, dsc1 = 0, dsc2 = 15.57),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)
#......................
# agebands
#......................
rawage <- readRDS("data/derived/DNK/DNK_agebands.RDS")
DNK_age_mod <- make_IFR_model_fit(num_mas = 5, maxMa = "ma5",
                                  groupvar = "ageband",  dat = rawage,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 192, init = 199, max = 206, dsc1 = 192, dsc2 = 206),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 15.57, dsc1 = 0, dsc2 = 15.57),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)


#............................................................
#---- ESP1-2 #-----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate"),
                                min =   c(0.83,    0.50,     10),
                                init =  c(0.85,    0.99,     15),
                                max =   c(1.00,    1.00,     30),
                                dsc1 =  c(123.5,   156.5,    2.8),
                                dsc2 =  c(30.5,    0.5,      0.1))
# https://www.thelancet.com/cms/10.1016/S0140-6736(20)31483-5/attachment/25c80941-a8c5-470e-a6a8-fde7397b9547/mmc1.pdf
# based on supp table 3
#......................
# regions
#......................
# NB in ecdc data, DNK was one day "ahead" of ESP
rawrgn <- readRDS("data/derived/ESP/ESP_regions.RDS")
ESP_rgn_mod <- make_IFR_model_fit(num_mas = 17, maxMa = "ma17", # La Rioja Region
                                  groupvar = "region",  dat = rawrgn,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 193, init = 200, max = 205, dsc1 = 193, dsc2 = 205),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 17.66, dsc1 = 0, dsc2 = 17.66),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)

#......................
# agebands
#......................
rawage <- readRDS("data/derived/ESP/ESP_agebands.RDS")
ESP_age_mod <- make_IFR_model_fit(num_mas = 10, maxMa = "ma10",
                                  groupvar = "ageband",  dat = rawage,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 193, init = 200, max = 205, dsc1 = 193, dsc2 = 205),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 17.66, dsc1 = 0, dsc2 = 17.66),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)


#...........................................................
#---- GBR2 #-----
#...........................................................
# sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1"),
#                                 min =   c(0.50,    0.50,     0.5,        94),
#                                 init =  c(0.99,    0.99,     1,          107),
#                                 max =   c(1.00,    1.00,     1.5,         120),
#                                 dsc1 =  c(990.5,   990.5,    0.5,        94),
#                                 dsc2 =  c(10.5,   10.5,      1.5,         120))
# #......................
# # regions
# #......................
# rawrgn <- readRDS("data/derived/UK/GBR_regions.RDS")
# GBR_rgn_mod <- make_IFR_model_fit(num_mas = 7, maxMa = "ma4",
#                                   groupvar = "region",  dat = rawrgn,
#                                   num_xs = 4, max_xveclist = list("name" = "x4", min = 174, init = 181, max = 188, dsc1 = 174, dsc2 = 188),
#                                   num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 17.82, dsc1 = 0, dsc2 = 17.82),
#                                   sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
#                                   serodayparams = "sero_day1")
#
# #......................
# # agebands
# #......................
# rawage <- readRDS("data/derived/UK/GBR_agebands.RDS")
# GBR_age_mod <- make_IFR_model_fit(num_mas = 10, maxMa = "ma10",
#                                   groupvar = "ageband",  dat = rawage,
#                                   num_xs = 4, max_xveclist = list("name" = "x4", min = 174, init = 181, max = 188, dsc1 = 174, dsc2 = 188),
#                                   num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 17.82, dsc1 = 0, dsc2 = 17.82),
#                                   sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
#                                   serodayparams = "sero_day1")
#


#............................................................
#---- NLD1 #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate"),
                                min =   c(0.50,    0.50,     10),
                                init =  c(0.85,    0.99,     15),
                                max =   c(1.00,    1.00,     30),
                                dsc1 =  c(171.5,   281.5,    2.8),
                                dsc2 =  c(3.5,     1.5,      0.1))
#......................
# regions
#......................


#......................
# agebands
#......................
rawage <- readRDS("data/derived/NLD/NLD_agebands.RDS")
NLD_age_mod <- make_IFR_model_fit(num_mas = 6, maxMa = "ma6",
                                  groupvar = "ageband",  dat = rawage,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 192, init = 199, max = 206, dsc1 = 192, dsc2 = 206),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 16.67, dsc1 = 0, dsc2 = 16.67),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)



# #............................................................
# #---- ITA1 #----
# #...........................................................
# sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate"),
#                                 min =   c(0.50,    0.50,     10),
#                                 init =  c(0.85,    0.99,     15),
#                                 max =   c(1.00,    1.00,     30),
#                                 dsc1 =  c(90.5,   95.5,      2.8),
#                                 dsc2 =  c(10.5,     5.5,      0.1))
# #......................
# #regions
# #......................
# rawrgn <- readRDS("data/derived/ITA/ITA_regions.RDS")
# ITA_rgn_mod <- make_IFR_model_fit(num_mas = 21, maxMa = "ma10", # Lombardia (Lombardy) Region
#                                   groupvar = "region",  dat = rawrgn,
#                                   num_xs = 4, max_xveclist = list("name" = "x4", min = 192, init = 199, max = 206, dsc1 = 192, dsc2 = 206),
#                                   num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 17.91, dsc1 = 0, dsc2 = 17.91),
#                                   sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)
#
# #......................
# # agebands
# #......................
# rawage <- readRDS("data/derived/ITA/ITA_agebands.RDS")
# ITA_age_mod <- make_IFR_model_fit(num_mas = 10, maxMa = "ma10",
#                                   groupvar = "ageband",  dat = rawage,
#                                   num_xs = 4, max_xveclist = list("name" = "x4", min = 192, init = 199, max = 206, dsc1 = 192, dsc2 = 206),
#                                   num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 17.91, dsc1 = 0, dsc2 = 17.91),
#                                   sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)


#............................................................
#---- LUX1 #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate"),
                                min =   c(0.50,    0.50,    10),
                                init =  c(0.85,    0.99,    15),
                                max =   c(1.00,    1.00,    30),
                                dsc1 =  c(12.5,   181.5,    2.8),
                                dsc2 =  c(14.5,     4.5,    0.1))

#......................
# agebands
#......................
rawage <- readRDS("data/derived/LUX/LUX_agebands.RDS")
LUX_age_mod <- make_IFR_model_fit(num_mas = 7, maxMa = "ma7",
                                  groupvar = "ageband",  dat = rawage,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 192, init = 199, max = 206, dsc1 = 192, dsc2 = 206),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 13.34, dsc1 = 0, dsc2 = 13.34),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)
# fixed inside code too?
LUX_age_mod$data$obs_deaths$Deaths <- ifelse(is.na(LUX_age_mod$data$obs_deaths$Deaths), -1, LUX_age_mod$data$obs_deaths$Deaths)

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
#---- CHN1 #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate"),
                                min =   c(0.50,    0.50,     10),
                                init =  c(0.85,    0.97,     15),
                                max =   c(1.00,    1.00,     30),
                                dsc1 =  c(91.5,    198.5,    2.8),
                                dsc2 =  c(14.5,    202.5,    0.1))

#......................
# agebands
#......................
rawage <- readRDS("data/derived/CHN/CHN_agebands.RDS")
CHN_age_mod <- make_IFR_model_fit(num_mas = 9, maxMa = "ma9",
                                  groupvar = "ageband",  dat = rawage,
                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 199, init = 206, max = 213, dsc1 = 199, dsc2 = 213),
                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 16.20, dsc1 = 0, dsc2 = 16.20),
                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)





#............................................................
#---- Come Together #----
#...........................................................
bvec <- seq(5, 2.5, length.out = 50)

fit_map <- tibble::tibble(
  name = c("BRA1_age", "BRA1_rgn",
           "CHE1_age",
           "DNK1_age", "DNK1_rgn",
           "ESP1-2_age", "ESP1-2_rgn",
           #"GBR_age", "GBR_rgn",
           "NLD1_age",
           #"ITA1_age", "ITA1_rgn",
           "LUX1_age",
           "CHN1_age"),
  modelobj = list(BRA_age_mod, BRA_rgn_mod,
                  CHE_age_mod,
                  DNK_age_mod, DNK_rgn_mod,
                  ESP_age_mod, ESP_rgn_mod,
                  #GBR_age_mod, GBR_rgn_mod,
                  NLD_age_mod,
                  #ITA_age_mod, ITA_rgn_mod,
                  LUX_age_mod,
                  CHN_age_mod),
  rungs = 50,
  GTI_pow = list(bvec),
  burnin = 1e4,
  samples = 1e4,
  thinning = 10)



# select what we need for fits and make outpaths
dir.create("data/param_map/Modfits/", recursive = T)
lapply(split(fit_map, 1:nrow(fit_map)), function(x){
  saveRDS(x, paste0("data/param_map/Modfits/",
                    x$name, "_rung", x$rungs, "_burn", x$burnin, "_smpl", x$samples, ".RDS"))
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
                                      reparamSeros = FALSE,
                                      reparamNe = TRUE,
                                      chains = n_chains,
                                      burnin = mod$burnin,
                                      samples = mod$samples,
                                      rungs = mod$rungs,
                                      GTI_pow = mod$GTI_pow[[1]],
                                      cluster = cl,
                                      thinning = 10)
  parallel::stopCluster(cl)
  gc()

  # out
  dir.create("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/Modfits/", recursive = TRUE)
  outpath = paste0("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/Modfits/",
                   mod$name, "_rung", mod$rungs, "_burn", mod$burnin, "_smpl", mod$samples, ".RDS")
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




