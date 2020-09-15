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
#......................
# onset to deaths
#......................
tod_paramsdf <- tibble::tibble(name = c("mod", "sod", "sero_con_rate"),
                               min  = c(18,     0,     16),
                               init = c(19,     0.90,  18),
                               max =  c(20,     1,     21),
                               dsc1 = c(19.66,  2700,  18.3),
                               dsc2 = c(0.1,    300,   0.1))

#............................................................
#---- BRA1 #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.85,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(850.5,   990.5),
                                dsc2 =  c(150.5,    10.5))

#......................
# agebands
#......................
rawage <- readRDS("data/derived/BRA1/BRA1_agebands.RDS")
BRA1_age_mod <- make_noSeroRev_IFR_model_fit(num_mas = 10, maxMa = "ma10",
                                             groupvar = "ageband",  dat = rawage,
                                             num_xs = 4, max_xveclist = list("name" = "x4", min = 175, init = 181, max = 189, dsc1 = 175, dsc2 = 189),
                                             num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 19.17, dsc1 = 0, dsc2 = 19.17),
                                             sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)



#............................................................
#---- BRA5 #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.85,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(447.5,   515.5),
                                dsc2 =  c(80.5,    5.5))

#......................
# basic
#......................
rawbase <- readRDS("data/derived/BRA5/BRA5_regions.RDS")
BRA5_rgn_mod <- make_noSeroRev_IFR_model_fit(num_mas = 1, maxMa = "ma1",
                                             groupvar = "region",  dat = rawbase,
                                             num_xs = 4, max_xveclist = list("name" = "x4", min = 175, init = 181, max = 189, dsc1 = 175, dsc2 = 189),
                                             num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 19.17, dsc1 = 0, dsc2 = 19.17),
                                             sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)


#............................................................
#---- CHE1 #-----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.85,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(156.5,   176.5),
                                dsc2 =  c(25.5,    0.5))

#......................
# agebands
#......................
rawage <- readRDS("data/derived/CHE1/CHE1_agebands.RDS")
CHE1_age_mod <- make_noSeroRev_IFR_model_fit(num_mas = 9, maxMa = "ma9",
                                             groupvar = "ageband",  dat = rawage,
                                             num_xs = 4, max_xveclist = list("name" = "x4", min = 134, init = 140, max = 146, dsc1 = 134, dsc2 = 146),
                                             num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 13.12, dsc1 = 0, dsc2 = 13.12),
                                             sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)


#............................................................
#---- CHE2 #-----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.85,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(49.5,   5497.5),
                                dsc2 =  c(5.5,    6.5))

#......................
# agebands
#......................
rawage <- readRDS("data/derived/CHE2/CHE2_agebands.RDS")
CHE2_age_mod <- make_noSeroRev_IFR_model_fit(num_mas = 9, maxMa = "ma9",
                                             groupvar = "ageband",  dat = rawage,
                                             num_xs = 4, max_xveclist = list("name" = "x4", min = 199, init = 206, max = 213, dsc1 = 199, dsc2 = 213),
                                             num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 14.24, dsc1 = 0, dsc2 = 14.24),
                                             sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)


#............................................................
#---- DNK1 #-----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.85,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(128.5,   647.5),
                                dsc2 =  c(27.5,    4.5))

#......................
# agebands
#......................
rawbase <- readRDS("data/derived/DNK1/DNK1_regions.RDS")
DNK_basic_mod <- make_noSeroRev_IFR_model_fit(num_mas = 1, maxMa = "ma1",
                                              groupvar = "region",  dat = rawbase,
                                              num_xs = 4, max_xveclist = list("name" = "x4", min = 219, init = 226, max = 233, dsc1 = 219, dsc2 = 233),
                                              num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 15.57, dsc1 = 0, dsc2 = 15.57),
                                              sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)


#............................................................
#---- ESP1-2 #-----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.83,    0.50),
                                init =  c(0.85,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(123.5,   156.5),
                                dsc2 =  c(30.5,    0.5))
# https://www.thelancet.com/cms/10.1016/S0140-6736(20)31483-5/attachment/25c80941-a8c5-470e-a6a8-fde7397b9547/mmc1.pdf
# based on supp table 3


#......................
# agebands
#......................
rawage <- readRDS("data/derived/ESP1-2/ESP1-2_agebands.RDS")
ESP_age_mod <- make_noSeroRev_IFR_model_fit(num_mas = 10, maxMa = "ma10",
                                            groupvar = "ageband",  dat = rawage,
                                            num_xs = 4, max_xveclist = list("name" = "x4", min = 219, init = 226, max = 233, dsc1 = 219, dsc2 = 233),
                                            num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 17.66, dsc1 = 0, dsc2 = 17.66),
                                            sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)


#...........................................................
#---- GBR3 #-----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.85,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(38.5,   493.5),
                                dsc2 =  c(10.5,    7.5))

#......................
# agebands
#......................
rawage <- readRDS("data/derived/GBR3/GBR3_agebands.RDS")
GBR3_age_mod <- make_noSeroRev_IFR_model_fit(num_mas = 4, maxMa = "ma4",
                                             groupvar = "ageband",  dat = rawage,
                                             num_xs = 4, max_xveclist = list("name" = "x4", min = 206, init = 210, max = 213, dsc1 = 199, dsc2 = 213),
                                             num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 17.857, dsc1 = 0, dsc2 = 17.857),
                                             sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)


#............................................................
#---- SWE1 #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.99,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(156.5,   267.5),
                                dsc2 =  c(1.5,     3.5))
#......................
# agebands
#......................
rawbase <- readRDS("data/derived/SWE1/SWE1_regions.RDS")
SWE_basic_mod <- make_noSeroRev_IFR_model_fit(num_mas = 1, maxMa = "ma1",
                                              groupvar = "region",  dat = rawbase,
                                              num_xs = 4, max_xveclist = list("name" = "x4", min = 219, init = 226, max = 233, dsc1 = 219, dsc2 = 233),
                                              num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 16.15, dsc1 = 0, dsc2 = 16.15),
                                              sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)


#............................................................
#---- NLD1 #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.85,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(171.5,   281.5),
                                dsc2 =  c(3.5,     1.5))

#......................
# agebands
#......................
rawage <- readRDS("data/derived/NLD1/NLD1_agebands.RDS")
NLD_age_mod <- make_noSeroRev_IFR_model_fit(num_mas = 5, maxMa = "ma5",
                                            groupvar = "ageband",  dat = rawage,
                                            num_xs = 4, max_xveclist = list("name" = "x4", min = 219, init = 226, max = 233, dsc1 = 219, dsc2 = 233),
                                            num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 16.67, dsc1 = 0, dsc2 = 16.67),
                                            sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)



#............................................................
#---- ITA1 #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.85,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(90.5,   95.5),
                                dsc2 =  c(10.5,     5.5))

#......................
# agebands
#......................
rawage <- readRDS("data/derived/ITA1/ITA1_agebands.RDS")
ITA_age_mod <- make_noSeroRev_IFR_model_fit(num_mas = 10, maxMa = "ma10",
                                            groupvar = "ageband",  dat = rawage,
                                            num_xs = 4, max_xveclist = list("name" = "x4", min = 192, init = 199, max = 206, dsc1 = 192, dsc2 = 206),
                                            num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 17.91, dsc1 = 0, dsc2 = 17.91),
                                            sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)


#............................................................
#---- LUX1 #----
#...........................................................
# using GBR numbers for sens bc LUX study sens so small
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.85,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(56.5,    181.5),
                                dsc2 =  c(19.5,     4.5))

#......................
# agebands
#......................
rawage <- readRDS("data/derived/LUX1/LUX1_agebands.RDS")
LUX_age_mod <- make_noSeroRev_IFR_model_fit(num_mas = 7, maxMa = "ma7",
                                            groupvar = "ageband",  dat = rawage,
                                            num_xs = 4, max_xveclist = list("name" = "x4", min = 219, init = 226, max = 233, dsc1 = 219, dsc2 = 233),
                                            num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 13.34, dsc1 = 0, dsc2 = 13.34),
                                            sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)

#............................................................
#---- Los Angeles County #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.73,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(27.5,    30.5),
                                dsc2 =  c(10.5,    0.5))

#......................
# agebands
#......................
rawbase <- readRDS("data/derived/USA/LA_CA1_regions.RDS")
LACA_basic_mod <- make_noSeroRev_IFR_model_fit(num_mas = 1, maxMa = "ma1",
                                               groupvar = "region",  dat = rawbase,
                                               num_xs = 4, max_xveclist = list("name" = "x4", min = 219, init = 226, max = 233, dsc1 = 219, dsc2 = 233),
                                               num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 16.12, dsc1 = 0, dsc2 = 16.12),
                                               sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)


#............................................................
#---- New York City #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.85,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(204.5,   92.5),
                                dsc2 =  c(30.5,    0.5))

#......................
# agebands
#......................
rawage <- readRDS("data/derived/USA/NYC_NY_1agebands.RDS")
NYC_age_mod <- make_noSeroRev_IFR_model_fit(num_mas = 5, maxMa = "ma5",
                                            groupvar = "ageband",  dat = rawage,
                                            num_xs = 4, max_xveclist = list("name" = "x4", min = 219, init = 226, max = 233, dsc1 = 219, dsc2 = 233),
                                            num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 15.94, dsc1 = 0, dsc2 = 15.94),
                                            sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)


#............................................................
#---- New York State #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.85,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(204.5,   92.5),
                                dsc2 =  c(30.5,    0.5))

#......................
# basic
#......................
rawbase <- readRDS("data/derived/USA/NYS1_regions.RDS")
NYS_basic_mod <- make_noSeroRev_IFR_model_fit(num_mas = 1, maxMa = "ma1",
                                              groupvar = "region",  dat = rawbase,
                                              num_xs = 4, max_xveclist = list("name" = "x4", min = 219, init = 226, max = 233, dsc1 = 219, dsc2 = 233),
                                              num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 15.94, dsc1 = 0, dsc2 = 15.94),
                                              sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)



#............................................................
#---- CHN1 #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.85,    0.97),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(91.5,    198.5),
                                dsc2 =  c(14.5,    4.5))

#......................
# agebands
#......................
rawage <- readRDS("data/derived/CHN1/CHN1_agebands.RDS")
CHN_age_mod <- make_noSeroRev_IFR_model_fit(num_mas = 9, maxMa = "ma9",
                                            groupvar = "ageband",  dat = rawage,
                                            num_xs = 4, max_xveclist = list("name" = "x4", min = 199, init = 206, max = 213, dsc1 = 199, dsc2 = 213),
                                            num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 16.20, dsc1 = 0, dsc2 = 16.20),
                                            sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)


#............................................................
#---- KEN1 #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.83,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(15.5,    901.5),
                                dsc2 =  c(3.4,     9.5))

#......................
# agebands
#......................
rawage <- readRDS("data/derived/KEN1/KEN1_agebands.RDS")
KEN1_age_mod <- make_noSeroRev_IFR_model_fit(num_mas = 7, maxMa = "ma7",
                                             groupvar = "ageband",  dat = rawage,
                                             num_xs = 4, max_xveclist = list("name" = "x4", min = 219, init = 226, max = 233, dsc1 = 219, dsc2 = 233),
                                             num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 17.11, dsc1 = 0, dsc2 = 17.11),
                                             sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)



#............................................................
#---- Come Together #----
#...........................................................
bvec <- seq(5, 2.5, length.out = 50)

fit_map <- tibble::tibble(
  name = c("BRA1_age",
           "BRA5_rgn",
           "CHE1_age",
           "CHE2_age",
           "DNK1_rgn", # basic but use region for splitting
           "ESP1-2_age",
           "GBR3_age",
           "SWE_rgn", # basic but use region for splitting
           "NLD1_age",
           "ITA1_age",
           "LUX1_age",
           "LA_CA1_rgn", # basic but use region for splitting
           "NYC_NY_1_age",
           "NYS_rgn", # basic but use region for splitting
           "CHN1_age",
           "KEN1_age"),
  modelobj = list(BRA1_age_mod,
                  BRA5_rgn_mod,
                  CHE1_age_mod,
                  CHE2_age_mod,
                  DNK_basic_mod,
                  ESP_age_mod,
                  GBR3_age_mod,
                  SWE_basic_mod,
                  NLD_age_mod,
                  ITA_age_mod,
                  LUX_age_mod,
                  LACA_basic_mod,
                  NYC_age_mod,
                  NYS_basic_mod,
                  CHN_age_mod,
                  KEN1_age_mod),
  rungs = 50,
  GTI_pow = list(bvec),
  burnin = 1e4,
  samples = 1e4,
  thinning = 10)

#......................
# manual adjustments to fit map
#......................
fit_map$GTI_pow[grepl("BRA1", fit_map$name)] <- list(seq(6, 4, length.out = 50))

#......................
# fitmap out
#......................
# select what we need for fits and make outpaths
dir.create("data/param_map/Modfits_noserorev/", recursive = T)
lapply(split(fit_map, 1:nrow(fit_map)), function(x){
  saveRDS(x, paste0("data/param_map/Modfits_noserorev/",
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

  if (grepl("ITA", basename(path))) {
    # logit case for ITA with multiple age groups
    fit <- COVIDCurve::run_IFRmodel_age(IFRmodel = mod$modelobj[[1]],
                                        reparamIFR = TRUE,
                                        reparamInfxn = TRUE,
                                        reparamKnots = TRUE,
                                        binomial_likelihood = FALSE,
                                        chains = n_chains,
                                        burnin = mod$burnin,
                                        samples = mod$samples,
                                        rungs = mod$rungs,
                                        GTI_pow = mod$GTI_pow[[1]],
                                        cluster = cl,
                                        thinning = 10)
  } else if (grepl("SWE|DNK", basename(path))) {
    # logit case where only one age group
    fit <- COVIDCurve::run_IFRmodel_age(IFRmodel = mod$modelobj[[1]],
                                        reparamIFR = FALSE,
                                        reparamInfxn = TRUE,
                                        reparamKnots = TRUE,
                                        binomial_likelihood = FALSE,
                                        chains = n_chains,
                                        burnin = mod$burnin,
                                        samples = mod$samples,
                                        rungs = mod$rungs,
                                        GTI_pow = mod$GTI_pow[[1]],
                                        cluster = cl,
                                        thinning = 10)

  } else if (nrow(mod$modelobj[[1]]$IFRdictkey == 1)) {
    # basic case
    fit <- COVIDCurve::run_IFRmodel_age(IFRmodel = mod$modelobj[[1]],
                                        reparamIFR = FALSE,
                                        reparamInfxn = TRUE,
                                        reparamKnots = TRUE,
                                        binomial_likelihood = TRUE,
                                        chains = n_chains,
                                        burnin = mod$burnin,
                                        samples = mod$samples,
                                        rungs = mod$rungs,
                                        GTI_pow = mod$GTI_pow[[1]],
                                        cluster = cl,
                                        thinning = 10)
  } else {
    # normal binomial case
    fit <- COVIDCurve::run_IFRmodel_age(IFRmodel = mod$modelobj[[1]],
                                        reparamIFR = TRUE,
                                        reparamInfxn = TRUE,
                                        reparamKnots = TRUE,
                                        binomial_likelihood = TRUE,
                                        chains = n_chains,
                                        burnin = mod$burnin,
                                        samples = mod$samples,
                                        rungs = mod$rungs,
                                        GTI_pow = mod$GTI_pow[[1]],
                                        cluster = cl,
                                        thinning = 10)
  }
  parallel::stopCluster(cl)
  gc()

  # out
  dir.create("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/Modfits_serorev/", recursive = TRUE)
  outpath = paste0("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/Modfits_serorev/",
                   mod$name, "_rung", mod$rungs, "_burn", mod$burnin, "_smpl", mod$samples, "_SeroRev.RDS")
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
file_param_map <- list.files(path = "data/param_map/Modfits_noserorev/",
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
        clustermq.template = "drake_clst/slurm_clustermq_LL.tmpl")
make(plan, parallelism = "clustermq", jobs = nrow(file_param_map),
     log_make = "Modfits_drake_noserorev.log", verbose = 2,
     log_progress = TRUE,
     log_build_times = FALSE,
     recoverable = FALSE,
     history = FALSE,
     session_info = FALSE,
     lock_envir = FALSE, # unlock environment so parallel::clusterApplyLB in drjacoby can work
     lock_cache = FALSE)



cat("************** Drake Finished **************************")




