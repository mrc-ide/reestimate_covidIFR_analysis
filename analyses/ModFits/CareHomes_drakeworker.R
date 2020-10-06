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
                               init = c(19,     0.85,  18),
                               max =  c(20,     1,     21),
                               dsc1 = c(19.8,   2550,  18.3),
                               dsc2 = c(0.1,    450,   0.1))



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
# carehomes
#......................
rawch <- readRDS("data/derived/carehomes/CHE1_agebands_noCH.RDS")
CHE1_carehomes_mod <- make_noSeroRev_IFR_model_fit(num_mas = 7, maxMa = "ma7",
                                                   groupvar = "ageband",  dat = rawch,
                                                   num_xs = 4, max_xveclist = list("name" = "x4", min = 216, init = 223, max = 230, dsc1 = 216, dsc2 = 230),
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
# carehomes
#......................
rawch <- readRDS("data/derived/carehomes/CHE2_agebands_noCH.RDS")
CHE2_carehomes_mod <- make_noSeroRev_IFR_model_fit(num_mas = 7, maxMa = "ma7",
                                                   groupvar = "ageband",  dat = rawch,
                                                   num_xs = 4, max_xveclist = list("name" = "x4", min = 216, init = 223, max = 230, dsc1 = 216, dsc2 = 230),
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
# get fits from stan model
#......................
dnksens <- readr::read_csv("results/prior_inputs/DNK1_sens_reg_age.csv")
sens <- fitdistrplus::fitdist(unlist(dnksens), distr = "beta", method = "mme")
dnkspec <- readr::read_csv("results/prior_inputs/DNK1_spec_reg_age.csv")
spec <- fitdistrplus::fitdist(unlist(dnkspec), distr = "beta", method = "mme")
sens_spec_tbl$dsc1[sens_spec_tbl$name == "sens"] <- sens$estimate[["shape1"]]
sens_spec_tbl$dsc2[sens_spec_tbl$name == "sens"] <- sens$estimate[["shape2"]]
sens_spec_tbl$dsc1[sens_spec_tbl$name == "spec"] <- spec$estimate[["shape1"]]
sens_spec_tbl$dsc2[sens_spec_tbl$name == "spec"] <- spec$estimate[["shape2"]]

#......................
# agebands
#......................
rawch <- readRDS("data/derived/carehomes/DNK1_agebands_noCH.RDS")
DNK_carehomes_mod <- make_noSeroRev_IFR_model_fit(num_mas = 2, maxMa = "ma2",
                                                  groupvar = "ageband",  dat = rawch,
                                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 216, init = 223, max = 230, dsc1 = 216, dsc2 = 230),
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
#......................
# get fits from stan model
#......................
espsens <- readr::read_csv("results/prior_inputs/ESP1-2_sens_reg_age.csv")
sens <- fitdistrplus::fitdist(unlist(espsens), distr = "beta", method = "mme")
espspec <- readr::read_csv("results/prior_inputs/ESP1-2_spec_reg_age.csv")
spec <- fitdistrplus::fitdist(unlist(espspec), distr = "beta", method = "mme")
sens_spec_tbl$dsc1[sens_spec_tbl$name == "sens"] <- sens$estimate[["shape1"]]
sens_spec_tbl$dsc2[sens_spec_tbl$name == "sens"] <- sens$estimate[["shape2"]]
sens_spec_tbl$dsc1[sens_spec_tbl$name == "spec"] <- spec$estimate[["shape1"]]
sens_spec_tbl$dsc2[sens_spec_tbl$name == "spec"] <- spec$estimate[["shape2"]]


#......................
# agebands
#......................
rawch <- readRDS("data/derived/carehomes/ESP1-2_agebands_noCH.RDS")
ESP_carehomes_mod <- make_noSeroRev_IFR_model_fit(num_mas = 7, maxMa = "ma7",
                                                  groupvar = "ageband",  dat = rawch,
                                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 216, init = 223, max = 230, dsc1 = 216, dsc2 = 230),
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
# get fits from stan model
#......................
gbrsens <- readr::read_csv("results/prior_inputs/GBR3_sens_reg_age.csv")
sens <- fitdistrplus::fitdist(unlist(gbrsens), distr = "beta", method = "mme")
gbrspec <- readr::read_csv("results/prior_inputs/GBR3_spec_reg_age.csv")
spec <- fitdistrplus::fitdist(unlist(gbrspec), distr = "beta", method = "mme")
sens_spec_tbl$dsc1[sens_spec_tbl$name == "sens"] <- sens$estimate[["shape1"]]
sens_spec_tbl$dsc2[sens_spec_tbl$name == "sens"] <- sens$estimate[["shape2"]]
sens_spec_tbl$dsc1[sens_spec_tbl$name == "spec"] <- spec$estimate[["shape1"]]
sens_spec_tbl$dsc2[sens_spec_tbl$name == "spec"] <- spec$estimate[["shape2"]]


#......................
# agebands
#......................
rawch <- readRDS("data/derived/carehomes/GBR3_agebands_noCH.RDS")
GBR3_carehomes_mod <- make_noSeroRev_IFR_model_fit(num_mas = 3, maxMa = "ma3",
                                                   groupvar = "ageband",  dat = rawch,
                                                   num_xs = 4, max_xveclist = list("name" = "x4", min = 216, init = 223, max = 230, dsc1 = 216, dsc2 = 230),
                                                   num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 17.857, dsc1 = 0, dsc2 = 17.857),
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
# get fits from stan model
#......................
nyssens <- readr::read_csv("results/prior_inputs/NYS_sens_reg_age.csv")
sens <- fitdistrplus::fitdist(unlist(nyssens), distr = "beta", method = "mme")
nysspec <- readr::read_csv("results/prior_inputs/NYS_spec_reg_age.csv")
spec <- fitdistrplus::fitdist(unlist(nysspec), distr = "beta", method = "mme")
sens_spec_tbl$dsc1[sens_spec_tbl$name == "sens"] <- sens$estimate[["shape1"]]
sens_spec_tbl$dsc2[sens_spec_tbl$name == "sens"] <- sens$estimate[["shape2"]]
sens_spec_tbl$dsc1[sens_spec_tbl$name == "spec"] <- spec$estimate[["shape1"]]
sens_spec_tbl$dsc2[sens_spec_tbl$name == "spec"] <- spec$estimate[["shape2"]]


#......................
# agebands
#......................
rawch <- readRDS("data/derived/carehomes/NYS1_agebands_noCH.RDS")
NYS_carehomes_mod <- make_noSeroRev_IFR_model_fit(num_mas = 7, maxMa = "ma7",
                                                  groupvar = "ageband",  dat = rawch,
                                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 216, init = 223, max = 230, dsc1 = 216, dsc2 = 230),
                                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 16.77, dsc1 = 0, dsc2 = 16.77),
                                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)


#............................................................
#---- Come Together #----
#...........................................................
bvec <- seq(5, 2.5, length.out = 50)

fit_map <- tibble::tibble(
  name = c("CHE1_carehomes",
           "CHE2_carehomes",
           "DNK1_carehomes",
           "ESP1-2_carehomes",
           "GBR3_carehomes",
           "NYS1_carehomes"),
  modelobj = list(CHE1_carehomes_mod,
                  CHE2_carehomes_mod,
                  DNK_carehomes_mod,
                  ESP_carehomes_mod,
                  GBR3_carehomes_mod,
                  NYS_carehomes_mod),
  rungs = 50,
  GTI_pow = list(bvec),
  burnin = 1e4,
  samples = 1e4,
  thinning = 10)

#......................
# fitmap out
#......................
# select what we need for fits and make outpaths
dir.create("data/param_map/Modfits_CareHomes/", recursive = T)
lapply(split(fit_map, 1:nrow(fit_map)), function(x){
  saveRDS(x, paste0("data/param_map/Modfits_CareHomes/",
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

  if (grepl("DNK", basename(path))) {
    # logit case for DNK with multiple age groups
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
  dir.create("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/Modfits_carehomes/", recursive = TRUE)
  outpath = paste0("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/Modfits_carehomes/",
                   mod$name, "_rung", mod$rungs, "_burn", mod$burnin, "_smpl", mod$samples, "_carehomes.RDS")
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
file_param_map <- list.files(path = "data/param_map/Modfits_CareHomes/",
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
     log_make = "Modfits_drake_carehomes.log", verbose = 2,
     log_progress = TRUE,
     log_build_times = FALSE,
     recoverable = FALSE,
     history = FALSE,
     session_info = FALSE,
     lock_envir = FALSE, # unlock environment so parallel::clusterApplyLB in drjacoby can work
     lock_cache = FALSE)



cat("************** Drake Finished **************************")

