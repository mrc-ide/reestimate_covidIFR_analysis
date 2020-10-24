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
#---- time series CHN1 reg spec, reg Ma #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.85,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(91.5,    26.5),
                                dsc2 =  c(14.5,    0.5))

#......................
# agebands
#......................
rawage <- readRDS("data/derived/CHN1/CHN1_agebands.RDS")
timeseriesiggCHN_age_mod <- make_noSeroRev_IFR_model_fit(num_mas = 9, maxMa = "ma9",
                                            groupvar = "ageband",  dat = rawage,
                                            num_xs = 4, max_xveclist = list("name" = "x4", min = 247, init = 254, max = 261, dsc1 = 247, dsc2 = 261),
                                            num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 16.20, dsc1 = 0, dsc2 = 16.20),
                                            sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)


#............................................................
#---- time series CHN1 reg spec, higher Ma #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.85,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(91.5,    26.5),
                                dsc2 =  c(14.5,    0.5))

#......................
# agebands
#......................
rawage <- readRDS("data/derived/CHN1/CHN1_agebands.RDS")
timeseriesiggCHN_age_mod_higherma <- make_noSeroRev_IFR_model_fit(num_mas = 9, maxMa = "ma9",
                                                     groupvar = "ageband",  dat = rawage,
                                                     num_xs = 4, max_xveclist = list("name" = "x4", min = 247, init = 254, max = 261, dsc1 = 247, dsc2 = 261),
                                                     num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 16.20, dsc1 = 0, dsc2 = 16.20),
                                                     sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                                     upperMa = 0.6)



#............................................................
#---- time series CHN1 higher spec, reg Ma #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.85,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(91.5,    260.5),
                                dsc2 =  c(14.5,    0.5))

#......................
# agebands
#......................
rawage <- readRDS("data/derived/CHN1/CHN1_agebands.RDS")
timeseriesiggCHN_age_mod_higherspec <- make_noSeroRev_IFR_model_fit(num_mas = 9, maxMa = "ma9",
                                            groupvar = "ageband",  dat = rawage,
                                            num_xs = 4, max_xveclist = list("name" = "x4", min = 247, init = 254, max = 261, dsc1 = 247, dsc2 = 261),
                                            num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 16.20, dsc1 = 0, dsc2 = 16.20),
                                            sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)


#............................................................
#---- igm/igg CHN1 reg spec, reg Ma #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.85,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(94.5,    200.5),
                                dsc2 =  c(11.5,    2.5))

#......................
# agebands
#......................
rawage <- readRDS("gti_tune_wuhan/igm_igg_CHN1_agebands.RDS")
igmiggCHN_age_mod <- make_noSeroRev_IFR_model_fit(num_mas = 9, maxMa = "ma9",
                                                  groupvar = "ageband",  dat = rawage,
                                                  num_xs = 4, max_xveclist = list("name" = "x4", min = 247, init = 254, max = 261, dsc1 = 247, dsc2 = 261),
                                                  num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 16.20, dsc1 = 0, dsc2 = 16.20),
                                                  sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)


#............................................................
#---- igm/igg CHN1 reg spec, higher Ma #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.85,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(94.5,    26.5),
                                dsc2 =  c(11.5,    0.5))

#......................
# agebands
#......................
rawage <- readRDS("data/derived/CHN1/CHN1_agebands.RDS")
igmiggCHN_age_mod_higherma <- make_noSeroRev_IFR_model_fit(num_mas = 9, maxMa = "ma9",
                                                           groupvar = "ageband",  dat = rawage,
                                                           num_xs = 4, max_xveclist = list("name" = "x4", min = 247, init = 254, max = 261, dsc1 = 247, dsc2 = 261),
                                                           num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 16.20, dsc1 = 0, dsc2 = 16.20),
                                                           sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf,
                                                           upperMa = 0.6)



#............................................................
#---- igg, igm CHN1 higher spec, reg Ma #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.85,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(94.5,    400.5),
                                dsc2 =  c(11.5,    4.5))

#......................
# agebands
#......................
rawage <- readRDS("gti_tune_wuhan/igm_igg_CHN1_agebands.RDS")
igmiggCHN_age_mod_higherspec <- make_noSeroRev_IFR_model_fit(num_mas = 9, maxMa = "ma9",
                                                             groupvar = "ageband",  dat = rawage,
                                                             num_xs = 4, max_xveclist = list("name" = "x4", min = 247, init = 254, max = 261, dsc1 = 247, dsc2 = 261),
                                                             num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 16.20, dsc1 = 0, dsc2 = 16.20),
                                                             sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)




#............................................................
#---- Come Together #----
#...........................................................
bvec1 <- seq(5, 3, length.out = 50)
bvec2 <- seq(5, 3.5, length.out = 50)
bvec3 <- seq(5, 4, length.out = 50)
bvec4 <- seq(5, 4.5, length.out = 50)
bveclist <- list(bvec1, bvec2, bvec3, bvec4,
                 bvec1, bvec2, bvec3, bvec4,
                 bvec1, bvec2, bvec3, bvec4,
                 bvec1, bvec2, bvec3, bvec4,
                 bvec1, bvec2, bvec3, bvec4,
                 bvec1, bvec2, bvec3, bvec4)

fit_map <- tibble::tibble(
  name = c("timeseriesiggCHN1_reg_bvec1",
           "timeseriesiggCHN1_reg_bvec2",
           "timeseriesiggCHN1_reg_bvec3",
           "timeseriesiggCHN1_reg_bvec4",
           "timeseriesiggCHN1_higherma_bvec1",
           "timeseriesiggCHN1_higherma_bvec2",
           "timeseriesiggCHN1_higherma_bvec3",
           "timeseriesiggCHN1_higherma_bvec4",
           "timeseriesiggCHN1_higherspec_bvec1",
           "timeseriesiggCHN1_higherspec_bvec2",
           "timeseriesiggCHN1_higherspec_bvec3",
           "timeseriesiggCHN1_higherspec_bvec4",
           "igmiggCHN1_reg_bvec1",
           "igmiggCHN1_reg_bvec2",
           "igmiggCHN1_reg_bvec3",
           "igmiggCHN1_reg_bvec4",
           "igmiggCHN1_higherma_bvec1",
           "igmiggCHN1_higherma_bvec2",
           "igmiggCHN1_higherma_bvec3",
           "igmiggCHN1_higherma_bvec4",
           "igmiggCHN1_higherspec_bvec1",
           "igmiggCHN1_higherspec_bvec2",
           "igmiggCHN1_higherspec_bvec3",
           "igmiggCHN1_higherspec_bvec4"),
  modelobj = list(timeseriesiggCHN_age_mod,
                  timeseriesiggCHN_age_mod,
                  timeseriesiggCHN_age_mod,
                  timeseriesiggCHN_age_mod,
                  timeseriesiggCHN_age_mod_higherma,
                  timeseriesiggCHN_age_mod_higherma,
                  timeseriesiggCHN_age_mod_higherma,
                  timeseriesiggCHN_age_mod_higherma,
                  timeseriesiggCHN_age_mod_higherspec,
                  timeseriesiggCHN_age_mod_higherspec,
                  timeseriesiggCHN_age_mod_higherspec,
                  timeseriesiggCHN_age_mod_higherspec,
                  igmiggCHN_age_mod,
                  igmiggCHN_age_mod,
                  igmiggCHN_age_mod,
                  igmiggCHN_age_mod,
                  igmiggCHN_age_mod_higherma,
                  igmiggCHN_age_mod_higherma,
                  igmiggCHN_age_mod_higherma,
                  igmiggCHN_age_mod_higherma,
                  igmiggCHN_age_mod_higherspec,
                  igmiggCHN_age_mod_higherspec,
                  igmiggCHN_age_mod_higherspec,
                  igmiggCHN_age_mod_higherspec),
  rungs = 50,
  GTI_pow = bveclist,
  burnin = 1e4,
  samples = 1e4,
  thinning = 10)

#......................
# fitmap out
#......................
# select what we need for fits and make outpaths
dir.create("data/param_map/tune_wuhan/", recursive = T)
lapply(split(fit_map, 1:nrow(fit_map)), function(x){
  saveRDS(x, paste0("data/param_map/tune_wuhan/",
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

  if (grepl("ITA|DNK|SWE", basename(path))) {
    # logit cases
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
                                        thinning = mod$thinning)

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
                                        thinning = mod$thinning)
  }
  parallel::stopCluster(cl)
  gc()

  # out
  dir.create("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/tune_wuhan/", recursive = TRUE)
  outpath = paste0("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/tune_wuhan/",
                   mod$name, "_rung", mod$rungs, "_burn", mod$burnin, "_smpl", mod$samples, "_tune_wuhan.RDS")
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
file_param_map <- list.files(path = "data/param_map/tune_wuhan/",
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
     log_make = "Modfits_drake_tune_wuhan.log", verbose = 2,
     log_progress = TRUE,
     log_build_times = FALSE,
     recoverable = FALSE,
     history = FALSE,
     session_info = FALSE,
     lock_envir = FALSE, # unlock environment so parallel::clusterApplyLB in drjacoby can work
     lock_cache = FALSE)



cat("************** Drake Finished **************************")




