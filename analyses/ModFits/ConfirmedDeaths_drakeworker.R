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
                               min  = c(17,     0,     16),
                               init = c(19,     0.79,  18),
                               max =  c(22,     1,     21),
                               dsc1 = c(19.26,  2370,  18.3),
                               dsc2 = c(0.5,    630,   0.5))
#......................
# seroreversion weibull scale/shape for various assay
# (for later rbinds)
#......................
empty <- tibble::tibble(name = c("sero_rev_shape", "sero_rev_scale"),
                        min  = c(NA,                 NA),
                        init = c(NA,                 NA),
                        max =  c(NA,                 NA),
                        dsc1 = c(NA,                 NA),
                        dsc2 = c(NA,                 NA))


#...........................................................
#---- GBR3 #-----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.85,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(38.5,   493.5),
                                dsc2 =  c(10.5,    7.5))
sens_spec_tbl_noserorev <- rbind(sens_spec_tbl, empty)


#......................
# agebands
#......................
rawage <- readRDS("data/derived/confirmeddeaths/GBR3_confirmed_deaths.rds")
GBR3_confirmeddeaths_mod <- make_IFR_model_fit(num_mas = 4, maxMa = "ma4",
                                         groupvar = "ageband",  dat = rawage,
                                         num_xs = 4, max_xveclist = list("name" = "x4", min = 206, init = 210, max = 213, dsc1 = 199, dsc2 = 213),
                                         num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 17.857, dsc1 = 0, dsc2 = 17.857),
                                         sens_spec_tbl = sens_spec_tbl_noserorev, tod_paramsdf = tod_paramsdf)



#............................................................
#---- New York City #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.85,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(204.5,   92.5),
                                dsc2 =  c(30.5,    0.5))
sens_spec_tbl_noserorev <- rbind(sens_spec_tbl, empty)

#......................
# agebands
#......................
rawage <- readRDS("data/derived/confirmeddeaths/NYC_confirmed_deaths.rds")
NYC_confirmeddeaths_mod <- make_IFR_model_fit(num_mas = 5, maxMa = "ma5",
                                        groupvar = "ageband",  dat = rawage,
                                        num_xs = 4, max_xveclist = list("name" = "x4", min = 219, init = 226, max = 233, dsc1 = 219, dsc2 = 233),
                                        num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 15.94, dsc1 = 0, dsc2 = 15.94),
                                        sens_spec_tbl = sens_spec_tbl_noserorev, tod_paramsdf = tod_paramsdf)


#............................................................
#---- Come Together #----
#...........................................................
bvec <- seq(5, 2.5, length.out = 50)

fit_map <- tibble::tibble(
  name = c("GBR3_confirmeddeaths",
           "NYC_NY_1_confirmeddeaths"),
  modelobj = list(GBR3_confirmeddeaths_mod,
                  NYC_confirmeddeaths_mod),
  rungs = 50,
  GTI_pow = list(bvec),
  burnin = 1e4,
  samples = 1e4,
  thinning = 10)

#......................
# fitmap out
#......................
# select what we need for fits and make outpaths
dir.create("data/param_map/Modfits_ConfirmedDeaths/", recursive = T)
lapply(split(fit_map, 1:nrow(fit_map)), function(x){
  saveRDS(x, paste0("data/param_map/Modfits_ConfirmedDeaths/",
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

  fit <- COVIDCurve::run_IFRmodel_age(IFRmodel = mod$modelobj[[1]],
                                      reparamIFR = TRUE,
                                      reparamInfxn = TRUE,
                                      reparamKnots = TRUE,
                                      reparamDelays = FALSE,
                                      reparamNe = FALSE,
                                      binomial_likelihood = TRUE,
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
  dir.create("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/Modfits_ConfirmedDeaths/", recursive = TRUE)
  outpath = paste0("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/Modfits_ConfirmedDeaths/",
                   mod$name, "_rung", mod$rungs, "_burn", mod$burnin, "_smpl", mod$samples, "_ConfirmedDeaths.RDS")
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
file_param_map <- list.files(path = "data/param_map/Modfits_ConfirmedDeaths/",
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
     log_make = "Modfits_drake_confirmeddeaths.log", verbose = 2,
     log_progress = TRUE,
     log_build_times = FALSE,
     recoverable = FALSE,
     history = FALSE,
     session_info = FALSE,
     lock_envir = FALSE, # unlock environment so parallel::clusterApplyLB in drjacoby can work
     lock_cache = FALSE)



cat("************** Drake Finished **************************")




