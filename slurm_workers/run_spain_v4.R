####################################################################################
## Purpose: This is a script that will wrap a DRAKE make process to
##          run multiple ESP data iterations to fine to the MCoupling
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
make_IFR_model_spain <- function(num_mas, maxMa, groupvar, dat) {
  # make dfs
  ifr_paramsdf <- make_ma_reparamdf(num_mas = num_mas)
  knot_paramsdf <- make_splinex_reparamdf(max_xvec = list("name" = "x4", min = 125, init = 131, max = 137, dsc1 = 125, dsc2 = 137),
                                          num_xs = 4)

  infxn_paramsdf <- make_spliney_reparamdf(max_yvec = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                           num_ys = 5)

  sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day1"),
                                  min =   c(0.83,     0.50,   10,         117),
                                  init =  c(0.85,     0.99,   10,         125),
                                  max =   c(0.87,     1.00,   10,         131),
                                  dsc1 =  c(850,      10,     5,         117),
                                  dsc2 =  c(150,      3,      15,        131))
  noise_paramsdf <- make_noiseeff_reparamdf(num_Nes = num_mas, min = 0, init = 5, max = 10)

  # bring together
  df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, knot_paramsdf, sens_spec_tbl, noise_paramsdf)
  #......................
  # format data
  #......................
  dictkey <- tibble::tibble(groupvar = unlist(dat$seroprev_group[, groupvar]), "Strata" = paste0("ma", 1:num_mas))
  colnames(dictkey) <- c(paste(groupvar), "Strata")
  # deaths
  dat$deaths <- dplyr::left_join(dat$deaths, dictkey) %>%
    dplyr::select(c("ObsDay", "Strata", "Deaths"))
  # seroprev
  dat$obs_serology <- dplyr::left_join(dat$seroprev_group, dictkey) %>%
    dplyr::mutate(SeroDay = "sero_day1") %>%
    dplyr::rename(SeroPrev = seroprevalence) %>%
    dplyr::select(c("SeroDay", "Strata", "SeroPrev"))

  inputdata <- list(obs_deaths = dat$deaths,
                    obs_serology = dat$obs_serology)

  demog <- dat$prop_pop %>%
    dplyr::left_join(., dictkey) %>%
    dplyr::select(c("Strata", "popN")) %>%
    dplyr::mutate(popN = round(popN))


  # make mod
  mod1 <- make_IFRmodel_agg$new()
  mod1$set_MeanOnset(18.8)
  mod1$set_CoefVarOnset(0.45)
  mod1$set_level("Time-Series")
  mod1$set_IFRparams(paste0("ma", 1:num_mas))
  mod1$set_maxMa(maxMa)
  mod1$set_Knotparams(paste0("x", 1:4))
  mod1$set_relKnot("x4")
  mod1$set_Infxnparams(paste0("y", 1:5))
  mod1$set_relInfxn("y5")
  mod1$set_Serotestparams(c("sens", "spec", "sero_rate"))
  mod1$set_Serodayparams(c("sero_day1"))
  mod1$set_Noiseparams(paste0("Ne", 1:num_mas))
  mod1$set_data(inputdata)
  mod1$set_demog(demog)
  mod1$set_paramdf(df_params)
  mod1$set_rho(rep(1, num_mas))
  mod1$set_rcensor_day(.Machine$integer.max)
  # out
  mod1
}

# read raw data and make param maps
rawage <- readRDS("data/derived/ESP/ESP_agebands.RDS")
rawrgn <- readRDS("data/derived/ESP/ESP_regions.RDS")
agemap <- tibble::as_tibble(expand.grid(rungs = c(10, 25, 50),
                                        GTI_pow = c(2, 2.5, 3, 3.5, 4.0, 4.5, 5.0, 5.5, 6),
                                        burnin = 1e3,
                                        samples = 1e3)) %>%
  dplyr::mutate(lvl = "age",
                num_mas = 10,
                maxMa = "ma10",
                groupvar = "ageband",
                dat = list(rawage))

rgnmap <- tibble::as_tibble(expand.grid(rungs = c(10, 25, 50),
                                        GTI_pow = c(2, 2.5, 3, 3.5, 4.0, 4.5, 5.0, 5.5, 6),
                                        burnin = 1e3,
                                        samples = 1e3)) %>%
  dplyr::mutate(lvl = "region",
                num_mas = 17,
                maxMa = "ma14",
                groupvar = "region",
                dat = list(rawrgn))

param_map <- dplyr::bind_rows(agemap, rgnmap)
param_map$modelobj <- purrr::pmap(param_map[, c("num_mas", "maxMa", "groupvar", "dat")], make_IFR_model_spain)
# select what we need for fits and make outpaths
dir.create("data/param_map")
param_map.fit <- param_map %>%
  dplyr::select(c("lvl", "modelobj", "rungs", "GTI_pow", "burnin", "samples"))
lapply(split(param_map.fit, 1:nrow(param_map.fit)), function(x){
  saveRDS(x, paste0("data/param_map/",
                    x$lvl, "_GTI", x$GTI_pow, "_rung", x$rungs, "_burn", x$burnin, "_smpl", x$samples, ".RDS"))
})

#......................
# for snake
#......................
param_map_snake <- gsub(".RDS", "", list.files(path = "data/param_map/", pattern = ".RDS"))
param_map_snake <- data.frame(parampath = param_map_snake)
colnames(param_map_snake) <- c("#parampath")
write.table(x = param_map_snake,
            file = "data/param_map/snake_map.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)


#............................................................
# MCMC Object
#...........................................................
run_MCMC <- function(path) {
  mod <- readRDS(path)
  #......................
  # make cluster object to parallelize chains
  #......................
  n_chains <- 10
  cl <- parallel::makeCluster(n_chains)
  fit <- COVIDCurve::run_IFRmodel_agg(IFRmodel = mod$modelobj[[1]],
                                      reparamIFR = TRUE,
                                      reparamInfxn = TRUE,
                                      reparamKnots = TRUE,
                                      chains = n_chains,
                                      burnin = mod$burnin,
                                      samples = mod$samples,
                                      rungs = mod$rungs,
                                      GTI_pow = mod$GTI_pow,
                                      cluster = cl,
                                      silent = FALSE)
  mc_accept_mean <- mean(fit$mcmcout$diagnostics$mc_accept$value)
  mc_accept_min <- min(fit$mcmcout$diagnostics$mc_accept$value)
  # out
  dir.create("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/", recursive = TRUE)
  outpath = paste0("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/", lvl, "_GTI", GTI_pow, "_rung", rungs, "_burn", burnin, "_smpl", samples, ".RDS")
  out <- list(fit = fit,
              mc_accept_mean = mc_accept_mean,
              mc_accept_min = mc_accept_min)

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
file_param_map <- list.files(path = "data/param_map/",
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
        clustermq.template = "slurm_workers/slurm_clustermq_LL.tmpl")
make(plan, parallelism = "clustermq", jobs = nrow(file_param_map),
     console_log_file = "drake.log", verbose = 2,
     log_progress = FALSE,
     log_build_times = FALSE,
     recoverable = FALSE,
     history = FALSE,
     session_info = FALSE,
     lock_envir = FALSE)






