####################################################################################
## Purpose: Run out simulations
##
## Notes: Assumes SLURM cluster
####################################################################################
setwd("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis")
library(drake)
library(tidyverse)
library(COVIDCurve)
source("R/covidcurve_helper_functions.R")
source("R/simple_seir_model.R")
set.seed(48)


#............................................................
# Run Various Scenarios for Incidence Curves
#...........................................................
nsims <- 100
popN <- 3e6 # pop slightly smaller, so infxns don't "run" away
#......................
# exponential growth
#.......... ............
expgrowth <- lapply(1:nsims, function(x){
  run_simple_seir(N = popN,
                  E0 = 50,
                  R0 = 0,
                  betas = 0.3,
                  beta_changes = 1,
                  sigma = 0.2,
                  gamma = 0.2,
                  time = 200)})
expgrowth <- expgrowth %>%
  dplyr::bind_rows(.) %>%
  dplyr::group_by(step) %>%
  dplyr::summarise(
    infxns = mean(I)
  ) %>%
  dplyr::mutate_if(is.numeric, round, 0) %>%
  dplyr::rename(time = step) %>%
  dplyr::mutate(nm = "expgrowth")


#......................
# exponential growth with intervetions
#......................
intervene <- lapply(1:nsims, function(x){
  run_simple_seir(N = popN,
                  E0 = 50,
                  R0 = 0,
                  betas = c(0.33, 0.13, 0.12, 0.11),
                  beta_changes = c(1, 130, 140, 150),
                  sigma = 0.2,
                  gamma = 0.2,
                  time = 200)
})
intervene <- intervene %>%
  dplyr::bind_rows(.) %>%
  dplyr::group_by(step) %>%
  dplyr::summarise(
    infxns = mean(I)
  ) %>%
  dplyr::mutate_if(is.numeric, round, 0) %>%
  dplyr::rename(time = step) %>%
  dplyr::mutate(nm = "intervene")



#......................
# exponential growth with implementations and then second wave
#......................
secondwave <- lapply(1:nsims, function(x){
  run_simple_seir(N = popN,
                  E0 = 50,
                  R0 = 0,
                  betas = c(0.37, 0.13, 0.12, 0.18, 0.2, 0.21),
                  beta_changes = c(1, 100, 110, 120, 130, 140),
                  sigma = 0.2,
                  gamma = 0.2,
                  time = 200)
})
secondwave <- secondwave %>%
  dplyr::bind_rows(.) %>%
  dplyr::group_by(step) %>%
  dplyr::summarise(
    infxns = mean(I)
  ) %>%
  dplyr::mutate_if(is.numeric, round, 0) %>%
  dplyr::rename(time = step) %>%
  dplyr::mutate(nm = "secondwave")



#............................................................
# setup fatality data
#............................................................
# make up fatality data
fatalitydata <- tibble::tibble(Strata = c("ma1", "ma2", "ma3", "ma4", "ma5"),
                               IFR = c(1e-3, 0.15, 0.2, 0.299, 0.6),
                               Rho = 1,
                               Ne = 1)
demog <- tibble::tibble(Strata = c("ma1", "ma2", "ma3", "ma4", "ma5"),
                        popN = c(5e5, 5e5, 5e5, 2250000, 1250000))


#............................................................
# Simulate Under Model
#...........................................................
map <- expand.grid(curve = list(expgrowth, intervene, secondwave),
                   sens = c(0.85, 0.90),
                   spec = c(0.95, 0.99))
map <- tibble::as_tibble(map) %>%
  dplyr::mutate(fatalitydata = list(fatalitydata),
                demog = list(demog))
#......................
# rung covidcurve simulator
#......................
wrap_sim <- function(curve, sens, spec, fatalitydata, demog, sero_day) {

  dat <- COVIDCurve::Aggsim_infxn_2_death(
    fatalitydata = fatalitydata,
    m_od = 14.26,
    s_od = 0.79,
    curr_day = 200,
    infections = curve$infxns,
    simulate_seroprevalence = TRUE,
    sens = sens,
    spec = spec,
    sero_delay_rate = 13.3,
    demog = demog)

  # sero tidy up
  obs_serology <- dat$AggSeroPrev %>%
    dplyr::group_by(Strata) %>%
    dplyr::filter(event_obs_day == sero_day) %>%
    dplyr::rename(
      SeroDay = event_obs_day,
      SeroPrev = ObsPrev) %>%
    dplyr::select(c("SeroDay", "Strata", "SeroPrev")) %>%
    dplyr::mutate(SeroDay = "sero_day")

  # death dat
  dat$AggDeath <- dat$AggDeath %>%
    dplyr::mutate(Strata = as.character(Strata))

  # make out
  inputdata <- list(obs_deaths = dat$AggDeath,
                    obs_serology = obs_serology)
  out <- list(simdat = dat,
              inputdata = inputdata)
  return(out)

}

# run simdat and extract results into separate pieces
map$simdat <- purrr::pmap(map, wrap_sim, sero_day = 150)
map$inputdata <- purrr::map(map$simdat, "inputdata")
map$simdat <- purrr::map(map$simdat, "simdat", sero_day = 150)

#......................
# make IFR model
#......................
# sens/spec
get_sens_spec <- function(sens, spec) {
  tibble::tibble(name =  c("sens",          "spec",         "sero_rate",  "sero_day"),
                 min =   c(0.5,              0.5,            0,            140),
                 init =  c(0.8,              0.8,            1,            150),
                 max =   c(1,                1,              2,            160),
                 dsc1 =  c(sens*1e3,        spec*1e3,        0,            140),
                 dsc2 =  c((1e3-sens*1e3),  (1e3-spec*1e3),  2,            160))
}
map$sens_spec_tbl <- purrr::map2(map$sens, map$spec, get_sens_spec)

# onset to deaths
tod_paramsdf <- tibble::tibble(name = c("mod", "sod"),
                               min  = c(10,     0.01),
                               init = c(14,     0.7),
                               max =  c(30,     3.00),
                               dsc1 = c(2.657,  -0.236),
                               dsc2 = c(0.05,   0.05))

# everything else for region
wrap_make_IFR_model <- function(curve, inputdata, sens_spec_tbl, demog) {
  ifr_paramsdf <- make_ma_reparamdf(num_mas = 5)
  knot_paramsdf <- make_splinex_reparamdf(max_xvec = list("name" = "x4", min = 180, init = 190, max = 200, dsc1 = 180, dsc2 = 200),
                                          num_xs = 4)

  if (curve$nm == "expgrowth") {
    infxn_paramsdf <- make_spliney_reparamdf(max_yvec = list("name" = "y5", min = 0, init = 9, max = 15.42, dsc1 = 0, dsc2 = 15.42),
                                             num_ys = 5)
  } else {
    infxn_paramsdf <- make_spliney_reparamdf(max_yvec = list("name" = "y3", min = 0, init = 9, max = 15.42, dsc1 = 0, dsc2 = 15.42),
                                             num_ys = 5)
  }
  noise_paramsdf <- make_noiseeff_reparamdf(num_Nes = 5, min = 0, init = 5, max = 10)

  # bring together
  df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, knot_paramsdf, sens_spec_tbl, noise_paramsdf, tod_paramsdf)

  # make mod
  mod1 <- COVIDCurve::make_IFRmodel_agg$new()
  mod1$set_MeanTODparam("mod")
  mod1$set_CoefVarOnsetTODparam("sod")
  mod1$set_IFRparams(paste0("ma", 1:5))
  mod1$set_maxMa("ma5")
  mod1$set_Knotparams(paste0("x", 1:4))
  mod1$set_relKnot("x4")
  mod1$set_Infxnparams(paste0("y", 1:5))
  mod1$set_relInfxn("y3")
  mod1$set_Noiseparams(c(paste0("Ne", 1:5)))
  mod1$set_Serotestparams(c("sens", "spec", "sero_rate"))
  mod1$set_Serodayparams("sero_day")
  mod1$set_data(inputdata)
  mod1$set_demog(demog)
  mod1$set_paramdf(df_params)
  mod1$set_rho(rep(1, 5))
  mod1$set_rcensor_day(.Machine$integer.max)
  # out
  mod1
}

map$modelobj <-  purrr::pmap(map[,c("curve", "inputdata", "sens_spec_tbl", "demog")], wrap_make_IFR_model)

#......................
# names
#......................
fit_map <- map %>%
  dplyr::mutate(sim = paste0("sim", 1:nrow(.)),
                lvl = purrr::map_chr(curve, function(x){unique(x$nm)})) %>%
  dplyr::select(c("sim", "lvl", dplyr::everything()))


#............................................................
# Come Together
#...........................................................
# save out full for later manips
dir.create("data/param_map/SimCurves/", recursive = TRUE)
saveRDS(fit_map, "data/param_map/SimCurves/simfit_param_map.RDS")

# select what we need for fits and make outpaths
fit_map_modelobj <- fit_map %>%
  dplyr::select(c("sim", "modelobj"))
lapply(split(fit_map_modelobj, 1:nrow(fit_map_modelobj)), function(x){
  saveRDS(x, paste0("data/param_map/SimCurves/",
                    x$sim, ".RDS"))
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
                                      burnin = 1e4,
                                      samples = 1e4,
                                      rungs = 50,
                                      GTI_pow = 3,
                                      thinning = 10,
                                      cluster = cl)
  parallel::stopCluster(cl)
  gc()

  # out
  dir.create("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/SimCurves/", recursive = TRUE)
  outpath = paste0("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/SimCurves/",
                   mod$sim, ".RDS")
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
file_param_map <- list.files(path = "data/param_map/SimCurves/",
                             pattern = "*.RDS",
                             full.names = TRUE)
file_param_map <- tibble::tibble(path = file_param_map)
# remove non-fit items that are for carrying forward simulations
file_param_map <- file_param_map[!grepl("small_param_map.RDS", file_param_map$path),]
file_param_map <- file_param_map[!grepl("simfit_param_map.RDS", file_param_map$path),]


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
make(plan, parallelism = "clustermq",
     jobs = nrow(file_param_map),
     log_make = "SimCurves.log", verbose = 2,
     log_progress = TRUE,
     log_build_times = FALSE,
     recoverable = FALSE,
     history = FALSE,
     session_info = FALSE,
     lock_envir = FALSE, # unlock environment so parallel::clusterApplyLB in drjacoby can work
     lock_cache = FALSE)



cat("************** Drake Finished **************************")






