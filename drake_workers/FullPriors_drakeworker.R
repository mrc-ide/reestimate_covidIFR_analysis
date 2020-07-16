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
# setup fatality data
#............................................................
fatalitydata <- tibble::tibble(Strata = paste0("ma", 1:8),
                               IFR = c(0, 0, 0, 0.02, 0.03, 0.05, 0.1, 0.2),
                               Rho = 1,
                               Ne = 1/8)

#............................................................
# Basic Incidence Curve
#...........................................................
nsims <- 100
popN <- 3e6
infxns <- lapply(1:nsims, function(x){
  run_simple_seir(N = popN,
                  E0 = 50,
                  R0 = 0,
                  betas = c(0.25, 0.13, 0.12, 0.11),
                  beta_changes = c(1, 130, 140, 150),
                  sigma = 0.2,
                  gamma = 0.2,
                  time = 200)
})
intervene <- infxns %>%
  dplyr::bind_rows(.) %>%
  dplyr::group_by(step) %>%
  dplyr::summarise(
    infxns = mean(I)
  ) %>%
  dplyr::mutate_if(is.numeric, round, 0) %>%
  dplyr::rename(time = step)

# plot(intervene$infxns)
popNvec <- c(1.5, 2.5, 5, 10, 50, 100)
poplist <- lapply(popNvec, function(x){
  totN <- round( sum(intervene$infxns) * x )
  popvec <- totN * fatalitydata$Ne
  demog <- tibble::tibble(Strata = paste0("ma", 1:8),
                          popN = as.integer(popvec))
  return(demog)
})

#............................................................
# Simulate Under Model
#...........................................................
map <- expand.grid(curve = list(intervene),
                   sens = c(0.85, 0.90, 0.95, 0.99),
                   spec = c(0.85, 0.90, 0.95, 0.99),
                   demog = poplist)
map <- tibble::as_tibble(map)

#......................
# rung covidcurve simulator
#......................
wrap_sim <- function(curve, sens, spec, demog, sero_day) {

  dat <- COVIDCurve::Aggsim_infxn_2_death(
    fatalitydata = fatalitydata,
    m_od = 14.26,
    s_od = 0.79,
    curr_day = 200,
    infections = curve$infxns,
    simulate_seroprevalence = TRUE,
    sens = sens,
    spec = spec,
    sero_delay_rate = 10,
    demog = demog)

  # sero tidy up
  obs_serology <- dat$seroprev %>%
    dplyr::group_by(Strata) %>%
    dplyr::filter(event_obs_day == sero_day) %>%
    dplyr::rename(
      SeroDay = event_obs_day,
      SeroPrev = ObsPrev) %>%
    dplyr::select(c("SeroDay", "Strata", "SeroPrev")) %>%
    dplyr::mutate(SeroDay = "sero_day")

  # death dat
  dat$AggDat <- dat$AggDat %>%
    dplyr::mutate(Strata = as.character(Strata))


  datinput <- list(obs_deaths = dat$AggDat,
                   obs_serology = obs_serology)
  return(datinput)

}

map$inputdata <- purrr::pmap(map, wrap_sim, sero_day = 150)

#......................
# make IFR model
#......................
# sens/spec
get_sens_spec <- function(sens, spec) {
  tibble::tibble(name =  c("sens",          "spec",         "sero_rate",  "sero_day"),
                 min =   c(0.5,              0.5,            0,            140),
                 init =  c(0.8,              0.8,            0.5,          150),
                 max =   c(1,       1,                       1,            160),
                 dsc1 =  c(sens*1e3,        spec*1e2,        50,           140),
                 dsc2 =  c((1e3-sens*1e3),  (1e2-spec*1e2),  50,           160))
}
map$sens_spec_tbl <- purrr::map2(map$sens, map$spec, get_sens_spec)

# onset to deaths
tod_paramsdf <- tibble::tibble(name = c("mod", "sod"),
                               min  = c(10,    0.01),
                               init = c(14,    0.7),
                               max =  c(20,    1.00),
                               dsc1 = c(2.7,   -0.23),
                               dsc2 = c(0.05,   0.05))

# everything else
wrap_make_IFR_model <- function(inputdata, sens_spec_tbl, demog) {
  ifr_paramsdf <- make_ma_reparamdf(num_mas = 8)
  knot_paramsdf <- make_splinex_reparamdf(max_xvec = list("name" = "x4", min = 180, init = 190, max = 200, dsc1 = 180, dsc2 = 200),
                                          num_xs = 4)
  infxn_paramsdf <- make_spliney_reparamdf(max_yvec = list("name" = "y3", min = 0, init = 9, max = 12, dsc1 = 0, dsc2 = 12),
                                           num_ys = 5)
  noise_paramsdf <- make_noiseeff_reparamdf(num_Nes = 8, min = 0, init = 5, max = 10)

  # bring together
  df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, knot_paramsdf, sens_spec_tbl, noise_paramsdf, tod_paramsdf)

  # make mod
  mod1 <- COVIDCurve::make_IFRmodel_agg$new()
  mod1$set_MeanTODparam("mod")
  mod1$set_CoefVarOnsetTODparam("sod")
  mod1$set_IFRparams(paste0("ma", 1:8))
  mod1$set_maxMa("ma8")
  mod1$set_Knotparams(paste0("x", 1:4))
  mod1$set_relKnot("x4")
  mod1$set_Infxnparams(paste0("y", 1:5))
  mod1$set_relInfxn("y3")
  mod1$set_Noiseparams(c(paste0("Ne", 1:8)))
  mod1$set_Serotestparams(c("sens", "spec", "sero_rate"))
  mod1$set_Serodayparams("sero_day")
  mod1$set_data(inputdata)
  mod1$set_demog(demog)
  mod1$set_paramdf(df_params)
  mod1$set_rho(rep(1, 8))
  mod1$set_rcensor_day(.Machine$integer.max)
  # out
  mod1
}

map$modelobj <-  purrr::pmap(map[,c("inputdata", "sens_spec_tbl", "demog")], wrap_make_IFR_model)

#......................
# names
#......................
fit_map <- map %>%
  dplyr::mutate(sim = paste0("sim", 1:nrow(.))) %>%
  dplyr::select(c("sim", dplyr::everything()))

fit_map_sm <- fit_map %>%
  dplyr::select(c("sim", "sens", "spec", "demog"))

#............................................................
# Come Together
#...........................................................

# select what we need for fits and make outpaths
dir.create("data/param_map/full_prior_sims/", recursive = T)
lapply(split(fit_map, 1:nrow(fit_map)), function(x){
  saveRDS(x, paste0("data/param_map/full_prior_sims/",
                    x$sim, ".RDS"))
})
saveRDS(fit_map_sm, "data/param_map/full_prior_sims/small_param_map.RDS")



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
                                      burnin = 1e3,
                                      samples = 1e3,
                                      rungs = 25,
                                      GTI_pow = 3,
                                      cluster = cl)
  # out
  dir.create("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/full_prior_sims/", recursive = TRUE)
  outpath = paste0("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/full_prior_sims/",
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
file_param_map <- list.files(path = "data/param_map/full_prior_sims/",
                             pattern = "*.RDS",
                             full.names = TRUE)
file_param_map <- tibble::tibble(path = file_param_map)
file_param_map <- file_param_map[!grepl("small_param_map.RDS", file_param_map$path),]


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
     log_make = "FullPriorSims_drake.log", verbose = 2,
     log_progress = FALSE,
     log_build_times = FALSE,
     recoverable = FALSE,
     history = FALSE,
     session_info = FALSE,
     lock_envir = FALSE) # unlock environment so parallel::clusterApplyLB in drjacoby can work





