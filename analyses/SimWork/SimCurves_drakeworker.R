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
set.seed(48)


#............................................................
# Read in Various Scenarios for Incidence Curves
#...........................................................
infxn_shapes <- readr::read_csv("data/simdat/infxn_curve_shapes.csv")

#............................................................
# setup fatality data
#............................................................
# make up fatality data
fatalitydata <- tibble::tibble(Strata = c("ma1", "ma2", "ma3", "ma4", "ma5"),
                               IFR = c(1e-3, 1e-3, 0.05, 0.1, 0.2),
                               Rho = 1)
demog <- tibble::tibble(Strata = c("ma1", "ma2", "ma3", "ma4", "ma5"),
                        popN = c(1.2e6, 1.1e6, 1e6, 9e5, 8e5))


#............................................................
# Simulate Under Model
#...........................................................
map <- expand.grid(nm = c("expgrowth", "intervene", "secondwave"),
                   curve = list(infxn_shapes$expgrowth,
                                infxn_shapes$intervene,
                                infxn_shapes$secondwave),
                   sens = c(0.85, 0.90),
                   spec = c(0.95, 0.99))


map <- tibble::as_tibble(map) %>%
  dplyr::mutate(fatalitydata = list(fatalitydata),
                demog = list(demog))
#......................
# run covidcurve simulator
#......................
wrap_sim <- function(nm, curve, sens, spec, mod, sero_rate, fatalitydata, demog, sero_days) {

  dat <- COVIDCurve::Agesim_infxn_2_death(
    fatalitydata = fatalitydata,
    m_od = 19.26,
    s_od = 0.79,
    curr_day = 200,
    infections = curve,
    simulate_seroreversion = FALSE,
    sens = sens,
    spec = spec,
    sero_delay_rate = 18.3,
    demog = demog,
    smplfrac = 0.25,
    return_linelist = FALSE)

  # liftover proprtion deaths
  totdeaths <- sum(dat$StrataAgg_TimeSeries_Death$Deaths)
  prop_strata_obs_deaths <- dat$StrataAgg_TimeSeries_Death %>%
    dplyr::group_by(Strata) %>%
    dplyr::summarise(Deaths = sum(Deaths),
                     PropDeaths = Deaths/totdeaths) %>%
    dplyr::select(c("Strata", "PropDeaths"))

  # liftover obs serology
  obs_serology <- dat$StrataAgg_Seroprev %>%
    dplyr::group_by(Strata) %>%
    dplyr::filter(ObsDay %in% sero_days) %>%
    dplyr::mutate(
      SeroPos = round(ObsPrev * testedN),
      SeroN = testedN,
      SeroLCI = NA,
      SeroUCI = NA) %>%
    dplyr::rename(
      SeroPrev = ObsPrev) %>%
    dplyr::mutate(SeroStartSurvey = sero_days - 5,
                  SeroEndSurvey = sero_days + 5) %>%
    dplyr::select(c("SeroStartSurvey", "SeroEndSurvey", "Strata", "SeroPos", "SeroN", "SeroPrev", "SeroLCI", "SeroUCI")) %>%
    dplyr::ungroup(.) %>%
    dplyr::arrange(SeroStartSurvey, Strata)

  inputdata <- list(obs_deaths = dat$Agg_TimeSeries_Death,
                   prop_deaths = prop_strata_obs_deaths,
                   obs_serology = obs_serology)
  out <- list(simdat = dat,
              inputdata = inputdata)
  return(out)

}

# run simdat and extract results into separate pieces
map$simdat <- purrr::pmap(map, wrap_sim, sero_days = c(140, 160))
map$inputdata <- purrr::map(map$simdat, "inputdata")
map$simdat <- purrr::map(map$simdat, "simdat", sero_days = c(140, 160))

#......................
# make IFR model
#......................
# sens/spec
get_sens_spec_tbl <- function(sens, spec) {
  tibble::tibble(name =  c("sens",          "spec"),
                 min =   c(0.5,              0.5),
                 init =  c(0.9,              0.99),
                 max =   c(1,                1),
                 dsc1 =  c(sens*1e3,        spec*1e3),
                 dsc2 =  c((1e3-sens*1e3),  (1e3-spec*1e3)))
}
map$sens_spec_tbl <- purrr::map2(map$sens, map$spec, get_sens_spec_tbl)

# delay priors
tod_paramsdf <- tibble::tibble(name = c("mod", "sod", "sero_con_rate"),
                               min  = c(18,     0,     16),
                               init = c(19,     0.79,  18),
                               max =  c(20,     1,     21),
                               dsc1 = c(19.26,  2370,  18.3),
                               dsc2 = c(0.1,    630,   0.1))

# everything else for region
wrap_make_IFR_model <- function(nm, curve, inputdata, sens_spec_tbl, demog) {
  ifr_paramsdf <- make_ma_reparamdf(num_mas = 5, upperMa = 0.4)
  knot_paramsdf <- make_splinex_reparamdf(max_xvec = list("name" = "x4", min = 180, init = 190, max = 200, dsc1 = 180, dsc2 = 200),
                                          num_xs = 4)

  if (nm == "expgrowth") {
    infxn_paramsdf <- make_spliney_reparamdf(max_yvec = list("name" = "y5", min = 0, init = 9, max = 15.42, dsc1 = 0, dsc2 = 15.42),
                                             num_ys = 5)
  } else {
    infxn_paramsdf <- make_spliney_reparamdf(max_yvec = list("name" = "y3", min = 0, init = 9, max = 15.42, dsc1 = 0, dsc2 = 15.42),
                                             num_ys = 5)
  }
  noise_paramsdf <- make_noiseeff_reparamdf(num_Nes = 5, min = 0.5, init = 1, max = 1.5)

  # bring together
  df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, knot_paramsdf, sens_spec_tbl, noise_paramsdf, tod_paramsdf)

  # make mod
  mod1 <- COVIDCurve::make_IFRmodel_age$new()
  mod1$set_MeanTODparam("mod")
  mod1$set_CoefVarOnsetTODparam("sod")
  mod1$set_IFRparams(paste0("ma", 1:5))
  mod1$set_maxMa("ma5")
  mod1$set_Knotparams(paste0("x", 1:4))
  mod1$set_relKnot("x4")
  mod1$set_Infxnparams(paste0("y", 1:5))
  mod1$set_relInfxn("y3")
  mod1$set_Noiseparams(c(paste0("Ne", 1:5)))
  mod1$set_Serotestparams(c("sens", "spec", "sero_con_rate"))
  mod1$set_data(inputdata)
  mod1$set_demog(demog)
  mod1$set_paramdf(df_params)
  mod1$set_rcensor_day(.Machine$integer.max)
  # out
  mod1
}

map$modelobj <-  purrr::pmap(map[,c("nm", "curve", "inputdata", "sens_spec_tbl", "demog")], wrap_make_IFR_model)

#......................
# names
#......................
fit_map <- map %>%
  dplyr::mutate(sim = paste0("sim", 1:nrow(.))) %>%
  dplyr::select(c("sim", "nm", dplyr::everything()))


#............................................................
# Come Together
#...........................................................
# save out full for later manips
dir.create("data/param_map/SimCurves_noserorev/", recursive = TRUE)
saveRDS(fit_map, "data/param_map/SimCurves_noserorev/simfit_param_map.RDS")

# select what we need for fits and make outpaths
fit_map_modelobj <- fit_map %>%
  dplyr::select(c("sim", "modelobj"))
lapply(split(fit_map_modelobj, 1:nrow(fit_map_modelobj)), function(x){
  saveRDS(x, paste0("data/param_map/SimCurves_noserorev/",
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
  fit <- COVIDCurve::run_IFRmodel_age(IFRmodel = mod$modelobj[[1]],
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
  dir.create("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/SimCurves_noserorev/", recursive = TRUE)
  outpath = paste0("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/SimCurves_noserorev/",
                   mod$sim, "_NoSeroRev.RDS")
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
file_param_map <- list.files(path = "data/param_map/SimCurves_noserorev/",
                             pattern = "*.RDS",
                             full.names = TRUE)
file_param_map <- tibble::tibble(path = file_param_map)
# remove non-fit items that are for carrying forward simulations
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
        clustermq.template = "drake_clst/slurm_clustermq_LL.tmpl")
make(plan, parallelism = "clustermq",
     jobs = nrow(file_param_map),
     log_make = "SimCurves_noserorev.log", verbose = 2,
     log_progress = TRUE,
     log_build_times = FALSE,
     recoverable = FALSE,
     history = FALSE,
     session_info = FALSE,
     lock_envir = FALSE, # unlock environment so parallel::clusterApplyLB in drjacoby can work
     lock_cache = FALSE)



cat("************** Drake Finished **************************")






