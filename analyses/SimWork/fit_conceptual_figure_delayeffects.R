####################################################################################
## Purpose: Plot for Figure Showing Delays and Inference Framework
##
## Notes:
####################################################################################
set.seed(1234)
library(COVIDCurve)
library(tidyverse)
library(drake)
source("R/covidcurve_helper_functions.R")
source("R/my_themes.R")

#............................................................
# Read in Various Scenarios for Incidence Curves
#...........................................................
infxn_shapes <- readr::read_csv("data/simdat/infxn_curve_shapes.csv")
interveneflat <- infxn_shapes$intervene
# note need more infxns for sensitivity to be apparent on conceptual diagrams
interveneflat <- interveneflat * 1.5
interveneflat <- c(interveneflat, round(seq(from = interveneflat[200],
                                            to = 10, length.out = 100)))


# read in fitted rate of seroreversion parameter
weibullparams <- readRDS("results/prior_inputs/weibull_params.RDS")
weibullparams$wscale <- weibullparams$wscale - 13.3 # account for delay in onset of symptoms to seroconversion

#............................................................
# setup fatality data
#............................................................
# make up fatality data
fatalitydata <- tibble::tibble(Strata = c("ma1", "ma2", "ma3"),
                               IFR = c(1e-3, 0.05, 0.1),
                               Rho = 1)
demog <- tibble::tibble(Strata = c("ma1", "ma2", "ma3"),
                        popN = c(1.3e6, 9e5, 8e5))

# run COVIDCurve sims for no seroreversion and seroreversion
dat <- COVIDCurve::Agesim_infxn_2_death(
  fatalitydata = fatalitydata,
  demog = demog,
  m_od = 19.8,
  s_od = 0.85,
  curr_day = 300,
  infections = interveneflat,
  simulate_seroreversion = FALSE,
  smplfrac = 1e-3,
  sens = 0.85,
  spec = 0.95,
  sero_delay_rate = 18.3,
  return_linelist = FALSE)

serorev_dat <- COVIDCurve::Agesim_infxn_2_death(
  fatalitydata = fatalitydata,
  demog = demog,
  m_od = 19.8,
  s_od = 0.85,
  curr_day = 300,
  infections = interveneflat,
  simulate_seroreversion = TRUE,
  sero_rev_shape = weibullparams$wshape,
  sero_rev_scale = weibullparams$wscale,
  smplfrac = 1e-3,
  sens = 0.85,
  spec = 0.95,
  sero_delay_rate = 18.3,
  return_linelist = FALSE)



#............................................................
#----- Model & Fit #-----
#...........................................................
#......................
# wrangle input data from non-seroreversion fit
#......................
# liftover obs serology
sero_days <- c(150, 200)
sero_days <- lapply(sero_days, function(x){seq(from = (x-5), to = (x+5), by = 1)})
obs_serology <- dat$StrataAgg_Seroprev %>%
  dplyr::group_by(Strata) %>%
  dplyr::filter(ObsDay %in% unlist(sero_days)) %>%
  dplyr::mutate(serodaynum = sort(rep(1:length(sero_days), 11))) %>%
  dplyr::mutate(
    SeroPos = ObsPrev * testedN,
    SeroN = testedN ) %>%
  dplyr::group_by(Strata, serodaynum) %>%
  dplyr::summarise(SeroPos = mean(SeroPos),
                   SeroN = mean(SeroN)) %>% # seroN doesn't change
  dplyr::mutate(SeroStartSurvey = sapply(sero_days, median) - 5,
                SeroEndSurvey = sapply(sero_days, median) + 5,
                SeroPos = round(SeroPos),
                SeroPrev = SeroPos/SeroN,
                SeroLCI = NA,
                SeroUCI = NA) %>%
  dplyr::select(c("SeroStartSurvey", "SeroEndSurvey", "Strata", "SeroPos", "SeroN", "SeroPrev", "SeroLCI", "SeroUCI")) %>%
  dplyr::ungroup(.) %>%
  dplyr::arrange(SeroStartSurvey, Strata)


# proportion deaths
prop_deaths <- dat$StrataAgg_TimeSeries_Death %>%
  dplyr::group_by(Strata) %>%
  dplyr::summarise(deaths = sum(Deaths)) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate(PropDeaths = deaths/sum(dat$Agg_TimeSeries_Death$Deaths)) %>%
  dplyr::select(-c("deaths"))

# make data out
reginputdata <- list(obs_deaths = dat$Agg_TimeSeries_Death,
                     prop_deaths = prop_deaths,
                     obs_serology = obs_serology)

#......................
# wrangle input data from seroreversion fit
#......................
# sero tidy up
sero_days <- c(150, 200)
sero_days <- lapply(sero_days, function(x){seq(from = (x-5), to = (x+5), by = 1)})
obs_serology <- serorev_dat$StrataAgg_Seroprev %>%
  dplyr::group_by(Strata) %>%
  dplyr::filter(ObsDay %in% unlist(sero_days)) %>%
  dplyr::mutate(serodaynum = sort(rep(1:length(sero_days), 11))) %>%
  dplyr::mutate(
    SeroPos = ObsPrev * testedN,
    SeroN = testedN ) %>%
  dplyr::group_by(Strata, serodaynum) %>%
  dplyr::summarise(SeroPos = mean(SeroPos),
                   SeroN = mean(SeroN)) %>% # seroN doesn't change
  dplyr::mutate(SeroStartSurvey = sapply(sero_days, median) - 5,
                SeroEndSurvey = sapply(sero_days, median) + 5,
                SeroPos = round(SeroPos),
                SeroPrev = SeroPos/SeroN) %>%
  dplyr::select(c("SeroStartSurvey", "SeroEndSurvey", "Strata", "SeroPos", "SeroN", "SeroPrev")) %>%
  dplyr::ungroup(.) %>%
  dplyr::arrange(SeroStartSurvey, Strata) %>%
  dplyr::mutate(SeroLCI = NA,
                SeroUCI = NA) # just add these in for catch


# proportion deaths
prop_deaths <- serorev_dat$StrataAgg_TimeSeries_Death %>%
  dplyr::group_by(Strata) %>%
  dplyr::summarise(deaths = sum(Deaths)) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate(PropDeaths = deaths/sum(serorev_dat$Agg_TimeSeries_Death$Deaths)) %>%
  dplyr::select(-c("deaths"))

# make data out
serorev_inputdata <- list(obs_deaths = dat$Agg_TimeSeries_Death,
                          prop_deaths = prop_deaths,
                          obs_serology = obs_serology)



#......................
# make IFR model
#......................
# sens/spec
sens_spec_tbl <- tibble::tibble(name =  c("sens",  "spec"),
                                min =   c(0.5,      0.5),
                                init =  c(0.85,     0.95),
                                max =   c(1,        1),
                                dsc1 =  c(850.5,    950.5),
                                dsc2 =  c(150.5,    50.5))

# delay priors
tod_paramsdf <- tibble::tibble(name = c("mod", "sod", "sero_con_rate"),
                               min  = c(18,     0,     16),
                               init = c(19,     0.85,  18),
                               max =  c(20,     1,     21),
                               dsc1 = c(19.8,   2550,  18.3),
                               dsc2 = c(0.1,    450,   0.1))

serorev <- tibble::tibble(name = c("sero_rev_shape",     "sero_rev_scale"),
                         min  = c(1,                     199),
                         init = c(2,                     200),
                         max =  c(3.5,                   205),
                         dsc1 = c(weibullparams$wshape,  weibullparams$wscale),
                         dsc2 = c(0.1,                   0.1))

# combine
tod_paramsdf_serorev <- rbind(tod_paramsdf, serorev)



# make param dfs
ifr_paramsdf <- make_ma_reparamdf(num_mas = 3, upperMa = 0.4)
knot_paramsdf <- make_splinex_reparamdf(max_xval = 300,
                                        num_xs = 9)
infxn_paramsdf <- make_spliney_reparamdf(max_yval = 14.91,
                                         num_ys = 10)
noise_paramsdf <- make_noiseeff_reparamdf(num_Nes = 3, min = 0.5, init = 1, max = 1.5)
# bring together
df_params_reg <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, noise_paramsdf, knot_paramsdf, sens_spec_tbl, tod_paramsdf)
df_params_serorev <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, noise_paramsdf, knot_paramsdf, sens_spec_tbl, tod_paramsdf_serorev)


#......................
# make model for serorev and regular
#......................
# reg
mod1_reg <- COVIDCurve::make_IFRmodel_age$new()
mod1_reg$set_MeanTODparam("mod")
mod1_reg$set_CoefVarOnsetTODparam("sod")
mod1_reg$set_IFRparams(paste0("ma", 1:3))
mod1_reg$set_maxMa("ma3")
mod1_reg$set_Knotparams(paste0("x", 1:9))
mod1_reg$set_Infxnparams(paste0("y", 1:10))
mod1_reg$set_Noiseparams(c(paste0("Ne", 1:3)))
mod1_reg$set_Serotestparams(c("sens", "spec", "sero_con_rate"))
mod1_reg$set_data(reginputdata)
mod1_reg$set_demog(demog)
mod1_reg$set_paramdf(df_params_reg)
mod1_reg$set_rcensor_day(.Machine$integer.max)
# serorev
mod1_serorev <- COVIDCurve::make_IFRmodel_age$new()
mod1_serorev$set_MeanTODparam("mod")
mod1_serorev$set_CoefVarOnsetTODparam("sod")
mod1_serorev$set_IFRparams(paste0("ma", 1:3))
mod1_serorev$set_maxMa("ma3")
mod1_serorev$set_Knotparams(paste0("x", 1:9))
mod1_serorev$set_Infxnparams(paste0("y", 1:10))
mod1_serorev$set_Noiseparams(c(paste0("Ne", 1:3)))
mod1_serorev$set_Serotestparams(c("sens", "spec", "sero_con_rate", "sero_rev_shape", "sero_rev_scale"))
mod1_serorev$set_data(serorev_inputdata)
mod1_serorev$set_demog(demog)
mod1_serorev$set_paramdf(df_params_serorev)
mod1_serorev$set_rcensor_day(.Machine$integer.max)

#............................................................
#---- Come Together #----
#...........................................................
fit_map <- tibble::tibble(
  name = c("reg_mod", "serorev_mod"),
  infxns = list(interveneflat, NULL), # Null since same infections
  simdat = list(dat, serorev_dat),
  modelobj = list(mod1_reg, mod1_serorev),
  rungs = 50,
  burnin = 1e4,
  samples = 1e4,
  thinning = 10)


#......................
# fitmap out
#......................
# select what we need for fits and make outpaths
dir.create("data/param_map/Fig_ConceptualFits/", recursive = T)
lapply(split(fit_map, 1:nrow(fit_map)), function(x){
  saveRDS(x, paste0("data/param_map/Fig_ConceptualFits/",
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
                                      reparamInfxn = FALSE,
                                      reparamKnots = FALSE,
                                      chains = n_chains,
                                      burnin = mod$burnin,
                                      samples = mod$samples,
                                      rungs = mod$rungs,
                                      GTI_pow = 3.0,
                                      cluster = cl,
                                      thinning = mod$thinning)
  parallel::stopCluster(cl)
  gc()

  # out
  dir.create("results/Fig_ConceptualFits/", recursive = TRUE)
  outpath = paste0("results/Fig_ConceptualFits/",
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
file_param_map <- list.files(path = "data/param_map/Fig_ConceptualFits/",
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
     log_make = "ConceptFig_drake.log", verbose = 2,
     log_progress = TRUE,
     log_build_times = FALSE,
     recoverable = FALSE,
     history = FALSE,
     session_info = FALSE,
     lock_envir = FALSE, # unlock environment so parallel::clusterApplyLB in drjacoby can work
     lock_cache = FALSE)



cat("************** Drake Finished **************************")




