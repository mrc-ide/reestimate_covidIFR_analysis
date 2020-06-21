####################################################################################
## Purpose: Run out simulations
##
## Notes: Assumes SLURM cluster
####################################################################################
setwd("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis")
library(rslurm)
library(tidyverse)
library(COVIDCurve)
source("R/covidcurve_helper_functions.R")
source("R/simple_seir_model.R")

#............................................................
# Basic Incidence Curve
#...........................................................
nsims <- 1e2

#......................
# exponential growth with intervetions
#......................
intervene <- lapply(1:nsims, function(x){
  run_simple_seir(N = 3e7,
                  E0 = 50,
                  R0 = 0,
                  betas = c(0.17, 0.13, 0.12, 0.1),
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
  dplyr::rename(time = step)

#............................................................
# Simulate Under Model
#...........................................................
map <- expand.grid(curve = intervene,
                   sens = c(0.85, 0.90, 0.95),
                   spec = c(0.85, 0.90, 0.95),
                   popN = c(8e7, 5e7, 3e7, 1e7, 5e6))
map <- tibble::as_tibble(map)
#......................
# setup fatality data
#......................
fatalitydata <- data.frame(strata = paste0("ma", 1:8),
                           ifr = c(0, 0, 0, 0.02, 0.03, 0.05, 0.1, 0.2),
                           rho = 1/8)

#......................
# rung covidcurve simulator
#......................
wrap_sim <- function(curve, sens, spec, popN) {
  COVIDCurve::Aggsim_infxn_2_death(
    fatalitydata = fatalitydata,
    m_od = 18.8,
    s_od = 0.45,
    curr_day = 200,
    level = "Time-Series",
    infections = curve$infxns,
    simulate_seroprevalence = TRUE,
    sens = sens,
    spec = spec,
    sero_delay_rate = 10,
    popN = popN)
}

map$inputdata <- purrr::pmap(map, wrap_sim)

#......................
# make IFR model
#......................
get_sens_spec <- function(sens, spec) {
  tibble::tibble(name =  c("sens",          "spec",         "sero_rate", "sero_day"),
                 min =   c(sens-0.02,       0,              1,           140),
                 init =  c(sens,            spec,           10,           150),
                 max =   c(sens+0.02,       1,              20,           160),
                 dsc1 =  c(sens*1e4,        spec*1e2,       1,            140),
                 dsc2 =  c((1e4-sens*1e4),  (1e2-spec*1e2), 20,           160))
}

map$sens_spec_tbl <- purrr::map2(map$sens, map$spec, get_sens_spec)
wrap_make_IFR_model <- function(inputdata, sens_spec_tbl, popN) {
  ifr_paramsdf <- make_ma_reparamdf(num_mas = 5)
  knot_paramsdf <- make_splinex_reparamdf(max_xvec = list("name" = "x4", min = 180, init = 190, max = 200, dsc1 = 180, dsc2 = 200),
                                          num_xs = 4)
    infxn_paramsdf <- make_spliney_reparamdf(max_yvec = list("name" = "y3", min = 0, init = 9, max = 12, dsc1 = 0, dsc2 = 12),
                                             num_ys = 5)

  # bring together
  df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, knot_paramsdf, sens_spec_tbl)
  # format data
  inputdata <- list(obs_deaths = inputdata$AggDat,
                    obs_serologyrate = inputdata$seroprev$ObsPrev[150])
  # make mod
  mod1 <- COVIDCurve::make_IFRmodel_agg$new()
  mod1$set_MeanOnset(18.8)
  mod1$set_CoefVarOnset(0.45)
  mod1$set_level("Time-Series")
  mod1$set_data(inputdata)
  mod1$set_IFRparams(paste0("ma", 1:5))
  mod1$set_maxMa("ma5")
  mod1$set_Knotparams(paste0("x", 1:4))
  mod1$set_relKnot("x4")
  mod1$set_Infxnparams(paste0("y", 1:5))
  mod1$set_relInfxn("y3")
  mod1$set_Serotestparams(c("sens", "spec", "sero_rate"))
  mod1$set_Serodayparams("sero_day")
  mod1$set_popN(popN)
  mod1$set_paramdf(df_params)
  mod1$set_rho(rep(1/8, 8))
  mod1$set_rcensor_day(.Machine$integer.max)
  # out
  mod1
}

map$modelobj <-  purrr::pmap(map[,c("inputdata", "sens_spec_tbl", "popN")], wrap_make_IFR_model)

#............................................................
# run covidcurve
#...........................................................
parammap <- map %>%
  dplyr::select(c("curve", "sens", "spec", "modelobj")) # drop extras

wrap_run_covidcurve <- function(curve, sens, spec, modelobj) {
  COVIDCurve::run_IFRmodel_agg(IFRmodel = modelobj,
                               reparamIFR = TRUE,
                               reparamInfxn = TRUE,
                               reparamKnot = TRUE,
                               burnin = 1e4,
                               samples = 1e4,
                               chains = 10,
                               GTI_pow = 4,
                               rungs = 25)
}

#............................................................
# slurm
#...........................................................
# for slurm on LL
ntry <- 50
sjob <- rslurm::slurm_apply(f = wrap_run_covidcurve,
                            params = parammap,
                            jobname = 'sim_run_fullpriors',
                            nodes = ntry,
                            cpus_per_node = 1,
                            submit = T,
                            slurm_options = list(mem = "24g",
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "1-00:00:00"))


cat("************ Submission Successful ***********************")
