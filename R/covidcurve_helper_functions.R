source("R/assertions_v5.R")
library(tidyverse)
library(stringr)

#' @title Make IFR Model for MCMC Fitting
make_IFR_model_fit <- function(num_mas, maxMa,
                               groupvar, dat,
                               num_xs, max_xveclist,
                               num_ys, max_yveclist,
                               sens_spec_tbl, tod_paramsdf,
                               serodayparams) {

  # make dfs
  ifr_paramsdf <- make_ma_reparamdf(num_mas = num_mas)
  knot_paramsdf <- make_splinex_reparamdf(max_xvec = max_xveclist,
                                          num_xs = num_xs)
  infxn_paramsdf <- make_spliney_reparamdf(max_yvec = max_yveclist,
                                           num_ys = num_ys)
  noise_paramsdf <- make_noiseeff_reparamdf(num_Nes = num_mas, min = 1, init = 1, max = 1)


  # bring together
  df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, knot_paramsdf, sens_spec_tbl, noise_paramsdf, tod_paramsdf)
  #......................
  # format data
  #......................
  dictkey <- tibble::tibble(groupvar = as.character(unlist(unique(dat$seroprevMCMC[, groupvar]))), "Strata" = paste0("ma", 1:num_mas))
  colnames(dictkey) <- c(paste(groupvar), "Strata")

  # time series deaths
  obs_deaths <- dat$deaths_TSMCMC %>%
    dplyr::select(c("ObsDay", "deaths")) %>%
    dplyr::rename(Deaths = deaths)

  # prop deaths
  prop_deaths <- dplyr::left_join(dat$deaths_propMCMC, dictkey) %>%
    dplyr::select(c("Strata", "death_prop"))  %>%
    dplyr::rename(PropDeaths = death_prop) %>%
    dplyr::mutate(Strata = factor(Strata, levels = paste0("ma", 1:num_mas))) %>%
    dplyr::arrange(Strata) %>%
    dplyr::mutate(Strata = as.character(Strata)) # coerce back to char for backward compat

  # seroprev
  obs_serology <- dplyr::left_join(dat$seroprevMCMC, dictkey) %>%
    dplyr::rename(SeroPos = n_positive,
                  SeroN = n_tested,
                  SeroStartSurvey = ObsDaymin,
                  SeroEndSurvey = ObsDaymax) %>%
    dplyr::select(c("SeroStartSurvey", "SeroEndSurvey", "Strata", "SeroPos", "SeroN", "SeroPrev")) %>%
    dplyr::mutate(Strata = factor(Strata, levels = paste0("ma", 1:num_mas))) %>%
    dplyr::arrange(SeroStartSurvey, Strata) %>%
    dplyr::mutate(Strata = as.character(Strata)) # coerce back to char for backward compat


  inputdata <- list(obs_deaths = obs_deaths,
                    prop_deaths = prop_deaths,
                    obs_serology = obs_serology)

  demog <- dat$prop_pop %>%
    dplyr::left_join(., dictkey) %>%
    dplyr::select(c("Strata", "popN")) %>%
    dplyr::group_by(Strata) %>%
    dplyr::summarise(popN = round(sum(popN))) %>%
    dplyr::mutate(Strata = factor(Strata, levels = paste0("ma", 1:num_mas))) %>%
    dplyr::arrange(Strata) %>%
    dplyr::mutate(Strata = as.character(Strata)) # coerce back to char for backward compat


  # make mod
  mod1 <- make_IFRmodel_agg$new()
  mod1$set_MeanTODparam("mod")
  mod1$set_CoefVarOnsetTODparam("sod")
  mod1$set_IFRparams(paste0("ma", 1:num_mas))
  mod1$set_maxMa(maxMa)
  mod1$set_Knotparams(paste0("x", 1:num_xs))
  mod1$set_relKnot(max_xveclist[["name"]])
  mod1$set_Infxnparams(paste0("y", 1:num_ys))
  mod1$set_relInfxn(max_yveclist[["name"]])
  mod1$set_Serotestparams(c("sens", "spec", "sero_rate"))
  mod1$set_Noiseparams(paste0("Ne", 1:num_mas))
  mod1$set_data(inputdata)
  mod1$set_demog(demog)
  mod1$set_paramdf(df_params)
  mod1$set_rho(rep(1, num_mas))
  mod1$set_rcensor_day(.Machine$integer.max)
  mod1$set_IFRdictkey(dictkey)
  # out
  mod1
}



#' @title Make Simple Data Dictionary Key for IFR age-bands, regions, etc. to simple
#' @param strata string vector; Names of stata to simplify

make_ma_dict_key <- function(strata_names) {
  assert_string(strata_names)
  tibble::tibble(strata_name = strata_names,
                 param_name  = paste0("ma", 1:length(strata_name)))
}


#' @title Make IFR Uniform Distributed Reparameterized Param Df
#' @param num_mas positive interger; Number of IFR strata to infer

make_ma_reparamdf <- function(num_mas = 10) {
  assert_pos_int(num_mas)
  tibble::tibble(name = paste0("ma", 1:num_mas),
                 min  = rep(0, size = num_mas),
                 init = rep(0.1, size = num_mas),
                 max = rep(1, size = num_mas),
                 dsc1 = rep(0, size = num_mas),
                 dsc2 = rep(1, size = num_mas))
}


#' @title Make Spline Y Position Reparameterized Param Df
#' @param num_ys positive interger; Number of Spline Y positions to infer
#' @param max_yvec vector; Describes the ymax position. Rest will be uniform 0,1 for reparameterized positions

make_spliney_reparamdf <- function(max_yvec = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                   num_ys = 5) {
  assert_pos_int(num_ys)
  assert_in(names(max_yvec), c("name", "min", "init", "max", "dsc1", "dsc2"))
  assert_in(c("name", "min", "init", "max", "dsc1", "dsc2"), names(max_yvec))
  assert_string(max_yvec[["name"]])
  assert_numeric(max_yvec[["min"]])
  assert_numeric(max_yvec[["init"]])
  assert_numeric(max_yvec[["max"]])
  assert_numeric(max_yvec[["dsc1"]])
  assert_numeric(max_yvec[["dsc2"]])

  out <- tibble::tibble(name = paste0("y", 1:num_ys),
                        min  = rep(0, size = num_ys),
                        init = rep(0.1, size = num_ys),
                        max = rep(1, size = num_ys),
                        dsc1 = rep(0, size = num_ys),
                        dsc2 = rep(1, size = num_ys))
  out %>%
    dplyr::filter(name != max_yvec["name"]) %>%
    dplyr::bind_rows(., max_yvec) %>%
    dplyr::arrange(name)
}

#' @title Make Spline X (Knots) Position Reparameterized Param Df
#' @param num_xs positive interger; Number of Spline X positions to infer
#' @param max_xvec vector; Describes the xmax position. Rest will be uniform 0,1 for reparameterized positions

make_splinex_reparamdf <- function(max_xvec = list("name" = "x4", min = 180, init = 190, max = 200, dsc1 = 180, dsc2 = 200),
                                   num_xs = 4) {
  assert_pos_int(num_xs)
  assert_in(names(max_xvec), c("name", "min", "init", "max", "dsc1", "dsc2"))
  assert_in(c("name", "min", "init", "max", "dsc1", "dsc2"), names(max_xvec))
  assert_string(max_xvec[["name"]])
  assert_numeric(max_xvec[["min"]])
  assert_numeric(max_xvec[["init"]])
  assert_numeric(max_xvec[["max"]])
  assert_numeric(max_xvec[["dsc1"]])
  assert_numeric(max_xvec[["dsc2"]])

  out <- tibble::tibble(name = paste0("x", 1:num_xs),
                        min  = rep(0, num_xs),
                        init = seq(1e-3, 0.9, length.out = num_xs),
                        max =  rep(1, num_xs),
                        dsc1 = rep(0, num_xs),
                        dsc2 = rep(1, num_xs))
  out %>%
    dplyr::filter(name != max_xvec["name"]) %>%
    dplyr::bind_rows(., max_xvec) %>%
    dplyr::arrange(name)
}



#' @title Make Noise Effect Reparameterized Param Df
#' @param num_Ne positive interger; Number of Noise effect parameters to create
#' @details By default, the first noise effect parameter is used as scalar for the other
#'          noise effect parameters in a reparameterization framework. See \link{COVIDCurve}
#'          for further details.
make_noiseeff_reparamdf <- function(num_Nes = 4,
                                    min = 0,
                                    init = 5,
                                    max = 10) {
  assert_pos_int(num_Nes)
  assert_numeric(min)
  assert_numeric(init)
  assert_numeric(max)

  tibble::tibble(name = paste0("Ne", 1:num_Nes),
                 min  = rep(min, size = num_Nes),
                 init = rep(init, size = num_Nes),
                 max = rep(max, size = num_Nes),
                 dsc1 = rep(min, size = num_Nes),
                 dsc2 = rep(max, size = num_Nes))
}

