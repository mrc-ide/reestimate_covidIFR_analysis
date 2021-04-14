source("R/assertions_v5.R")
library(tidyverse)
library(stringr)


#' @title get mid age from agebands (factorized from cut)
get_mid_age <- function(ageband) {
  # character from factor
  ageband <- as.character(ageband)
  # extract and get mean
  age_mid <- purrr::map_dbl(ageband, function(x){
    lwr <- as.numeric(stringr::str_extract(ageband, "[0-9]+(?=\\,)"))
    lwr <- ifelse(lwr == 0, 0, lwr + 1) # treat as one-based
    upr <- as.numeric(stringr::str_extract(ageband, "[0-9]+?(?=])")) + 1 # treat as one-based
    # fix upper
    upr[upr == (999+1)] <- 100
    midages <- purrr::map2_dbl(lwr, upr, function(x, y) mean(c(x,y)))
    return(midages)})
  # out
  return(age_mid)
}

#' @title Log Transform IFR params
#' @details goal here is to be memory light
#' @importFrom magrittr %>%
#' @export

get_log10_transformed_IFR_cred_intervals <- function(path, by_chain = FALSE) {
  # read in
  IFRmodel_inf <- readRDS(path)
  # checks
  assert_custom_class(IFRmodel_inf$inputs$IFRmodel, "IFRmodel")
  assert_custom_class(IFRmodel_inf$mcmcout, "drjacoby_output")
  assert_custom_class(IFRmodel_inf, "IFRmodel_inf")
  assert_logical(by_chain)

  # grouping vars
  if (by_chain) {
    groupingvar <- c("chain", "param")
    params <- c("chain", IFRmodel_inf$inputs$IFRmodel$IFRparams)
  } else {
    groupingvar <- "param"
    params <- IFRmodel_inf$inputs$IFRmodel$IFRparams
  }

  IFRmodel_inf$mcmcout$output %>%
    dplyr::filter(stage == "sampling" & rung == "rung50") %>%
    dplyr::select_at(params) %>%
    tidyr::pivot_longer(., cols = params[!grepl("chain", params)], # if chain isn't included in vector, grepl won't do anything
                        names_to = "param", values_to = "est") %>%
    dplyr::mutate(est = log10(est)) %>%
    dplyr::group_by_at(groupingvar) %>%
    dplyr::summarise(
      min = min(est),
      LCI = quantile(est, 0.025),
      median = median(est),
      mean = mean(est),
      UCI = quantile(est, 0.975),
      max = max(est)
    )
}

#' @title extract data dictionary from IFRmodel_inf path
#' @details goal here is to be memory light

get_data_dict <- function(path) {
  modout <- readRDS(path)
  return(modout$inputs$IFRmodel$IFRdictkey)
}


#' @title Calculate seroprevalens from a path
#' @details goal here is to be memory light

get_strata_seroprevs <- function(path, dwnsmpl = 1e2) {
  modout <- readRDS(path)
  seroprevs <- COVIDCurve::draw_posterior_sero_curves(IFRmodel_inf = modout,
                                                      whichrung = "rung50",
                                                      dwnsmpl = dwnsmpl,
                                                      by_chain = FALSE)
  return(seroprevs)
}



#' @title Calculate stratified IFR from a path
#' @details goal here is to be memory light
get_strata_IFRs <- function(path) {
  modout <- readRDS(path)
  stratachar <- ifelse(grepl("age", basename(path)), "ageband",
                       ifelse(grepl("rgn", basename(path)), "region", NA))
  dictkey <- modout$inputs$IFRmodel$IFRdictkey %>%
    dplyr::rename(param = Strata)
  colnames(dictkey)[colnames(dictkey) == stratachar] <- "strata"
  # get ifrs
  ifrs <- COVIDCurve::get_cred_intervals(IFRmodel_inf = modout, whichrung = "rung50",
                                         what = "IFRparams", by_chain = FALSE)
  # out
  dplyr::left_join(dictkey, ifrs) %>%
    dplyr::mutate(param = factor(param, levels = paste0("ma", 1:nrow(dictkey))),
                  strata = forcats::fct_reorder(strata, as.numeric(param)))
}

#' @title Calculate overall IFR weighted for demography from a path
#' @details goal here is to be memory light
get_overall_IFRs <- function(path, whichstandard) {
  modout <- readRDS(path)
  out <- COVIDCurve::get_overall_IFR_cred_intervals(IFRmodel_inf = modout,
                                                    whichrung = "rung50",
                                                    whichstandard = whichstandard,
                                                    by_chain = FALSE)
  return(out)
}


#' @title Calculate stratified IFR precision from a path
#' @details goal here is to be memory light
get_strata_IFR_variance <- function(path, by_chain) {
  # read in
  IFRmodel_inf <- readRDS(path)
  # checks
  assert_custom_class(IFRmodel_inf$inputs$IFRmodel, "IFRmodel")
  assert_custom_class(IFRmodel_inf$mcmcout, "drjacoby_output")
  assert_custom_class(IFRmodel_inf, "IFRmodel_inf")
  assert_logical(by_chain)

  # grouping vars
  if (by_chain) {
    groupingvar <- c("chain", "param")
    params <- c("chain", IFRmodel_inf$inputs$IFRmodel$IFRparams)
  } else {
    groupingvar <- "param"
    params <- IFRmodel_inf$inputs$IFRmodel$IFRparams
  }

  IFRmodel_inf$mcmcout$output %>%
    dplyr::filter(stage == "sampling" & rung == "rung50") %>%
    dplyr::select_at(params) %>%
    tidyr::pivot_longer(., cols = params[!grepl("chain", params)], # if chain isn't included in vector, grepl won't do anything
                        names_to = "param", values_to = "est") %>%
    dplyr::group_by_at(groupingvar) %>%
    dplyr::summarise(var = var(est))
}

#' @title Calculate overall IFR from a path
#' @details goal here is to be memory light
get_sens_spec <- function(path) {
  modout <- readRDS(path)
  out <- COVIDCurve::get_cred_intervals(IFRmodel_inf = modout,
                                        what = "Serotestparams",
                                        whichrung = "rung50",
                                        by_chain = FALSE) %>%
    dplyr::filter(param %in% c("sens", "spec"))
  return(out)
}



#' @title Make IFR Model for MCMC Fitting
make_noSeroRev_IFR_model_fit <- function(num_mas, maxMa,
                                         groupvar, dat,
                                         num_xs, max_xveclist,
                                         num_ys, max_yveclist,
                                         sens_spec_tbl, tod_paramsdf,
                                         upperMa = 0.4) {


  ifr_paramsdf <- make_ma_reparamdf(num_mas = num_mas, upperMa = upperMa)

  knot_paramsdf <- make_splinex_reparamdf(max_xvec = max_xveclist,
                                          num_xs = num_xs)

  infxn_paramsdf <- make_spliney_reparamdf(max_yvec = max_yveclist,
                                           num_ys = num_ys)

  if (num_mas > 1) {
    noise_paramsdf <- make_noiseeff_reparamdf(num_Nes = num_mas, min = 0.5, init = 1, max = 1.5)
    df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, knot_paramsdf, sens_spec_tbl, noise_paramsdf, tod_paramsdf)
  } else {
    df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, knot_paramsdf, sens_spec_tbl, tod_paramsdf)
  }

  #......................
  # format data
  #......................
  dictkey <- tibble::tibble(groupvar = as.character(unlist(unique(dat$seroprevMCMC[, groupvar]))),
                            "Strata" = paste0("ma", 1:num_mas))
  colnames(dictkey) <- c(paste(groupvar), "Strata")

  # time series deaths
  obs_deaths <- dat$deaths_TSMCMC %>%
    dplyr::select(c("ObsDay", "deaths")) %>%
    dplyr::rename(Deaths = deaths)
  # sanity check
  if( any(duplicated(obs_deaths$ObsDay)) ) {
    stop("Time series has duplicated observed days")
  }

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
    dplyr::mutate(SeroLCI = NA,
                  SeroUCI = NA) %>%
    dplyr::select(c("SeroStartSurvey", "SeroEndSurvey", "Strata", "SeroPos", "SeroN", "SeroPrev", "SeroLCI", "SeroUCI")) %>%
    dplyr::mutate(Strata = factor(Strata, levels = paste0("ma", 1:num_mas))) %>%
    dplyr::arrange(SeroStartSurvey, Strata) %>%
    dplyr::mutate(Strata = as.character(Strata)) # coerce back to char for backward compat

  if ( all(c("SeroLCI", "SeroUCI") %in% colnames(dat$seroprevMCMC))) {
    if ( all( !is.na(c(dat$seroprevMCMC$SeroLCI, dat$seroprevMCMC$SeroUCI))) ) {
      obs_serology$SeroLCI <- dat$seroprevMCMC$SeroLCI
      obs_serology$SeroUCI <- dat$seroprevMCMC$SeroUCI
    } else {
      stop("You have LCI and UCI columns with incomplete information")
    }
  }

  inputdata <- list(obs_deaths = obs_deaths,
                    prop_deaths = prop_deaths,
                    obs_serology = obs_serology)

  demog <- dat$prop_pop %>%
    dplyr::left_join(., dictkey) %>%
    dplyr::select(c("Strata", "popN")) %>%
    dplyr::mutate(Strata = factor(Strata, levels = paste0("ma", 1:num_mas))) %>%
    dplyr::arrange(Strata) %>%
    dplyr::mutate(Strata = as.character(Strata)) # coerce back to char for backward compat


  # make mod
  if (num_mas > 1) {
    mod1 <- make_IFRmodel_age$new()
    mod1$set_MeanTODparam("mod")
    mod1$set_CoefVarOnsetTODparam("sod")
    mod1$set_IFRparams(paste0("ma", 1:num_mas))
    mod1$set_maxMa(maxMa)
    mod1$set_Knotparams(paste0("x", 1:num_xs))
    mod1$set_relKnot(max_xveclist[["name"]])
    mod1$set_Infxnparams(paste0("y", 1:num_ys))
    mod1$set_relInfxn(max_yveclist[["name"]])
    mod1$set_Serotestparams(c("sens", "spec", "sero_con_rate"))
    mod1$set_Noiseparams(paste0("Ne", 1:num_mas))
    mod1$set_data(inputdata)
    mod1$set_demog(demog)
    mod1$set_paramdf(df_params)
    mod1$set_rcensor_day(.Machine$integer.max)
    mod1$set_IFRdictkey(dictkey)
    # out
    mod1
  } else {
    mod1 <- make_IFRmodel_age$new()
    mod1$set_MeanTODparam("mod")
    mod1$set_CoefVarOnsetTODparam("sod")
    mod1$set_IFRparams("ma1")
    mod1$set_Knotparams(paste0("x", 1:num_xs))
    mod1$set_relKnot(max_xveclist[["name"]])
    mod1$set_Infxnparams(paste0("y", 1:num_ys))
    mod1$set_relInfxn(max_yveclist[["name"]])
    mod1$set_Serotestparams(c("sens", "spec", "sero_con_rate"))
    mod1$set_data(inputdata)
    mod1$set_demog(demog)
    mod1$set_paramdf(df_params)
    mod1$set_rcensor_day(.Machine$integer.max)
    mod1$set_IFRdictkey(dictkey)
    # out
    mod1
  }

}




#' @title Make IFR Model for MCMC Fitting
make_SeroRev_IFR_model_fit <- function(num_mas, maxMa,
                                       groupvar, dat,
                                       num_xs, max_xveclist,
                                       num_ys, max_yveclist,
                                       sens_spec_tbl, tod_paramsdf,
                                       upperMa = 0.4) {


  ifr_paramsdf <- make_ma_reparamdf(num_mas = num_mas, upperMa = upperMa)

  knot_paramsdf <- make_splinex_reparamdf(max_xvec = max_xveclist,
                                          num_xs = num_xs)

  infxn_paramsdf <- make_spliney_reparamdf(max_yvec = max_yveclist,
                                           num_ys = num_ys)

  if (num_mas > 1) {
    noise_paramsdf <- make_noiseeff_reparamdf(num_Nes = num_mas, min = 0.5, init = 1, max = 1.5)
    df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, knot_paramsdf, sens_spec_tbl, noise_paramsdf, tod_paramsdf)
  } else {
    df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, knot_paramsdf, sens_spec_tbl, tod_paramsdf)
  }

  #......................
  # format data
  #......................
  dictkey <- tibble::tibble(groupvar = as.character(unlist(unique(dat$seroprevMCMC[, groupvar]))), "Strata" = paste0("ma", 1:num_mas))
  colnames(dictkey) <- c(paste(groupvar), "Strata")

  # time series deaths
  obs_deaths <- dat$deaths_TSMCMC %>%
    dplyr::select(c("ObsDay", "deaths")) %>%
    dplyr::rename(Deaths = deaths)
  # sanity check
  if( any(duplicated(obs_deaths$ObsDay)) ) {
    stop("Time series has duplicated observed days")
  }

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
    dplyr::mutate(SeroLCI = NA,
                  SeroUCI = NA) %>%
    dplyr::select(c("SeroStartSurvey", "SeroEndSurvey", "Strata", "SeroPos", "SeroN", "SeroPrev", "SeroLCI", "SeroUCI")) %>%
    dplyr::mutate(Strata = factor(Strata, levels = paste0("ma", 1:num_mas))) %>%
    dplyr::arrange(SeroStartSurvey, Strata) %>%
    dplyr::mutate(Strata = as.character(Strata)) # coerce back to char for backward compat

  if ( all(c("SeroLCI", "SeroUCI") %in% colnames(dat$seroprevMCMC))) {
    if ( all( !is.na(c(dat$seroprevMCMC$SeroLCI, dat$seroprevMCMC$SeroUCI))) ) {
      obs_serology$SeroLCI <- dat$seroprevMCMC$SeroLCI
      obs_serology$SeroUCI <- dat$seroprevMCMC$SeroUCI
    } else {
      stop("You have LCI and UCI columns with incomplete information")
    }
  }

  inputdata <- list(obs_deaths = obs_deaths,
                    prop_deaths = prop_deaths,
                    obs_serology = obs_serology)

  demog <- dat$prop_pop %>%
    dplyr::left_join(., dictkey) %>%
    dplyr::select(c("Strata", "popN")) %>%
    dplyr::mutate(Strata = factor(Strata, levels = paste0("ma", 1:num_mas))) %>%
    dplyr::arrange(Strata) %>%
    dplyr::mutate(Strata = as.character(Strata)) # coerce back to char for backward compat


  # make mod
  if (num_mas > 1) {
    mod1 <- make_IFRmodel_age$new()
    mod1$set_MeanTODparam("mod")
    mod1$set_CoefVarOnsetTODparam("sod")
    mod1$set_IFRparams(paste0("ma", 1:num_mas))
    mod1$set_maxMa(maxMa)
    mod1$set_Knotparams(paste0("x", 1:num_xs))
    mod1$set_relKnot(max_xveclist[["name"]])
    mod1$set_Infxnparams(paste0("y", 1:num_ys))
    mod1$set_relInfxn(max_yveclist[["name"]])
    mod1$set_Serotestparams(c("sens", "spec", "sero_con_rate", "sero_rev_shape", "sero_rev_scale"))
    mod1$set_Noiseparams(paste0("Ne", 1:num_mas))
    mod1$set_data(inputdata)
    mod1$set_demog(demog)
    mod1$set_paramdf(df_params)
    mod1$set_rcensor_day(.Machine$integer.max)
    mod1$set_IFRdictkey(dictkey)
    # out
    mod1
  } else {
    mod1 <- make_IFRmodel_age$new()
    mod1$set_MeanTODparam("mod")
    mod1$set_CoefVarOnsetTODparam("sod")
    mod1$set_IFRparams("ma1")
    mod1$set_maxMa(maxMa)
    mod1$set_Knotparams(paste0("x", 1:num_xs))
    mod1$set_relKnot(max_xveclist[["name"]])
    mod1$set_Infxnparams(paste0("y", 1:num_ys))
    mod1$set_relInfxn(max_yveclist[["name"]])
    mod1$set_Serotestparams(c("sens", "spec", "sero_con_rate", "sero_rev_shape", "sero_rev_scale"))
    mod1$set_data(inputdata)
    mod1$set_demog(demog)
    mod1$set_paramdf(df_params)
    mod1$set_rcensor_day(.Machine$integer.max)
    mod1$set_IFRdictkey(dictkey)
    # out
    mod1
  }

}




#' @title Make IFR Uniform Distributed Reparameterized Param Df
#' @param num_mas positive interger; Number of IFR strata to infer

make_ma_reparamdf <- function(num_mas = 10, upperMa) {
  assert_pos_int(num_mas)
  tibble::tibble(name = paste0("ma", 1:num_mas),
                 min  = rep(0, size = num_mas),
                 init = rep(0.1, size = num_mas),
                 max = c(rep(1, times = num_mas-1), upperMa),
                 dsc1 = rep(0, size = num_mas),
                 dsc2 = c(rep(1, times = num_mas-1), upperMa))
}


#' @title Make Spline Y Position Reparameterized Param Df
#' @param num_ys positive interger; Number of Spline Y positions to infer
#' @param max_yval integer; Describes the maximum number of people that can be infected in the population -- the ymax position.

make_spliney_reparamdf <- function(max_yval,
                                   num_ys = 5) {
  assert_pos_int(num_ys)
  assert_pos_int(max_yval)

  out <- tibble::tibble(name = paste0("y", 1:num_ys),
                        min  = rep(0, size = num_ys),
                        init = rep(0.1, size = num_ys),
                        max = rep(max_yval, size = num_ys),
                        dsc1 = rep(0, size = num_ys),
                        dsc2 = rep(max_yval, size = num_ys))
  out %>%
    dplyr::arrange(name)
}

#' @title Make Spline X (Knots) Position Reparameterized Param Df
#' @param num_xs positive interger; Number of Spline X positions to infer
#' @param max_xval vector; Describes the maximum observed time -- the xmax position.
#' @description Will break up "X-knots" into "n" non-overlapping, contiguous discrete groups as specified by the user
make_splinex_reparamdf <- function(max_xval = 200,
                                   num_xs = 4) {
  assert_pos_int(num_xs)
  assert_pos_int(max_xval)
  day_dur <- max_xval/num_xs
  x_start <- c(1, cumsum(rep(day_dur, times = num_xs-1)))
  x_end <- c(x_start[-1]-1, max_xval) # not the first value, and then -1 for day so contiguous nonoverlap


  out %>%
    dplyr::arrange(name)
}



#' @title Make Noise Effect Reparameterized Param Df
#' @param num_Ne positive interger; Number of Noise effect parameters to create
make_noiseeff_reparamdf <- function(num_Nes = 4,  min = 0.5, init = 1, max = 1.5) {
  assert_pos_int(num_Nes)
  assert_numeric(min)
  assert_numeric(init)
  assert_numeric(max)

  # normal(1, 0.05) for all Nes
  tibble::tibble(name = paste0("Ne", 1:num_Nes),
                 min  = rep(min, size = num_Nes),
                 init = rep(init, size = num_Nes),
                 max = rep(max, size = num_Nes),
                 dsc1 = rep(1, size = num_Nes),
                 dsc2 = rep(0.05, size = num_Nes))
}

