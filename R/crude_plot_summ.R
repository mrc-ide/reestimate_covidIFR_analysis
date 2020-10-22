source("R/assertions_v5.R")



#' @title get mid age from agebands (factorized from cut)
#' @details function is copied over from covidcurve helper just for convenience
get_mid_age <- function(ageband) {
  # character from factor
  ageband <- as.character(ageband)
  # extract and get mean
  age_mid <- purrr::map_dbl(ageband, function(x){
    lwr <- as.numeric(stringr::str_extract(ageband, "[0-9]+(?=\\,)")) + 1 # treat as one-based
    upr <- as.numeric(stringr::str_extract(ageband, "[0-9]+?(?=])")) + 1 # treat as one-based
    # fix upper
    upr[upr == (999+1)] <- 100
    midages <- purrr::map2_dbl(lwr, upr, function(x, y) mean(c(x,y)))
    return(midages)})
  # out
  return(age_mid)
}

#' @title Rogan-Gladen Estimator for Correct Prevalence Observations
#' @return Rogan-Gladen Correct Sens/Spec
rogan_gladen <- function(obs_prev, sens, spec){
  ret <- (obs_prev + spec - 1)/(spec + sens - 1)
  ret <- ifelse(ret <= 0, NA, ret) # occurs when obs_prev + spec are less than 1
  return(ret)
}

#' @title Adjust Seroprevalence Studies for Sens/Spec
#' @param seroprevdat dataframe; dataframe of strata specific seroprevalences. Must contain SeroPrev column
#' @param sens numeric; sensitivity
#' @param spec numeric; specificity
adjust_seroprev <- function(seroprevdat, sens, spec) {
  assert_in("SeroPrev", colnames(seroprevdat))
  assert_bounded(sens, left = 0, right = 1)
  assert_bounded(spec, left = 0, right = 1)

  seroprevdat <- seroprevdat %>%
    magrittr::set_colnames(tolower(colnames(.))) %>%
    dplyr::mutate(seroprevadj = purrr::map_dbl(seroprev, rogan_gladen, sens = sens, spec = spec),
                  seromidpt = round( (obsdaymin + obsdaymax)/2 ) )
  return(seroprevdat)
}

#' @title Standardize Death Data
#' @param deathdat_long dataframe; Death data in long format, for every observation day
#' @param popdata dataframe;
#' @param groupingvar character; grouping variable
#' @param Nstandardization numeric; value to standardize deaths due; default is deaths per million
standardize_deathdat <- function(deathdat_long, popdat, groupingvar, Nstandardization = 1e6) {
  assert_in(c("ObsDay", "Deaths"), colnames(deathdat_long))
  assert_in("popN", colnames(popdat))
  assert_in(groupingvar, colnames(deathdat_long))
  assert_in(groupingvar, colnames(popdat))

  # tidy popdat
  popdat <- popdat %>%
    dplyr::group_by_at(groupingvar) %>%
    dplyr::summarise(
      popN = sum(popN)
    )
  # tidy death data
  deathdat_long <- deathdat_long %>%
    dplyr::filter(Deaths != -1) %>%  # missing value for cpp
    dplyr::group_by_at(groupingvar) %>%
    dplyr::mutate(
      cumdeaths = cumsum(Deaths)
    )
  # combine
  ret <- dplyr::left_join(deathdat_long, popdat, by = groupingvar) %>%
    dplyr::filter(Deaths != -1) %>%  # missing value for cpp
    dplyr::mutate(std_deaths = (Deaths/popN) * Nstandardization,
                  std_cum_deaths = (cumdeaths/popN) * Nstandardization) %>%
    magrittr::set_colnames(tolower(colnames(.)))

  # for plot
  if (groupingvar == "ageband") {
    ret <- ret %>%
      dplyr::mutate(age_mid = purrr::map_dbl(ageband, get_mid_age))
  } else {
    ret$age_mid = NA
  }

  # make sure factor order preserved which can be overwritten in the group_by_at, so for safety
  if (groupingvar == "ageband") {
    ret <- ret %>%
      dplyr::mutate(age_low = as.numeric(stringr::str_extract(ageband, "[0-9]+(?=\\,)"))) %>%
      dplyr::arrange(age_low) %>%
      dplyr::mutate(ageband = forcats::fct_reorder(.f = ageband, .x = age_low)) %>%
      dplyr::select(-c("age_low"))
  }

  # out
  return(ret)
}


