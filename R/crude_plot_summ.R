source("R/assertions_v5.R")
#' @title Quick Jpeg
jpgsnapshot <- function(outpath, plot, type = "wide") {
  assert_in(type, c("long", "wide"))
  if (type == "long") {
    jpeg(outpath, width = 8, height = 11, units = "in", res = 500)
    plot(plot)
    graphics.off()
  } else if (type == "wide") {
    jpeg(filename = outpath, width = 11, height = 8, units = "in", res = 500)
    plot(plot)
    graphics.off()
  }
}

#' @title Rogan-Gladen Estimator for Correct Prevalence Obersvations
rogan_gladen <- function(obs_prev, sens, spec){
  ret <- (obs_prev + spec - 1)/(spec + sens - 1)
  ret <- ifelse(ret <= 0, NA, ret) # occurs when obs_prev and spec are less than 1
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
                  seromidpt = round(mean(c(obsdaymin, obsdaymax))))
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
    dplyr::mutate(std_deaths = (Deaths/popN) * Nstandardization,
                  std_cum_deaths = (cumdeaths/popN) * Nstandardization) %>%
    magrittr::set_colnames(tolower(colnames(.)))

  # for plot
  if (groupingvar == "ageband") {
    ret <- ret %>%
      dplyr::mutate(age_mid = purrr::map_dbl(ageband, function(x){
        nums <- as.numeric(stringr::str_split_fixed(x, "-", n = 2))
        nums[nums == 999] <- 95
        return(mean(nums))})
      )
  } else {
    ret$age_mid = NA
  }

  # out
  return(ret)
}

#' @title
#'
get_crude_summarydf <- function(IFRmodinput, groupingvar) {
  rogan_gladen <- function(obs_prev, sens, spec){(obs_prev + spec -1)/(spec+sens-1)}
  #......................
  # setup new df
  #......................
  if (length(IFRmodinput$seroprev_group_adj) != 0) {
    # out
    ret <- dplyr::left_join(IFRmodinput$seroprev_group_adj,
                            IFRmodinput$prop_pop, by = groupingvar) %>%
      dplyr::left_join(., IFRmodinput$deaths_group, by = groupingvar) %>%
      dplyr::mutate(
        age_mid = purrr::map_dbl(ageband, function(x){
          nums <- as.numeric(stringr::str_split_fixed(x, "-", n = 2))
          nums[nums == 999] <- 95
          return(mean(nums))
        }),
        inf_pop_crude = seroprevalence_adj * popN,
        ifr_age_crude = deaths_at_sero/inf_pop_crude,
        deaths_per_pop = deaths_at_sero/popN,
        tot_deaths_per_pop = sum(deaths_per_pop),
        prop_deaths_per_pop = deaths_per_pop / tot_deaths_per_pop,
        deaths_per_million = (1e6*deaths_at_sero)/popN
      ) %>%
      dplyr::rename(seroprev = seroprevalence_adj) %>%
      dplyr::select(c("age_mid", "seroprev", "inf_pop_crude", "ifr_age_crude",
                      "deaths_per_pop", "prop_deaths_per_pop", "deaths_per_million"))
    if (length(IFRmodinput$sero_spec$specificity) > 0 & length(IFRmodinput$sero_sens$sensitivity) > 0) {
      ret <- ret %>%
        dplyr::mutate(seroprev_adj_ss = purrr::map_dbl(seroprev, rogan_gladen,
                                                       spec = IFRmodinput$sero_spec$specificity,
                                                       sens = IFRmodinput$sero_sens$sensitivity))
    }

  } else {
    ret <- dplyr::left_join(IFRmodinput$seroprev_group,
                            IFRmodinput$prop_pop, by = groupingvar) %>%
      dplyr::left_join(., IFRmodinput$deaths_group, by = groupingvar) %>%
      dplyr::mutate(
        age_mid = purrr::map_dbl(ageband, function(x){
          nums <- as.numeric(stringr::str_split_fixed(x, "-", n = 2))
          nums[nums == 999] <- 95
          return(mean(nums))
        }),
        inf_pop_crude = seroprevalence * popN,
        ifr_age_crude = deaths_at_sero/inf_pop_crude,
        deaths_per_pop = deaths_at_sero/popN,
        tot_deaths_per_pop = sum(deaths_per_pop),
        prop_deaths_per_pop = deaths_per_pop / tot_deaths_per_pop,
        deaths_per_million = (1e6*deaths_at_sero)/popN
      ) %>%
      dplyr::rename(seroprev = seroprevalence) %>%
      dplyr::select(c("age_mid", "seroprev", "inf_pop_crude", "ifr_age_crude",
                      "deaths_per_pop", "prop_deaths_per_pop", "deaths_per_million"))

    if (length(IFRmodinput$sero_spec$specificity) > 0 & length(IFRmodinput$sero_sens$sensitivity) > 0) {
      ret <- ret %>%
        dplyr::mutate(seroprev_adj_ss = purrr::map_dbl(seroprev, rogan_gladen,
                                                       spec = IFRmodinput$sero_spec$specificity,
                                                       sens = IFRmodinput$sero_sens$sensitivity))
    }
  }
  return(ret)
}



