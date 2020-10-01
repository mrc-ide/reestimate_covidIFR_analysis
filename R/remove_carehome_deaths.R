
#' @title  Remove care home deaths from older age groups, any containing >65
#' @param agebands.dat object made from process_data4
#' @param study_id character

remove_ch_deaths <-  function(ageband_dat, carehomesdf, studyid) {
  #......................
  # first figure out proportion of deaths
  #......................
  # we have cumulative deaths up to a given day for all age bands
  deaths_propMCMC_adj <- ageband_dat$deaths_propMCMC

  # we then have cumulative care home deaths for a given day (likely a different day than our given day above for age bands)
  # but we also have the proportion of all deaths that were in care homes in the population
  deaths_frac_ch <- carehomesdf %>%
    dplyr::filter(study_id == studyid) %>%
    dplyr::pull(percent_deaths)/100
  # because we assume a constant proportion through time, if we assume carehome deahts are only in 65+,
  # we can recast the proportions:
  deaths_propMCMC_adj <- deaths_propMCMC_adj %>%
    dplyr::mutate(
      age_high = as.numeric(stringr::str_split_fixed(ageband, "-", n = 2)[,2]),
      ageband = ifelse(age_high >= 65, "65-999", ageband)) %>%
    dplyr::group_by(ageband) %>%
    dplyr::summarise(death_num = sum(death_num))
  # deaths in care homes
  chdeaths_num <- round(deaths_propMCMC_adj$death_num[deaths_propMCMC_adj$ageband == "65-999"] * deaths_frac_ch)

  # carehome deaths
  chdeaths <- tibble::tibble(ageband = "carehomes",
                            death_num = chdeaths_num)

  # bring together
  deaths_propMCMC_adj <- dplyr::bind_rows(deaths_propMCMC_adj, chdeaths) %>%
      dplyr::mutate(death_num = ifelse(ageband == "65-999", death_num - chdeaths_num, death_num)) %>%
      dplyr::mutate(death_denom = sum(death_num),
                    death_prop = death_num/sum(death_num))

  #......................
  # now recast deaths from timeseries
  #......................
  # now recast proportions across days equally
  recast_deaths <- as.data.frame(matrix(NA, nrow = nrow(deaths_propMCMC_adj),
                                        ncol = max(ageband_dat$deaths_TSMCMC$ObsDay)))
  # know deaths are aligned and no missing days in middle of time-series (resolved on first pass)
  for (i in 1:ncol(recast_deaths)) {
    # need to round to nearest person
    recast_deaths[,i] <- round(ageband_dat$deaths_TSMCMC$deaths[i] * deaths_propMCMC_adj$death_prop)
  }
  # drop carehomes deaths from time series
  agenames <- deaths_propMCMC_adj$ageband[1:(length(deaths_propMCMC_adj$ageband)-1)]
  recast_deaths <- cbind.data.frame(agenames, recast_deaths[1:(nrow(recast_deaths)-1), ])
  # tidy out to long format
  recast_deaths <- recast_deaths %>%
    tidyr::pivot_longer(., cols = -c("agenames"), names_to = "ObsDay", values_to = "Deaths") %>%
    dplyr::mutate(ObsDay = as.numeric(gsub("V", "", ObsDay))) %>%
    dplyr::group_by(ObsDay) %>%
    dplyr::summarise(Deaths = sum(Deaths))

  #............................................................
  # bring together
  #...........................................................
  # tidy up prop deaths
  deaths_propMCMC_adj <- deaths_propMCMC_adj %>%
    dplyr::filter(ageband != "carehomes") %>%
    dplyr::mutate(death_denom = sum(death_num),
                  death_prop = death_num/death_denom) %>%
    dplyr::mutate(age_high = as.numeric(stringr::str_split_fixed(ageband, "-", n=2)[,2])) %>%
    dplyr::arrange(age_high) %>%
    dplyr::select(-c("age_high"))

  # tidy up population df
  pop_adj <- ageband_dat$prop_pop %>%
    dplyr::mutate(age_high = as.numeric(stringr::str_split_fixed(ageband, "-", n=2)[,2]),
                  ageband = ifelse(age_high >= 65, "65-999", ageband)) %>%
    dplyr::group_by(ageband) %>%
    dplyr::summarise(popN = sum(popN)) %>%
    dplyr::mutate(pop_prop = popN/sum(popN)) %>%
    dplyr::mutate(age_high = as.numeric(stringr::str_split_fixed(ageband, "-", n=2)[,2])) %>%
    dplyr::arrange(age_high) %>%
    dplyr::select(-c("age_high"))

  # tidy up seroprevalence df
  if (any(is.na(ageband_dat$seroprevMCMC$n_positive))) {
    seroprev_adj <- ageband_dat$seroprevMCMC %>%
      dplyr::mutate(age_high = as.numeric(stringr::str_split_fixed(ageband, "-", n=2)[,2]),
                    ageband = ifelse(age_high >= 65, "65-999", ageband)) %>%
      dplyr::group_by(ObsDaymin, ObsDaymax, ageband) %>%
      dplyr::summarise(SeroPrev = mean(SeroPrev),
                       SeroLCI = mean(SeroLCI),
                       SeroUCI = mean(SeroUCI)) %>%
      dplyr::mutate(n_positive = NA,
                    n_tested = NA) %>%
      dplyr::ungroup(.) %>%
      dplyr::arrange(ObsDaymin, ObsDaymax, ageband)
  } else {
    seroprev_adj <- ageband_dat$seroprevMCMC %>%
      dplyr::mutate(age_high = as.numeric(stringr::str_split_fixed(ageband, "-", n=2)[,2]),
                    ageband = ifelse(age_high >= 65, "65-999", ageband)) %>%
      dplyr::group_by(ObsDaymin, ObsDaymax, ageband) %>%
      dplyr::summarise(n_positive = sum(n_positive),
                       n_tested = sum(n_tested)) %>%
      dplyr::mutate(SeroPrev = n_positive/n_tested,
                    SeroLCI = NA,
                    SeroUCI = NA) %>%
      dplyr::ungroup(.) %>%
      dplyr::arrange(ObsDaymin, ObsDaymax, ageband)
  }



  # overwrite
  agebands_noCH_dat <- ageband_dat
  agebands_noCH_dat$deaths_TSMCMC <- recast_deaths
  agebands_noCH_dat$deaths_propMCMC <- deaths_propMCMC_adj
  agebands_noCH_dat$seroprevMCMC <- seroprev_adj
  agebands_noCH_dat$prop_pop <- pop_adj
  agebands_noCH_dat$deaths_group <- NULL
  agebands_noCH_dat$seroprev_group <- NULL
  # out
  return(agebands_noCH_dat)
}
