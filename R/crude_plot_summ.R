#' @title
get_orig_seroprev <- function(IFRmodinput, groupingvar) {
  rogan_gladen <- function(obs_prev, sens, spec){(obs_prev + spec -1)/(spec+sens-1)}
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



