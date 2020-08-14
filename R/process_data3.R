


#' @title Process Data for IFR Inference
#' @param deaths dataframe; death counts by strata
#' @param population dataframe; population counts by strata
#' @param seroval dataframe; serovalidation data read from the serovalidation file
#' @param seroprev ddataframe; seroprevalence data read from the seroprevalence file
#' @param cumulative logical; Are the deaths cumulative to a given data or time-series
#' @param recast_deaths_df ddataframe; Daily deaths counts used for recasting cumulative death proportions -- only evaluated if \code{cumulative} is TRUE. Columns names are date, georegion, and deaths, where deaths are an incident count and not a cumulative count
#' @param recast_deaths_geocode character;Geocode to be subset the recast_deaths_df (filtered on the georegion column in the recast data dataframe)
#' @param groupingvar character; name of strata(s) being considered
#' @param study_ids character; Study id to be considered
#' @param filtRegion character; region levels to keep
#' @param filtGender character; biological sex levels to keep
#' @param filtAgeBand character; age groups to keep -- note this will be a factor of the age_low and age_high concatenated together
#' @param death_agebreaks character; potential to customize the break points for the factorization of ages from the death data. Default NULL will use data to set breaks
#' @param sero_agebreaks character; potential to customize the break points for the factorization of ages from the death data. Default NULL will use data to set breaks
#' @param origin date; Day 1 to process all data relative to
#' @import tidyverse
#' NB, this isn't a package, so ^^ is just a reminder to users to have tidyverse loaded in order to allow embracing and piping to work as expected

source("R/assertions_v5.R")
process_data3 <- function(deaths = NULL, population = NULL, sero_val = NULL, seroprev = NULL,
                          cumulative = FALSE, recast_deaths_df = NULL,
                          groupingvar, study_ids, recast_deaths_geocode,
                          filtRegions = NULL, filtGender = NULL, filtAgeBand = NULL, death_agebreaks = NULL,
                          sero_agebreaks = NULL, origin = lubridate::ymd("2020-01-01")) {
  #......................
  # assertions and checks
  #......................
  assert_dataframe(deaths)
  assert_dataframe(population)
  assert_dataframe(sero_val)
  assert_dataframe(seroprev)
  assert_logical(cumulative)
  assert_string(groupingvar)
  assert_in(groupingvar, c("region", "ageband", "gender"))
  assert_string(study_ids)
  if (cumulative){
    assert_dataframe(recast_deaths_df)
    assert_string(recast_deaths_geocode)
  }
  if(!is.null(filtRegions)){
    assert_string(filtRegions)
  }
  if(!is.null(filtGender)){
    assert_string(filtGender)
  }
  if(!is.null(filtAgeBand)){
    assert_string(filtAgeBand)
  }

  # check columns for population df
  assert_in(c("country", "study_id", "age_low", "age_high", "region", "gender", "population", "age_breakdown", "for_regional_analysis", "gender_breakdown"),
            colnames(population))
  # check columns for deaths df and dates
  assert_in(c("country", "study_id", "age_low", "age_high", "region", "gender", "n_deaths", "date_start_survey", "date_end_survey", "age_breakdown", "for_regional_analysis", "gender_breakdown"),
            colnames(deaths))
  chckdf <- deaths %>%
    dplyr::filter(study_id %in% study_ids)
  assert_date(chckdf$date_start_survey, message = "Deaths date_start_survey is not in lubridate format \n")
  assert_date(chckdf$date_end_survey, message = "Deaths date_start_survey is not in lubridate format \n")

  # check columns for seroval
  assert_in(c("study_id", "sensitivity", "specificity"),
            colnames(sero_val))
  # check columns for seroprev df and dates
  assert_in(c("country", "study_id", "age_low", "age_high", "region", "gender", "seroprevalence_unadjusted",
              "seroprevalence_weighted", "n_tested", "n_positive", "date_start_survey", "date_end_survey",
              "age_breakdown", "for_regional_analysis", "gender_breakdown"),
            colnames(seroprev))
  chckdf <- seroprev %>%
    dplyr::filter(study_id %in% study_ids)
  assert_date(chckdf$date_start_survey, message = "Seroprevalence date_start_survey is not in lubridate format \n")
  assert_date(chckdf$date_end_survey, message = "Seroprevalence date_end_survey is not in lubridate format \n")

  # check ecdc
  if (cumulative){
    assert_in(c("date", "georegion", "deaths"),
              colnames(recast_deaths_df))
  }

  # check daily time steps are daily -- date start survey in deaths purely for sanity check
  if (!cumulative) {
    assert_eq(deaths$date_start_survey, deaths$date_end_survey)
  }

  #............................................................
  # process death data
  #...........................................................
  deaths_summ <- process_death_data(deaths, recast_deaths_df, recast_deaths_geocode,
                                    cumulative, study_ids, death_agebreaks, groupingvar,
                                    filtRegions, filtGender, filtAgeBand, origin = origin)

  #............................................................
  # process seroprev data
  #...........................................................
  seroprev_list <- process_seroprev_data(seroprev,
                                         study_ids, sero_agebreaks, groupingvar,
                                         filtRegions, filtGender, filtAgeBand, origin = origin)

  #...........................................................
  # process population
  #...........................................................
  pop_prop.summ <- process_population_data(population, deaths, death_agebreaks, study_ids, groupingvar,
                                           filtRegions, filtGender, filtAgeBand)

  #............................................................
  # process test sens/spec
  #...........................................................
  sero_val <- sero_val %>%
    dplyr::filter(study_id %in% study_ids)
  sensitivity <- sero_val %>%
    dplyr::select(c("n_test_pos_out_of_true_pos", "n_samples_true_pos", "sensitivity")) %>%
    magrittr::set_colnames(c("npos", "ntest", "sensitivity"))
  specificity <- sero_val %>%
    dplyr::select(c("n_test_neg_out_of_true_neg", "n_samples_true_neg", "specificity")) %>%
    magrittr::set_colnames(c("npos", "ntest", "specificity"))



  #...........................................................
  # out
  #...........................................................
  ret <- list(
    deathsMCMC = deaths_summ,
    seroprevMCMC = seroprev_list$seroprevMCMC,
    prop_pop = pop_prop.summ,
    sero_sens = sensitivity,
    sero_spec = specificity,
    seroprev_group = seroprev_list$seroprevFull
  )
  return(ret)

}



#' @title Internal Function for Processing Seroprevalence Data from our Systematic review
process_seroprev_data <- function(seroprev, origin, study_ids, sero_agebreaks, groupingvar,
                                  filtRegions, filtGender, filtAgeBand) {
  #......................
  # tidy up seropred data
  #......................
  seroprev <- seroprev %>%
    dplyr::mutate(date_start_survey = lubridate::ymd(date_start_survey),
                  date_end_survey = lubridate::ymd(date_end_survey),
                  ObsDaymin = as.numeric(date_start_survey - origin) + 1,
                  ObsDaymax = as.numeric(date_end_survey - origin) + 1,
                  n_positive = ifelse(is.na(n_positive) & !is.na(n_tested) & !is.na(seroprevalence_unadjusted),
                                      round(n_tested * seroprevalence_unadjusted), n_positive)
    ) %>%
    dplyr::filter(study_id %in% study_ids)

  #......................
  # Age Process if necessary
  #......................
  if (groupingvar == "ageband") {
    # do any neccesary age filtering
    seroprev <- seroprev %>%
      dplyr::filter(age_breakdown == 1)
    if(!is.null(filtAgeBand)) {
      seroprev <- seroprev %>%
        dplyr::filter(ageband %in% filtAgeBand)
    }
    if (!is.null(sero_agebreaks)) {
      # user provided age breaks
      assert_vector(sero_agebreaks)
      assert_greq(length(sero_agebreaks), 2)
      agebrks_sero <- sero_agebreaks
    } else {
      # if none are provided, factor age to one level
      agebrks_sero <- c(min(seroprev$age_low), sort(unique(seroprev$age_high)))
    }
    # make new age band categories
    seroprev <- seroprev %>%
      dplyr::mutate(
        ageband = cut(age_high,
                      breaks = agebrks_sero,
                      labels = c(paste0(agebrks_sero[1:(length(agebrks_sero)-1)], "-", lead(agebrks_sero)[1:(length(agebrks_sero)-1)]))),
        ageband = as.character(ageband),
        ageband = ifelse(age_low == 0 & age_high == 999, "all", ageband),
        ageband = forcats::fct_reorder(.f = ageband, .x = age_low)
      )
  }

  #......................
  # Region Process (if necessary)
  #......................
  if (groupingvar == "region"){
    seroprev <- seroprev %>%
      dplyr::filter(for_regional_analysis == 1)

    if(!is.null(filtRegions)) {
      seroprev <- seroprev %>%
        dplyr::filter(region %in% filtRegions)
    }
  }

  #......................
  # Gender Process (if necessary)
  #......................
  if (groupingvar == "gender") {
    seroprev <- seroprev %>%
      dplyr::filter(gender_breakdown == 1)

    if(!is.null(filtGender)) {
      seroprev <- seroprev %>%
        dplyr::filter(gender %in% filtGender)
    }
  }

  # sanity checks
  assert_gr(nrow(seroprev), 0, message = "There was a mismatch in filtering seroprevalence observations. Returned no observations")

  # split pieces for MCMC model vs. Descriptive plot
  # NB this group_by and mean/sum if there is only one day will just return the same value
  seroprevMCMC <- seroprev %>%
    dplyr::group_by_at(c("ObsDaymin", "ObsDaymax", groupingvar)) %>%
    dplyr::summarise(n_positive = sum(n_positive),
                     n_tested = sum(n_tested),
                     seroprevalence_unadjusted = n_positive/n_tested) %>%
    dplyr::rename(SeroPrev = seroprevalence_unadjusted) %>%
    dplyr::ungroup()
  # out
  out <- list(seroprevMCMC = seroprevMCMC,
              seroprevFull = seroprev)
  return(out)
}

#' @title Internal Function for Processing Population Data from Demographic Surveys
process_population_data <- function(population, deaths, death_agebreaks, study_ids, groupingvar,
                                    filtRegions, filtGender, filtAgeBand) {

  #......................
  # tidy population up
  #......................
  # subset to study id
  population <- population %>%
    dplyr::filter(study_id %in% study_ids)

  #......................
  # Process Age
  # NB, we always consider age stratification - use same age breaks as deaths.
  #......................
  #......................
  # Age Process
  #......................
  if (!is.null(death_agebreaks)) {
    # user provided age breaks
    assert_vector(death_agebreaks)
    assert_greq(length(death_agebreaks), 2)
    agebrks <- death_agebreaks
  } else {
    # if none are provided, factor age to one level
    agebrks <- c(0, sort(unique(deaths$age_high)))
  }
  # make new age band categories
  population <- population %>%
    dplyr::mutate(
      ageband = cut(age_high,
                    breaks = agebrks,
                    labels = c(paste0(agebrks[1:(length(agebrks)-1)], "-", lead(agebrks)[1:(length(agebrks)-1)]))),
      ageband = as.character(ageband),
      ageband = ifelse(age_low == 0 & age_high == 999, "all", ageband),
      ageband = forcats::fct_reorder(.f = ageband, .x = age_low)
    )
  # do any neccesary age filtering
  if (groupingvar == "ageband") {
    population <- population %>%
      dplyr::filter(age_breakdown == 1)

    if(!is.null(filtAgeBand)) {
      population <- population %>%
        dplyr::filter(ageband %in% filtAgeBand)
    }
  }

  #......................
  # Process Regions (if necessary)
  #......................
  if (groupingvar == "region") { # but always account age stratification
    population <- population %>%
      dplyr::filter(for_regional_analysis == 1)
    if (!is.null(filtRegions)) {
      population <- population %>%
        dplyr::filter(region %in% filtRegions)
    }
  }

  #......................
  # Process Gender (if necessary)
  #......................
  if (groupingvar == "gender") {
    population <- population %>%
      dplyr::filter(gender_breakdown == 1)
    if(!is.null(filtGender)) {
      population <- population %>%
        dplyr::filter(gender %in% filtGender)
    }
  }


  # SANITY CHECKs
  assert_gr(nrow(population), 0,
            message = "There was a mismatch in filtering population observations. Returned no observations")

  #......................
  # liftover to pop demographics
  #......................
  totalpop <- sum(population$population)
  if (groupingvar == "ageband") {
    pop_prop.summ <- population %>%
      dplyr::group_by_at(groupingvar) %>%
      dplyr::summarise(
        popN = sum(population),
        pop_prop = sum(population)/totalpop
      ) %>%
      dplyr::arrange_at(groupingvar) %>%
      dplyr::ungroup()
  } else {
    pop_prop.summ <- population %>%
      dplyr::group_by_at(c(groupingvar, "ageband")) %>%
      dplyr::summarise(
        popN = sum(population),
        pop_prop = sum(population)/totalpop
      ) %>%
      dplyr::arrange_at(c(groupingvar, "ageband")) %>%
      dplyr::ungroup()
  }
  return(pop_prop.summ)
}



#' @title Internal Function for Processing Death Data from COVID Surveys
process_death_data <- function(deaths, recast_deaths_df, recast_deaths_geocode,
                               cumulative, study_ids, death_agebreaks, groupingvar,
                               filtRegions, filtGender, filtAgeBand, origin) {
  #......................
  # tidy up death data
  #......................
  deaths <- deaths %>%
    dplyr::filter(study_id %in% study_ids)


  #......................
  # Age Process
  #......................
  if (!is.null(death_agebreaks)) {
    # user provided age breaks
    assert_vector(death_agebreaks)
    assert_greq(length(death_agebreaks), 2)
    agebrks <- death_agebreaks
  } else {
    # if none are provided, factor age to one level
    agebrks <- c(0, sort(unique(deaths$age_high)))
  }

  # make new age band categories
  deaths <- deaths %>%
    dplyr::mutate(
      ageband = cut(age_high,
                    breaks = agebrks,
                    labels = c(paste0(agebrks[1:(length(agebrks)-1)], "-", lead(agebrks)[1:(length(agebrks)-1)]))),
      ageband = as.character(ageband),
      ageband = ifelse(age_low == 0 & age_high == 999, "all", ageband),
      ageband = forcats::fct_reorder(.f = ageband, .x = age_low)
    )

  # do any neccesary age filtering
  if (groupingvar == "ageband") {
    deaths <- deaths %>%
      dplyr::filter(age_breakdown == 1)
    if(!is.null(filtAgeBand)) {
      deaths <- deaths %>%
        dplyr::filter(ageband %in% filtAgeBand)
    }
  }

  #......................
  # Region Process (if necessary)
  #......................
  if (groupingvar == "region") {
    deaths <- deaths %>%
      dplyr::filter(for_regional_analysis == 1)
    if (!is.null(filtRegions)) {
      deaths <- deaths %>%
        dplyr::filter(region %in% filtRegions)
    }
  }

  #......................
  # Gender Process (if necessary)
  #......................
  if (groupingvar == "gender") {
    deaths <- deaths %>%
      dplyr::filter(gender_breakdown == 1)
    if(!is.null(filtGender)) {
      deaths <- deaths %>%
        dplyr::filter(gender %in% filtGender)
    }
  }


  # SANITY CHECK
  assert_gr(nrow(deaths), 0,
            message = "There was a mismatch in filtering death observations. Returned no observations")

  #......................
  # RECAST DEATHS or TIME SERIES
  # if cumulative, recast from recast_deaths_df
  #......................
  if (cumulative){
    recast_deaths_df <- recast_deaths_df %>%
      dplyr::filter(georegion %in% recast_deaths_geocode) %>%
      dplyr::mutate(ObsDay = as.numeric(lubridate::ymd(date) - origin) + 1) %>%
      dplyr::filter(ObsDay >= 1) %>% # consistent origin
      dplyr::arrange(ObsDay)

    # now multiple proportions to get time series
    deaths.prop <- deaths %>%
      dplyr::mutate(date_start_survey = origin,
                    date_end_survey = lubridate::ymd(date_end_survey),
                    ObsDay = max(date_end_survey)) %>%
      dplyr::group_by_at(c(groupingvar, "ObsDay")) %>% # to get deaths at end of shorter survey period
      dplyr::summarise(death_num = sum(n_deaths),
                       age_low = mean(age_low),
                       age_high = mean(age_high)) %>%
      dplyr::arrange(age_low)


    # store group names
    groupvarnames <- group_keys(deaths.prop)

    # get proportion
    deaths.prop <- deaths.prop %>%
      dplyr::ungroup(.) %>%
      dplyr::select(-c("ObsDay")) %>%
      dplyr::mutate(death_denom = sum(death_num),
                    death_prop = death_num/death_denom) # protect against double counting of same person in multiple groups

    # now recast proportions across days equally
    deaths.summ <- as.data.frame(matrix(NA, nrow = nrow(deaths.prop), ncol = max(recast_deaths_df$ObsDay)))
    # start counter after we found study start date
    iter <- 1
    for (i in 1:ncol(deaths.summ)) {
      if (i < min(recast_deaths_df$ObsDay)) { # fill in for common origin
        deaths.summ[,i] <- 0
      } else {
        if (i %in% recast_deaths_df$ObsDay) {
          deaths.summ[,i] <- deaths.prop$death_prop * recast_deaths_df$deaths[iter]
          iter <- iter + 1
        } else { # catch instance where missing days in the middle of the series
          deaths.summ[,i] <- -1
          iter <- iter + 1
        }
      }
    }
    # need to round to nearest person
    deaths.summ <- round(deaths.summ)
    # tidy out to long format
    deaths.summ <- deaths.summ %>%
      cbind.data.frame(groupvarnames, .)
    first_non_group_col <- which(colnames(deaths.summ) == "V1")
    deaths.summ <- deaths.summ %>%
      tidyr::gather(., key = "ObsDay", value = "Deaths", first_non_group_col:ncol(.)) %>%
      dplyr::mutate(ObsDay = as.numeric(gsub("V", "", ObsDay)))

  } else {
    #......................
    # deaths summary out for time series
    #......................

    deaths.summ <- deaths %>%
      dplyr::mutate(date_start_survey = origin,
                    date_end_survey = lubridate::ymd(date_end_survey),
                    ObsDay = as.numeric(date_end_survey - origin) + 1)  %>%
      dplyr::group_by_at(c(groupingvar, "ObsDay")) %>%
      dplyr::summarise( Deaths = sum(n_deaths) ) %>%
      dplyr::arrange_at(c("ObsDay", groupingvar)) %>%
      dplyr::ungroup(.)
  }
  # out
  return(deaths.summ)
}

#' @title Utility function for matching age groups that aren't "exact"
#' @param x query age band
#' @param y target age band
#' @param wiggle how much "wiggle" room on matching
wiggle_age_matchfun <- function(x, y, wiggle = 2) {
  mtchlow <- c(as.numeric(x["age_low"]) + wiggle >= y$age_low | as.numeric(x["age_low"]) - wiggle >= y$age_low)
  mtchhigh <- c(as.numeric(x["age_high"]) + wiggle <= y$age_high | as.numeric(x["age_high"]) - wiggle <= y$age_high)
  # rows
  mtchrows <- mtchlow &  mtchhigh
  if (sum(mtchrows) == 0) {
    return(as.numeric(x["seroprevalence"]))
  } else if (sum(mtchrows) == 1) {
    return(as.numeric(y[mtchrows, "seroprevalence"]))
  } else {
    stop("Error")
  }
}


