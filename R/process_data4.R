
#' @title Process Data for IFR Inference
#' @param cum_tp_deaths dataframe; cumulative deaths at a particular time point in the epidemic
#' @param population dataframe; population counts by strata
#' @param seroval dataframe; serovalidation data read from the serovalidation file
#' @param seroprev dataframe; seroprevalence data read from the seroprevalence file
#' @param time_series_totdeaths_df dataframe; Daily deaths counts used for MCMC fit. Will also be used for recasting cumulative death proportions -- only evaluated if \code{get_descriptive_dat} is TRUE -- for descriptive plotting. Columns names are date, georegion, and deaths, where deaths are an incident count and not a cumulative count
#' @param time_series_totdeaths_geocode character; Geocode to be subset the time series deaths df (filtered on the georegion column in the recast data dataframe)
#' @param get_descriptive_dat logical; Whether to recast deaths for descriptive plotting
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
process_data4 <- function(cum_tp_deaths = NULL, population = NULL, sero_val = NULL, seroprev = NULL,
                          get_descriptive_dat = TRUE, time_series_totdeaths_df = NULL,
                          groupingvar, study_ids, time_series_totdeaths_geocode,
                          filtRegions = NULL, filtGender = NULL, filtAgeBand = NULL,
                          agebreaks = NULL,
                          origin = lubridate::ymd("2020-01-01")) {
  #......................
  # assertions and checks
  #......................
  assert_dataframe(cum_tp_deaths)
  assert_dataframe(population)
  assert_dataframe(sero_val)
  assert_dataframe(seroprev)
  assert_logical(get_descriptive_dat)
  assert_string(groupingvar)
  assert_in(groupingvar, c("region", "ageband", "gender"))
  assert_string(study_ids)
  assert_dataframe(time_series_totdeaths_df)
  assert_string(time_series_totdeaths_geocode)
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
  assert_in(c("country", "study_id", "age_low", "age_high", "region", "gender", "population", "age_breakdown", "for_regional_analysis"),
            colnames(population))
  # check columns for deaths df and dates
  assert_in(c("country", "study_id", "age_low", "age_high", "region", "gender", "n_deaths", "age_breakdown", "for_regional_analysis"),
            colnames(cum_tp_deaths))
  chckdf <- cum_tp_deaths %>%
    dplyr::filter(study_id %in% study_ids)
  assert_greq(nrow(chckdf), 1, message = "Study ID not found in the Cumulative TimePoint Dataframe")

  # check columns for seroval
  assert_in(c("study_id", "n_samples_true_neg", "n_test_neg_out_of_true_neg", "n_samples_true_pos", "n_test_pos_out_of_true_pos"),
            colnames(sero_val))
  # check columns for seroprev df and dates
  assert_in(c("country", "study_id", "age_low", "age_high", "region", "gender", "seroprevalence_unadjusted",
              "seroprevalence_weighted", "n_tested", "n_positive", "date_start_survey", "date_end_survey",
              "age_breakdown", "for_regional_analysis"),
            colnames(seroprev))
  chckdf <- seroprev %>%
    dplyr::filter(study_id %in% study_ids)
  assert_date(chckdf$date_start_survey, message = "Seroprevalence date_start_survey is not in lubridate format \n")
  assert_date(chckdf$date_end_survey, message = "Seroprevalence date_end_survey is not in lubridate format \n")

  # check time series
  assert_eq(c("date", "georegion", "deaths"), colnames(time_series_totdeaths_df))
  assert_date(time_series_totdeaths_df$date, message = "Deaths Time Series Date is not in lubridate format \n")
  checkdf <- time_series_totdeaths_df %>%
    dplyr::filter(georegion == time_series_totdeaths_geocode)
  assert_increasing(as.numeric(checkdf$date), message = "Dates not monotonically increasing")

  #............................................................
  # process death data
  #...........................................................
  deaths_list <- process_death_data(cum_tp_deaths, time_series_totdeaths_df, time_series_totdeaths_geocode,
                                    get_descriptive_dat, study_ids, agebreaks, groupingvar,
                                    filtRegions, filtGender, filtAgeBand, origin = origin)

  #............................................................
  # process seroprev data
  #...........................................................
  seroprev_list <- process_seroprev_data(seroprev,
                                         study_ids, agebreaks, groupingvar,
                                         filtRegions, filtGender, filtAgeBand, origin = origin)

  #...........................................................
  # process population
  #...........................................................
  pop_prop.summ <- process_population_data(population, cum_tp_deaths, agebreaks, study_ids, groupingvar,
                                           filtRegions, filtGender, filtAgeBand)

  #............................................................
  # process test sens/spec
  #...........................................................
  sero_val <- sero_val %>%
    dplyr::filter(study_id %in% study_ids)
  sensitivity <- sero_val %>%
    dplyr::select(c("n_test_pos_out_of_true_pos", "n_samples_true_pos")) %>%
    magrittr::set_colnames(c("npos", "ntest")) %>%
    dplyr::mutate(sensitivity = npos/ntest)
  specificity <- sero_val %>%
    dplyr::select(c("n_test_neg_out_of_true_neg", "n_samples_true_neg")) %>%
    magrittr::set_colnames(c("npos", "ntest")) %>%
    dplyr::rename(nneg = npos) %>%  # test negatives (positive result but for clarity)
    dplyr::mutate(specificity = nneg/ntest)


  #...........................................................
  # out
  #...........................................................
  ret <- list(
    deaths_TSMCMC = deaths_list$deaths_TSMCMC,
    deaths_propMCMC = deaths_list$deaths_propMCMC,
    seroprevMCMC = seroprev_list$seroprevMCMC,
    prop_pop = pop_prop.summ,
    sero_sens = sensitivity,
    sero_spec = specificity,
    deaths_group = deaths_list$deaths_group,
    seroprev_group = seroprev_list$seroprevFull
  )
  return(ret)

}



#' @title Internal Function for Processing Seroprevalence Data from our Systematic review
#' @inheritParams process_data4

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
      # if none are provided, factor based on data
      agebrks_sero <- c(min(seroprev$age_low), sort(unique(seroprev$age_high)))
    }
    # make new age band categories
    seroprev <- seroprev %>%
      dplyr::mutate(
        ageband = cut(age_high,
                      breaks = agebrks_sero,
                      labels = c(paste0(agebrks_sero[1:(length(agebrks_sero)-1)], "-", lead(agebrks_sero)[1:(length(agebrks_sero)-1)])),
                      ageband = as.character(ageband))

      ) %>%
      dplyr::arrange(age_low)
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
    dplyr::ungroup() %>%
    dplyr::select(c(groupingvar, dplyr::everything()))

  # make sure age order perserved which can be overwritten in the group_by_at, so for safety
  if (groupingvar == "ageband") {
    seroprevMCMC <- seroprevMCMC %>%
      dplyr::mutate(age_low = as.numeric(stringr::str_extract(ageband, "[0-9]+?(?=-)"))) %>%
      dplyr::arrange(age_low) %>%
      dplyr::select(-c("age_low"))
  }


  # out
  out <- list(seroprevMCMC = seroprevMCMC,
              seroprevFull = seroprev)
  return(out)
}

#' @title Internal Function for Processing Population Data from Demographic Surveys
#' @inheritParams process_data4

process_population_data <- function(population, cum_tp_deaths, death_agebreaks, study_ids, groupingvar,
                                    filtRegions, filtGender, filtAgeBand) {

  #......................
  # tidy population up
  #......................
  # subset to study id
  population <- population %>%
    dplyr::filter(study_id %in% study_ids)
  cum_tp_deaths <- cum_tp_deaths %>%
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
    # if none are provided, factor age based on data
    agebrks <- c(0, sort(unique(cum_tp_deaths$age_high)))
  }
  # make new age band categories
  population <- population %>%
    dplyr::mutate(
      ageband = cut(age_high,
                    breaks = agebrks,
                    labels = c(paste0(agebrks[1:(length(agebrks)-1)], "-", lead(agebrks)[1:(length(agebrks)-1)]))),
      ageband = as.character(ageband))

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

  # make sure  order perserved which can be overwritten in the group_by_at, so for safety
  if (groupingvar == "ageband") {
    pop_prop.summ <- pop_prop.summ %>%
      dplyr::mutate(age_low = as.numeric(stringr::str_extract(ageband, "[0-9]+?(?=-)"))) %>%
      dplyr::arrange(age_low) %>%
      dplyr::select(-c("age_low"))
  }

  # out
  return(pop_prop.summ)
}



#' @title Internal Function for Processing Death Data from COVID Surveys
#' @inheritParams process_data4

process_death_data <- function(cum_tp_deaths, time_series_totdeaths_df, time_series_totdeaths_geocode,
                               get_descriptive_dat, study_ids, death_agebreaks, groupingvar,
                               filtRegions, filtGender, filtAgeBand, origin) {

  #......................
  # Data Piece 1
  # get time series total deaths as MCMC data input
  #......................
  time_series_totdeaths_df <- time_series_totdeaths_df %>%
    dplyr::filter(georegion %in% time_series_totdeaths_geocode) %>%
    dplyr::mutate(ObsDay = as.numeric(lubridate::ymd(date) - origin) + 1) %>%
    dplyr::filter(ObsDay >= 1) %>% # consistent origin
    dplyr::arrange(ObsDay) %>%
    dplyr::select(-c("date"))

  # catch
  assert_greq(nrow(time_series_totdeaths_df), 1, message = "There was a mismatch in filtering the Time-Series data")


  if (min(time_series_totdeaths_df$ObsDay) != 1) {
    backfilldf <- tibble::tibble(ObsDay = 1:(min(time_series_totdeaths_df$ObsDay)-1),
                                 georegion = time_series_totdeaths_geocode,
                                 deaths = 0,
    )
    time_series_totdeaths_df <- dplyr::bind_rows(backfilldf, time_series_totdeaths_df)
  }


  #......................
  # tidy up cumulative death data for proportions as MCMC data input (data piece 2)
  #......................
  cum_tp_deaths <- cum_tp_deaths %>%
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
    # if none are provided, factor based on age in data
    agebrks <- c(0, sort(unique(cum_tp_deaths$age_high)))

  }

  # make new age band categories
  cum_tp_deaths <- cum_tp_deaths %>%
    dplyr::mutate(
      ageband = cut(age_high,
                    breaks = agebrks,
                    labels = c(paste0(agebrks[1:(length(agebrks)-1)], "-", lead(agebrks)[1:(length(agebrks)-1)]))),
      ageband = as.character(ageband)
    )

  # do any neccesary age filtering
  if (groupingvar == "ageband") {
    cum_tp_deaths <- cum_tp_deaths %>%
      dplyr::filter(age_breakdown == 1)
    if(!is.null(filtAgeBand)) {
      cum_tp_deaths <- cum_tp_deaths %>%
        dplyr::filter(ageband %in% filtAgeBand)
    }
  }

  #......................
  # Region Process (if necessary)
  #......................
  if (groupingvar == "region") {
    cum_tp_deaths <- cum_tp_deaths %>%
      dplyr::filter(for_regional_analysis == 1)
    if (!is.null(filtRegions)) {
      cum_tp_deaths <- cum_tp_deaths %>%
        dplyr::filter(region %in% filtRegions)
    }
  }

  #......................
  # Gender Process (if necessary)
  #......................
  if (groupingvar == "gender") {
    cum_tp_deaths <- cum_tp_deaths %>%
      dplyr::filter(gender_breakdown == 1)
    if(!is.null(filtGender)) {
      cum_tp_deaths <- cum_tp_deaths %>%
        dplyr::filter(gender %in% filtGender)
    }
  }

  # SANITY CHECK
  assert_gr(nrow(cum_tp_deaths), 0,
            message = "There was a mismatch in filtering death observations. Returned no observations")

  #......................
  # Data Piece 2
  # Extract death data for proportions as MCMC data input
  #......................
  cum_tp_deaths.prop <- cum_tp_deaths %>%
    dplyr::group_by_at(c(groupingvar)) %>%
    dplyr::summarise(death_num = sum(n_deaths)) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(death_denom = sum(death_num),
                  death_prop = death_num/death_denom) # protect against double counting of same person in multiple groups

  # make sure age order perserved which can be overwritten in the group_by_at, so for safety
  if (groupingvar == "ageband") {
    cum_tp_deaths.prop <- cum_tp_deaths.prop %>%
      dplyr::mutate(age_low = as.numeric(stringr::str_extract(ageband, "[0-9]+?(?=-)"))) %>%
      dplyr::arrange(age_low) %>%
      dplyr::select(-c("age_low"))
  }
  # store group names
  groupvarnames <- cum_tp_deaths.prop %>%
    dplyr::group_by_at(c(groupingvar))  %>%
    dplyr::select(groupingvar) %>%
    dplyr::ungroup()

  if (get_descriptive_dat){
    # consider whether to get multiple proportions to get time series for descriptive data
    # now recast proportions across days equally
    cum_tp_deaths.summ <- as.data.frame(matrix(NA, nrow = nrow(cum_tp_deaths.prop), ncol = max(time_series_totdeaths_df$ObsDay)))
    # start counter after we found study start date
    iter <- 1
    for (i in 1:ncol(cum_tp_deaths.summ)) {
      if (i %in% time_series_totdeaths_df$ObsDay) {
        cum_tp_deaths.summ[,i] <- cum_tp_deaths.prop$death_prop * time_series_totdeaths_df$deaths[iter]
        iter <- iter + 1
      } else { # catch instance where missing days in the middle of the series
        cum_tp_deaths.summ[,i] <- -1
        iter <- iter + 1
      }
    }
    # need to round to nearest person
    #cum_tp_deaths.summ <- round(cum_tp_deaths.summ)
    # tidy out to long format
    cum_tp_deaths.summ <- cum_tp_deaths.summ %>%
      cbind.data.frame(groupvarnames, .)
    first_non_group_col <- which(colnames(cum_tp_deaths.summ) == "V1")
    cum_tp_deaths.summ <- cum_tp_deaths.summ %>%
      tidyr::gather(., key = "ObsDay", value = "Deaths", first_non_group_col:ncol(.)) %>%
      dplyr::mutate(ObsDay = as.numeric(gsub("V", "", ObsDay)))

  } else {
    cum_tp_deaths.summ <- NULL
  }

  # out
  ret <- list(
    deaths_group = cum_tp_deaths.summ,
    deaths_propMCMC = cum_tp_deaths.prop,
    deaths_TSMCMC = time_series_totdeaths_df)
  return(ret)
}


#' @title  Remove care home deaths from older age groups, any containing >65
#' @param agebands.dat object made from process_data4
#' @param study_id character

remove_ch_deaths <- function(agebands.dat, study_id) {
  # get new propMCMC or deaths
  deaths_propMCMC_adj <- agebands.dat$deaths_propMCMC %>%
    dplyr::mutate(age_low = as.numeric(stringr::str_extract(ageband, "[0-9]+?(?=-)")),
                  age_high = as.numeric(gsub( "(.*)-(.*)", "\\2",  ageband)))
  # find individuals over 65
  inds <- which(deaths_propMCMC_adj$age_high > 65)
  lower_age_old <- min(deaths_propMCMC_adj$age_low[inds])
  deaths_propMCMC_adj$ageband2 <- deaths_propMCMC_adj$ageband
  deaths_propMCMC_adj$ageband2[inds] <- paste0(lower_age_old,"-999")
  new_age <- deaths_propMCMC_adj %>%
    select(ageband, ageband2)  %>%
    dplyr::mutate(age_low = as.numeric(stringr::str_split_fixed(ageband, "-", n = 2)[,1])) %>% # need ages in correct order
    dplyr::arrange(age_low) %>%
    dplyr::select(-c("age_low"))

  ## save new agebands and adjust
  inds <- which(deaths_ch$study_id==study_id)
  deaths_propMCMC_adj <- deaths_propMCMC_adj %>%
    dplyr::group_by(ageband2) %>%
    dplyr::summarise(death_num = sum(death_num),
                     death_denom = mean(death_denom)) %>%
    dplyr::mutate(ch_deaths = round(death_denom*deaths_ch$percent_deaths[inds]/100),
                  death_prop=NA,
                  prop_chd_oldest=ifelse(ageband2 == new_age$ageband2[nrow(new_age)],ch_deaths/death_num,0),
                  death_num=ifelse(ageband2 == new_age$ageband2[nrow(new_age)],death_num-ch_deaths,death_num),
                  death_denom_noch=death_denom - ch_deaths,
                  death_prop=death_num/death_denom_noch) %>%
    dplyr::rename(ageband=ageband2)  %>%
    dplyr::mutate(age_low = as.numeric(stringr::str_split_fixed(ageband, "-", n = 2)[,1])) %>% # need ages in correct order
    dplyr::arrange(age_low) %>%
    dplyr::select(-c("age_low"))

  # overwrite w/ new age groups
  agebands_noCH.dat <- agebands.dat
  agebands_noCH.dat$deaths_propMCMC <- deaths_propMCMC_adj

  # adjust group object for desc data now as well
  deaths_group_adj <- full_join(agebands.dat$deaths_group, new_age, by="ageband")
  high_ageband<-new_age$ageband2[nrow(new_age)]
  deaths_group_adj <- deaths_group_adj %>%
    dplyr::group_by(ageband2,ObsDay) %>%
    dplyr::summarise(Deaths=sum(Deaths)) %>%
    dplyr::mutate(Deaths=ifelse(ageband2 == high_ageband,
                                Deaths * (1-deaths_propMCMC_adj$prop_chd_oldest[nrow(deaths_propMCMC_adj)]),
                                Deaths)) %>%
    dplyr::rename(ageband=ageband2) %>%
    dplyr::mutate(age_low = as.numeric(stringr::str_split_fixed(ageband, "-", n = 2)[,1])) %>% # need ages in correct order
    dplyr::arrange(age_low) %>%
    dplyr::select(-c("age_low"))
  agebands_noCH.dat$deaths_group <- deaths_group_adj

  # adjust population
  pop_adj <- full_join(agebands.dat$prop_pop, new_age,by="ageband")
  pop_adj <- pop_adj %>%
    dplyr::group_by(ageband2) %>%
    dplyr::summarise(popN = sum(popN),
                     pop_prop = sum(pop_prop)) %>%
    dplyr::rename(ageband = ageband2) %>%
    dplyr::mutate(age_low = as.numeric(stringr::str_split_fixed(ageband, "-", n = 2)[,1])) %>% # need ages in correct order
    dplyr::arrange(age_low) %>%
    dplyr::select(-c("age_low"))
  agebands_noCH.dat$prop_pop<-pop_adj

  # seroprevalence adjustment
  seroprevMCMC_adj <- full_join(agebands.dat$seroprevMCMC, new_age, by="ageband")
  seroprevMCMC_adj <- seroprevMCMC_adj %>%
    dplyr::filter(!duplicated(cbind(ObsDaymin, ObsDaymax, SeroPrev, ageband2)))  %>%
    dplyr::group_by(ObsDaymin, ObsDaymax, ageband2) %>%
    dplyr::summarise(n_positive = sum(n_positive),
                     n_tested = sum(n_tested)) %>%
    dplyr::mutate(SeroPrev = n_positive/n_tested) %>%
    dplyr::rename(ageband = ageband2) %>%
    dplyr::mutate(age_low = as.numeric(stringr::str_split_fixed(ageband, "-", n = 2)[,1])) %>% # need ages in correct order
    dplyr::arrange(ObsDaymin, ObsDaymax, age_low) %>%
    dplyr::select(-c("age_low"))
  agebands_noCH.dat$seroprevMCMC <- seroprevMCMC_adj

  # out
  return(agebands_noCH.dat)
}
