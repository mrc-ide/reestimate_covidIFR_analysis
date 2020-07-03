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
#' @import tidyverse
#' NB, this isn't a package, so ^^ is just a reminder to users to have tidyverse loaded in order to allow embracing and piping to work as expected

source("R/assertions_v5.R")
process_data2 <- function(deaths = NULL, population = NULL, sero_val = NULL, seroprev = NULL,
                          cumulative = FALSE, recast_deaths_df = NULL,
                          groupingvar, study_ids, recast_deaths_geocode,
                          filtRegions = NULL, filtGender = NULL, filtAgeBand = NULL, death_agebreaks = NULL,
                          sero_agebreaks = NULL) {
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
  assert_in(c("country", "study_id", "age_low", "age_high", "region", "gender", "population"),
            colnames(population))
  # check columns for deaths df and dates
  assert_in(c("country", "study_id", "age_low", "age_high", "region", "gender", "n_deaths", "date_start_survey", "date_end_survey"),
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
              "seroprevalence_weighted", "n_tested", "date_start_survey", "date_end_survey"),
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

  # check daily time steps are daily
  if (!cumulative) {
    assert_eq(deaths$date_start_survey, deaths$date_end_survey)
  }

  #............................................................
  # process death data
  #...........................................................
  if (cumulative) {
    deaths <- deaths %>%
      dplyr::mutate(date_start_survey = lubridate::ymd(date_start_survey),
                    date_end_survey = lubridate::ymd(date_end_survey),
                    start_date = min(date_start_survey),
                    ObsDay = as.numeric(date_end_survey - start_date)) %>%
      dplyr::filter(study_id %in% study_ids)

  } else {
    deaths <- deaths %>%
      dplyr::mutate(start_date = min(lubridate::ymd(deaths$date_start_survey)),
                    ObsDay = as.numeric(lubridate::ymd(date_end_survey) - start_date)) %>%
      dplyr::filter(study_id %in% study_ids)
    warning("Assuming you are using European date format")
  }


  # handle age
  if (!is.null(death_agebreaks)) {
    assert_vector(death_agebreaks)
    assert_greq(length(death_agebreaks), 2)
    agebrks <- death_agebreaks
  } else {
    agebrks <- c(0, sort(unique(deaths$age_high)))
  }
  deaths <- deaths %>%
    dplyr::mutate(
      ageband = cut(age_high,
                    breaks = agebrks,
                    labels = c(paste0(agebrks[1:(length(agebrks)-1)], "-", lead(agebrks)[1:(length(agebrks)-1)]))),
      ageband = as.character(ageband),
      ageband = ifelse(age_low == 0 & age_high == 999, "all", ageband)
    )

  # save study start date for later -- this is our index time 0
  start_date <- unique(deaths$start_date)

  ### Iran data processing. Need to rescale recast_deaths_df as we only have one region.
  if (study_ids == "IRN1") {
    ### extract total deaths in Guilan for the time point where we have region specific deaths.
    ## (the age specific deaths are not complete and wrong date so cannot use them for total deaths absolute numbers)
    tot_deaths_iran<-deaths %>%
      dplyr::filter(study_id == "IRN1" & age_low == 0 & age_high == 999 & gender == "both")
  }

  # various filters for death data
  if (groupingvar == "region") {
    deaths <- deaths %>%
      dplyr::filter(for_regional_analysis == 1)
    if (!is.null(filtRegions)) {
      deaths <- deaths %>%
        dplyr::filter(region %in% filtRegions)
    }
  }

  if (groupingvar == "gender") {
    deaths <- deaths %>%
      dplyr::filter(gender_breakdown == 1)
    if(!is.null(filtGender)) {
      deaths <- deaths %>%
        dplyr::filter(gender %in% filtGender)
    }
  }

  if (groupingvar == "ageband") {
    deaths <- deaths %>%
      dplyr::filter(age_breakdown == 1)
    if(!is.null(filtAgeBand)) {
      deaths <- deaths %>%
        dplyr::filter(ageband %in% filtAgeBand)
    }
  }
  # check filters
  assert_gr(nrow(deaths), 0,
            message = "There was a mismatch in filtering death observations. Returned no observations")


  if (cumulative){
    if(length(unique(deaths$ObsDay)) > 1) {
      stop("Cumulative data has multiple end dates")
    }
    #......................
    # if cumulative, recast from recast_deaths_df
    #......................
    recast_deaths_df <- recast_deaths_df %>%
      dplyr::filter(georegion %in% recast_deaths_geocode) %>%
      dplyr::mutate(ObsDay = as.numeric(lubridate::ymd(date) - start_date)) %>%
      dplyr::filter(ObsDay >= 1) %>% # consistent origin
      dplyr::arrange(ObsDay)

    ### For Iran, scale deaths down to represent Guilan region (assuming it has a similar time course)
    if (study_ids == "IRN1") {
      # obtain scaling factor
      recast_tot_deaths_iran <- recast_deaths_df %>%
        dplyr::filter(ObsDay <= tot_deaths_iran$ObsDay)
      scale_iran <- tot_deaths_iran$n_deaths / sum(recast_tot_deaths_iran$deaths)
      ## now scale all recast_deaths_df deaths for this region before we use them
      recast_deaths_df$deaths<- recast_deaths_df$deaths * scale_iran
    }

    # now multiple proportions to get time series
    deaths.prop <- deaths %>%
      dplyr::mutate(ObsDay = max(date_end_survey)) %>%
      dplyr::group_by_at(c(groupingvar, "ObsDay")) %>%
      dplyr::summarise(death_num = sum(n_deaths),
                       age_low = mean(age_low),
                       age_high = mean(age_high)) %>%
      dplyr::arrange(age_low)
    # store group names
    groupvarnames <- group_keys(deaths.prop)
    # get proportion
    deaths.prop <- deaths.prop %>%
      dplyr::ungroup(.) %>%
      dplyr::mutate(death_denom = sum(death_num),
                    death_prop = death_num/death_denom) # protect against double counting of same person in multiple groups

    # now recast proportions across days equally
    deaths.summ <- as.data.frame(matrix(NA, nrow = nrow(deaths.prop), ncol = max(recast_deaths_df$ObsDay)))
    for (i in 1:ncol(deaths.summ)) {
      deaths.summ[,i] <- deaths.prop$death_prop * recast_deaths_df$deaths[i]
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
      dplyr::mutate(ObsDay = as.numeric(lubridate::ymd(date_start_survey) - start_date)) %>%
      dplyr::group_by_at(c(groupingvar, "ObsDay")) %>%
      dplyr::summarise( Deaths = sum(n_deaths) )
  }



  #............................................................
  # process seroprev data
  #...........................................................
  seroprev <- seroprev %>%
    dplyr::mutate(date_start_survey = lubridate::ymd(date_start_survey),
                  date_end_survey = lubridate::ymd(date_end_survey),
                  ObsDaymin = as.numeric(date_start_survey - start_date),
                  ObsDaymax = as.numeric(date_end_survey - start_date)) %>%
    dplyr::filter(study_id %in% study_ids)

  # various filters for serology data
  if (groupingvar == "region"){
    seroprev <- seroprev %>%
      dplyr::filter(for_regional_analysis == 1)
  }

  if (groupingvar == "gender") {
    seroprev <- seroprev %>%
      dplyr::filter(gender_breakdown == 1)
  }

  if (groupingvar == "ageband") {
    seroprev <- seroprev %>%
      dplyr::filter(age_breakdown == 1)
  }

  if (length(unique(seroprev$ObsDaymin)) > 1 | length(unique(seroprev$ObsDaymax)) > 1) {
    stop("Serology data has multiple start or end dates")
  }

  if (groupingvar == "ageband") {
    # handle age
    if (!is.null(sero_agebreaks)) {
      assert_vector(sero_agebreaks)
      assert_greq(length(sero_agebreaks), 2)
      agebrks_sero <- sero_agebreaks
    } else {
      agebrks_sero <- c(min(seroprev$age_low), sort(unique(seroprev$age_high)))
    }
    seroprev <- seroprev %>%
      dplyr::mutate(
        ageband = cut(age_high,
                      breaks = agebrks_sero,
                      labels = c(paste0(agebrks_sero[1:(length(agebrks_sero)-1)], "-", lead(agebrks_sero)[1:(length(agebrks_sero)-1)]))),
        ageband = as.character(ageband),
        ageband = ifelse(age_low == 0 & age_high == 999, "all", ageband)
      )
  }

  # check filters
  assert_gr(nrow(seroprev), 0,
            message = "There was a mismatch in filtering seroprevalence observations. Returned no observations")

  # COMPUTE OVERALL AVERAGE SEROPREVALENCE
  ### fill in gaps for studies which only have number positive, and those which only have seroprevalence.
  seroprev$seroprevalence <- NA
  if (!is.na(seroprev$seroprevalence_weighted[1])) {   ## check if we have weighted/adjusted data for this study.
    seroprev$seroprevalence <- seroprev$seroprevalence_weighted
  } else {
    seroprev$seroprevalence <- seroprev$seroprevalence_unadjusted
  }
  inds <- which(is.na(seroprev$n_positive))
  if (length(inds) >= 1) {
    seroprev$n_positive[inds] <- seroprev$n_tested[inds]*seroprev$seroprevalence[inds]
    inds <- which(is.na(seroprev$seroprevalence))
    seroprev$seroprevalence[inds] <- seroprev$n_positive[inds]/seroprev$n_tested[inds]
  }

  if (is.na(seroprev$seroprevalence_weighted[1]) & is.na(seroprev$seroprevalence_unadjusted[1])) {
    seroprev$seroprevalence<-rowMeans(cbind(seroprev$range_sero_low,seroprev$range_sero_high))
  }
  inds <- which(is.na(seroprev$n_positive))
  seroprev$n_positive[inds] <- seroprev$n_tested[inds]*seroprev$seroprevalence[inds]
  inds <- which(is.na(seroprev$seroprevalence))
  seroprev$seroprevalence[inds] <- seroprev$n_positive[inds]/seroprev$n_tested[inds]


  seroprev.summ <- seroprev %>%
    dplyr::group_by_at(c("ObsDaymin", "ObsDaymax")) %>%
    dplyr::summarise(n_tested = sum(n_tested),
                     n_positive = sum(n_positive)) %>%
    dplyr::mutate(seroprev = n_positive/n_tested) %>%
    dplyr::select(ObsDaymin, ObsDaymax, seroprev) %>%
    dplyr::ungroup()

  ### summarise over grouping variable
  if (all(is.na(seroprev$n_tested)) | all(is.na(seroprev$n_positive))) { ## if no information on sample size to compute weighted average, output current data.
    seroprev.summ.group <- seroprev %>%
      select(c("ObsDaymin", "ObsDaymax", groupingvar, "age_low", "age_high", "n_tested", "n_positive", "seroprevalence"))
  } else {
    seroprev.summ.group <- seroprev %>%
      dplyr::group_by_at(c("ObsDaymin", "ObsDaymax",groupingvar)) %>%
      #dplyr::summarise(seroprev = mean(seroprevalence)) %>%
      dplyr::summarise(age_low = mean(age_low),
                       age_high = mean(age_high),
                       n_tested = sum(n_tested),
                       n_positive = sum(n_positive)) %>%
      dplyr::mutate(seroprevalence = n_positive/n_tested) %>%
      dplyr::arrange(age_low) %>%
      dplyr::ungroup()
  }

  #......................
  # get deaths at midpoint of survey
  #......................
  seromidpt <- ceiling((0.5*(seroprev$ObsDaymax[1] + seroprev$ObsDaymin[1])))
  if (cumulative){
    recast_deaths_at_sero <- recast_deaths_df %>%
      dplyr::filter(ObsDay <= seromidpt)
    deaths.prop <- deaths.prop %>%
      mutate(deaths_at_sero = sum(recast_deaths_at_sero$deaths)*death_prop,
             deaths_denom_at_sero = sum(recast_deaths_at_sero$deaths))
  } else {
    deaths.prop <- deaths %>%
      dplyr::mutate(ObsDay = as.numeric(lubridate::ymd(date_start_survey) - start_date)) %>%
      dplyr::group_by_at(c(groupingvar)) %>%
      dplyr::summarise( deaths_at_sero = sum(n_deaths) ) %>%
      dplyr::filter(ObsDay <= seromidpt) %>%
      dplyr::mutate(deaths_denom_at_sero = sum(deaths_at_sero))

  }

  #...........................................................
  # process population
  #...........................................................
  # NB, we always consider age stratification - use same age breaks as deaths.
  population <- population %>%
    dplyr::mutate(
      ageband = cut(age_high,
                    breaks = agebrks,
                    labels = c(paste0(agebrks[1:(length(agebrks)-1)], "-", lead(agebrks)[1:(length(agebrks)-1)]))),
      ageband = as.character(ageband),
      ageband = ifelse(age_low == 0 & age_high == 999, "all", ageband)
    )

  if (groupingvar == "ageband") {
    population <- population %>%
      dplyr::filter(age_breakdown == 1)
  }

  # various filters for population data
  if (groupingvar == "region") { # but always account age stratification
    population <- population %>%
      dplyr::filter(for_regional_analysis == 1)
  }

  if (groupingvar == "gender") {
    population <- population %>%
      dplyr::filter(gender_breakdown == 1)
  }

  # subset to study id
  population <- population %>%
    dplyr::filter(study_id %in% study_ids)

  # check filters
  assert_gr(nrow(population), 0,
            message = "There was a mismatch in filtering population observations. Returned no observations")

  # get pop group demographics
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
    deaths = deaths.summ,
    seroprev = seroprev.summ,
    prop_pop = pop_prop.summ,
    sero_sens = sensitivity,
    sero_spec = specificity,
    seroprev_group = seroprev.summ.group,
    deaths_group = deaths.prop
  )
  return(ret)

}



#............................................................
#' @title Separate Process Data function for US
#' @param deaths character path; path to death counts by strata
#' @param population character path; path to population counts by strata
#' @param sero_val character path; path to serovalidation file
#' @param seroprev character path; path to seroprevalence file
#' @param cumulative logical; Are the deaths cumulative to a given data or time-series
#' @param timeSeriesFile character path; path to timeSeriesFile aggregate death counts file -- only evaluated if \code{cumulative} is TRUE
#' @param groupingvar character; name of strata(s) being considered
#' @param study_ids character; Study id to be considered
#' @param state character; State to be considered
#' @param county character; County to be considered
#' @import tidyverse
#' NB, this isn't a package, so ^^ is just a reminder to users to have tidyverse loaded in order to allow embracing and piping to work as expected

process_data_usa_facts <- function(deaths = NULL, population = NULL, sero_val = NULL, seroprev = NULL,
                                   cumulative = FALSE, timeSeriesFile = NULL,
                                   groupingvar = NULL, study_ids, state, county) {
  #......................
  # assertions and checks
  #......................
  assert_string(deaths)
  assert_string(population)
  assert_string(sero_val)
  assert_string(seroprev)
  assert_logical(cumulative)
  assert_string(groupingvar)
  assert_in(groupingvar, c("region", "ageband", "gender"))
  assert_string(study_ids)

  #......................
  # read in
  #......................
  deaths <- readr::read_csv(deaths)
  population <- readr::read_csv(population)
  sero_val <- readr::read_csv(sero_val)
  seroprev <- readr::read_csv(seroprev)
  timeSeries <- readr::read_csv(timeSeriesFile)

  #............................................................
  # process time series death data
  #...........................................................
  # save study start date for later -- this is our index time 0
  start_date <- min(lubridate::mdy(timeSeries$Date))
  timeSeries <- timeSeries %>%
    dplyr::mutate(start_date = min(lubridate::ymd(timeSeries$Date)),
                  ObsDay = as.numeric(lubridate::ymd(Date) - start_date)) %>%
    dplyr::select(Date,ObsDay,study_ids) %>%
    dplyr::arrange(ObsDay) %>%
    dplyr::rename(deaths=study_ids)
  warning("Assuming you are using USA date format in time series data")

  #............................................................
  # process death data
  #...........................................................
  deaths <- deaths %>%
    filter(study_id == study_ids)
  # various filters for death data
  if (groupingvar == "region") {
    deaths <- deaths %>%
      dplyr::filter(for_regional_analysis == 1)
  }

  if (groupingvar == "gender") {
    deaths <- deaths %>%
      dplyr::filter(gender_breakdown == 1)
  }

  if (groupingvar == "ageband") {
    deaths <- deaths %>%
      dplyr::filter(age_breakdown == 1)
  }
  # check filtering
  assert_gr(nrow(deaths), 0,
            message = "There was a mismatch in filtering death observations. Returned no observations")

  ## dates
  deaths <- deaths %>%
    dplyr::mutate(ObsDay = as.numeric(lubridate::ymd(date_end_survey) - start_date))

  # handle age
  agebrks <- c(0, sort(unique(deaths$age_high)))

  deaths <- deaths %>%
    dplyr::mutate(
      ageband = cut(age_high,
                    breaks = agebrks,
                    labels = c(paste0(agebrks[1:(length(agebrks)-1)], "-", lead(agebrks)[1:(length(agebrks)-1)]))),
      ageband = as.character(ageband),
      ageband = ifelse(age_low == 0 & age_high == 999, "all", ageband)
    )


  if(length(unique(deaths$ObsDay)) > 1) {
    stop("Cumulative data has multiple end dates")
  }

  #......................
  # recast from time series
  #......................
  if(!is.null(groupingvar)) {
    deaths.prop <- deaths %>%
      dplyr::group_by_at(c(groupingvar, "ObsDay")) %>%
      dplyr::summarise(death_num = sum(n_deaths),
                       age_low=mean(age_low),
                       age_high=mean(age_high)) %>%
      dplyr::arrange(age_low)
    # store group names
    groupvarnames <- group_keys(deaths.prop)
    # get proportion
    deaths.prop <- deaths.prop %>%
      dplyr::ungroup(.) %>%
      dplyr::mutate(death_denom = sum(death_num),
                    death_prop = death_num/death_denom) # protect against double counting of same person in multiple groups
    # now recast proportions across days equally
    deaths.summ <- as.data.frame(matrix(NA, nrow = nrow(deaths.prop), ncol = max(timeSeries$ObsDay)))
    for (i in 1:ncol(deaths.summ)) {
      deaths.summ[,i] <- deaths.prop$death_prop * timeSeries$deaths[i]
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

  }

  #............................................................
  # process seroprev data
  #...........................................................

  # handle sero dates and subset
  seroprev <- seroprev %>%
    dplyr::filter(study_id %in% study_ids) %>%
    dplyr::mutate(date_start_survey = lubridate::ymd(date_start_survey),
                  date_end_survey = lubridate::ymd(date_end_survey),
                  ObsDaymin = as.numeric(date_start_survey - start_date),
                  ObsDaymax = as.numeric(date_end_survey - start_date))
  warning("Assuming you are using European date format for seroprevalence")
  # check filtering
  assert_gr(nrow(seroprev), 0,
            message = "There was a mismatch in filtering seroprevalence observations. Returned no observations")


  # various filters for serology data
  if (groupingvar == "region"){
    seroprev <- seroprev %>%
      dplyr::filter(for_regional_analysis == 1)
  }

  if (groupingvar == "gender") {
    seroprev <- seroprev %>%
      dplyr::filter(gender_breakdown == 1)
  }

  if (groupingvar == "ageband") {
    seroprev <- seroprev %>%
      dplyr::filter(age_breakdown == 1)
  }

  if (groupingvar == "ageband") {
    # handle age

    agebrks_sero <- c(min(seroprev$age_low), sort(unique(seroprev$age_high)))
    seroprev <- seroprev %>%
      dplyr::mutate(
        ageband = cut(age_high,
                      breaks = agebrks_sero,
                      labels = c(paste0(agebrks_sero[1:(length(agebrks_sero)-1)], "-", lead(agebrks_sero)[1:(length(agebrks_sero)-1)]))),
        ageband = as.character(ageband),
        ageband = ifelse(age_low == 0 & age_high == 999, "all", ageband)
      )
  }

  if (length(unique(seroprev$ObsDaymin)) > 1 | length(unique(seroprev$ObsDaymax)) > 1) {
    stop("Serology data has multiple start or end dates")
  }

  # COMPUTE OVERALL AVERAGE SEROPREVALENCE
  ### fill in gaps for studies which only have number positive, and those which only have seroprevalence.
  seroprev$seroprevalence <- NA
  if (!is.na(seroprev$seroprevalence_weighted[1])) {   ## check if we have weighted/adjusted data for this study.
    seroprev$seroprevalence <- seroprev$seroprevalence_weighted
  } else {
    seroprev$seroprevalence <- seroprev$seroprevalence_unadjusted
  }
  inds <- which(is.na(seroprev$n_positive))
  seroprev$n_positive[inds] <- seroprev$n_tested[inds]*seroprev$seroprevalence[inds]
  inds <- which(is.na(seroprev$seroprevalence))
  seroprev$seroprevalence[inds] <- seroprev$n_positive[inds]/seroprev$n_tested[inds]

  seroprev.summ <- seroprev %>%
    dplyr::group_by_at(c("ObsDaymin", "ObsDaymax")) %>%
    dplyr::summarise(n_tested = sum(n_tested),
                     n_positive = sum(n_positive)) %>%
    dplyr::mutate(seroprev = n_positive/n_tested) %>%
    dplyr::select(ObsDaymin, ObsDaymax, seroprev) %>%
    dplyr::ungroup()

  ## deaths at midpoint of survey
  timeSeries_deaths_at_sero <- timeSeries %>%
    dplyr::filter(ObsDay <= ceiling((0.5*(seroprev$ObsDaymax[1] + seroprev$ObsDaymin[1]))))
  if(!is.null(groupingvar)) {
    deaths.prop<-deaths.prop %>%
      mutate(deaths_at_sero=sum(timeSeries_deaths_at_sero$deaths)*death_prop,
             deaths_denom_at_sero=sum(timeSeries_deaths_at_sero$deaths))
  } else {
    deaths.prop<-data.frame(deaths_at_sero=sum(timeSeries_deaths_at_sero$deaths))
  }


  ### summarise over grouping variable
  if (all(is.na(seroprev$n_tested)) | all(is.na(seroprev$n_positive))) { ## if no information on sample size to compute weighted average, output current data.
    seroprev.summ.group <- seroprev %>%
      select(c("ObsDaymin", "ObsDaymax", groupingvar, "age_low", "age_high", "n_tested", "n_positive", "seroprevalence"))
  } else {
    seroprev.summ.group <- seroprev %>%
      dplyr::group_by_at(c("ObsDaymin", "ObsDaymax",groupingvar)) %>%
      #dplyr::summarise(seroprev = mean(seroprevalence)) %>%
      dplyr::summarise(age_low = mean(age_low),
                       age_high = mean(age_high),
                       n_tested = sum(n_tested),
                       n_positive = sum(n_positive)) %>%
      dplyr::mutate(seroprevalence = n_positive/n_tested) %>%
      dplyr::arrange(age_low) %>%
      dplyr::ungroup()
  }

  #...........................................................
  # process population
  #...........................................................
  # subset to study id
  population <- population %>%
    dplyr::filter(State == state & County %in% county)
  # check filtering
  # assert_leq(nrow(population), 1,    ### allow more than one county - for NYC
  #            message = "Cannot select more than one county population")
  assert_gr(nrow(population), 0,
            message = "There was a mismatch in filtering population observations. Returned no observations")
  pop_cols <- grep("Both",names(population))
  pop_cols <- pop_cols[which(pop_cols!=3)]
  population<-population[,pop_cols]
  ## usually only one column but sum over rows in case (e.g. New York City)
  population <- data.frame(age_low=seq(0,85,5),age_high=seq(5,90,5),pop=colSums(population[,1:ncol(population)]))



  if (groupingvar == "ageband") {  ## use same age breaks as deaths.
    # handle age
    population <- population %>%
      dplyr::mutate(
        ageband = cut(age_high,
                      breaks = agebrks,
                      labels = c(paste0(agebrks[1:(length(agebrks)-1)], "-", lead(agebrks)[1:(length(agebrks)-1)]))),
        ageband = as.character(ageband),
        ageband = ifelse(age_low == 0 & age_high == 999, "all", ageband)
      )
  }

  # get pop group demographics
  popN <- sum(population$pop)
  if(!is.null(groupingvar)) {
    pop_prop.summ <- population %>%
      dplyr::group_by_at(groupingvar) %>%
      dplyr::summarise(
        pop_prop = sum(pop)/popN,
        age_low = mean(age_low),
        age_high = mean(age_high)
      ) %>%
      dplyr::arrange(age_low)
  }

  #............................................................
  # process test sens/spec
  #...........................................................
  sero_val <- sero_val %>%
    dplyr::filter(study_id %in% study_ids)

  #...........................................................
  # out
  #...........................................................
  ret <- list(
    deaths = deaths.summ,
    seroprev = seroprev.summ,
    prop_pop = pop_prop.summ,
    popN = popN,
    sero_sens = sero_val$sensitivity,
    sero_spec = sero_val$specificity,
    seroprev_group = seroprev.summ.group,
    deaths_group = deaths.prop
  )
  return(ret)

}


#...........................................................
#' @title Separate Process Data function for US where we have a time series of deaths but no other info
#' @param population character path; path to population counts by strata
#' @param sero_val character path; path to serovalidation file
#' @param seroprev character path; path to seroprevalence file
#' @param timeSeriesFile character path; path to timeSeriesFile aggregate death counts file
#' @param study_ids character; Study id to be considered
#' @param state character; State to be considered
#' @param county character; County to be considered
#' @import tidyverse
#' NB, this isn't a package, so ^^ is just a reminder to users to have tidyverse loaded in order to allow embracing and piping to work as expected

process_usa_basic_data_timeseries <- function(population = NULL, sero_val = NULL, seroprev = NULL,
                                              timeSeriesFile = NULL,
                                              study_ids, state, county) {
  population <- readr::read_csv(population)
  sero_val <- readr::read_csv(sero_val)
  seroprev <- readr::read_csv(seroprev)
  timeSeries<-readr::read_csv(timeSeriesFile)

  # save study start date for later -- this is our index time 0
  start_date <- min(lubridate::mdy(timeSeries$Date))

  #......................
  # process time series
  #......................
  timeSeries <- timeSeries %>%
    dplyr::mutate(start_date = min(lubridate::mdy(timeSeries$Date)),
                  ObsDay = as.numeric(lubridate::mdy(Date) - start_date)) %>%
    dplyr::select(ObsDay,study_ids) %>%
    dplyr::arrange(ObsDay) %>%
    dplyr::rename(deaths=study_ids)

  warning("Assuming you are using USA date format in time series data")

  #......................
  #  process seroprevalence
  #......................
  # handle sero dates and subset
  seroprev <- seroprev %>%
    dplyr::filter(study_id %in% study_ids & for_regional_analysis == 1)

  # check filters
  assert_leq(nrow(seroprev), 1,
             message = "Cannot return more than one seroprevalence observations")
  assert_gr(nrow(seroprev), 0,
            message = "There was a mismatch in filtering seroprevalence observations. Returned no observations")
  # fix dates
  seroprev <- seroprev %>%
    dplyr::mutate(date_start_survey = lubridate::ymd(date_start_survey),
                  date_end_survey = lubridate::ymd(date_end_survey),
                  ObsDaymin = as.numeric(date_start_survey - start_date),
                  ObsDaymax = as.numeric(date_end_survey - start_date))
  warning("Assuming you are using European date format for seroprevalence")

  if (length(unique(seroprev$ObsDaymin)) > 1 | length(unique(seroprev$ObsDaymax)) > 1) {
    stop("Serology data has multiple start or end dates")
  }

  # COMPUTE OVERALL AVERAGE SEROPREVALENCE
  ### fill in gaps for studies which only have number positive, and those which only have seroprevalence.
  seroprev$seroprevalence <- NA
  if (!is.na(seroprev$seroprevalence_weighted[1])) {   ## check if we have weighted/adjusted data for this study.
    seroprev$seroprevalence <- seroprev$seroprevalence_weighted
  } else {
    seroprev$seroprevalence <- seroprev$seroprevalence_unadjusted
  }
  inds <- which(is.na(seroprev$n_positive))
  seroprev$n_positive[inds] <- seroprev$n_tested[inds]*seroprev$seroprevalence[inds]
  inds <- which(is.na(seroprev$seroprevalence))
  seroprev$seroprevalence[inds] <- seroprev$n_positive[inds]/seroprev$n_tested[inds]

  seroprev.summ <- seroprev %>%
    dplyr::select(ObsDaymin, ObsDaymax, seroprevalence)

  ## deaths at midpoint of survey
  timeSeries_deaths_at_sero <- timeSeries %>%
    dplyr::filter(ObsDay <= ceiling((0.5*(seroprev$ObsDaymax[1] + seroprev$ObsDaymin[1]))))
  deaths.prop <- data.frame(ObsDay = ceiling(0.5*(seroprev$ObsDaymax[1] + seroprev$ObsDaymin[1])),
                            deaths_at_sero = sum(timeSeries_deaths_at_sero$deaths))

  # get pop group demographics
  # subset
  population <- population %>%
    dplyr::filter(State == state & County %in% county)
  #if(nrow(population)>1) stop("more than one county population selected")
  assert_gr(nrow(population), 0,
            message = "There was a mismatch in filtering population observations. Returned no observations")

  pop_cols<-grep("Both_Total",names(population))
  popN<-as.numeric(population[1,pop_cols])

  ret <- list(
    deaths = timeSeries,
    seroprev = seroprev.summ,
    popN = popN,
    # sero_sens = sero_val$sensitivity,
    # sero_spec = sero_val$specificity,
    #seroprev_group = seroprev.summ.group,
    deaths_group = deaths.prop
  )
  return(ret)

}


#...........................................................
#' @title Separate Process Data function for US data curated by JHU
#' @param jhufile character path; path to JHU time-series deaths
#' @param study_ids character; Study id to be considered
#' @param state character; State to be considered
#' @param county character; County to be considered -- corresponds to JHU Admin2
process_jhu_timeseries_long <- function(jhufile, county, state) {
  assert_string(jhufile)
  assert_string(county)
  assert_string(state)
  # catch state abbrev
  assert_gr(nchar(state), 2, message = "State abbreviations are not allowed. Full state names need")
  #......................
  # read and wrangle
  #......................
  jhu <- readr::read_csv(jhufile) %>%
    magrittr::set_colnames(tolower(colnames(.)))
  warning("Assuming data is in \"wide format\" and the first date of deaths is the thirteenth column" )

  # apply filters
  jhu.filt <- jhu %>%
    dplyr::filter(admin2 == county & province_state == state)
  # check
  assert_gr(ncol(jhu.filt), 0, message = "There was a mismatch in filtering. Returned no observations")

  # take to long
  jhu.filt <- jhu.filt %>%
    tidyr::gather(., key = date_survey, value = Deaths, 13:ncol(.)) %>%
    dplyr::rename(state = province_state,
                  county = admin2) %>%
    dplyr::mutate(date_survey = lubridate::mdy(date_survey),
                  date_start_survey = min(date_survey),
                  ObsDay = as.numeric(date_survey - date_start_survey)) %>%
    dplyr::select(c("state", "county", "ObsDay", "Deaths"))
  warning("Assuming you are using USA date format")

  # out
  return(jhu.filt)
}



