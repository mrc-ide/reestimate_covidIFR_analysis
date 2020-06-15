#' @title Process Data for IFR Inference
#' @param deaths character path; path to death counts by strata
#' @param population character path; path to population counts by strata
#' @param sero_val character path; path to serovalidation file
#' @param seroprev character path; path to seroprevalence file
#' @param cumulative logical; Are the deaths cumulative to a given data or time-series
#' @param ECDC character path; path to ECDC aggregate death counts file -- only evaluated if \code{cumulative} is TRUE and USAdata is FALSE
#' @param USAdata logical; Is this USA data (which is not included in ECDC)
#' @param JH* character path; path to ECDC aggregate death counts file -- only evaluated if \code{cumulative} is TRUE and USAdata is TRUE
#' @param groupingvar character; name of strata(s) being considered
#' @param filtRegion character; region levels to keep
#' @param filtGender character; biological sex levels to keep
#' @param filtAgeBand character; age groups to keep -- note this will be a factor of the age_low and age_high concatenated together
#' @param death_agebreaks character; potential to customize the break points for the factorization of ages from the death data. Default NULL will use data to set breaks
#' @param sero_agebreaks character; potential to customize the break points for the factorization of ages from the death data. Default NULL will use data to set breaks
#' @import tidyverse
#' NB, this isn't a package, so ^^ is just a reminder to users to have tidyverse loaded in order to allow embracing and piping to work as expected

source("R/assertions_v5.R")
process_data2 <- function(deaths = NULL, population = NULL, sero_val = NULL, seroprev = NULL,
                          cumulative = FALSE, ECDC = NULL, USAdata = FALSE, JHU = NULL,
                          groupingvar, study_ids, geocode,
                          filtRegions = NULL, filtGender = NULL, filtAgeBand = NULL, death_agebreaks = NULL, sero_agebreaks = NULL) {
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
  if (cumulative & !USAdata){
    assert_string(ECDC)
    assert_string(geocode)
  }
  if (cumulative & USAdata){
    assert_string(JHU)
    assert_string(geocode)
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

  #......................
  # read in
  #......................
  deaths <- readr::read_csv(deaths)
  population <- readr::read_csv(population)
  sero_val <- readr::read_csv(sero_val)
  seroprev <- readr::read_csv(seroprev)
  if (cumulative & !USAdata){
    ECDC <- readr::read_csv(ECDC)
  }

  if (cumulative & USAdata){
    JHU <- readr::read_csv(JHU)
  }

  #......................
  # more assertions and checks
  #......................
  assert_in(c("country", "study_id", "age_low", "age_high", "region", "gender", "n_deaths", "date_start_survey", "date_end_survey"),
            colnames(deaths))
  assert_in(c("country", "study_id", "age_low", "age_high", "region", "gender", "population", "date"),
            colnames(population))
  assert_in(c("study_id", "sensitivity", "specificity"),
            colnames(sero_val))
  assert_in(c("country", "study_id", "age_low", "age_high", "region", "gender", "seroprevalence_unadjusted", "date_start_survey", "date_end_survey"),
            colnames(seroprev))
  # check ecdc
  if (cumulative & !USAdata){
    assert_in(c("dateRep", "deaths", "countryterritoryCode"),
              colnames(ECDC))
  }
  # check JHU
  if (cumulative & USAdata){
    assert_in(c("admin2", "date", "deaths"),
              colnames(JHU))
  }

  # check daily time steps are daily
  if (!cumulative) {
    assert_eq(deaths$date_start_survey, deaths$date_end_survey)
  }

  #............................................................
  # process death data
  #...........................................................
  if (cumulative) {
    if (USAdata) { # protect dates
      deaths <- deaths %>%
        dplyr::mutate(date_start_survey = lubridate::mdy(date_start_survey),
                      date_end_survey = lubridate::mdy(date_end_survey),
                      start_date = min(date_start_survey),
                      ObsDay = as.numeric(date_end_survey - start_date)) %>%
        dplyr::filter(study_id %in% study_ids)
      warning("Assuming you are using USA date format")

    } else {
      deaths <- deaths %>%
        dplyr::mutate(date_start_survey = lubridate::dmy(date_start_survey),
                      date_end_survey = lubridate::dmy(date_end_survey),
                      start_date = min(date_start_survey),
                      ObsDay = as.numeric(date_end_survey - start_date)) %>%
        dplyr::filter(study_id %in% study_ids)
      warning("Assuming you are using European date format")

    }


  } else {
    if (USAdata) { # protect dates
      start_date <- min(lubridate::mdy(deaths$date_start_survey))
      deaths <- deaths %>%
        dplyr::mutate(start_date = min(lubridate::mdy(deaths$date_start_survey)),
                      ObsDay = as.numeric(lubridate::mdy(date_end_survey) - start_date)) %>%
        dplyr::filter(study_id %in% study_ids)

      warning("Assuming you are using USA date format")
    } else {

      deaths <- deaths %>%
        dplyr::mutate(start_date = min(lubridate::dmy(deaths$date_start_survey)),
                      ObsDay = as.numeric(lubridate::dmy(date_end_survey) - start_date)) %>%
        dplyr::filter(study_id %in% study_ids)
      warning("Assuming you are using European date format")
    }


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

  ### Iran data processing. Need to rescale ECDC deaths as we only have one region.
  if(study_ids=="IRN1") {
    ### extract total deaths in Guilan for the time point where we have region specific deaths.
    ## (the age specific deaths are not complete and wrong date so cannot use them for total deaths absolute numbers)
    tot_deaths_iran<-deaths %>%
      dplyr::filter(study_id=="IRN1" & age_low==0 & age_high==999 & gender=="both")
  }

  # various filters for death data
  if (groupingvar == "region") {
    deaths <- deaths %>%
      dplyr::filter(for_regional_analysis == 1)
    if(!is.null(filtRegions)) {
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


  if (cumulative & !USAdata){
    if(length(unique(deaths$ObsDay)) > 1) {
      stop("Cumulative data has multiple end dates")
    }
    #......................
    # if cumulative, recast from ECDC
    #......................
    #upperlim <- unique(deaths$ObsDay)
    ECDC <- ECDC %>%
      dplyr::filter(countryterritoryCode %in% geocode) %>%
      dplyr::mutate(ObsDay = as.numeric(lubridate::dmy(dateRep) - start_date)) %>%
      dplyr::filter(ObsDay >= 1) %>%
      dplyr::arrange(ObsDay)
    #dplyr::filter(ObsDay <= upperlim & ObsDay >= 1) %>%  # cut off days greater than study period in ECDC and before study period
    # now multiple proportions to get time series

    ### For Iran, scale deaths down to represent Guilan region (assuming it has a similar time course)
    if(study_ids=="IRN1") {
      # obtain scaling factor
      ecdc_tot_deaths_iran<-ECDC %>%
        dplyr::filter(ObsDay<=tot_deaths_iran$ObsDay)
      scale_iran<-tot_deaths_iran$n_deaths / sum(ecdc_tot_deaths_iran$deaths)
      ## now scale all ECDC deaths for this region before we use them
      ECDC$deaths<- ECDC$deaths * scale_iran
    }

    deaths.prop <- deaths %>%
      dplyr::mutate(ObsDay = max(date_end_survey)) %>%
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
    deaths.summ <- as.data.frame(matrix(NA, nrow = nrow(deaths.prop), ncol = max(ECDC$ObsDay)))
    for (i in 1:ncol(deaths.summ)) {
      deaths.summ[,i] <- deaths.prop$death_prop * ECDC$deaths[i]
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

  } else if (cumulative & USAdata) {
    if(length(unique(deaths$ObsDay)) > 1) {
      stop("Cumulative data has multiple end dates")
    }
    #......................
    # if cumulative, recast from JHU
    #......................
    JHU <- JHU %>%
      dplyr::filter(province_state %in% geocode) %>%
      dplyr::mutate(ObsDay = as.numeric(lubridate::dmy(date) - start_date)) %>%
      dplyr::filter(ObsDay >= 1) %>%
      dplyr::arrange(ObsDay)

    deaths.prop <- deaths %>%
      dplyr::mutate(ObsDay = max(date_end_survey)) %>%
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
    deaths.summ <- as.data.frame(matrix(NA, nrow = nrow(deaths.prop), ncol = max(ECDC$ObsDay)))
    for (i in 1:ncol(deaths.summ)) {
      deaths.summ[,i] <- deaths.prop$death_prop * ECDC$deaths[i]
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
    if (USAdata) {
      deaths.summ <- deaths %>%
        dplyr::mutate(ObsDay = as.numeric(lubridate::mdy(date_start_survey) - start_date)) %>%
        dplyr::group_by_at(c(groupingvar, "ObsDay")) %>%
        dplyr::summarise( Deaths = sum(n_deaths) )
      warning("Assuming you are using USA date format")
    } else {
      deaths.summ <- deaths %>%
        dplyr::mutate(ObsDay = as.numeric(lubridate::dmy(date_start_survey) - start_date)) %>%
        dplyr::group_by_at(c(groupingvar, "ObsDay")) %>%
        dplyr::summarise( Deaths = sum(n_deaths) )
      warning("Assuming you are using European date format")
    }

  }


  #............................................................
  # process seroprev data
  #...........................................................

  # handle sero dates and subset
  if (USAdata) {
    seroprev <- seroprev %>%
      dplyr::mutate(date_start_survey = lubridate::mdy(date_start_survey),
                    date_end_survey = lubridate::mdy(date_end_survey),
                    ObsDaymin = as.numeric(date_start_survey - start_date),
                    ObsDaymax = as.numeric(date_end_survey - start_date)) %>%
      dplyr::filter(study_id %in% study_ids)
    warning("Assuming you are using US date format")
  } else {
    seroprev <- seroprev %>%
      dplyr::mutate(date_start_survey = lubridate::dmy(date_start_survey),
                    date_end_survey = lubridate::dmy(date_end_survey),
                    ObsDaymin = as.numeric(date_start_survey - start_date),
                    ObsDaymax = as.numeric(date_end_survey - start_date)) %>%
      dplyr::filter(study_id %in% study_ids)
    warning("Assuming you are using European date format")
  }


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
    seroprev.summ.group<-seroprev %>%
            select(ObsDaymin, ObsDaymax, groupingvar,age_low, age_high,n_tested,n_positive,seroprevalence)
  } else {
    seroprev.summ.group <- seroprev %>%
      dplyr::group_by_at(c("ObsDaymin", "ObsDaymax",groupingvar)) %>%
      #dplyr::summarise(seroprev = mean(seroprevalence)) %>%
      dplyr::summarise(age_low=mean(age_low),
                       age_high=mean(age_high),
                       n_tested = sum(n_tested),
                       n_positive = sum(n_positive)) %>%
      dplyr::mutate(seroprevalence = n_positive/n_tested) %>%
      dplyr::arrange(age_low) %>%
      dplyr::ungroup()
  }

  ## deaths at midpoint of survey
  ecdc_deaths_at_sero <- ECDC %>%
    dplyr::filter(ObsDay <= (0.5*(seroprev$ObsDaymax[1] + seroprev$ObsDaymin[1])))
  deaths.prop<-deaths.prop %>%
    mutate(deaths_at_sero=sum(ecdc_deaths_at_sero$deaths)*death_prop,
           deaths_denom_at_sero=sum(ecdc_deaths_at_sero$deaths))

  #...........................................................
  # process population
  #...........................................................
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

  # various filters for population data
  if (groupingvar == "region") {
    population <- population %>%
      dplyr::filter(for_regional_analysis == 1)
  }

  if (groupingvar == "gender") {
    population <- population %>%
      dplyr::filter(gender_breakdown == 1)
  }

  if (groupingvar == "ageband") {
    population <- population %>%
      dplyr::filter(age_breakdown == 1)
  }

  # subset to study id
  population <- population %>%
    dplyr::filter(study_id %in% study_ids)

  # get pop group demographics
  popN <- sum(population$population)
  pop_prop.summ <- population %>%
    dplyr::group_by_at(groupingvar) %>%
    dplyr::summarise(
      pop_prop = sum(population)/popN,
      age_low = mean(age_low),
      age_high = mean(age_high)
    ) %>%
    dplyr::arrange(age_low)

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



###########################################
# Separate function for US.
process_data_usa <- function(deaths = NULL, population = NULL, sero_val = NULL, seroprev = NULL,
                          cumulative = FALSE, USAdata = FALSE, JHU = NULL,
                          groupingvar, study_ids, geocode,
                          filtRegions = NULL, filtGender = NULL, filtAgeBand = NULL, death_agebreaks = NULL, sero_agebreaks = NULL) {
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
  if (cumulative & USAdata){
    assert_string(JHU)
    assert_string(geocode)
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

  #......................
  # read in
  #......................
  deaths <- readr::read_csv(deaths)
  population <- readr::read_csv(population)
  sero_val <- readr::read_csv(sero_val)
  seroprev <- readr::read_csv(seroprev)

  if (cumulative & USAdata){
    JHU <- readr::read_csv(JHU)
  }

  #......................
  # more assertions and checks
  #......................
  assert_in(c("country", "study_id", "age_low", "age_high", "region", "gender", "n_deaths", "date_start_survey", "date_end_survey"),
            colnames(deaths))
  assert_in(c("country", "study_id", "age_low", "age_high", "region", "gender", "population", "date"),
            colnames(population))
  assert_in(c("study_id", "sensitivity", "specificity"),
            colnames(sero_val))
  assert_in(c("country", "study_id", "age_low", "age_high", "region", "gender", "seroprevalence_unadjusted", "date_start_survey", "date_end_survey"),
            colnames(seroprev))

  # check JHU
  if (cumulative & USAdata){
    assert_in(c("admin2", "date", "deaths"),
              colnames(JHU))
  }

  # check daily time steps are daily
  if (!cumulative) {
    assert_eq(deaths$date_start_survey, deaths$date_end_survey)
  }

  #............................................................
  # process death data
  #...........................................................
  if (cumulative) {
    if (USAdata) { # protect dates
      deaths <- deaths %>%
        dplyr::mutate(date_start_survey = lubridate::mdy(date_start_survey),
                      date_end_survey = lubridate::mdy(date_end_survey),
                      start_date = min(date_start_survey),
                      ObsDay = as.numeric(date_end_survey - start_date)) %>%
        dplyr::filter(study_id %in% study_ids)
      warning("Assuming you are using USA date format")

    } else {
      deaths <- deaths %>%
        dplyr::mutate(date_start_survey = lubridate::dmy(date_start_survey),
                      date_end_survey = lubridate::dmy(date_end_survey),
                      start_date = min(date_start_survey),
                      ObsDay = as.numeric(date_end_survey - start_date)) %>%
        dplyr::filter(study_id %in% study_ids)
      warning("Assuming you are using European date format")

    }


  } else {
    if (USAdata) { # protect dates
      start_date <- min(lubridate::mdy(deaths$date_start_survey))
      deaths <- deaths %>%
        dplyr::mutate(start_date = min(lubridate::mdy(deaths$date_start_survey)),
                      ObsDay = as.numeric(lubridate::mdy(date_end_survey) - start_date)) %>%
        dplyr::filter(study_id %in% study_ids)

      warning("Assuming you are using USA date format")
    } else {

      deaths <- deaths %>%
        dplyr::mutate(start_date = min(lubridate::dmy(deaths$date_start_survey)),
                      ObsDay = as.numeric(lubridate::dmy(date_end_survey) - start_date)) %>%
        dplyr::filter(study_id %in% study_ids)
      warning("Assuming you are using European date format")
    }


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

  # various filters for death data
  if (groupingvar == "region") {
    deaths <- deaths %>%
      dplyr::filter(for_regional_analysis == 1)
    if(!is.null(filtRegions)) {
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


  if (cumulative & !USAdata){
    if(length(unique(deaths$ObsDay)) > 1) {
      stop("Cumulative data has multiple end dates")
    }

  if (cumulative & USAdata) {
    if(length(unique(deaths$ObsDay)) > 1) {
      stop("Cumulative data has multiple end dates")
    }
    #......................
    # if cumulative, recast from JHU
    #......................
    JHU <- JHU %>%
      dplyr::filter(province_state %in% geocode) %>%
      dplyr::mutate(ObsDay = as.numeric(lubridate::dmy(date) - start_date)) %>%
      dplyr::filter(ObsDay >= 1) %>%
      dplyr::arrange(ObsDay)

    deaths.prop <- deaths %>%
      dplyr::mutate(ObsDay = max(date_end_survey)) %>%
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
    deaths.summ <- as.data.frame(matrix(NA, nrow = nrow(deaths.prop), ncol = max(ECDC$ObsDay)))
    for (i in 1:ncol(deaths.summ)) {
      deaths.summ[,i] <- deaths.prop$death_prop * ECDC$deaths[i]
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
    if (USAdata) {
      deaths.summ <- deaths %>%
        dplyr::mutate(ObsDay = as.numeric(lubridate::mdy(date_start_survey) - start_date)) %>%
        dplyr::group_by_at(c(groupingvar, "ObsDay")) %>%
        dplyr::summarise( Deaths = sum(n_deaths) )
      warning("Assuming you are using USA date format")
    } else {
      deaths.summ <- deaths %>%
        dplyr::mutate(ObsDay = as.numeric(lubridate::dmy(date_start_survey) - start_date)) %>%
        dplyr::group_by_at(c(groupingvar, "ObsDay")) %>%
        dplyr::summarise( Deaths = sum(n_deaths) )
      warning("Assuming you are using European date format")
    }

  }


  #............................................................
  # process seroprev data
  #...........................................................

  # handle sero dates and subset
  if (USAdata) {
    seroprev <- seroprev %>%
      dplyr::mutate(date_start_survey = lubridate::mdy(date_start_survey),
                    date_end_survey = lubridate::mdy(date_end_survey),
                    ObsDaymin = as.numeric(date_start_survey - start_date),
                    ObsDaymax = as.numeric(date_end_survey - start_date)) %>%
      dplyr::filter(study_id %in% study_ids)
    warning("Assuming you are using US date format")
  } else {
    seroprev <- seroprev %>%
      dplyr::mutate(date_start_survey = lubridate::dmy(date_start_survey),
                    date_end_survey = lubridate::dmy(date_end_survey),
                    ObsDaymin = as.numeric(date_start_survey - start_date),
                    ObsDaymax = as.numeric(date_end_survey - start_date)) %>%
      dplyr::filter(study_id %in% study_ids)
    warning("Assuming you are using European date format")
  }


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
    seroprev.summ.group<-seroprev %>%
      select(ObsDaymin, ObsDaymax, groupingvar,age_low, age_high,n_tested,n_positive,seroprevalence)
  } else {
    seroprev.summ.group <- seroprev %>%
      dplyr::group_by_at(c("ObsDaymin", "ObsDaymax",groupingvar)) %>%
      #dplyr::summarise(seroprev = mean(seroprevalence)) %>%
      dplyr::summarise(age_low=mean(age_low),
                       age_high=mean(age_high),
                       n_tested = sum(n_tested),
                       n_positive = sum(n_positive)) %>%
      dplyr::mutate(seroprevalence = n_positive/n_tested) %>%
      dplyr::arrange(age_low) %>%
      dplyr::ungroup()
  }

  ## deaths at midpoint of survey
  ecdc_deaths_at_sero <- ECDC %>%
    dplyr::filter(ObsDay <= (0.5*(seroprev$ObsDaymax[1] + seroprev$ObsDaymin[1])))
  deaths.prop<-deaths.prop %>%
    mutate(deaths_at_sero=sum(ecdc_deaths_at_sero$deaths)*death_prop,
           deaths_denom_at_sero=sum(ecdc_deaths_at_sero$deaths))

  #...........................................................
  # process population
  #...........................................................
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

  # various filters for population data
  if (groupingvar == "region") {
    population <- population %>%
      dplyr::filter(for_regional_analysis == 1)
  }

  if (groupingvar == "gender") {
    population <- population %>%
      dplyr::filter(gender_breakdown == 1)
  }

  if (groupingvar == "ageband") {
    population <- population %>%
      dplyr::filter(age_breakdown == 1)
  }

  # subset to study id
  population <- population %>%
    dplyr::filter(study_id %in% study_ids)

  # get pop group demographics
  popN <- sum(population$population)
  pop_prop.summ <- population %>%
    dplyr::group_by_at(groupingvar) %>%
    dplyr::summarise(
      pop_prop = sum(population)/popN,
      age_low = mean(age_low),
      age_high = mean(age_high)
    ) %>%
    dplyr::arrange(age_low)

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

