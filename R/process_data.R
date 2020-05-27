#' @title Process Data for IFR Inference
#' @param deaths character path; path to death counts by strata
#' @param population character path; path to population counts by strata
#' @param sero_val character path; path to serovalidation file
#' @param seroprev character path; path to seroprevalence file
#' @param cumulative logical; Are the deaths cumulative to a given data or time-series
#' @param ECDC character path; path to ECDC aggregate death counts file -- only evaluated if \code{cumulative} is TRUE
#' @param groupingvar character; name of strata(s) being considered
#' @param filtRegion character; region levels to keep
#' @param filtGender character; biological sex levels to keep
#' @param filtAgeBand character; age groups to keep -- note this will be a factor of the age_low and age_high concatenated together
#' @import tidyverse
#' NB, this isn't a package, so ^^ is just a reminder to users to have tidyverse loaded in order to allow embracing and piping to work as expected

source("R/assertions_v5.R")
process_data <- function(deaths = NULL, population = NULL, sero_val = NULL, seroprev = NULL, cumulative = FALSE, ECDC = NULL,
                         groupingvar, study_ids, ecdc_countrycode,
                         filtRegions = NULL, filtGender = NULL, filtAgeBand = NULL) {
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
  if (cumulative){
    assert_string(ECDC)
    assert_string(ecdc_countrycode)
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
  if (cumulative){
    ECDC <- readr::read_csv(ECDC)
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
  assert_in(c("country", "study_id", "age_low", "age_high", "region", "gender", "seroprevalence", "date_start_survey", "date_end_survey"),
            colnames(seroprev))
  # check ecdc
  if (cumulative){
    assert_in(c("dateRep", "deaths", "countryterritoryCode"),
              colnames(ECDC))
  }
  # check daily time steps are daily
  if (!cumulative) {
    assert_eq(deaths$date_start_survey, deaths$date_end_survey)
  }

  #............................................................
  # process death data
  #...........................................................
  deaths <- deaths %>%
    dplyr::mutate(date_start_survey = lubridate::dmy(date_start_survey),
                  date_end_survey = lubridate::dmy(date_end_survey),
                  start_date = min(date_start_survey),
                  ObsDay = as.numeric(date_end_survey - start_date)) %>%
    dplyr::filter(study_id %in% study_ids)
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

  # save study start date for later -- this is our index time 0
  start_date <- unique(deaths$start_date)

  # variuos filters for death data
  if (!is.null(filtRegions)) {
    deaths <- deaths %>%
      dplyr::filter(region %in% filtRegions)
  }

  if (!is.null(filtGender)) {
    deaths <- deaths %>%
      dplyr::filter(gender %in% filtGender)
  }

  if (!is.null(filtAgeBand)) {
    deaths <- deaths %>%
      dplyr::filter(ageband %in% filtAgeBand)
  }


  if (cumulative){
    if(length(unique(deaths$ObsDay)) > 1) {
      stop("Cumulative data has multiple end dates")
    }
    #......................
    # if cumulative, recast from ECDC
    #......................
    upperlim <- unique(deaths$ObsDay)
    ECDC <- ECDC %>%
      dplyr::filter(countryterritoryCode %in% ecdc_countrycode) %>%
      dplyr::mutate(ObsDay = as.numeric(lubridate::dmy(dateRep) - start_date)) %>%
      dplyr::filter(ObsDay <= upperlim & ObsDay >= 1) %>%  # cut off days greater than study period in ECDC and before study period
      dplyr::arrange(ObsDay)
    # now multiple proportions to get time series

    deaths.prop <- deaths %>%
      dplyr::mutate(ObsDay = max(date_end_survey)) %>%
      dplyr::group_by_at(c(groupingvar, "ObsDay")) %>%
      dplyr::summarise(death_num = sum(n_deaths) )
    # store group names
    groupvarnames <- group_keys(deaths.prop)
    # get proportion
    deaths.prop <- deaths.prop %>%
      dplyr::ungroup(.) %>%
      dplyr::mutate(death_denom = sum(death_num),
                    death_prop = death_num/death_denom) # protect against double counting of same person in multiple groups
    # now recast proportions across days equally
    deaths.summ <- as.data.frame(matrix(NA, nrow = nrow(deaths.prop), ncol = upperlim))
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
    deaths.summ <- deaths %>%
      dplyr::mutate(ObsDay = as.numeric(lubridate::dmy(date_start_survey) - start_date)) %>%
      dplyr::group_by_at(c(groupvar, "ObsDay")) %>%
      dplyr::summarise( Deaths = sum(n_deaths) )
  }


  #............................................................
  # process seroprev data
  #...........................................................
  if (groupingvar == "ageband") {
    # handle age
    seroprev <- seroprev %>%
      dplyr::mutate(
        ageband = cut(age_high,
                      breaks = agebrks,
                      labels = c(paste0(agebrks[1:(length(agebrks)-1)], "-", lead(agebrks)[1:(length(agebrks)-1)]))),
        ageband = as.character(ageband),
        ageband = ifelse(age_low == 0 & age_high == 999, "all", ageband)
      )
  }
  # handle sero dates and subset
  seroprev <- seroprev %>%
    dplyr::mutate(date_start_survey = lubridate::dmy(date_start_survey),
                  date_end_survey = lubridate::dmy(date_end_survey),
                  ObsDaymin = as.numeric(date_start_survey - start_date),
                  ObsDaymax = as.numeric(date_end_survey - start_date)) %>%
    dplyr::filter(study_id %in% study_ids)

  if (length(unique(seroprev$ObsDaymin)) > 1 | length(unique(seroprev$ObsDaymax)) > 1) {
    stop("Serology data has multiple start or end dates")
  }

  # variuos filters for serology data
  if (!is.null(filtRegions)) {
    seroprev <- seroprev %>%
      dplyr::filter(region %in% filtRegions)
  }

  if (!is.null(filtGender)) {
    seroprev <- seroprev %>%
      dplyr::filter(gender %in% filtGender)
  }

  if (!is.null(filtAgeBand)) {
    seroprev <- seroprev %>%
      dplyr::filter(ageband %in% filtAgeBand)
  }

  # TODO apply sampling weights?
  # TODO this mean is not really what we want
  seroprev.summ <- seroprev %>%
    dplyr::group_by_at(c("ObsDaymin", "ObsDaymax")) %>%
    dplyr::summarise(seroprev = mean(seroprevalence)) %>%
    dplyr::ungroup()

  #...........................................................
  # process population
  #...........................................................
  if (groupingvar == "ageband") {
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

  # variuos filters for population data
  if (!is.null(filtRegions)) {
    population <- population %>%
      dplyr::filter(region %in% filtRegions)
  }

  if (!is.null(filtGender)) {
    population <- population %>%
      dplyr::filter(gender %in% filtGender)
  }

  if (!is.null(filtAgeBand)) {
    population <- population %>%
      dplyr::filter(ageband %in% filtAgeBand)
  }

  # subset to study id
  population <- population %>%
    dplyr::filter(study_id %in% study_ids)

  # get pa
  popN <- sum(population$population)
  pa <- population %>%
    dplyr::group_by_at(groupingvar) %>%
    dplyr::summarise(
      pop_prop = sum(population)/popN
    )



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
      pa = pa,
      popN = popN,
      sero_sens = sero_val$sensitivity,
      sero_spec = sero_val$specificity
    )
    return(ret)
}
