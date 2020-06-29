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


  ## Serology data not broken down by any of these aspects, so remove filters.
  # various filters for serology data

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
  if(!is.null(groupingvar)) {
    deaths.prop<-deaths.prop %>%
      mutate(deaths_at_sero=sum(timeSeries_deaths_at_sero$deaths)*death_prop,
             deaths_denom_at_sero=sum(timeSeries_deaths_at_sero$deaths))
  } else {
    deaths.prop<-data.frame(deaths_at_sero=sum(timeSeries_deaths_at_sero$deaths))
  }

  #...........................................................
  # process population
  #...........................................................
  # subset to study id
  population <- population %>%
    dplyr::filter(State == state & County == county)
  # check filtering
  assert_leq(nrow(population), 1,
             message = "Cannot select more than one county population")
  assert_gr(nrow(population), 0,
            message = "There was a mismatch in filtering population observations. Returned no observations")
  pop_cols <- grep("Both",names(population))
  pop_cols <- pop_cols[which(pop_cols!=3)]
  population <- data.frame(age_low=seq(0,85,5),age_high=seq(5,90,5),pop=as.numeric(population[1,pop_cols]))

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
    #seroprev_group = seroprev.summ.group,
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
