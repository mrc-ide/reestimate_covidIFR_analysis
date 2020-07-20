#................................................................................................
## Purpose: Process Country datasets for analysis
##
## Notes: Process data functions and files are specific to regions
#................................................................................................
library(tidyverse)
source("R/process_data3.R")
#..................................................................................
#---- Preprocess Eurasia Data #----
#..................................................................................
# deaths
deathsdf <- readr::read_csv("data/raw/deaths.csv") %>%
  dplyr::select(-c("ref", "notes")) %>%
  dplyr::mutate(date_start_survey = lubridate::dmy(date_start_survey), # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format
                date_end_survey = lubridate::dmy(date_end_survey))
ECDCdf <- readr::read_csv("data/raw/daily_deaths_ECDC20200518.csv") %>%
  dplyr::select(c("dateRep", "countryterritoryCode", "deaths")) %>%
  dplyr::rename(date = dateRep,
                georegion = countryterritoryCode) %>%
  dplyr::mutate(date = lubridate::dmy(date), # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format
                deaths = ifelse(deaths < 0, 0, deaths)) # remove deaths typos
# demography
populationdf <- readr::read_csv("data/raw/population.csv") %>%
  dplyr::select(-c("reference")) %>%
  dplyr::mutate(age_low = ifelse(age_low == 0 & age_high == 0, 1, age_low),
                age_high = ifelse(age_low == 1 & age_high == 0, 1, age_high))  # liftover "zero" year olds to be 1, 1 as well

# seroprev
sero_valdf <- readr::read_csv("data/raw/seroassay_validation.csv")
sero_prevdf <- readr::read_csv("data/raw/seroprevalence.csv") %>%
  dplyr::select(-c("ref", "notes")) %>%
  dplyr::mutate(date_start_survey = lubridate::dmy(date_start_survey), # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format
                date_end_survey = lubridate::dmy(date_end_survey),
                seroprevalence_unadjusted = ifelse(is.na(seroprevalence_unadjusted), n_positive/n_tested, seroprevalence_unadjusted)
  )



#..................................................................................
#----- Process Eurasia Data #-----
#..................................................................................
#............................................................
# Spain
#...........................................................
#......................
# regions
#......................
ESP.regions.dat <- process_data3(deaths = deathsdf,
                                 population = populationdf,
                                 sero_val = sero_valdf,
                                 seroprev = sero_prevdf,
                                 cumulative = TRUE,
                                 recast_deaths_df = ECDCdf,
                                 groupingvar = "region",
                                 study_ids = "ESP1",
                                 recast_deaths_geocode = "ESP",
                                 filtRegions = c("Andalucia", "Aragon",    "Asturias",  "Baleares",  "C Valenciana", "Canarias",
                                                 "Cantabria", "Castilla La Mancha", "Castilla y Leon", "Cataluna", "Extremadura",
                                                 "Galicia",   "La Rioja",  "Madrid", "Murcia", "Navarra", "Pais Vasco"), # limit to mainland Spain
                                 filtGender = NULL,
                                 filtAgeBand = NULL)


#......................
# agebands
#......................
ESP.agebands.dat <- process_data3(deaths = deathsdf,
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  cumulative = TRUE,
                                  recast_deaths_df = ECDCdf,
                                  groupingvar = "ageband",
                                  study_ids = "ESP1",
                                  recast_deaths_geocode = "ESP",
                                  death_agebreaks = c(0, 10, 20, 30, 40,
                                                      50, 60, 70, 80, 90, 999),
                                  sero_agebreaks = c(0, 10, 20, 30, 40,
                                                     50, 60, 70, 80, 90, 999))
#......................
# MANUAL ADJUSTMENTS
#......................
# suspicious 0 followed by very large increase
# Monday, April 27 appears suspicious, such that all deaths from Apr 27 got pushed to Apr 28
ESP.agebands.dat$deathsMCMC$Deaths[ESP.agebands.dat$deathsMCMC$ObsDay == 117] <- -1
ESP.agebands.dat$deathsMCMC$Deaths[ESP.agebands.dat$deathsMCMC$ObsDay == 118] <- -1
ESP.regions.dat$deathsMCMC$Deaths[ESP.regions.dat$deathsMCMC$ObsDay == 117] <- -1
ESP.regions.dat$deathsMCMC$Deaths[ESP.regions.dat$deathsMCMC$ObsDay == 118] <- -1

# No adjustments to serology as there is alignment of deaths and seroprevalence age groups.

#......................
# get rho
#......................
ESP.agebands.dat$rho <- rep(1, length(unique(ESP.agebands.dat$deathsMCMC$ageband)))
# multiple through demog and age-standardize for region
ESP.regions.dat$rho <- rep(1, length(unique(ESP.regions.dat$deathsMCMC$region)))


#......................
# save out
#......................
dir.create("data/derived/ESP/ESP_agebands.RDS")
saveRDS(ESP.agebands.dat, "data/derived/ESP/ESP_agebands.RDS")
saveRDS(ESP.regions.dat, "data/derived/ESP/ESP_regions.RDS")

#............................................................
# Netherlands
#...........................................................
#......................
# regions
#......................
## NB one or two regions have missing seroprevalence, as map regions did not match to current regions.
NLD.regions.dat <- process_data3(deaths = deathsdf,
                                 population = populationdf,
                                 sero_val = sero_valdf,
                                 seroprev = sero_prevdf,
                                 cumulative = TRUE,
                                 recast_deaths_df = ECDCdf,
                                 groupingvar = "region",
                                 study_ids = "NLD1",
                                 recast_deaths_geocode = "NLD")

#......................
# ages
#......................
NLD.agebands.dat <- process_data3(deaths = deathsdf,
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  cumulative = TRUE,
                                  recast_deaths_df = ECDCdf,
                                  groupingvar = "ageband",
                                  study_ids = "NLD1",
                                  recast_deaths_geocode = "NLD",
                                  death_agebreaks = c(0, 9, 19, 29, 39, 49, 59, 69, 79, 89, 999),
                                  sero_agebreaks = c(0, 9, 19, 29, 39, 49, 59, 69, 79, 89, 999))


#......................
# MANUAL ADJUSTMENTS
#......................
# Netherlands seroprevalence and deaths not perfectly aligned.
# Assumptions.
# 1) 31-40 seroprevalence will be equivalent to the 30-39 age group, same for other 10 year age bands 30-70.
# 2) 18-30 seroprevalence = 20-29 seroprevalence
# 3) 60-72 seroprevalence = 60-69 seroprevalence
# 4) all other seroprevalence = national average.
agebands <- unique(NLD.agebands.dat$deathsMCMC$ageband)
nld_adj_seroprev <- tibble::tibble(
  ObsDaymin = unique(NLD.agebands.dat$seroprevMCMC$ObsDaymin),
  ObsDaymax = unique(NLD.agebands.dat$seroprevMCMC$ObsDaymax),
  ageband = agebands,
  age_low = as.numeric(stringr::str_split_fixed(agebands, "-[0-9]+", n=2)[,1]),
  age_high= as.numeric(stringr::str_split_fixed(agebands, "[0-9]+-", n=2)[,2]),
  seroprevalence = NA) %>%
  dplyr::arrange(age_low)

nld_org_seroprev <- NLD.agebands.dat$seroprev_group %>%
  dplyr::select(c("age_low", "age_high", "seroprevalence"))

nld_adj_seroprev$seroprevalence <- apply(nld_adj_seroprev, 1, wiggle_age_matchfun, wiggle = 2, y = nld_org_seroprev)
# assuming missing is mean
nldmean <- mean(nld_adj_seroprev$seroprevalence, na.rm = T)
nld_adj_seroprev <- nld_adj_seroprev %>%
  dplyr::mutate(seroprevalence = ifelse(is.na(seroprevalence), nldmean, seroprevalence)) %>%
  dplyr::rename(SeroPrev = seroprevalence) %>%
  dplyr::select(c("ObsDaymin", "ObsDaymax", "ageband", "SeroPrev"))

# write new serology df
NLD.agebands.dat$seroprevMCMC <- nld_adj_seroprev

# Netherlands seroprevalence missing in some regions
# assume that this missing values can be imputed as the mean of the other regions
NLD.regions.dat$seroprevMCMC$SeroPrev <- NLD.regions.dat$seroprev_group$seroprevalence
nldmean <- mean(NLD.regions.dat$seroprev_group$seroprevalence, na.rm = T)
NLD.regions.dat$seroprevMCMC <- NLD.regions.dat$seroprevMCMC %>%
  dplyr::mutate(SeroPrev = ifelse(is.na(SeroPrev), nldmean, SeroPrev))



#......................
# get rho
#......................
NLD.agebands.dat$rho <- rep(1, length(unique(NLD.agebands.dat$deathsMCMC$ageband)))
# multiple through demog and age-standardize for region
NLD.regions.dat$rho <- rep(1, length(unique(NLD.regions.dat$deathsMCMC$region)))

#......................
# save out
#......................
dir.create("data/derived/NLD", recursive = T)
saveRDS(NLD.regions.dat, "data/derived/NLD/NLD_regions.RDS")
saveRDS(NLD.agebands.dat, "data/derived/NLD/NLD_agebands.RDS")



#............................................................
# Denmark
#...........................................................
#......................
# regions
#......................
DNK.regions.dat <- process_data3(deaths = deathsdf,
                                 population = populationdf,
                                 sero_val = sero_valdf,
                                 seroprev = sero_prevdf,
                                 cumulative = TRUE,
                                 recast_deaths_df = ECDCdf,
                                 groupingvar = "region",
                                 study_ids = "DNK1",
                                 recast_deaths_geocode = "DNK")

#......................
# ages
#......................
DNK.agebands.dat <- process_data3(deaths = deathsdf,
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  cumulative = TRUE,
                                  recast_deaths_df = ECDCdf,
                                  groupingvar = "ageband",
                                  study_ids = "DNK1",
                                  recast_deaths_geocode = "DNK",
                                  death_agebreaks = c(0, 59, 69, 79, 89, 999),
                                  sero_agebreaks = c(0, 59, 69, 79, 89, 999),
                                  filtGender = "both")


#......................
# MANUAL ADJUSTMENTS
#......................
# Denmark deaths agebands do not overlap well -- potentially have one age group in 59-69
#  but instead assume ALL national average
agebands <- unique(DNK.agebands.dat$deathsMCMC$ageband)
dnk_adj_seroprev <- tibble::tibble(
  ObsDaymin = unique(DNK.agebands.dat$seroprevMCMC$ObsDaymin),
  ObsDaymax = unique(DNK.agebands.dat$seroprevMCMC$ObsDaymax),
  ageband = agebands,
  SeroPrev = mean(DNK.agebands.dat$seroprev_group$seroprevalence))
DNK.agebands.dat$seroprevMCMC <- dnk_adj_seroprev

#......................
# get rho
#......................
DNK.agebands.dat$rho <- rep(1, length(unique(DNK.agebands.dat$deathsMCMC$ageband)))
# multiple through demog and age-standardize for region
DNK.regions.dat$rho <- rep(1, length(unique(DNK.regions.dat$deathsMCMC$region)))

#......................
# save out
#......................
dir.create("data/derived/DNK", recursive = T)
saveRDS(DNK.regions.dat, "data/derived/DNK/DNK_regions.RDS")
saveRDS(DNK.agebands.dat, "data/derived/DNK/DNK_agebands.RDS")

#............................................................
# Switzerland
#...........................................................
#......................
# ages
#......................
## CHE1. TODO include change in seroprevalence over 3 weeks. (Currently output average over all of them)
## Has specific Geneva deaths time series, do not need ECDC
## rename vars so will work like ECDC file
CHE1TimeSeries <- readr::read_csv("data/raw/deaths_time_series.csv") %>%
  dplyr::filter(study_id == "CHE1")
# sanity check
identical(CHE1TimeSeries$date_start_survey, CHE1TimeSeries$date_end_survey)
CHE1TimeSeries <- CHE1TimeSeries %>%
  dplyr::rename(date = date_end_survey,
                deaths = n_deaths,
                georegion = region) %>%
  dplyr::mutate(date = lubridate::dmy(date)) # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format

#......................
# ages
#......................
CHE.agebands.dat <- process_data3(deaths = deathsdf,
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  cumulative = TRUE,
                                  recast_deaths_df = CHE1TimeSeries,
                                  groupingvar = "ageband",
                                  study_ids = "CHE1",
                                  recast_deaths_geocode = "Geneva")   ## use study id in case we get more studies later.

#......................
# MANUAL ADJUSTMENTS
#......................
# Assume:
# (1) seroprevalence 5-19 is representative of deaths for 0-10 and 10-20, respectively
# (1) seroprevalence 19-49 is representative of deaths for 20-30, 30-40, and 40-50, respectively
# (1) seroprevalence 49+ is representative of deaths for 49+ age bands
agebands <- unique(CHE.agebands.dat$deathsMCMC$ageband)
che_adj_seroprev <- tibble::tibble(
  ObsDaymin = unique(CHE.agebands.dat$seroprevMCMC$ObsDaymin),
  ObsDaymax = unique(CHE.agebands.dat$seroprevMCMC$ObsDaymax),
  ageband = agebands,
  age_low = as.numeric(stringr::str_split_fixed(agebands, "-[0-9]+", n=2)[,1]),
  age_high= as.numeric(stringr::str_split_fixed(agebands, "[0-9]+-", n=2)[,2]),
  seroprevalence = NA) %>%
  dplyr::arrange(age_low)

che_org_seroprev <- CHE.agebands.dat$seroprev_group %>%
  dplyr::select(c("age_low", "age_high", "seroprevalence"))

che_adj_seroprev$seroprevalence <- apply(che_adj_seroprev, 1, wiggle_age_matchfun, wiggle = 2,
                                         y = che_org_seroprev)
# imputer missing age group write new serology df
meanche <- mean(che_adj_seroprev$seroprevalence, na.rm = T)
che_adj_seroprev <- che_adj_seroprev %>%
  dplyr::mutate(seroprevalence = ifelse(is.na(seroprevalence), meanche, seroprevalence)) %>%
  dplyr::rename(SeroPrev = seroprevalence) %>%
  dplyr::select(-c("age_low", "age_high"))
CHE.agebands.dat$seroprevMCMC <- che_adj_seroprev

#......................
# get rho
#......................
CHE.agebands.dat$rho <- rep(1, length(unique(che_adj_seroprev$ageband)))

#......................
# save out
#......................
dir.create("data/derived/CHE", recursive = T)
saveRDS(CHE.agebands.dat, "data/derived/CHE/CHE_agebands.RDS")


#............................................................
# Iran
#...........................................................
# extract total deaths in Guilan for the time point where we have region specific deaths.
# (the age specific deaths are not complete and wrong date so cannot use them for total deaths absolute
IRNdeaths <- deathsdf %>%
  dplyr::filter(study_id == "IRN1") %>%
  dplyr::filter(study_id == "IRN1" & age_low == 0 & age_high == 999 & gender == "both") %>%
  dplyr::mutate(age_breakdown = 1)

IRN_recast_df <- ECDCdf %>%
  dplyr::filter(georegion == "IRN") %>%
  dplyr::filter(date <= IRNdeaths$date_end_survey)
# now scale all recast_deaths_df deaths for this region before we use them
scale_iran <- IRNdeaths %>%
  dplyr::mutate(n_deaths = n_deaths/sum(IRN_recast_df$deaths)) %>%
  .$n_deaths
IRN_recast_df <- IRN_recast_df %>%
  dplyr::mutate(deaths = round(deaths * scale_iran))


IRN.agebands.dat <- process_data3(deaths = IRNdeaths,
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  cumulative = TRUE,
                                  recast_deaths_df = IRN_recast_df,
                                  groupingvar = "ageband",
                                  study_ids = "IRN1",
                                  recast_deaths_geocode = "IRN")

#......................
# MANUAL ADJUSTMENTS
#......................
# see above


#......................
# get rho
#......................
# basic study, so rho is 1
IRN.agebands.dat$rho <- 1

#......................
# save out
#......................
dir.create("data/derived/IRN", recursive = T)
saveRDS(IRN.agebands.dat, "data/derived/IRN/IRN_agebands.RDS")

#............................................................
# Sweden
#...........................................................
#......................
# regions
#......................
SWE.regions.dat <- process_data3(deaths = deathsdf,
                                 population = populationdf,
                                 sero_val = sero_valdf,
                                 seroprev = sero_prevdf,
                                 cumulative = TRUE,
                                 recast_deaths_df = ECDCdf,
                                 groupingvar = "region",
                                 study_ids = "SWE1",
                                 recast_deaths_geocode = "SWE")
#......................
# age bands
#......................
### For age analysis, assume data from the 9 regions, about 70% of the regions, is representative.
SWE.agebands.dat <- process_data3(deaths = deathsdf,
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  cumulative = TRUE,
                                  recast_deaths_df = ECDCdf,
                                  groupingvar = "ageband",
                                  study_ids = "SWE1",
                                  recast_deaths_geocode = "SWE")

#......................
# MANUAL ADJUSTMENTS
#......................
# none
# TODO decide if we should try and infer sens and spec from test listing


#......................
# get rho
#......................
SWE.agebands.dat$rho <- rep(1, length(unique(SWE.agebands.dat$deathsMCMC$region)))
# standarize for age
SWE.regions.dat$rho <-  rep(1, length(unique(SWE.regions.dat$deathsMCMC$region)))


#......................
# save out
#......................
dir.create("data/derived/SWE", recursive = T)
saveRDS(SWE.regions.dat, "data/derived/SWE/SWE_regions.RDS")
saveRDS(SWE.agebands.dat, "data/derived/SWE/SWE_agebands.RDS")

#............................................................
# Great Britian
#...........................................................
# TODO temp location/data for GB
# Process GB death data
gbrdeaths <- readRDS("~/Desktop/GBR_deaths_timeseries.rds") %>%
  dplyr::filter(!is.na(ageband)) # remove missing obs in age bands
# TODO check w/ orderly/Lucy why these are missing
# TODO check why serology start dates for GRB2 have a range?
GBRsero_prevdf <- sero_prevdf %>%
  dplyr::filter(study_id == "GBR2") %>%
  dplyr::mutate(date_start_survey = lubridate::ymd("2020-05-10"),
                date_end_survey = lubridate::ymd("2020-06-05"))


GBRTimeSeries <- gbrdeaths %>%
  dplyr::mutate(ObsDay = factor(ObsDay, levels = c(0:max(ObsDay))), # fix ObsDay/dates lapses
                region = factor(region),
                gender = factor(gender)) %>%
  dplyr::group_by_at(c("ObsDay", "ageband", "region", "gender"), .drop = F) %>%
  dplyr::summarise(
    n_deaths = sum(n_deaths)
  ) %>%
  dplyr::mutate( # now lift back over
    ObsDay = as.integer(as.character(ObsDay)),
    date = min(lubridate::ymd(gbrdeaths$dod)) + ObsDay,
    region = as.character(region),
    gender = as.character(gender),
    age_low = as.numeric( stringr::str_extract(ageband, "[0-9]+(?=\\,)") ),
    age_high = as.numeric( stringr::str_extract(ageband, "(?<=\\,)[0-9]+") ),
    age_high = ifelse(age_high == 120, 999, age_high)
  ) %>%
  dplyr::ungroup(.) %>%
  dplyr::select(c("date", "age_low", "age_high", "region", "gender", "n_deaths")) %>%
  dplyr::arrange(date, age_low, region, gender)

# now liftover for process data function
GBRTimeSeries <- GBRTimeSeries %>%
  dplyr::mutate(
    country = "GBR",
    study_id = "GBR2",
    age_breakdown = 1,
    gender_breakdown = 1,
    for_regional_analysis = 1) %>%
  dplyr::rename(date_start_survey = date) %>%
  dplyr::mutate(date_end_survey = date_start_survey) %>%
  dplyr::select(c("country", "study_id", "age_low", "age_high", "region", "gender", "n_deaths", "date_start_survey", "date_end_survey", "age_breakdown", "for_regional_analysis", "gender_breakdown"))

# population data
GBRpopdf <- readr::read_csv("data/raw/UK_ONS_2016_Population_Data.csv") %>%
  dplyr::filter(study_id == "GBR2") %>%
  dplyr::mutate(region = factor(region, # to match Marc, Lilith, and Co.
                                levels = c("East",
                                           "East Midlands",
                                           "London",
                                           "North East",
                                           "North West",
                                           "South East",
                                           "South West",
                                           "West Midlands",
                                           "Yorkshire and The Humber"),
                                labels = c("East of England",
                                           "Midlands",
                                           "London",
                                           "North East and Yorkshire",
                                           "North West",
                                           "South East",
                                           "South West",
                                           "Midlands",
                                           "North East and Yorkshire"
                                )),
                region = as.character(region),
                country = "GBR",
                study_id = "GBR2",
                age_breakdown = 1,
                gender_breakdown = 1,
                for_regional_analysis = 1) %>%
  dplyr::select(c("country", "study_id", "age_low", "age_high", "region", "gender", "population", "age_breakdown", "for_regional_analysis", "gender_breakdown"))

#......................
# regions
#......................
GBR.regions.dat <- process_data3(deaths = GBRTimeSeries,
                                 population = GBRpopdf,
                                 sero_val = sero_valdf,
                                 seroprev = GBRsero_prevdf,
                                 cumulative = FALSE,
                                 groupingvar = "region",
                                 study_ids = "GBR2")

#......................
# agebands
#......................
GBR.agebands.dat <- process_data3(deaths = GBRTimeSeries,
                                  population = GBRpopdf,
                                  sero_val = sero_valdf,
                                  seroprev = GBRsero_prevdf,
                                  cumulative = FALSE,
                                  groupingvar = "ageband",
                                  study_ids = "GBR2",
                                  death_agebreaks = c(0, 10, 20, 30, 40,
                                                      50, 60, 70, 80, 90, 999))
#......................
# MANUAL ADJUSTMENTS
#......................
gbr_org_seroprev <- GBR.agebands.dat$seroprev_group %>%
  dplyr::select(c("age_low", "age_high", "seroprevalence"))
agebands <- unique(GBR.agebands.dat$deathsMCMC$ageband)
GBR_adj_seroprev <- tibble::tibble(
  ObsDaymin = unique(GBR.agebands.dat$seroprevMCMC$ObsDaymin),
  ObsDaymax = unique(GBR.agebands.dat$seroprevMCMC$ObsDaymax),
  ageband = agebands,
  age_low = as.numeric(stringr::str_split_fixed(agebands, "-[0-9]+", n=2)[,1]),
  age_high = as.numeric(stringr::str_split_fixed(agebands, "[0-9]+-", n=2)[,2]),
  seroprevalence = NA) %>%
  dplyr::arrange(age_low)

# wiggle
GBR_adj_seroprev$seroprevalence <- apply(GBR_adj_seroprev, 1, wiggle_age_matchfun, wiggle = 2, y = gbr_org_seroprev)

# Great Britian seroprevalence and deaths not perfectly aligned.
# Assumptions.
# 1) 0-10 and 10-20 seroprevalence will be equivalent to the 20-30 age group
# 2) 70-80, 80-90, 90+ seroprevalence will be equivalent to the 60-70 age group
GBR_adj_seroprev$seroprevalence[1:2] <- gbr_org_seroprev$seroprevalence[1]
GBR_adj_seroprev$seroprevalence[8:10] <- gbr_org_seroprev$seroprevalence[6]
GBR_adj_seroprev <- GBR_adj_seroprev %>%
  dplyr::rename(SeroPrev = seroprevalence)
# write over
GBR.agebands.dat$seroprevMCMC <- GBR_adj_seroprev

#......................
# get rho
#......................
GBR.agebands.dat$rho <- rep(1, length(unique(GBR.agebands.dat$deathsMCMC$ageband)))
GBR.regions.dat$rho <-  rep(1, length(unique(GBR.regions.dat$deathsMCMC$region)))


#......................
# save out
#......................
dir.create("data/derived/UK/")
saveRDS(GBR.agebands.dat, "data/derived/UK/GBR_agebands.RDS")
saveRDS(GBR.regions.dat, "data/derived/UK/GBR_regions.RDS")


#..................................................................................
#---- Preprocess USA Data  #-----
populationdf <- readr::read_csv("data/raw/USA_County_Demographic_Data.csv") %>%
  tidyr::gather(., key = "strata", value = "population", 3:ncol(.)) %>%
  dplyr::filter(stringr::str_detect(strata, "Both_", negate = TRUE)) %>%
  dplyr::filter(stringr::str_detect(strata, "_Total", negate = TRUE)) %>%
  dplyr::mutate(
    country = "USA",
    Countysp = gsub(" County", "", County),
    Countysp = gsub(" ", "-", Countysp),
    region = paste0(State, "_", Countysp),
    ageband = stringr::str_split_fixed(strata, "[A-Za-z]_", n = 2)[,2],
    ageband = ifelse(stringr::str_detect(ageband, "\\+"),
                     paste0(stringr::str_extract_all(ageband, "[0-9]+", simplify = TRUE), "-", 999),
                     ageband),
    age_low = as.numeric( stringr::str_split_fixed(ageband, "-[0-9]+", n = 2)[,1] ),
    age_high = as.numeric( stringr::str_split_fixed(ageband, "[0-9]-", n = 2)[,2] ),
    gender = stringr::str_extract_all(strata, "[A-Za-z]+", simplify = TRUE)[,1],
    age_breakdown = 1,
    for_regional_analysis = 1,
    gender_breakdown = 1
  ) %>%
  dplyr::select(c("country", "age_low", "age_high", "region", "gender", "population", "age_breakdown", "for_regional_analysis", "gender_breakdown")) %>%
  dplyr::left_join(., readr::read_csv("data/raw/usa_study_id_county_key.csv"), by = "region")

# seroprevalence
sero_valdf <-  readr::read_csv("data/raw/seroassay_validation.csv")
sero_prevdf <- readr::read_csv("data/raw/seroprevalence.csv") %>%
  dplyr::select(-c("ref", "notes")) %>%
  dplyr::mutate(date_start_survey = lubridate::dmy(date_start_survey), # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format
                date_end_survey = lubridate::dmy(date_end_survey))


# deaths
deathsdf <- readr::read_csv("data/raw/deaths.csv") %>%
  dplyr::select(-c("ref", "notes")) %>%
  dplyr::mutate(date_start_survey = lubridate::dmy(date_start_survey), # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format
                date_end_survey = lubridate::dmy(date_end_survey))
JHUdf <- readr::read_csv("data/raw/JHU_time_series_covid19_deaths_US_july72020.csv") %>%
  tidyr::gather(., key = "date", value = "deaths", 13:ncol(.)) %>%
  dplyr::filter(!is.na(Admin2)) %>%
  dplyr::mutate(Admin2sp = sub(" ", "-", Admin2),
                georegion = paste0(Province_State, "_", Admin2sp)) %>%
  dplyr::mutate(date = lubridate::mdy(date)) %>% # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format
  dplyr::group_by(georegion) %>% # group by region for daily deaths
  dplyr::rename(cumdeaths = deaths) %>%
  dplyr::mutate(deaths = cumdeaths - dplyr::lag(cumdeaths),
                deaths = ifelse(is.na(deaths), 0, deaths)) %>% # take care of first value
  dplyr::select(c("date", "georegion", "deaths"))

# fill in from origin
georegions <- dplyr::group_keys(dplyr::group_by(JHUdf, georegion)) %>% dplyr::pull(.)
origindf <- tibble::as_tibble(expand.grid(date = seq(lubridate::ymd("2020-01-01"), min(JHUdf$date), by = "1 day"),
                  georegion = list(georegions),
                  deaths = 0)) %>%
  tidyr::unnest(cols = georegion)
# combine now
JHUdf <- dplyr::bind_rows(JHUdf, origindf) %>%
  dplyr::arrange(georegion, date)


#..................................................................................
#---- Process USA Data #-----
#..................................................................................
#............................................................
# New York City, NY
#...........................................................
NYCJHU <- JHUdf %>%
  dplyr::filter(georegion == "New York_New-York")

#......................
# regions
#......................
# TODO buroughs? -- state level ... ?


#......................
# agebands
#......................
NYC_NY_1.agebands.dat <- process_data3(deaths = deathsdf,
                                       population = populationdf,
                                       sero_val = sero_valdf,
                                       seroprev = sero_prevdf,
                                       cumulative = TRUE,
                                       recast_deaths_df = JHUdf,
                                       groupingvar = "ageband",
                                       study_ids = "NYC_NY_1",
                                       recast_deaths_geocode = "New York_New-York",
                                       death_agebreaks = c(0, 18, 45, 65, 75, 999)) # need this for pop liftover
#......................
# MANUAL ADJUSTMENTS
#......................
# NYC seroprevalence and deaths not perfectly aligned because blood donor data
# Assumptions.
# 1) 18-34 and 34-44 seroprevalence will be averaged for the 0-45 age group
# 2) Seroprev in the 0-45 age group will be equivalent to the 0-18 age group
# 2) Seroprev in the 44-54 age group will be equivalent to the 45-65 age group
# 3) Seroprev in the 54+ age group will be equivalent to the 65-75 and 75+ age group

nyc_adj_seroprev <- tibble::tibble(
  ObsDaymin = unique(NYC_NY_1.agebands.dat$seroprevMCMC$ObsDaymin),
  ObsDaymax = unique(NYC_NY_1.agebands.dat$seroprevMCMC$ObsDaymax),
  ageband = unique(NYC_NY_1.agebands.dat$deathsMCMC$ageband),
  age_low = unique(as.numeric(stringr::str_split_fixed(NYC_NY_1.agebands.dat$deathsMCMC$ageband, "-[0-9]+", n=2)[,1])),
  age_high= unique(as.numeric(stringr::str_split_fixed(NYC_NY_1.agebands.dat$deathsMCMC$ageband, "[0-9]+-", n=2)[,2])),
  seroprevalence = NA) %>%
  dplyr::arrange(age_low)
# lift over
nyc_adj_seroprev$seroprevalence[1:2] <- NYC_NY_1.agebands.dat$seroprevMCMC %>%
  dplyr::filter(ageband %in% c("18-34", "34-44")) %>%
  dplyr::summarise(
    n_positive = sum(n_positive),
    n_tested = sum(n_tested),
    SeroPrev = n_positive/n_tested
  ) %>%
  dplyr::select(c("SeroPrev")) %>%
  unlist(.) %>%
  unname(.)
nyc_adj_seroprev$seroprevalence[3] <- NYC_NY_1.agebands.dat$seroprevMCMC$SeroPrev[3]
nyc_adj_seroprev$seroprevalence[4:5] <- NYC_NY_1.agebands.dat$seroprevMCMC$SeroPrev[4]
nyc_adj_seroprev <- nyc_adj_seroprev %>%
  dplyr::rename(SeroPrev = seroprevalence)


# write over
NYC_NY_1.agebands.dat$seroprevMCMC <- nyc_adj_seroprev

#......................
# get rho
#......................
NYC_NY_1.agebands.dat$rho <- rep(1, length(unique(NYC_NY_1.agebands.dat$deathsMCMC$ageband)))

#......................
# save out
#......................
dir.create("data/derived/USA", recursive = T)
saveRDS(NYC_NY_1.agebands.dat, "data/derived/USA/NYC_NY_1_agebands.RDS")





#..................................................................................
#---- _BASIC_ USA data  #----
# For LA_CA, SC_CA,  CH_MA, MD_FL -
#..................................................................................
#............................................................
# Los Angeles, CA Regional (Basic)
#...........................................................
LACAdeathsdf <- JHUdf %>%
  dplyr::filter(georegion == "California_Los-Angeles") %>%
  dplyr::mutate(
    country = "USA",
    study_id = "LA_CA",
    age_low = 0,
    age_high = 999,
    region = "California_Los-Angeles",
    gender = "both",
    age_breakdown = 0,
    gender_breakdown = 0,
    for_regional_analysis = 1) %>%
  dplyr::rename(date_start_survey = date,
                n_deaths = deaths) %>%
  dplyr::mutate(date_end_survey = date_start_survey)



LA_CA.regions.dat <- process_data3(deaths = LACAdeathsdf,
                                   population = populationdf,
                                   sero_val = sero_valdf,
                                   seroprev = sero_prevdf,
                                   cumulative = FALSE,
                                   groupingvar = "region",
                                   study_ids = "LA_CA",
                                   filtRegions = NULL, # some regions combined in serosurvey
                                   filtGender = NULL,
                                   filtAgeBand = NULL)
#......................
# MANUAL ADJUSTMENTS
#......................
# assume blood group donors are representative
LA_CA.regions.dat$seroprev_group$region <- "California_Los-Angeles"
#......................
# get rho
#......................
# one because basic
LA_CA.regions.dat$rho <- 1

#......................
# save out
#......................
saveRDS(LA_CA.regions.dat, "data/derived/USA/LA_CA_regions.RDS")


#............................................................
# Santa Clara, CA Regional (Basic)
#...........................................................
SCCAdeathsdf <- JHUdf %>%
  dplyr::filter(georegion == "California_Santa-Clara") %>%
  dplyr::mutate(
    country = "USA",
    study_id = "SC_CA",
    age_low = 0,
    age_high = 999,
    region = "California_Santa-Clara",
    gender = "both",
    age_breakdown = 0,
    gender_breakdown = 0,
    for_regional_analysis = 1) %>%
  dplyr::rename(date_start_survey = date,
                n_deaths = deaths) %>%
  dplyr::mutate(date_end_survey = date_start_survey)



SC_CA.regions.dat <- process_data3(deaths = SCCAdeathsdf,
                                   population = populationdf,
                                   sero_val = sero_valdf,
                                   seroprev = sero_prevdf,
                                   cumulative = FALSE,
                                   groupingvar = "region",
                                   study_ids = "SC_CA",
                                   filtRegions = NULL, # some regions combined in serosurvey
                                   filtGender = NULL,
                                   filtAgeBand = NULL)
#......................
# MANUAL ADJUSTMENTS
#......................
SC_CA.regions.dat$seroprev_group$region <- "California_Santa-Clara"
#......................
# get rho
#......................
# one because basic
SC_CA.regions.dat$rho <- 1
#......................
# save out
#......................
saveRDS(SC_CA.regions.dat, "data/derived/USA/SC_CA_regions.RDS")

#............................................................
# Chelsea, MA Regional (Basic)
#...........................................................
### NB matching to Suffolk county may not be quite right (Chelsea is a city)
CHMAdeathsdf <- JHUdf %>%
  dplyr::filter(georegion == "Massachusetts_Suffolk") %>%
  dplyr::mutate(
    country = "USA",
    study_id = "CH_MA",
    age_low = 0,
    age_high = 999,
    region = "Massachusetts_Suffolk",
    gender = "both",
    age_breakdown = 0,
    gender_breakdown = 0,
    for_regional_analysis = 1) %>%
  dplyr::rename(date_start_survey = date,
                n_deaths = deaths) %>%
  dplyr::mutate(date_end_survey = date_start_survey)



CH_MA.regions.dat <- process_data3(deaths = CHMAdeathsdf,
                                   population = populationdf,
                                   sero_val = sero_valdf,
                                   seroprev = sero_prevdf,
                                   cumulative = FALSE,
                                   groupingvar = "region",
                                   study_ids = "CH_MA",
                                   filtRegions = NULL, # some regions combined in serosurvey
                                   filtGender = NULL,
                                   filtAgeBand = NULL)
#......................
# MANUAL ADJUSTMENTS
#......................
CH_MA.regions.dat$seroprev_group$region <- "Massachusetts_Suffolk"

#......................
# get rho
#......................
# one because basic
CH_MA.regions.dat$rho <- 1

#......................
# save out
#......................
saveRDS(CH_MA.regions.dat, "data/derived/USA/CH_MA_regions.RDS")


#............................................................
# Miama, FL Regional (Basic)
#...........................................................
MDFLdeathsdf <- JHUdf %>%
  dplyr::filter(georegion == "Florida_Miami-Dade") %>%
  dplyr::mutate(
    country = "USA",
    study_id = "MD_FL",
    age_low = 0,
    age_high = 999,
    region = "Florida_Miami-Dade",
    gender = "both",
    age_breakdown = 0,
    gender_breakdown = 0,
    for_regional_analysis = 1) %>%
  dplyr::rename(date_start_survey = date,
                n_deaths = deaths) %>%
  dplyr::mutate(date_end_survey = date_start_survey)



MD_FL.regions.dat <- process_data3(deaths = MDFLdeathsdf,
                                   population = populationdf,
                                   sero_val = sero_valdf,
                                   seroprev = sero_prevdf,
                                   cumulative = FALSE,
                                   groupingvar = "region",
                                   study_ids = "MD_FL",
                                   filtRegions = NULL, # some regions combined in serosurvey
                                   filtGender = NULL,
                                   filtAgeBand = NULL)
#......................
# MANUAL ADJUSTMENTS
#......................
# assume this group is representative
MD_FL.regions.dat$seroprev_group$region <- "Florida_Miami-Dade"

#......................
# get rho
#......................
# one because basic
MD_FL.regions.dat$rho <- 1

#......................
# save out
#......................
saveRDS(MD_FL.regions.dat, "data/derived/USA/MD_FL_regions.RDS")


#..................................................................................
#----- Process Latin America Data #-----
#..................................................................................
#............................................................
# Brazil
#...........................................................
#......................
# pre-process Brazil
#......................
bradeaths <- readRDS("data/raw/Brazil_state_age_sex_deaths.rds") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::rename(gender = sex) %>%
  dplyr::group_by(date, region, age, gender) %>%
  dplyr::summarise(
    deaths = sum(count)
  ) %>%
  dplyr::mutate(
    country = "BRA",
    study_id = "BRA1",
    age_low = age,
    age_high = age + 1, # make 1-based for cuts
    age_breakdown = 1,
    gender_breakdown = 1,
    for_regional_analysis = 1) %>%
  dplyr::select(-c("age")) %>%
  dplyr::rename(date_start_survey = date,
                n_deaths = deaths) %>%
  dplyr::mutate(date_end_survey = date_start_survey) %>%
  dplyr::ungroup(.) %>%
  dplyr::arrange(date_start_survey, region)

# seroprevalence
sero_valdf <-  readr::read_csv("data/raw/seroassay_validation.csv")

#TODO eventually put regional seroprev in our seroprev csv
rgnsero_prevdf <- readr::read_csv("data/raw/Brazil_seroprevalence_first_survey.csv") %>%
  dplyr::select(-c(dplyr::starts_with("X"))) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(!is.na(positive)) %>% # CHARLIE, note this drop
  dplyr::group_by(region) %>%
  dplyr::summarise(
    n_tested = sum(tested),
    n_positive = sum(positive)
  ) %>%
  dplyr::mutate(
    seroprevalence_unadjusted = n_positive/n_tested,
    country = "BRA",
    study_id = "BRA1",
    age_low = 0,
    age_high = 999,
    gender = "both",
    seroprevalence_weighted = NA,
    date_start_survey = lubridate::ymd("2020-05-15"),
    date_end_survey = lubridate::ymd("2020-05-22"),
    age_breakdown = 0,
    for_regional_analysis = 1,
    gender_breakdown = 0)

agesero_prevdf <- readr::read_csv("data/raw/seroprevalence.csv") %>%
  dplyr::select(-c("ref", "notes")) %>%
  dplyr::mutate(date_start_survey = lubridate::dmy(date_start_survey), # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format
                date_end_survey = lubridate::dmy(date_end_survey))

bra_populationdf <- readr::read_csv("data/raw/Brazil_2020_Population_Data.csv") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::mutate(
    age_group = stringr::str_replace(age_group, "\\+", "-1000"),
    age_low = as.numeric(stringr::str_split_fixed(age_group, "-[0-9]+", n = 2)[,1]),
    age_high = as.numeric(stringr::str_split_fixed(age_group, "[0-9]-", n = 2)[,2]) - 1 # make 0-based like seroprev
  ) %>%
  dplyr::group_by(region, age_low, age_high, sex) %>%
  dplyr::summarise(population = sum(population)) %>%
  dplyr::rename(gender =sex) %>%
  dplyr::mutate(
    country = "BRA",
    study_id = "BRA1",
    age_breakdown = 1,
    for_regional_analysis = 1,
    gender_breakdown = 1)

#......................
# regions
#......................
BRA.regions.dat <- process_data3(deaths = bradeaths,
                                 population = bra_populationdf,
                                 sero_val = sero_valdf,
                                 seroprev = rgnsero_prevdf,
                                 cumulative = FALSE,
                                 groupingvar = "region",
                                 study_ids = "BRA1")


#......................
# ages
#......................
BRA.agebands.dat <- process_data3(deaths = bradeaths,
                                  population = bra_populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = agesero_prevdf,
                                  cumulative = FALSE,
                                  groupingvar = "ageband",
                                  study_ids = "BRA1",
                                  filtRegions = NULL,
                                  filtGender = NULL,
                                  filtAgeBand = NULL,
                                  death_agebreaks = c(0, 4, 9,
                                                      19, 29, 39,
                                                      49, 59, 69,
                                                      79, 999))



#......................
# MANUAL ADJUSTMENTS
#......................
# none because we had line list

#......................
# get rho
#......................
BRA.regions.dat$rho <- rep(1, length(unique(BRA.regions.dat$deathsMCMC$region)))
BRA.agebands.dat$rho <- rep(1, length(unique(BRA.agebands.dat$deathsMCMC$ageband)))

#......................
# save out
#......................
dir.create("data/derived/BRA/", recursive = T)
saveRDS(BRA.regions.dat, "data/derived/BRA/BRA_regions.RDS")
saveRDS(BRA.agebands.dat, "data/derived/BRA/BRA_agebands.RDS")


















#
