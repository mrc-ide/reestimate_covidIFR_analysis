#................................................................................................
## Purpose: Process Country datasets for analysis
##
## Notes: Process data functions and files are specific to regions
#................................................................................................
library(tidyverse)
source("R/process_data3.R")

#......................
# global data
#......................
# serovalidation
sero_valdf <-  readr::read_tsv("data/raw/serovalidation_final_raw.tsv")
# seroprevalence
sero_prevdf <- readr::read_tsv("data/raw/seroprevalence_final_raw.tsv") %>%
  dplyr::select(-c("ref", "notes")) %>%
  dplyr::mutate(date_start_survey = lubridate::ymd(date_start_survey), # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format
                date_end_survey = lubridate::ymd(date_end_survey),
                seroprevalence_unadjusted = ifelse(is.na(seroprevalence_unadjusted), n_positive/n_tested, seroprevalence_unadjusted))


#..................................................................................
#---- Preprocess Latin America Data  #-----
#..................................................................................

#............................................................
#---- BRA1 #----
#...........................................................
#......................
# pre-process Brazil
#......................
bradeaths <- readRDS("data/raw/Brazil_state_age_sex_deaths.rds") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::rename(gender = sex) %>%
  dplyr::mutate(age = ifelse(age >= 100, 100, age), # liftover some very old ages -- capping at 100
                age = factor(age, levels = c(0:100), labels = c(0:100)),
                tempday = as.numeric(lubridate::ymd(date) - lubridate::ymd("2020-01-01")) + 1,
                tempday = factor(tempday, levels = c(1:max(tempday)))) %>% # ugly code to fill in dates
  dplyr::select(-c("date")) %>%
  dplyr::group_by(tempday, region, age, gender, .drop = FALSE) %>%
  dplyr::summarise(
    deaths = sum(count)
  ) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate(
    date = lubridate::ymd("2020-01-01") + as.numeric(as.character(tempday)) - 1,
    age = as.numeric(as.character(age)),
    country = "BRA",
    study_id = "BRA1",
    age_low = age,
    age_high = age + 1, # make 1-based for cuts
    age_breakdown = 1,
    gender_breakdown = 1,
    for_regional_analysis = 1) %>%
  dplyr::select(-c("age", "tempday")) %>%
  dplyr::rename(date_start_survey = date,
                n_deaths = deaths) %>%
  dplyr::mutate(date_end_survey = date_start_survey) %>%
  dplyr::ungroup(.) %>%
  dplyr::arrange(date_start_survey, region)

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
    gender_breakdown = 1) %>%
  dplyr::ungroup(.)

#......................
# regions
#......................
BRA.regions.dat <- process_data3(deaths = bradeaths,
                                 population = bra_populationdf,
                                 sero_val = sero_valdf,
                                 seroprev = sero_prevdf,
                                 cumulative = FALSE,
                                 groupingvar = "region",
                                 study_ids = "BRA1",
                                 origin = lubridate::ymd("2020-01-01"))


#......................
# ages
#......................
BRA.agebands.dat <- process_data3(deaths = bradeaths,
                                  population = bra_populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  cumulative = FALSE,
                                  groupingvar = "ageband",
                                  study_ids = "BRA1",
                                  origin = lubridate::ymd("2020-01-01"),
                                  death_agebreaks = c(0, 4, 9,
                                                      19, 29, 39,
                                                      49, 59, 69,
                                                      79, 999))
#......................
# MANUAL ADJUSTMENTS
#......................
# seroprevalence not absolutely 0
BRA.regions.dat$seroprevMCMC <- BRA.regions.dat$seroprevMCMC %>%
  dplyr::mutate(SeroPrev = ifelse(SeroPrev == 0, 1e-10, SeroPrev))

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



#..................................................................................
#---- Preprocess European Data #----
#..................................................................................
# deaths
deathsdf <- readr::read_csv("data/raw/cumulative_deaths.csv") %>%
  dplyr::select(-c("ref", "notes")) %>%
  dplyr::mutate(date_start_survey = lubridate::ymd(date_start_survey), # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format
                date_end_survey = lubridate::ymd(date_end_survey))
ECDCdf <- readr::read_csv("data/raw/daily_deaths_ECDC20200724.csv") %>%
  dplyr::select(c("dateRep", "countryterritoryCode", "deaths")) %>%
  dplyr::rename(date = dateRep,
                georegion = countryterritoryCode) %>%
  dplyr::mutate(date = lubridate::mdy(date), # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format
                deaths = ifelse(deaths < 0, 0, deaths)) # remove deaths typos
# demography
populationdf <- readr::read_tsv("data/raw/population.tsv") %>%
  dplyr::select(-c("reference")) %>%
  dplyr::mutate(age_low = ifelse(age_low == 0 & age_high == 0, 1, age_low),
                age_high = ifelse(age_low == 1 & age_high == 0, 1, age_high))  # liftover "zero" year olds to be 1, 1 as well

#............................................................
#---- CHE1 #----
#...........................................................
#......................
# ages
#......................
# TODO check on this CHE timeseries
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
CHE.region.dat <- process_data3(deaths = deathsdf,
                                population = populationdf,
                                sero_val = sero_valdf,
                                seroprev = sero_prevdf,
                                cumulative = TRUE,
                                recast_deaths_df = CHE1TimeSeries,
                                groupingvar = "region",
                                study_ids = "CHE1",
                                death_agebreaks = c(0, 999), # for pop splits
                                recast_deaths_geocode = "Geneva")   ## use study id in case we get more studies later.


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
                                  death_agebreaks = c(0, 10, 20, 30,
                                                      40, 50, 60, 70,
                                                      80, 999), # for pop splits
                                  recast_deaths_geocode = "Geneva")   ## use study id in case we get more studies later.

#......................
# MANUAL ADJUSTMENTS
#......................
# Assume:
# (1) seroprevalence 5-9 is representative of ages 0-10
# (2) seroprevalence 10-19 is representative of ages 10-20
# (3) seroprevalence 20-49 is representative of ages 20-30 and 30-40 and 40-50
# (4) seroprevalence 50-64 is representative of ages 50-60 and 60-70
# (5) seroprevalence 65+ is representative of ages 70-80 and 80+
ageband <- unique(CHE.agebands.dat$deathsMCMC$ageband)
che_adj_seroprev <- tibble::tibble(
  ObsDaymin = unique(CHE.agebands.dat$seroprevMCMC$ObsDaymin),
  ObsDaymax = unique(CHE.agebands.dat$seroprevMCMC$ObsDaymax))
che_adj_seroprev <- tidyr::expand_grid(che_adj_seroprev, ageband) %>%
  dplyr::mutate(age_low = as.numeric(stringr::str_split_fixed(ageband, "-[0-9]+", n=2)[,1]),
                age_high= as.numeric(stringr::str_split_fixed(ageband, "[0-9]+-", n=2)[,2])) %>%
  dplyr::arrange(age_low)

# pull out original
che_org_seroprev <- CHE.agebands.dat$seroprev_group %>%
  dplyr::select(c("ObsDaymin", "ObsDaymax", "age_low", "age_high", "seroprevalence_unadjusted"))

# because of multiple dates, this is a slightly harder wiggle than below
# assumption 1
che_org_seroprev_assum1 <- che_org_seroprev %>%
  dplyr::filter(age_low == 5 & age_high == 9) %>%
  dplyr::mutate(age_low = 0,
                age_high = 10)
# assumption 2
che_org_seroprev_assum2 <- che_org_seroprev %>%
  dplyr::filter(age_low == 10 & age_high == 19) %>%
  dplyr::mutate(age_low = 10,
                age_high = 20)
# assumption 3
che_org_seroprev_assum3.1 <- che_org_seroprev %>%
  dplyr::filter(age_low == 20 & age_high == 49) %>%
  dplyr::mutate(age_low = 20,
                age_high = 30)
che_org_seroprev_assum3.2 <- che_org_seroprev %>%
  dplyr::filter(age_low == 20 & age_high == 49) %>%
  dplyr::mutate(age_low = 30,
                age_high = 40)
che_org_seroprev_assum3.3 <- che_org_seroprev %>%
  dplyr::filter(age_low == 20 & age_high == 49) %>%
  dplyr::mutate(age_low = 40,
                age_high = 50)
che_org_seroprev_assum3 <- dplyr::bind_rows(che_org_seroprev_assum3.1, che_org_seroprev_assum3.2, che_org_seroprev_assum3.3)
# assumption 4
che_org_seroprev_assum4.1 <- che_org_seroprev %>%
  dplyr::filter(age_low == 50 & age_high == 64) %>%
  dplyr::mutate(age_low = 50,
                age_high = 60)
che_org_seroprev_assum4.2 <- che_org_seroprev %>%
  dplyr::filter(age_low == 50 & age_high == 64) %>%
  dplyr::mutate(age_low = 60,
                age_high = 70)
che_org_seroprev_assum4 <- dplyr::bind_rows(che_org_seroprev_assum4.1, che_org_seroprev_assum4.2)
# assumption 5
che_org_seroprev_assum5.1 <- che_org_seroprev %>%
  dplyr::filter(age_low == 65 & age_high == 999) %>%
  dplyr::mutate(age_low = 70,
                age_high = 80)
che_org_seroprev_assum5.2 <- che_org_seroprev %>%
  dplyr::filter(age_low == 65 & age_high == 999) %>%
  dplyr::mutate(age_low = 80,
                age_high = 999)
che_org_seroprev_assum5 <- dplyr::bind_rows(che_org_seroprev_assum5.1, che_org_seroprev_assum5.2)
# bring toghether this monstrosity
che_org_seroprev <- dplyr::bind_rows(che_org_seroprev_assum1,
                                     che_org_seroprev_assum2,
                                     che_org_seroprev_assum3,
                                     che_org_seroprev_assum4,
                                     che_org_seroprev_assum5)

che_adj_seroprev <- dplyr::left_join(che_adj_seroprev, che_org_seroprev,
                                     by = c("ObsDaymin", "ObsDaymax", "age_low", "age_high")) %>%
  dplyr::rename(SeroPrev = seroprevalence_unadjusted) %>%
  dplyr::select(c("ObsDaymin", "ObsDaymax", "ageband", "SeroPrev")) %>%
  dplyr::arrange(ObsDaymin, ObsDaymax, ageband)

# bring together
CHE.agebands.dat$seroprevMCMC <- che_adj_seroprev

# seroprevalence not absolutely 0
CHE.agebands.dat$seroprevMCMC <- CHE.agebands.dat$seroprevMCMC %>%
  dplyr::mutate(SeroPrev = ifelse(SeroPrev == 0, 1e-10, SeroPrev))


#......................
# get rho
#......................
CHE.agebands.dat$rho <- rep(1, length(unique(che_adj_seroprev$ageband)))

#......................
# save out
#......................
dir.create("data/derived/CHE", recursive = T)
saveRDS(CHE.region.dat, "data/derived/CHE/CHE_region.RDS")
saveRDS(CHE.agebands.dat, "data/derived/CHE/CHE_agebands.RDS")


#............................................................
#---- DNK1 #----
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
# Denmark deaths agebands do not overlap well
# will take mean for the 0-59 age group for blood donors less then 60
# will assume > 60 for rest
agebands <- unique(DNK.agebands.dat$deathsMCMC$ageband)
dnk_adj_seroprev <- tibble::tibble(
  ObsDaymin = unique(DNK.agebands.dat$seroprevMCMC$ObsDaymin),
  ObsDaymax = unique(DNK.agebands.dat$seroprevMCMC$ObsDaymax),
  ageband = agebands,
  SeroPrev = NA)
dnk_adj_seroprev$SeroPrev[1] <- mean(DNK.agebands.dat$seroprev_group$seroprevalence_unadjusted[1:4])
dnk_adj_seroprev$SeroPrev[2:5] <- DNK.agebands.dat$seroprev_group$seroprevalence_unadjusted[5]
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
#--- ESP1-2 #----
#...........................................................
sero_prevdfESP <- sero_prevdf %>%
  dplyr::mutate(age_high = ifelse(study_id == "ESP1-2" & age_low == 0 & age_high == 0, 0.99, # to make the cut easier
                                   age_high))

#......................
# regions
#......................
ESP.regions.dat <- process_data3(deaths = deathsdf,
                                 population = populationdf,
                                 sero_val = sero_valdf,
                                 seroprev = sero_prevdfESP,
                                 cumulative = TRUE,
                                 recast_deaths_df = ECDCdf,
                                 groupingvar = "region",
                                 study_ids = "ESP1-2",
                                 recast_deaths_geocode = "ESP",
                                 filtRegions = c("Andalucía", "Aragón", "Asturias, Principado de",  "Balears, Illes",  "Comunitat Valenciana",
                                                 "Canarias", "Cantabria", "Castilla La Mancha", "Castilla y León", "Cataluña", "Extremadura",
                                                 "Galicia",   "Rioja, La",  "Madrid, Comunidad de",
                                                 "Murcia, Región de", "Navarra, Comunidad Foral de", "País Vasco")) # limit to mainland Spain



#......................
# agebands
#......................
ESP.agebands.dat <- process_data3(deaths = deathsdf,
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdfESP,
                                  cumulative = TRUE,
                                  recast_deaths_df = ECDCdf,
                                  groupingvar = "ageband",
                                  study_ids = "ESP1-2",
                                  recast_deaths_geocode = "ESP",
                                  death_agebreaks = c(0, 10, 20, 30, 40,
                                                      50, 60, 70, 80, 90, 999),
                                  sero_agebreaks = c(0, 10, 20, 30, 40,
                                                     50, 60, 70, 80, 90, 999))
#......................
# MANUAL ADJUSTMENTS
#......................
# # suspicious 0 followed by very large increase
# # Monday, April 27 appears suspicious, such that all deaths from Apr 27 got pushed to Apr 28
# ESP.agebands.dat$deathsMCMC$Deaths[ESP.agebands.dat$deathsMCMC$ObsDay == 117] <- -1
# ESP.agebands.dat$deathsMCMC$Deaths[ESP.agebands.dat$deathsMCMC$ObsDay == 118] <- -1
# ESP.regions.dat$deathsMCMC$Deaths[ESP.regions.dat$deathsMCMC$ObsDay == 117] <- -1
# ESP.regions.dat$deathsMCMC$Deaths[ESP.regions.dat$deathsMCMC$ObsDay == 118] <- -1
#
# # typo for 1179 deaths on May 25
# ESP.agebands.dat$deathsMCMC$Deaths[ESP.agebands.dat$deathsMCMC$ObsDay == 171] <- -1
# ESP.regions.dat$deathsMCMC$Deaths[ESP.regions.dat$deathsMCMC$ObsDay == 171] <- -1


# No adjustments to serology as there is alignment of deaths and seroprevalence age groups.

#......................
# get rho
#......................
ESP.agebands.dat$rho <- rep(1, length(unique(ESP.agebands.dat$deathsMCMC$ageband)))
ESP.regions.dat$rho <- rep(1, length(unique(ESP.regions.dat$deathsMCMC$region)))


#......................
# save out
#......................
dir.create("data/derived/ESP/ESP_agebands.RDS")
saveRDS(ESP.agebands.dat, "data/derived/ESP/ESP_agebands.RDS")
saveRDS(ESP.regions.dat, "data/derived/ESP/ESP_regions.RDS")

#............................................................
#---- GBR2 #----
#...........................................................

# population data
GBR2popdf <- readr::read_csv("data/raw/UK_ONS_2016_Population_Data.csv") %>%
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

#......................
# agebands
#......................

#......................
# MANUAL ADJUSTMENTS
#......................


#......................
# save out
#......................


#............................................................
#---- NLD1 #-----
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
  age_high = as.numeric(stringr::str_split_fixed(agebands, "[0-9]+-", n=2)[,2]),
  seroprevalence = NA) %>%
  dplyr::arrange(age_low)

nld_org_seroprev <- NLD.agebands.dat$seroprev_group %>%
  dplyr::select(c("age_low", "age_high", "seroprevalence_unadjusted")) %>%
  dplyr::rename(seroprevalence = seroprevalence_unadjusted)

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
NLD.regions.dat$seroprevMCMC$SeroPrev <- (NLD.regions.dat$seroprev_group$range_sero_low + NLD.regions.dat$seroprev_group$range_sero_high)/2
nldmean <- mean(NLD.regions.dat$seroprevMCMC$SeroPrev, na.rm = T)
NLD.regions.dat$seroprevMCMC <- NLD.regions.dat$seroprevMCMC %>%
  dplyr::mutate(SeroPrev = ifelse(is.na(SeroPrev), nldmean, SeroPrev))


# seroprevalence not absolutely 0
NLD.regions.dat$seroprevMCMC <- NLD.regions.dat$seroprevMCMC %>%
  dplyr::mutate(SeroPrev = ifelse(SeroPrev == 0, 1e-10, SeroPrev))


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
#--- SWE1 #----
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
#---- ITA1 #----
#...........................................................
#......................
# regions
#......................
ITA.regions.dat <- process_data3(deaths = deathsdf,
                                 population = populationdf,
                                 sero_val = sero_valdf,
                                 seroprev = sero_prevdf,
                                 cumulative = TRUE,
                                 recast_deaths_df = ECDCdf,
                                 groupingvar = "region",
                                 study_ids = "ITA1",
                                 recast_deaths_geocode = "ITA")

#......................
# ages
#......................
ITA.agebands.dat <- process_data3(deaths = deathsdf,
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  cumulative = TRUE,
                                  recast_deaths_df = ECDCdf,
                                  groupingvar = "ageband",
                                  study_ids = "ITA1",
                                  recast_deaths_geocode = "ITA",
                                  death_agebreaks = c(0,9,19,29,39,49,59,69,79,89, 999),
                                  sero_agebreaks = c(0, 17,  34,  49,  59,  69, 999),
                                  filtGender = "both")

#......................
# MANUAL ADJUSTMENTS
#......................
## TODO - decide prior with Nick (sensitivity>90%, specificity>95%). Just chose some values for now.
ITA.agebands.dat$sero_sens$sensitivity<-0.9  ## unknown
ITA.agebands.dat$sero_spec$specificity <-0.975
ITA.regions.dat$sero_sens$sensitivity<-0.9  ## unknown
ITA.regions.dat$sero_spec$specificity <-0.975


# some lack of overlap between serology and deaths data.
# use sero 0-17 for 9-19 yr olds
# use the mean of 18-34 and 35-49 for the 30-39 group (they are v similar anyway)
#
agebands <- unique(ITA.agebands.dat$deathsMCMC$ageband)
ita_adj_seroprev <- tibble::tibble(
  ObsDaymin = unique(ITA.agebands.dat$seroprevMCMC$ObsDaymin),
  ObsDaymax = unique(ITA.agebands.dat$seroprevMCMC$ObsDaymax),
  ageband = agebands,
  SeroPrev = NA)
ita_adj_seroprev$SeroPrev[1:2] <- ITA.agebands.dat$seroprev_group$seroprevalence_unadjusted[1]
ita_adj_seroprev$SeroPrev[3] <- ITA.agebands.dat$seroprev_group$seroprevalence_unadjusted[2]
ita_adj_seroprev$SeroPrev[4] <- mean(ITA.agebands.dat$seroprev_group$seroprevalence_unadjusted[2:3])
ita_adj_seroprev$SeroPrev[5] <- ITA.agebands.dat$seroprev_group$seroprevalence_unadjusted[3]
ita_adj_seroprev$SeroPrev[6:7] <- ITA.agebands.dat$seroprev_group$seroprevalence_unadjusted[4:5]
ita_adj_seroprev$SeroPrev[8:10] <- ITA.agebands.dat$seroprev_group$seroprevalence_unadjusted[6]
ITA.agebands.dat$seroprevMCMC <- ita_adj_seroprev

#......................
# get rho
#......................
ITA.agebands.dat$rho <- rep(1, length(unique(ITA.agebands.dat$deathsMCMC$ageband)))
# multiple through demog and age-standardize for region
ITA.regions.dat$rho <- rep(1, length(unique(ITA.regions.dat$deathsMCMC$region)))

#......................
# save out
#......................
dir.create("data/derived/ITA", recursive = T)
saveRDS(ITA.regions.dat, "data/derived/ITA/ITA_regions.RDS")
saveRDS(ITA.agebands.dat, "data/derived/ITA/ITA_agebands.RDS")





#............................................................
#---- LUX1 #----
#...........................................................
#......................

#......................
# ages
#......................
LUX.agebands.dat <- process_data3(deaths = deathsdf,
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  cumulative = TRUE,
                                  recast_deaths_df = ECDCdf,
                                  groupingvar = "ageband",
                                  study_ids = "LUX1",
                                  recast_deaths_geocode = "LUX",
                                  death_agebreaks = c(0, 29,39,49,59,69,79, 999),
                                  sero_agebreaks = NULL,
                                  filtGender = "both")


#......................
# MANUAL ADJUSTMENTS
#......................
# Only one ageband - assume all the same.
agebands <- unique(LUX.agebands.dat$deathsMCMC$ageband)
lux_adj_seroprev <- tibble::tibble(
  ObsDaymin = unique(LUX.agebands.dat$seroprevMCMC$ObsDaymin),
  ObsDaymax = unique(LUX.agebands.dat$seroprevMCMC$ObsDaymax),
  ageband = agebands,
  SeroPrev = NA)
lux_adj_seroprev$SeroPrev <- LUX.agebands.dat$seroprev_group$seroprevalence_unadjusted
LUX.agebands.dat$seroprevMCMC <- lux_adj_seroprev

#......................
# get rho
#......................
LUX.agebands.dat$rho <- rep(1, length(unique(LUX.agebands.dat$deathsMCMC$ageband)))

#......................
# save out
#......................
dir.create("data/derived/LUX", recursive = T)
saveRDS(LUX.agebands.dat, "data/derived/LUX/LUX_agebands.RDS")


#..................................................................................
#---- Preprocess USA Data  #-----
#..................................................................................
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


#............................................................
#---- LA_CA1 #----
# Los Angeles, CA Regional (Basic)
#...........................................................
LACAdeathsdf <- JHUdf %>%
  dplyr::filter(georegion %in% "California_Los-Angeles") %>%
  dplyr::mutate(
    country = "USA",
    study_id = "LA_CA1",
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
#---- NYC_NY_1 #-----
#...........................................................
NYCJHU <- JHUdf %>%
  dplyr::filter(georegion == "New York_New-York")
NYpopdf <- populationdf %>%
  dplyr::filter(study_id == "NYC_NY_1")
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
# # NYC seroprevalence and deaths not perfectly aligned because blood donor data
# # Assumptions.
# # 1) 18-34 and 34-44 seroprevalence will be averaged for the 0-18 and 18-45 age group
# # 2) Seroprev in the 44-54 age group will be equivalent to the 45-65 age group
# # 3) Seroprev in the 54+ age group will be equivalent to the 65-75 and 75+ age group
#
# nyc_adj_seroprev <- tibble::tibble(
#   ObsDaymin = unique(NYC_NY_1.agebands.dat$seroprevMCMC$ObsDaymin),
#   ObsDaymax = unique(NYC_NY_1.agebands.dat$seroprevMCMC$ObsDaymax),
#   ageband = unique(NYC_NY_1.agebands.dat$deathsMCMC$ageband),
#   age_low = unique(as.numeric(stringr::str_split_fixed(NYC_NY_1.agebands.dat$deathsMCMC$ageband, "-[0-9]+", n=2)[,1])),
#   age_high= unique(as.numeric(stringr::str_split_fixed(NYC_NY_1.agebands.dat$deathsMCMC$ageband, "[0-9]+-", n=2)[,2])),
#   seroprevalence = NA) %>%
#   dplyr::arrange(age_low)
# # lift over
# nyc_adj_seroprev$seroprevalence[1:2] <- NYC_NY_1.agebands.dat$seroprevMCMC %>%
#   dplyr::filter(ageband %in% c("18-34", "34-44")) %>%
#   dplyr::summarise(
#     n_positive = sum(n_positive),
#     n_tested = sum(n_tested),
#     SeroPrev = n_positive/n_tested
#   ) %>%
#   dplyr::select(c("SeroPrev")) %>%
#   unlist(.) %>%
#   unname(.)
# nyc_adj_seroprev$seroprevalence[3] <- NYC_NY_1.agebands.dat$seroprevMCMC$SeroPrev[3]
# nyc_adj_seroprev$seroprevalence[4:5] <- NYC_NY_1.agebands.dat$seroprevMCMC$SeroPrev[4]
# nyc_adj_seroprev <- nyc_adj_seroprev %>%
#   dplyr::rename(SeroPrev = seroprevalence)
#

#......................
# adding in CDC_1 data as well
#......................
cdc1 <- sero_prevdf %>%
  dplyr::filter(study_id == "CDC_1" & stringr::str_detect(region, "New York") & gender == "both") %>%
  dplyr::filter(!(age_low == 0 & age_high == 999)) %>%
  dplyr::rename(ObsDaymin = date_start_survey,
                ObsDaymax = date_end_survey) %>%
  dplyr::mutate(ObsDaymin = as.numeric(ObsDaymin - lubridate::ymd("20200101")),
                ObsDaymax = as.numeric(ObsDaymax - lubridate::ymd("20200101"))) %>%
  dplyr::select(c("ObsDaymin", "ObsDaymax", "age_low", "age_high", "seroprevalence_unadjusted")) %>%
  dplyr::rename(seroprevalence = seroprevalence_unadjusted)

nyc_adj_seroprev2 <- tibble::tibble(
  ObsDaymin = unique(cdc1$ObsDaymin),
  ObsDaymax = unique(cdc1$ObsDaymax),
  ageband = unique(NYC_NY_1.agebands.dat$deathsMCMC$ageband),
  age_low = unique(as.numeric(stringr::str_split_fixed(NYC_NY_1.agebands.dat$deathsMCMC$ageband, "-[0-9]+", n=2)[,1])),
  age_high= unique(as.numeric(stringr::str_split_fixed(NYC_NY_1.agebands.dat$deathsMCMC$ageband, "[0-9]+-", n=2)[,2])),
  seroprevalence = NA) %>%
  dplyr::arrange(age_low)

# Assumptions
# 1) 0-18 matches 0-18
# 2) 18-45 matches 19-49
# 3) 50-64 matches 45-65
# 4) 65-75 and 75+ matches 65+

nyc_adj_seroprev2$seroprevalence <- apply(nyc_adj_seroprev2, 1, wiggle_age_matchfun, wiggle = 5, y = cdc1)
nyc_adj_seroprev2 <- nyc_adj_seroprev2 %>%
  dplyr::rename(SeroPrev = seroprevalence)


# write over
#NYC_NY_1.agebands.dat$seroprevMCMC <- dplyr::bind_rows(nyc_adj_seroprev2, nyc_adj_seroprev)
NYC_NY_1.agebands.dat$seroprevMCMC <- nyc_adj_seroprev2

#......................
# get rho
#......................
NYC_NY_1.agebands.dat$rho <- rep(1, length(unique(NYC_NY_1.agebands.dat$deathsMCMC$ageband)))

#......................
# save out
#......................
dir.create("data/derived/USA", recursive = T)
saveRDS(NYC_NY_1.agebands.dat, "data/derived/USA/NYC_NY_1_cdc1_agebands.RDS")


#............................................................
#--- SF_CA1 #----
# San Francisco Bay Area, CA (Basic)
#...........................................................
bay_area_vec <- c("California_Sonoma", "California_Marin", "California_Napa", "California_Contra-Costa",
                  "California_Alameda", "California_Santa-Clara", "California_San-Mateo", "California_Sacramento",
                  "California_San-Joaquin")
SF_CAdeathsdf <- JHUdf %>%
  dplyr::filter(georegion %in% bay_area_vec) %>%
  dplyr::mutate(
    country = "USA",
    study_id = "SF_CA1",
    age_low = 0,
    age_high = 999,
    region = "California_San-Francisco",
    gender = "both",
    age_breakdown = 0,
    gender_breakdown = 0,
    for_regional_analysis = 1) %>%
  dplyr::rename(date_start_survey = date,
                n_deaths = deaths) %>%
  dplyr::mutate(date_end_survey = date_start_survey)

SF_CA.regions.dat <- process_data3(deaths = SF_CAdeathsdf,
                                   population = populationdf,
                                   sero_val = sero_valdf,
                                   seroprev = sero_prevdf,
                                   cumulative = FALSE,
                                   groupingvar = "region",
                                   study_ids = "SF_CA1",
                                   filtRegions = NULL, # some regions combined in serosurvey
                                   filtGender = NULL,
                                   filtAgeBand = NULL)
#......................
# MANUAL ADJUSTMENTS
#......................
# assume blood group donors are representative
SF_CA.regions.dat$seroprev_group$region <- "California_San-Francisco"
#......................
# get rho
#......................
# one because basic
SF_CA.regions.dat$rho <- 1

#......................
# save out
#......................
saveRDS(SF_CA.regions.dat, "data/derived/USA/SF_CA_regions.RDS")



#..................................................................................
#---- Preprocess Asian Data  #-----
#..................................................................................
# demography
populationdf <- readr::read_tsv("data/raw/population.tsv") %>%
  dplyr::select(-c("reference")) %>%
  dplyr::mutate(age_low = ifelse(age_low == 0 & age_high == 0, 1, age_low),
                age_high = ifelse(age_low == 1 & age_high == 0, 1, age_high))  # liftover "zero" year olds to be 1, 1 as well

#............................................................
#---- CHN1 #----
#...........................................................
#......................
CHN1TimeSeries <- readr::read_csv("data/raw/deaths_time_series.csv") %>%
  dplyr::filter(study_id == "CHN1")
# sanity check
CHN1TimeSeries <- CHN1TimeSeries %>%
  dplyr::rename(date = date_end_survey,
                deaths = n_deaths) %>%
  dplyr::mutate(date = lubridate::dmy(date)) # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format
CHN1TimeSeries$georegion <- "CHN"
#......................
# ages
#......................
CHN.agebands.dat <- process_data3(deaths = deathsdf,
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  cumulative = TRUE,
                                  recast_deaths_df = CHN1TimeSeries,
                                  groupingvar = "ageband",
                                  study_ids = "CHN1",
                                  recast_deaths_geocode = "CHN",
                                  death_agebreaks = c(0,9,19, 29,39,49,59,69,79, 999),
                                  sero_agebreaks = c(0,19,  24,  29,  34,  39,  44,  49,  54,  59,999),
                                  filtGender = "both")


#......................
# MANUAL ADJUSTMENTS
#......................
#
agebands <- unique(CHN.agebands.dat$deathsMCMC$ageband)
chn_adj_seroprev <- tibble::tibble(
  ObsDaymin = unique(CHN.agebands.dat$seroprevMCMC$ObsDaymin),
  ObsDaymax = unique(CHN.agebands.dat$seroprevMCMC$ObsDaymax),
  ageband = agebands,
  SeroPrev = NA)
chn_adj_seroprev$SeroPrev[1:2] <- mean(CHN.agebands.dat$seroprev_group$seroprevalence_unadjusted)
chn_adj_seroprev$SeroPrev[3] <- sum(CHN.agebands.dat$seroprev_group$n_positive[2:3]) /sum(CHN.agebands.dat$seroprev_group$n_tested[2:3])
chn_adj_seroprev$SeroPrev[4] <- sum(CHN.agebands.dat$seroprev_group$n_positive[4:5]) /sum(CHN.agebands.dat$seroprev_group$n_tested[4:5])
chn_adj_seroprev$SeroPrev[5] <- sum(CHN.agebands.dat$seroprev_group$n_positive[6:7]) /sum(CHN.agebands.dat$seroprev_group$n_tested[6:7])
chn_adj_seroprev$SeroPrev[6] <- sum(CHN.agebands.dat$seroprev_group$n_positive[8:9]) /sum(CHN.agebands.dat$seroprev_group$n_tested[8:9])
chn_adj_seroprev$SeroPrev[7:9] <- sum(CHN.agebands.dat$seroprev_group$n_positive[10]) /sum(CHN.agebands.dat$seroprev_group$n_tested[10])
CHN.agebands.dat$seroprevMCMC <- chn_adj_seroprev

#......................
# get rho
#......................
CHN.agebands.dat$rho <- rep(1, length(unique(CHN.agebands.dat$deathsMCMC$ageband)))

#......................
# save out
#......................
dir.create("data/derived/CHN", recursive = T)
saveRDS(CHN.agebands.dat, "data/derived/CHN/CHN_agebands.RDS")














