#................................................................................................
## Purpose: Process Country datasets for analysis
##
## Notes: Process data functions and files are specific to regions
#................................................................................................
library(tidyverse)
source("R/process_data2.R")
source("R/contact_mat_helpers.R")
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
  dplyr::mutate(date = lubridate::dmy(date)) # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format

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
                date_end_survey = lubridate::dmy(date_end_survey))



#..................................................................................
#----- Process Eurasia Data #-----
#..................................................................................
#............................................................
# Spain
#...........................................................
#......................
# regions
#......................
ESP.regions.dat <- process_data2(deaths = deathsdf,
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
ESP.agebands.dat <- process_data2(deaths = deathsdf,
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  cumulative = TRUE,
                                  recast_deaths_df = ECDCdf,
                                  groupingvar = "ageband",
                                  study_ids = "ESP1",
                                  recast_deaths_geocode = "ESP",
                                  filtRegions = NULL,
                                  filtGender = NULL,
                                  filtAgeBand = c("0-10", "10-20", "20-30",
                                                  "30-40", "40-50", "50-60",
                                                  "60-70", "70-80", "80-90", "90-999"),
                                  sero_agebreaks = c(0, 10, 20, 30, 40,
                                                     50, 60, 70, 80, 90, 999))
#......................
# MANUAL ADJUSTMENTS
#......................
# suspicious 0 followed by very large increase
# Monday, April 27 appears suspicious, such that all deaths from Apr 27 got pushed to Apr 28
ESP.agebands.dat$deaths$Deaths[ESP.agebands.dat$deaths$ObsDay == 117] <- -1
ESP.agebands.dat$deaths$Deaths[ESP.agebands.dat$deaths$ObsDay == 118] <- -1
ESP.regions.dat$deaths$Deaths[ESP.regions.dat$deaths$ObsDay == 117] <- -1
ESP.regions.dat$deaths$Deaths[ESP.regions.dat$deaths$ObsDay == 118] <- -1

# No adjustments to serology as there is alignment of deaths and seroprevalence age groups.

#......................
# get rho
#......................
ESPcontact <- get_contact_mat(country = "Spain",
                              strict = FALSE,
                              nboots_extrapolate = 5,
                              surveyDOI = "https://doi.org/10.5281/zenodo.1043437",
                              agebands = c(seq(from = 0, to = 80, by = 10), 999))
# assume mixing matrix for 90+ (missing data) is same as 80-90
ESPcontact <- rbind.data.frame(ESPcontact, ESPcontact[nrow(ESPcontact), ])
ESPcontact <- cbind.data.frame(ESPcontact, ESPcontact[, ncol(ESPcontact)])
colnames(ESPcontact)[(length(ESPcontact)-1):length(ESPcontact)] <- c("[80, 90)", "90+")

# multiple through demog for age -- these essentially are age standardized "contact counts"
ESPrho.age <- ESP.agebands.dat$prop_pop$popN %*% as.matrix(ESPcontact)
# standardize for model stability
ESP.agebands.dat$rho <- ESPrho.age/sd(ESPrho.age)

# multiple through demog and age-standardize for region
ESPrho.region <- get_rgnal_contacts(rgndemog = ESP.regions.dat$prop_pop,
                                    contactmat = as.matrix(ESPcontact))
# standardize for model stability
ESP.regions.dat$rho <- ESPrho.region/sd(ESPrho.region)


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
NLD.regions.dat <- process_data2(deaths = deathsdf,
                                 population = populationdf,
                                 sero_val = sero_valdf,
                                 seroprev = sero_prevdf,
                                 cumulative = TRUE,
                                 recast_deaths_df = ECDCdf,
                                 groupingvar = "region",
                                 study_ids = "NLD1",
                                 recast_deaths_geocode = "NLD",
                                 filtRegions = NULL, # some regions combined in serosurvey
                                 filtGender = NULL,
                                 filtAgeBand = NULL)

#......................
# ages
#......................
NLD.agebands.dat <- process_data2(deaths = deathsdf,
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  cumulative = TRUE,
                                  recast_deaths_df = ECDCdf,
                                  groupingvar = "ageband",
                                  study_ids = "NLD1",
                                  recast_deaths_geocode = "NLD",
                                  filtRegions = NULL, # some regions combined in serosurvey
                                  filtGender = NULL,
                                  filtAgeBand = NULL,
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
agebands <- unique(NLD.agebands.dat$deaths$ageband)
nld_adj_seroprev <- tibble::tibble(
  ObsDaymin = NLD.agebands.dat$seroprev$ObsDaymin,
  ObsDaymax = NLD.agebands.dat$seroprev$ObsDaymax,
  ageband = agebands,
  age_low = as.numeric(stringr::str_split_fixed(agebands, "-[0-9]+", n=2)[,1]),
  age_high= as.numeric(stringr::str_split_fixed(agebands, "[0-9]+-", n=2)[,2]),
  seroprevalence_adj = NLD.agebands.dat$seroprev$seroprev) %>%
  dplyr::arrange(age_low)

nld_org_seroprev <- NLD.agebands.dat$seroprev_group %>%
  dplyr::select(c("age_low", "age_high", "seroprevalence"))

nld_adj_seroprev$seroprevalence <- apply(nld_adj_seroprev, 1, wiggle_age_matchfun, wiggle = 2, y = nld_org_seroprev)
# write new serology df
NLD.agebands.dat$seroprev_group_adj <- nld_adj_seroprev

# Netherlands seroprevalence missing in some regions
# assume that this missing values can be imputed as the mean of the other regions
NLD.regions.dat$seroprev_group_adj <- NLD.regions.dat$seroprev_group
NLD.regions.dat$seroprev_group_adj$seroprevalence_adj <- NLD.regions.dat$seroprev_group_adj$seroprevalence
NLD.regions.dat$seroprev_group_adj$seroprevalence_adj[is.na(NLD.regions.dat$seroprev_group_adj$seroprevalence_adj)] <- mean(NLD.regions.dat$seroprev_group_adj$seroprevalence_adj, na.rm = T)


#......................
# get rho
#......................
nld_contact_agebands <- unique(c(0, NLD.agebands.dat$deaths_group$age_high))
# don' have ends, so shorten
nld_contact_agebands <- nld_contact_agebands[1:(length(nld_contact_agebands)-4)]
NLDcontact <- get_contact_mat(country = "Netherlands",
                              strict = TRUE,
                              surveyDOI = "https://doi.org/10.5281/zenodo.1043437",
                              agebands = nld_contact_agebands)$matrix
# assume mixing matrix for 89-94, and 94+ (missing data) is same as 84-94
NLDcontact <- rbind.data.frame(NLDcontact, NLDcontact[nrow(NLDcontact), ], NLDcontact[nrow(NLDcontact), ], NLDcontact[nrow(NLDcontact), ], NLDcontact[nrow(NLDcontact), ])
NLDcontact <- cbind.data.frame(NLDcontact, NLDcontact[, ncol(NLDcontact)], NLDcontact[, ncol(NLDcontact)], NLDcontact[, ncol(NLDcontact)], NLDcontact[, ncol(NLDcontact)])
colnames(NLDcontact)[(length(NLDcontact)-3):length(NLDcontact)] <- c("[79, 84)", "[84, 89)", "[89, 94)", "94+")

# multiple through demog for age -- these essentially are age standardized "contact counts"
NLDrho.age <- NLD.agebands.dat$prop_pop$popN %*% as.matrix(NLDcontact)
# standardize for model stability
NLD.agebands.dat$rho <- NLDrho.age/sd(NLDrho.age)

#TODO FIX NLD MISSING ISSUE
# multiple through demog and age-standardize for region
NLDrho.region <- get_rgnal_contacts(rgndemog = NLD.regions.dat$prop_pop,
                                    contactmat = mean(unlist(as.matrix(NLDcontact))))
# standardize for model stability
NLD.regions.dat$rho <- NLDrho.region/sd(NLDrho.region)


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
DNK.regions.dat <- process_data2(deaths = deathsdf,
                                 population = populationdf,
                                 sero_val = sero_valdf,
                                 seroprev = sero_prevdf,
                                 cumulative = TRUE,
                                 recast_deaths_df = ECDCdf,
                                 groupingvar = "region",
                                 study_ids = "DNK1",
                                 recast_deaths_geocode = "DNK",
                                 filtRegions = NULL, # some regions combined in serosurvey
                                 filtGender = NULL,
                                 filtAgeBand = NULL)

#......................
# ages
#......................
DNK.agebands.dat <- process_data2(deaths = deathsdf,
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  cumulative = TRUE,
                                  recast_deaths_df = ECDCdf,
                                  groupingvar = "ageband",
                                  study_ids = "DNK1",
                                  recast_deaths_geocode = "DNK",
                                  filtRegions = NULL, # some regions combined in serosurvey
                                  filtGender = NULL,
                                  filtAgeBand = NULL)


#......................
# MANUAL ADJUSTMENTS
#......................
# Denmark deaths agebands do not overlap well -- potentially have one age group in 59-69
#  but instead assume ALL national average
agebands <- unique(DNK.agebands.dat$deaths$ageband)
dnk_adj_seroprev <- tibble::tibble(
  ObsDaymin = DNK.agebands.dat$seroprev$ObsDaymin,
  ObsDaymax = DNK.agebands.dat$seroprev$ObsDaymax,
  ageband = agebands,
  age_low = as.numeric(stringr::str_split_fixed(agebands, "-[0-9]+", n=2)[,1]),
  age_high= as.numeric(stringr::str_split_fixed(agebands, "[0-9]+-", n=2)[,2]),
  seroprevalence_adj = DNK.agebands.dat$seroprev$seroprev) %>%
  dplyr::arrange(age_low)
DNK.agebands.dat$seroprev_group_adj <- dnk_adj_seroprev

#......................
# get rho
#......................
dnk_contact_agebands <- unique(c(0, DNK.agebands.dat$deaths_group$age_high))
dnk_contact_agebands <- dnk_contact_agebands[1:(length(dnk_contact_agebands)-2)]
DNKcontact <- get_contact_mat(country = "Denmark",
                              strict = FALSE,
                              nboots_extrapolate = 5,
                              surveyDOI = "https://doi.org/10.5281/zenodo.1043437",
                              agebands = dnk_contact_agebands)

# assume mixing matrix for 79-89, and 89+ (missing data) is same as 79+
DNKcontact <- rbind.data.frame(DNKcontact, DNKcontact[nrow(DNKcontact), ])
DNKcontact <- cbind.data.frame(DNKcontact, DNKcontact[, ncol(DNKcontact)])
colnames(DNKcontact)[(length(DNKcontact)-3):length(DNKcontact)] <- c("[79, 89)", "89+")

# multiple through demog for age -- these essentially are age standardized "contact counts"
DNKrho.age <- DNK.agebands.dat$prop_pop$popN %*% as.matrix(DNKcontact)
# standardize for model stability
DNK.agebands.dat$rho <- DNKrho.age/sd(DNKrho.age)

# multiple through demog and age-standardize for region
DNKrho.region <- get_rgnal_contacts(rgndemog = DNK.regions.dat$prop_pop,
                                    contactmat = as.matrix(DNKcontact))
# standardize for model stability
DNK.regions.dat$rho <- DNKrho.region/sd(DNKrho.region)


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
CHE.agebands.dat <- process_data2(deaths = deathsdf,
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  cumulative = TRUE,
                                  recast_deaths_df = CHE1TimeSeries,
                                  groupingvar = "ageband",
                                  study_ids = "CHE1",
                                  recast_deaths_geocode = "Geneva",   ## use study id in case we get more studies later.
                                  filtRegions = NULL, # some regions combined in serosurvey
                                  filtGender = NULL,
                                  filtAgeBand = NULL)

#......................
# MANUAL ADJUSTMENTS
#......................
# Assume:
# (1) seroprevalence 5-19 is representative of deaths for 0-10 and 10-20, respectively
# (1) seroprevalence 19-49 is representative of deaths for 20-30, 30-40, and 40-50, respectively
# (1) seroprevalence 49+ is representative of deaths for 49+ age bands
agebands <- unique(CHE.agebands.dat$deaths$ageband)
che_adj_seroprev <- tibble::tibble(
  ObsDaymin = CHE.agebands.dat$seroprev$ObsDaymin,
  ObsDaymax = CHE.agebands.dat$seroprev$ObsDaymax,
  ageband = agebands,
  age_low = as.numeric(stringr::str_split_fixed(agebands, "-[0-9]+", n=2)[,1]),
  age_high= as.numeric(stringr::str_split_fixed(agebands, "[0-9]+-", n=2)[,2]),
  seroprevalence_adj = CHE.agebands.dat$seroprev$seroprev) %>%
  dplyr::arrange(age_low)

che_org_seroprev <- CHE.agebands.dat$seroprev_group %>%
  dplyr::select(c("age_low", "age_high", "seroprevalence"))

che_adj_seroprev$seroprevalence <- apply(che_adj_seroprev, 1, wiggle_age_matchfun, wiggle = 2,
                                         y = che_org_seroprev)
# write new serology df
CHE.agebands.dat$seroprev_group_adj <- che_adj_seroprev

#......................
# get rho
#......................
CHEcontact <- get_contact_mat(country = "Switzerland",
                              strict = FALSE,
                              nboots_extrapolate = 5,
                              surveyDOI = "https://doi.org/10.5281/zenodo.1043437",
                              agebands = c(seq(from = 0, to = 80, by = 10), 999))

# multiple through demog for age -- these essentially are age standardized "contact counts"
CHErho.age <- CHE.agebands.dat$prop_pop$popN %*% as.matrix(CHEcontact)
# standardize for model stability
CHE.agebands.dat$rho <- CHErho.age/sd(CHErho.age)

#......................
# save out
#......................
dir.create("data/derived/CHE", recursive = T)
saveRDS(CHE.agebands.dat, "data/derived/CHE/CHE_agebands.RDS")


#............................................................
# Iran
#...........................................................
## Dealt with specially within data processing function to study region of interest
## Do not have info for more than region.
IRN.agebands.dat<-process_data2(deaths = deathsdf,
                                population = populationdf,
                                sero_val = sero_valdf,
                                seroprev = sero_prevdf,
                                cumulative = TRUE,
                                recast_deaths_df = ECDCdf,
                                groupingvar = "ageband",
                                study_ids = "IRN1",
                                recast_deaths_geocode = "IRN",
                                filtRegions = NULL, # some regions combined in serosurvey
                                filtGender = NULL,
                                filtAgeBand = NULL)

#......................
# MANUAL ADJUSTMENTS
#......................
# none, taken care of in script


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
SWE.regions.dat <- process_data2(deaths = deathsdf,
                                 population = populationdf,
                                 sero_val = sero_valdf,
                                 seroprev = sero_prevdf,
                                 cumulative = TRUE,
                                 recast_deaths_df = ECDCdf,
                                 groupingvar = "region",
                                 study_ids = "SWE1",
                                 recast_deaths_geocode = "SWE",
                                 filtRegions = NULL, # some regions combined in serosurvey
                                 filtGender = NULL,
                                 filtAgeBand = NULL)
#......................
# age bands
#......................
### For age analysis, assume data from the 9 regions, about 70% of the regions, is representative.
SWE.agebands.dat <- process_data2(deaths = deathsdf,
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  cumulative = TRUE,
                                  recast_deaths_df = ECDCdf,
                                  groupingvar = "ageband",
                                  study_ids = "SWE1",
                                  recast_deaths_geocode = "SWE",
                                  filtRegions = NULL, # some regions combined in serosurvey
                                  filtGender = NULL,
                                  filtAgeBand = NULL)

#......................
# MANUAL ADJUSTMENTS
#......................
# none
# TODO decide if we should try and infer sens and spec from test listing


#......................
# get rho
#......................
swe_contact_agebands <- unique(c(0, SWE.agebands.dat$deaths_group$age_high))
swe_contact_agebands <- swe_contact_agebands[1:(length(swe_contact_agebands)-2)]
SWEcontact <- get_contact_mat(country = "Sweden",
                              strict = FALSE,
                              nboots_extrapolate = 5,
                              surveyDOI = "https://doi.org/10.5281/zenodo.1043437",
                              agebands = swe_contact_agebands)

# assume mixing matrix for 79-89, and 89+ (missing data) is same as 79+
SWEcontact <- rbind.data.frame(SWEcontact, SWEcontact[nrow(SWEcontact), ])
SWEcontact <- cbind.data.frame(SWEcontact, SWEcontact[, ncol(SWEcontact)])
colnames(SWEcontact)[(length(SWEcontact)-1):length(SWEcontact)] <- c("[79, 89)", "89+")

# multiple through demog for age -- these essentially are age standardized "contact counts"
SWErho.age <- SWE.agebands.dat$prop_pop$popN %*% as.matrix(SWEcontact)
# standardize for model stability
SWE.agebands.dat$rho <- SWErho.age/sd(SWErho.age)

# multiple through demog and age-standardize for region
SWErho.region <- get_rgnal_contacts(rgndemog = SWE.regions.dat$prop_pop,
                                    contactmat = mean(unlist(as.matrix(SWEcontact))))
# standardize for model stability
SWE.regions.dat$rho <- SWErho.region/sd(SWErho.region)


#......................
# save out
#......................
dir.create("data/derived/SWE", recursive = T)
saveRDS(SWE.regions.dat, "data/derived/SWE/SWE_regions.RDS")
saveRDS(SWE.agebands.dat, "data/derived/SWE/SWE_agebands.RDS")






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
  dplyr::select(c("date", "georegion", "deaths")) %>%
  dplyr::mutate(date = lubridate::mdy(date)) # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format


#..................................................................................
#---- Process USA Data #-----
#..................................................................................
#............................................................
# New York City, NY
#...........................................................
#......................
# regions
#......................
# TODO buroughs? -- state level ... ?


#......................
# agebands
#......................
NYC_NY_1.agebands.dat <- process_data2(deaths = deathsdf,
                                       population = populationdf,
                                       sero_val = sero_valdf,
                                       seroprev = sero_prevdf,
                                       cumulative = TRUE,
                                       recast_deaths_df = JHUdf,
                                       groupingvar = "ageband",
                                       study_ids = "NYC_NY_1",
                                       recast_deaths_geocode = "New York_New-York",
                                       filtRegions = NULL, # some regions combined in serosurvey
                                       filtGender = NULL,
                                       filtAgeBand = NULL)
#......................
# MANUAL ADJUSTMENTS
#......................
# TODO -- we have blood age donors sero but death in 17 yr age bands

#......................
# get rho
#......................
# TODO -- assume same as UK?

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



LA_CA.regions.dat <- process_data2(deaths = LACAdeathsdf,
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
# none because basic

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



SC_CA.regions.dat <- process_data2(deaths = SCCAdeathsdf,
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
# none because basic

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



CH_MA.regions.dat <- process_data2(deaths = CHMAdeathsdf,
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
# none because basic

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



MD_FL.regions.dat <- process_data2(deaths = MDFLdeathsdf,
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
# none because basic

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
rgn_key <- readr::read_csv("data/raw/BRA_state_region_key.csv") %>%
  dplyr::select(-c("state")) %>%
  dplyr::rename(state = state_abbr)
bradeaths <- readRDS("data/raw/Brazil_state_age_sex_deaths.rds") %>%
  dplyr::left_join(., rgn_key, by = "state") %>%
  dplyr::select(-c("state")) %>%
  dplyr::filter(!is.na(date)) %>%
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
BRA.regions.dat <- process_data2(deaths = bradeaths,
                                 population = bra_populationdf,
                                 sero_val = sero_valdf,
                                 seroprev = rgnsero_prevdf,
                                 cumulative = FALSE,
                                 groupingvar = "region",
                                 study_ids = "BRA1",
                                 filtRegions = NULL,
                                 filtGender = NULL,
                                 filtAgeBand = NULL)


#......................
# ages
#......................
BRA.agebands.dat <- process_data2(deaths = bradeaths,
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
BRA.regions.dat$rho <- rep(1, length(unique(BRA.regions.dat$deaths$region)))
BRA.agebands.dat$rho <- rep(1, length(unique(BRA.agebands.dat$deaths$ageband)))
# TODO talk to Brazil team if we use contact mats

#......................
# save out
#......................
dir.create("data/derived/BRA/", recursive = T)
saveRDS(BRA.regions.dat, "data/derived/BRA/BRA_regions.RDS")
saveRDS(BRA.agebands.dat, "data/derived/BRA/BRA_agebands.RDS")


















#
