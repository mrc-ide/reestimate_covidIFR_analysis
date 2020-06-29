####################################################################################
## Purpose: Process Country datasets for analysis
##
## Notes: Process data functions and files are specific to regions
####################################################################################
library(tidyverse)
source("R/process_data2.R")

#..................................................................................
# Preprocess Data for numerous Eurasia studies
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
  dplyr::select(-c("reference"))

# seroprev
sero_valdf <- readr::read_csv("data/raw/seroassay_validation.csv")
sero_prevdf <- readr::read_csv("data/raw/seroprevalence.csv") %>%
  dplyr::select(-c("ref", "notes", "X22", "X23")) # extra note cols



#............................................................
# Spain
#...........................................................
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
# manual adjustment for
# suspicious 0 followed by very large increase
# not sure why April 27 was an issue...a Monday? But looks like all deaths from Apr 27 got pushed to Apr 28
#......................
ESP.regions.dat$deaths$Deaths[ESP.regions.dat$deaths$ObsDay == 117] <- -1
ESP.regions.dat$deaths$Deaths[ESP.regions.dat$deaths$ObsDay == 118] <- -1

#......................
# save out
#......................
dir.create("data/derived/ESP", recursive = T)
saveRDS(ESP.regions.dat, "data/derived/ESP/ESP_regions.RDS")


#............................................................
# from dropbox, focus on age bands
#...........................................................
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
                                                  "60-70", "70-80", "80-90", "90-999"))
#......................
# manual adjustment for
# suspicious 0 followed by very large increase
# not sure why April 27 was an issue...a Monday? But looks like all deaths from Apr 27 got pushed to Apr 28
#......................
ESP.agebands.dat$deaths$Deaths[ESP.agebands.dat$deaths$ObsDay == 117] <- -1
ESP.agebands.dat$deaths$Deaths[ESP.agebands.dat$deaths$ObsDay == 118] <- -1
#......................
# save out
#......................
saveRDS(ESP.agebands.dat, "data/derived/ESP/ESP_agebands.RDS")



#............................................................
# Denmark
#...........................................................
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

dir.create("data/derived/DNK", recursive = T)
saveRDS(DNK.regions.dat, "data/derived/DNK/DNK_regions.RDS")
saveRDS(DNK.agebands.dat, "data/derived/DNK/DNK_agebands.RDS")

#............................................................
# Netherlands
#...........................................................
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
                                  filtAgeBand = NULL)

dir.create("data/derived/NLD", recursive = T)
saveRDS(NLD.regions.dat, "data/derived/NLD/NLD_regions.RDS")
saveRDS(NLD.agebands.dat, "data/derived/NLD/NLD_agebands.RDS")

#............................................................
#Iran
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

dir.create("data/derived/IRN", recursive = T)
saveRDS(IRN.agebands.dat, "data/derived/IRN/IRN_agebands.RDS")

#............................................................
# Switzerland
#...........................................................
# TODO fix here -- with DF format this should be easier -- also this should be able
# to go in as time series without the deathsdf?
## CHE1. TODO include change in seroprevalence over 3 weeks. (Currently output average over all of them)
## Has specific Geneva deaths time series, do not need ECDC
## rename vars so will work like ECDC file
deathsTimeSeries <- readr::read_csv("data/raw/deaths_time_series.csv")
deathsTimeSeries <- deathsTimeSeries %>%
  rename(dateRep=date_start_survey,
         deaths=n_deaths)
deathsTimeSeries$countryterritoryCode <- NA
deathsTimeSeries$countryterritoryCode[which(deathsTimeSeries$study_id=="CHE1")] <- "CHE1"
write.csv(deathsTimeSeries, file="data/raw/deathsTimeSeriesR.csv")

CHE.agebands.dat<-process_data2(deaths = deathsdf,
                                population = populationdf,
                                sero_val = sero_valdf,
                                seroprev = sero_prevdf,
                                cumulative = TRUE,
                                recast_deaths_df = "data/raw/deathsTimeSeriesR.csv",
                                groupingvar = "ageband",
                                study_ids = "CHE1",
                                recast_deaths_geocode = "CHE1",   ## use study id in case we get more studies later.
                                filtRegions = NULL, # some regions combined in serosurvey
                                filtGender = NULL,
                                filtAgeBand = NULL)

dir.create("data/derived/CHE", recursive = T)
saveRDS(CHE.agebands.dat, "data/derived/CHE/CHE_agebands.RDS")


#............................................................
# Sweden
#...........................................................
### For age analysis, assume data from the 9 regions, about 70% of the regions, is representative.
SWE.agebands.dat<-process_data2(deaths = deathsdf,
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

dir.create("data/derived/SWE", recursive = T)
saveRDS(SWE.agebands.dat, "data/derived/SWE/SWE_agebands.RDS")


#..................................................................................
# USA Data
#..................................................................................
source("R/misc_USA.R")
#######################
# USA data.
### For LA_CA, SC_CA,  CH_MA, MD_FL - use USA facts (?). And process_usa_basic_data_timeseries to get overall estimates.

deathsdf <- "data/raw/deaths.csv"
populationdf <- "data/raw/USA_County_Demographic_Data.csv"
sero_valdf <- "data/raw/seroassay_validation.csv"
sero_prevdf <- "data/raw/seroprevalence.csv"
timeSeriesFile <- "data/raw/covid_deaths_usafacts_study_countys.csv"


#............................................................
# Los Angeles, CA
#...........................................................
# LA_CA - process for age and for region.
LA_CA.regions.dat <- process_usa_basic_data_timeseries(population = populationdf,
                                                       sero_val = sero_valdf,
                                                       seroprev = sero_prevdf,
                                                       timeSeriesFile = timeSeriesFile,
                                                       study_ids = "LA_CA",
                                                       state = "California",
                                                       county = "Los Angeles County")



LA_CA.agebands.dat <- process_data_usa_facts(deaths = deathsdf,
                                             population = populationdf,
                                             sero_val = sero_valdf,
                                             seroprev = sero_prevdf,
                                             timeSeriesFile = timeSeriesFile,
                                             cumulative = FALSE,
                                             groupingvar = "ageband",
                                             study_ids = "LA_CA",
                                             state = "California",
                                             county = "Los Angeles County")

dir.create("data/derived/USA", recursive = T)
saveRDS(LA_CA.agebands.dat, "data/derived/USA/LA_CA_agebands.RDS")
saveRDS(LA_CA.regions.dat, "data/derived/USA/LA_CA_regions.RDS")

#............................................................
# Santa Clara, CA
#...........................................................
SC_CA.regions.dat <- process_usa_basic_data_timeseries(population = populationdf,
                                                       sero_val = sero_valdf,
                                                       seroprev = sero_prevdf,
                                                       timeSeriesFile=timeSeriesFile,
                                                       study_ids = "SC_CA",
                                                       state = "California",
                                                       county = "Santa Clara County")

saveRDS(SC_CA.regions.dat, "data/derived/USA/SC_CA_regions.RDS")

#............................................................
# Chelsea, MA
#...........................................................
### NB matching to Suffolk county may not be quite right (Chelsea is a city)
CH_MA.regions.dat <- process_usa_basic_data_timeseries(population = populationdf,
                                                       sero_val = sero_valdf,
                                                       seroprev = sero_prevdf,
                                                       timeSeriesFile = timeSeriesFile,
                                                       study_ids = "CH_MA",
                                                       state = "Massachusetts",
                                                       county = "Suffolk County")

saveRDS(CH_MA.regions.dat, "data/derived/USA/CH_MA_regions.RDS")

#............................................................
# Miama, FL
#...........................................................
MD_FL.regions.dat <- process_usa_basic_data_timeseries(population = populationdf,
                                                       sero_val = sero_valdf,
                                                       seroprev = sero_prevdf,
                                                       timeSeriesFile=timeSeriesFile,
                                                       study_ids = "MD_FL",
                                                       state = "Florida",
                                                       county = "Miami-Dade County")
saveRDS(MD_FL.regions.dat, "data/derived/USA/MD_FL_regions.RDS")


#............................................................
# New York City, NY
#...........................................................
NYC_NY_1.regions.dat<-process_usa_basic_data_timeseries(population = populationdf,
                                                        sero_val = sero_valdf,
                                                        seroprev = sero_prevdf,
                                                        timeSeriesFile=timeSeriesFile,
                                                        study_ids = "NYC_NY_1",
                                                        state = "New York",
                                                        county = c("New York County","Kings County","Bronx County",
                                                                   "Richmond County","Queens County"))
## fix population. TODO - fix properly in code.
NYC_NY_1.regions.dat$popN <- 8399000
saveRDS(NYC_NY_1.regions.dat, "data/derived/USA/NYC_NY_1_regions.RDS")


#..................................................................................
# Latin America Data
#..................................................................................
#............................................................
# Brazil
#...........................................................





















#
