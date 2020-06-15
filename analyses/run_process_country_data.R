####################################################################################
## Purpose: Process Country datasets for analysis
##
## Author: Nick Brazeau
##
## Date: 26 May, 2020
####################################################################################
library(tidyverse)
source("R/process_data2.R")

# deathsFile<-"https://www.dropbox.com/s/z42zfiyr9tr86qe/deaths_v2.csv?dl=1"
# populationFile<- "https://www.dropbox.com/s/hv4woy1zdlveg72/population.csv?dl=1"
# sero_valFile<-"https://www.dropbox.com/s/nu7ek2t2bwo9gxo/seroassay_validation.csv?dl=1"
# seroprevFile<-"https://www.dropbox.com/s/kc3mle86os2e6g1/seroprevalence.csv?dl=1"
# ECDCFile<-"https://www.dropbox.com/s/a2ds6orlpl5ashs/daily_deaths_ECDC20200518.csv?dl=1"

deathsFile<-"data/deaths.csv"
populationFile<- "data/population.csv"
sero_valFile<-"data/seroassay_validation.csv"
seroprevFile<-"data/seroprevalence.csv"
ECDCFile<-"data/daily_deaths_ECDC20200518.csv"

#............................................................
# from dropbox, focus on regions
#...........................................................
ESP.regions.dat <- process_data2(deaths = deathsFile,
                                population = populationFile,
                                sero_val = sero_valFile,
                                seroprev = seroprevFile,
                                cumulative = TRUE,
                                USAdata = FALSE,
                                ECDC = ECDCFile,
                                groupingvar = "region",
                                study_ids = "ESP1",
                                geocode = "ESP",
                                filtRegions = c("Andalucia", "Aragon",    "Asturias",  "Baleares",  "C Valenciana", "Canarias",
                                                "Cantabria", "Castilla La Mancha", "Castilla y Leon", "Cataluna", "Extremadura",
                                                "Galicia",   "La Rioja",  "Madrid", "Murcia", "Navarra", "Pais Vasco"), # limit to mainland Spain
                                filtGender = NULL,
                                filtAgeBand = NULL)
#......................
# manual adjustment for
# suspicious 0
#......................
ESP.regions.dat$deaths$Deaths[ESP.regions.dat$deaths$ObsDay == 117] <- -1

#......................
# save out
#......................
dir.create("data/derived/ESP", recursive = T)
saveRDS(ESP.regions.dat, "data/derived/ESP/ESP_regions.RDS")


#............................................................
# from dropbox, focus on age bands
#...........................................................
ESP.agebands.dat <- process_data2(deaths = deathsFile,
                                 population = populationFile,
                                 sero_val = sero_valFile,
                                 seroprev = seroprevFile,
                                 cumulative = TRUE,
                                 ECDC = ECDCFile,
                                 groupingvar = "ageband",
                                 study_ids = "ESP1",
                                 geocode = "ESP",
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



#########################
# Denmark

DNK.regions.dat <- process_data2(deaths = deathsFile,
                                 population = populationFile,
                                 sero_val = sero_valFile,
                                 seroprev = seroprevFile,
                                 cumulative = TRUE,
                                 ECDC = ECDCFile,
                                 groupingvar = "region",
                                 study_ids = "DNK1",
                                 geocode = "DNK",
                                 filtRegions = NULL, # some regions combined in serosurvey
                                 filtGender = NULL,
                                 filtAgeBand = NULL)


DNK.agebands.dat <- process_data2(deaths = deathsFile,
                                  population = populationFile,
                                  sero_val = sero_valFile,
                                  seroprev = seroprevFile,
                                  cumulative = TRUE,
                                  ECDC = ECDCFile,
                                  groupingvar = "ageband",
                                  study_ids = "DNK1",
                                  geocode = "DNK",
                                  filtRegions = NULL, # some regions combined in serosurvey
                                  filtGender = NULL,
                                  filtAgeBand = NULL)

dir.create("data/derived/DNK", recursive = T)
saveRDS(DNK.regions.dat, "data/derived/DNK/DNK_regions.RDS")
saveRDS(DNK.agebands.dat, "data/derived/DNK/DNK_agebands.RDS")

#########################
# Netherlands
## NB one or two regions have missing seroprevalence, as map regions did not match to current regions.
NLD.regions.dat <- process_data2(deaths = deathsFile,
                                 population = populationFile,
                                 sero_val = sero_valFile,
                                 seroprev = seroprevFile,
                                 cumulative = TRUE,
                                 ECDC = ECDCFile,
                                 groupingvar = "region",
                                 study_ids = "NLD1",
                                 geocode = "NLD",
                                 filtRegions = NULL, # some regions combined in serosurvey
                                 filtGender = NULL,
                                 filtAgeBand = NULL)

NLD.agebands.dat <- process_data2(deaths = deathsFile,
                                  population = populationFile,
                                  sero_val = sero_valFile,
                                  seroprev = seroprevFile,
                                  cumulative = TRUE,
                                  ECDC = ECDCFile,
                                  groupingvar = "ageband",
                                  study_ids = "NLD1",
                                  geocode = "NLD",
                                  filtRegions = NULL, # some regions combined in serosurvey
                                  filtGender = NULL,
                                  filtAgeBand = NULL)

dir.create("data/derived/NLD", recursive = T)
saveRDS(NLD.regions.dat, "data/derived/NLD/NLD_regions.RDS")
saveRDS(NLD.agebands.dat, "data/derived/NLD/NLD_agebands.RDS")

########
### Iran. (dealt with specially within data processing function to study region of interest)
## Do not have info for more than region.
IRN.agebands.dat<-process_data2(deaths = deathsFile,
                                population = populationFile,
                                sero_val = sero_valFile,
                                seroprev = seroprevFile,
                                cumulative = TRUE,
                                ECDC = ECDCFile,
                                groupingvar = "ageband",
                                study_ids = "IRN1",
                                geocode = "IRN",
                                filtRegions = NULL, # some regions combined in serosurvey
                                filtGender = NULL,
                                filtAgeBand = NULL)

dir.create("data/derived/IRN", recursive = T)
saveRDS(IRN.agebands.dat, "data/derived/IRN/IRN_agebands.RDS")

# Switzerland. CHE1. TODO include change in seroprevalence over 3 weeks. (Currently output average over all of them)
### Has specific Geneva deaths time series, do not need ECDC
## rename vars so will work like ECDC file
deathsTimeSeries<-readr::read_csv("data/deaths_time_series.csv")
deathsTimeSeries<- deathsTimeSeries %>%
  rename(dateRep=date_start_survey,
         deaths=n_deaths)
deathsTimeSeries$countryterritoryCode<-NA
deathsTimeSeries$countryterritoryCode[which(deathsTimeSeries$study_id=="CHE1")]<-"CHE1"
write.csv(deathsTimeSeries,file="data/deathsTimeSeriesR.csv")

CHE.agebands.dat<-process_data2(deaths = deathsFile,
                                population = populationFile,
                                sero_val = sero_valFile,
                                seroprev = seroprevFile,
                                cumulative = TRUE,
                                ECDC = "data/deathsTimeSeriesR.csv",
                                groupingvar = "ageband",
                                study_ids = "CHE1",
                                geocode = "CHE1",   ## use study id in case we get more studies later.
                                filtRegions = NULL, # some regions combined in serosurvey
                                filtGender = NULL,
                                filtAgeBand = NULL)

dir.create("data/derived/CHE", recursive = T)
saveRDS(CHE.agebands.dat, "data/derived/CHE/CHE_agebands.RDS")


### Sweden needs customised function as not national survey (similar to Iran).


#######################
# USA data.
deathsFile<-"data/deaths.csv"
populationFile<- "data/population.csv"
sero_valFile<-"data/seroassay_validation.csv"
seroprevFile<-"data/seroprevalence.csv"
JHUFile<-"data/daily_deaths_ECDC20200518.csv"

jhu<-read.csv("data/JHU.csv")

process_data_usa <- function(deaths = NULL, population = NULL, sero_val = NULL, seroprev = NULL,
                             cumulative = FALSE, USAdata = FALSE, JHU = NULL,
                             groupingvar, study_ids, geocode,
                             filtRegions = NULL, filtGender = NULL, filtAgeBand = NULL, death_agebreaks = NULL, sero_agebreaks = NULL)
