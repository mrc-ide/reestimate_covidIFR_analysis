####################################################################################
## Purpose: Process spanish datasets for analysis
##
## Author: Nick Brazeau
##
## Date: 26 May, 2020
####################################################################################
library(tidyverse)
source("R/process_data2.R")
#............................................................
# from dropbox, focus on regions
#...........................................................
ESP.regions.dat <- process_data2(deaths = "data/deaths.csv",
                                population = "data/population.csv",
                                sero_val = "data/seroassay_validation.csv",
                                seroprev = "data/seroprevalence.csv",
                                cumulative = TRUE,
                                ECDC = "data/daily_deaths_ECDC20200518.csv",
                                groupingvar = "region",
                                study_ids = "ESP1",
                                ecdc_countrycode = "ESP",
                                filtRegions = NULL, # limit to mainland Spain
                                # filtRegions = c("Andalucia", "Aragon",    "Asturias",  "Baleares",  "C Valenciana", "Canarias",
                                #                 "Cantabria", "Castilla La Mancha", "Castilla y Leon", "Cataluna", "Extremadura",
                                #                 "Galicia",   "La Rioja",  "Madrid", "Murcia", "Navarra", "Pais Vasco"), # limit to mainland Spain
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
ESP.agebands.dat <- process_data2(deaths = "data/deaths.csv",
                                 population = "data/population.csv",
                                 sero_val = "data/seroassay_validation.csv",
                                 seroprev = "data/seroprevalence.csv",
                                 cumulative = TRUE,
                                 ECDC = "data/daily_deaths_ECDC20200518.csv",
                                 groupingvar = "ageband",
                                 study_ids = "ESP1",
                                 ecdc_countrycode = "ESP",
                                 filtRegions = NULL,
                                 filtGender = NULL,
                                 filtAgeBand = NULL)
                                                                  # filtAgeBand = c("0-10", "10-20", "20-30",
                                 #                 "30-40", "40-50", "50-60",
                                 #                 "60-70", "70-80", "80-90", "90-999"))
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

DNK.regions.dat <- process_data2(deaths = "data/deaths.csv",
                                population = "data/population.csv",
                                sero_val = "data/seroassay_validation.csv",
                                seroprev = "data/seroprevalence.csv",
                                cumulative = TRUE,
                                ECDC = "data/daily_deaths_ECDC20200518.csv",
                                groupingvar = "region",
                                study_ids = "DNK1",
                                ecdc_countrycode = "DNK",
                                filtRegions = NULL, # some regions combined in serosurvey
                                filtGender = NULL,
                                filtAgeBand = NULL)

DNK.agebands.dat <- process_data2(deaths = "data/deaths.csv",
                                population = "data/population.csv",
                                sero_val = "data/seroassay_validation.csv",
                                seroprev = "data/seroprevalence.csv",
                                cumulative = TRUE,
                                ECDC = "data/daily_deaths_ECDC20200518.csv",
                                groupingvar = "ageband",
                                study_ids = "DNK1",
                                ecdc_countrycode = "DNK",
                                filtRegions = NULL, # some regions combined in serosurvey
                                filtGender = NULL,
                                filtAgeBand = NULL)


#########################
# Netherlands

NLD.regions.dat <- process_data2(deaths = "data/deaths.csv",
                                 population = "data/population.csv",
                                 sero_val = "data/seroassay_validation.csv",
                                 seroprev = "data/seroprevalence.csv",
                                 cumulative = TRUE,
                                 ECDC = "data/daily_deaths_ECDC20200518.csv",
                                 groupingvar = "region",
                                 study_ids = "NLD1",
                                 ecdc_countrycode = "NLD",
                                 filtRegions = NULL, # some regions combined in serosurvey
                                 filtGender = NULL,
                                 filtAgeBand = NULL)

NLD.agebands.dat <- process_data2(deaths = "data/deaths.csv",
                                  population = "data/population.csv",
                                  sero_val = "data/seroassay_validation.csv",
                                  seroprev = "data/seroprevalence.csv",
                                  cumulative = TRUE,
                                  ECDC = "data/daily_deaths_ECDC20200518.csv",
                                  groupingvar = "ageband",
                                  study_ids = "NLD1",
                                  ecdc_countrycode = "NLD",
                                  filtRegions = NULL, # some regions combined in serosurvey
                                  filtGender = NULL,
                                  filtAgeBand = NULL)

########
# Switzerland, Iran, Sweden need customised function as they are not national surveys.


if(length(unique(ESP.agebands.dat$seroprev_group$ObsDaymax))>1 |
   length(unique(ESP.agebands.dat$seroprev_group$ObsDaymin))>1) print("more than one unique serosurvey time window")
sero_day<-round(mean(c(ESP.agebands.dat$seroprev_group$ObsDaymax[1], ESP.agebands.dat$seroprev_group$ObsDaymin[1])))
ESP.agebands.dat$deaths %>%
  filter(ObsDay<=sero_day) %>%
  group_by(ageband) %>%
  summarise(sum(Deaths))

