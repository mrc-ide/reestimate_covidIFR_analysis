####################################################################################
## Purpose: Process spanish datasets for analysis
##
## Author: Nick Brazeau
##
## Date: 26 May, 2020
####################################################################################
library(tidyverse)
source("R/process_data.R")
#............................................................
# from dropbox, focus on regions
#...........................................................
ESP.regions.dat <- process_data(deaths = "https://www.dropbox.com/s/z42zfiyr9tr86qe/deaths_v2.csv?dl=1",
                                population = "https://www.dropbox.com/s/hv4woy1zdlveg72/population.csv?dl=1",
                                sero_val = "https://www.dropbox.com/s/nu7ek2t2bwo9gxo/seroassay_validation.csv?dl=1",
                                seroprev = "https://www.dropbox.com/s/kc3mle86os2e6g1/seroprevalence.csv?dl=1",
                                cumulative = TRUE,
                                ECDC = "https://www.dropbox.com/s/a2ds6orlpl5ashs/daily_deaths_ECDC20200518.csv?dl=1",
                                groupingvar = "region",
                                study_ids = "ESP1",
                                ecdc_countrycode = "ESP",
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
ESP.agebands.dat <- process_data(deaths = "https://www.dropbox.com/s/z42zfiyr9tr86qe/deaths_v2.csv?dl=1",
                                 population = "https://www.dropbox.com/s/hv4woy1zdlveg72/population.csv?dl=1",
                                 sero_val = "https://www.dropbox.com/s/nu7ek2t2bwo9gxo/seroassay_validation.csv?dl=1",
                                 seroprev = "https://www.dropbox.com/s/kc3mle86os2e6g1/seroprevalence.csv?dl=1",
                                 cumulative = TRUE,
                                 ECDC = "https://www.dropbox.com/s/a2ds6orlpl5ashs/daily_deaths_ECDC20200518.csv?dl=1",
                                 groupingvar = "ageband",
                                 study_ids = "ESP1",
                                 ecdc_countrycode = "ESP",
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


