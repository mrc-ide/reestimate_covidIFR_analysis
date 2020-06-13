##....................................................................
## Purpose: Process various datasets for analysis
##
## Notes: Countries are delimited with R-comment sectioning
##....................................................................
library(tidyverse)
source("R/process_data2.R")

#............................................................
# Spain ####
#............................................................
#......................
# from dropbox, focus on regions
#......................
ESP.regions.dat <- process_data2(deaths = "https://www.dropbox.com/s/r4rct1rsp8e5rne/deaths.csv?dl=1",
                                 population = "https://www.dropbox.com/s/hv4woy1zdlveg72/population.csv?dl=1",
                                 sero_val = "https://www.dropbox.com/s/nu7ek2t2bwo9gxo/seroassay_validation.csv?dl=1",
                                 seroprev = "https://www.dropbox.com/s/kc3mle86os2e6g1/seroprevalence.csv?dl=1",
                                 cumulative = TRUE,
                                 USAdata = FALSE,
                                 ECDC = "https://www.dropbox.com/s/a2ds6orlpl5ashs/daily_deaths_ECDC20200518.csv?dl=1",
                                 groupingvar = "region",
                                 study_ids = "ESP1",
                                 geocode = "ESP",
                                 filtRegions = NULL, # limit to mainland Spain
                                 # filtRegions = c("Andalucia", "Aragon",    "Asturias",  "Baleares",  "C Valenciana", "Canarias",
                                 #                 "Cantabria", "Castilla La Mancha", "Castilla y Leon", "Cataluna", "Extremadura",
                                 #                 "Galicia",   "La Rioja",  "Madrid", "Murcia", "Navarra", "Pais Vasco"), # limit to mainland Spain
                                 filtGender = NULL,
                                 filtAgeBand = NULL)

#......................
# from dropbox, focus on age bands
#......................
ESP.agebands.dat <- process_data2(deaths = "https://www.dropbox.com/s/r4rct1rsp8e5rne/deaths.csv?dl=1",
                                  population = "https://www.dropbox.com/s/hv4woy1zdlveg72/population.csv?dl=1",
                                  sero_val = "https://www.dropbox.com/s/nu7ek2t2bwo9gxo/seroassay_validation.csv?dl=1",
                                  seroprev = "https://www.dropbox.com/s/kc3mle86os2e6g1/seroprevalence.csv?dl=1",
                                  cumulative = TRUE,
                                  USAdata = FALSE,
                                  ECDC = "https://www.dropbox.com/s/a2ds6orlpl5ashs/daily_deaths_ECDC20200518.csv?dl=1",
                                  groupingvar = "ageband",
                                  study_ids = "ESP1",
                                  geocode = "ESP",
                                  filtRegions = NULL,
                                  filtGender = NULL,
                                  filtAgeBand = NULL,
                                  death_agebreaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 999),
                                  sero_agebreaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 999))
#......................
# manual adjustment for
# suspicious 0 followed by very large increase
# not sure why April 27 was an issue...a Monday? But looks like all deaths from Apr 27 got pushed to Apr 28
#......................
ESP.regions.dat$deaths$Deaths[ESP.regions.dat$deaths$ObsDay == 117] <- -1
ESP.regions.dat$deaths$Deaths[ESP.regions.dat$deaths$ObsDay == 118] <- -1

ESP.agebands.dat$deaths$Deaths[ESP.agebands.dat$deaths$ObsDay == 117] <- -1
ESP.agebands.dat$deaths$Deaths[ESP.agebands.dat$deaths$ObsDay == 118] <- -1
#......................
# save out
#......................
dir.create("data/derived/ESP", recursive = T)
saveRDS(ESP.regions.dat, "data/derived/ESP/ESP_regions.RDS")
saveRDS(ESP.agebands.dat, "data/derived/ESP/ESP_agebands.RDS")


#............................................................
# Denmark ####
#............................................................
DNK.regions.dat <- process_data2(deaths = "https://www.dropbox.com/s/r4rct1rsp8e5rne/deaths.csv?dl=1",
                                 population = "https://www.dropbox.com/s/hv4woy1zdlveg72/population.csv?dl=1",
                                 sero_val = "https://www.dropbox.com/s/nu7ek2t2bwo9gxo/seroassay_validation.csv?dl=1",
                                 seroprev = "https://www.dropbox.com/s/kc3mle86os2e6g1/seroprevalence.csv?dl=1",
                                 cumulative = TRUE,
                                 USAdata = FALSE,
                                 ECDC = "https://www.dropbox.com/s/a2ds6orlpl5ashs/daily_deaths_ECDC20200518.csv?dl=1",
                                 groupingvar = "region",
                                 study_ids = "DNK1",
                                 geocode = "DNK",
                                 filtRegions = NULL, # some regions combined in serosurvey
                                 filtGender = NULL,
                                 filtAgeBand = NULL)

DNK.agebands.dat <- process_data2(deaths = "https://www.dropbox.com/s/r4rct1rsp8e5rne/deaths.csv?dl=1",
                                  population = "https://www.dropbox.com/s/hv4woy1zdlveg72/population.csv?dl=1",
                                  sero_val = "https://www.dropbox.com/s/nu7ek2t2bwo9gxo/seroassay_validation.csv?dl=1",
                                  seroprev = "https://www.dropbox.com/s/kc3mle86os2e6g1/seroprevalence.csv?dl=1",
                                  cumulative = TRUE,
                                  USAdata = FALSE,
                                  ECDC = "https://www.dropbox.com/s/a2ds6orlpl5ashs/daily_deaths_ECDC20200518.csv?dl=1",
                                  groupingvar = "ageband",
                                  study_ids = "DNK1",
                                  geocode = "DNK",
                                  filtRegions = NULL, # some regions combined in serosurvey
                                  filtGender = NULL,
                                  filtAgeBand = NULL)
#......................
# manual adjustment for
# suspicious negative deaths
#......................
DNK.agebands.dat$deaths$Deaths[DNK.agebands.dat$deaths$ObsDay == 133] <- -1
DNK.agebands.dat$deaths$Deaths[DNK.agebands.dat$deaths$ObsDay == 133] <- -1

#......................
# save out
#......................
dir.create("data/derived/DNK", recursive = T)
saveRDS(DNK.regions.dat, "data/derived/DNK/DNK_regions.RDS")
saveRDS(DNK.agebands.dat, "data/derived/DNK/DNK_agebands.RDS")


#............................................................
# Netherlands ####
#............................................................
NLD.regions.dat <- process_data2(deaths = "https://www.dropbox.com/s/r4rct1rsp8e5rne/deaths.csv?dl=1",
                                 population = "https://www.dropbox.com/s/hv4woy1zdlveg72/population.csv?dl=1",
                                 sero_val = "https://www.dropbox.com/s/nu7ek2t2bwo9gxo/seroassay_validation.csv?dl=1",
                                 seroprev = "https://www.dropbox.com/s/kc3mle86os2e6g1/seroprevalence.csv?dl=1",
                                 cumulative = TRUE,
                                 USAdata = FALSE,
                                 ECDC = "https://www.dropbox.com/s/a2ds6orlpl5ashs/daily_deaths_ECDC20200518.csv?dl=1",
                                 groupingvar = "region",
                                 study_ids = "NLD1",
                                 geocode = "NLD",
                                 filtRegions = NULL, # some regions combined in serosurvey
                                 filtGender = NULL,
                                 filtAgeBand = NULL)

NLD.agebands.dat <- process_data2(deaths = "https://www.dropbox.com/s/r4rct1rsp8e5rne/deaths.csv?dl=1",
                                  population = "https://www.dropbox.com/s/hv4woy1zdlveg72/population.csv?dl=1",
                                  sero_val = "https://www.dropbox.com/s/nu7ek2t2bwo9gxo/seroassay_validation.csv?dl=1",
                                  seroprev = "https://www.dropbox.com/s/kc3mle86os2e6g1/seroprevalence.csv?dl=1",
                                  cumulative = TRUE,
                                  USAdata = FALSE,
                                  ECDC = "https://www.dropbox.com/s/a2ds6orlpl5ashs/daily_deaths_ECDC20200518.csv?dl=1",
                                  groupingvar = "ageband",
                                  study_ids = "NLD1",
                                  geocode = "NLD",
                                  filtRegions = NULL, # some regions combined in serosurvey
                                  filtGender = NULL,
                                  filtAgeBand = NULL,
                                  death_agebreaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 999))
#......................
# save out
#......................
dir.create("data/derived/NLD", recursive = T)
saveRDS(NLD.regions.dat, "data/derived/NLD/NLD_regions.RDS")
saveRDS(NLD.agebands.dat, "data/derived/NLD/NLD_agebands.RDS")




#........................................................
# New York #####
#........................................................
NYC.regions.dat <- process_data2(deaths = "https://www.dropbox.com/s/r4rct1rsp8e5rne/deaths.csv?dl=1",
                                 population = "https://www.dropbox.com/s/hv4woy1zdlveg72/population.csv?dl=1",
                                 sero_val = "https://www.dropbox.com/s/nu7ek2t2bwo9gxo/seroassay_validation.csv?dl=1",
                                 seroprev = "https://www.dropbox.com/s/7jsqb4l7hd5i4i2/seroprevalence_USA.xlsx?dl=1",
                                 cumulative = TRUE,
                                 USAdata = TRUE,
                                 JHU = "https://www.dropbox.com/s/a2ds6orlpl5ashs/daily_deaths_ECDC20200518.csv?dl=1",
                                 groupingvar = "region",
                                 study_ids = "NLD1",
                                 ecdc_countrycode = "NLD",
                                 filtRegions = NULL, # some regions combined in serosurvey
                                 filtGender = NULL,
                                 filtAgeBand = NULL)

NYC.agebands.dat <- process_data2(deaths = "https://www.dropbox.com/s/r4rct1rsp8e5rne/deaths.csv?dl=1",
                                  population = "https://www.dropbox.com/s/hv4woy1zdlveg72/population.csv?dl=1",
                                  sero_val = "https://www.dropbox.com/s/nu7ek2t2bwo9gxo/seroassay_validation.csv?dl=1",
                                  seroprev = "https://www.dropbox.com/s/kc3mle86os2e6g1/seroprevalence.csv?dl=1",
                                  cumulative = TRUE,
                                  USAdata = TRUE,
                                  JHU = "https://www.dropbox.com/s/a2ds6orlpl5ashs/daily_deaths_ECDC20200518.csv?dl=1",
                                  groupingvar = "ageband",
                                  study_ids = "NLD1",
                                  ecdc_countrycode = "NLD",
                                  filtRegions = NULL, # some regions combined in serosurvey
                                  filtGender = NULL,
                                  filtAgeBand = NULL,
                                  death_agebreaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 999))




#............................................................
# Customization ####
#............................................................
# Switzerland, Iran, Sweden need customised function as they are not national surveys.
#........................................................
# Switzerland #####
#........................................................


# Colorado, Massachusetts, and California are single serological estimates. Fit one curve?
#........................................................
# Massachusetts #####
#........................................................






#........................................................
##### parking lot #####
#........................................................

if(length(unique(ESP.agebands.dat$seroprev_group$ObsDaymax))>1 |
   length(unique(ESP.agebands.dat$seroprev_group$ObsDaymin))>1) print("more than one unique serosurvey time window")
sero_day<-round(mean(c(ESP.agebands.dat$seroprev_group$ObsDaymax[1], ESP.agebands.dat$seroprev_group$ObsDaymin[1])))
ESP.agebands.dat$deaths %>%
  filter(ObsDay<=sero_day) %>%
  group_by(ageband) %>%
  summarise(sum(Deaths))

