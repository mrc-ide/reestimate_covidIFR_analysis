#................................................................................................
## Purpose: Process Country datasets for analysis
##
## Notes: Process data functions and files are specific to regions
#................................................................................................
library(tidyverse)
source("R/process_data4.R")

#......................
# global data
#......................
# JHU data for new recast df
JHUdf <- readr::read_csv("data/raw/time_series_covid19_deaths_global.csv") %>%
  tidyr::pivot_longer(., cols = -c("Province/State", "Country/Region", "Lat", "Long"),
                      names_to = "date", values_to = "deaths") %>%
  magrittr::set_colnames(c("province", "country_region", "lat", "long", "date", "deaths")) %>%
  dplyr::filter(is.na(province)) %>% # only want by country
  dplyr::mutate(georegion = countrycode::countryname(country_region, destination = "iso3c"),
                date = lubridate::mdy(date)) %>% # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format
  dplyr::group_by(georegion) %>% # group by region for daily deaths
  dplyr::rename(cumdeaths = deaths) %>%
  dplyr::mutate(deaths = cumdeaths - dplyr::lag(cumdeaths),
                deaths = ifelse(is.na(deaths), 0, deaths), # take care of first value
                deaths = ifelse(deaths < 1, 0, deaths)) %>% # take care of cumulative death correction
  dplyr::select(c("date", "georegion", "deaths")) %>%
  dplyr::filter(date <= lubridate::mdy("08-17-2020")) %>%
  dplyr::ungroup(.)

# serovalidation
sero_valdf <-  readr::read_tsv("data/raw/serovalidation_final_raw.tsv")
# seroprevalence
sero_prevdf <- readr::read_tsv("data/raw/seroprevalence_final_raw.tsv") %>%
  dplyr::select(-c("ref", "notes")) %>%
  dplyr::mutate(date_start_survey = lubridate::ymd(date_start_survey), # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format
                date_end_survey = lubridate::ymd(date_end_survey),
                seroprevalence_unadjusted = ifelse(is.na(seroprevalence_unadjusted), n_positive/n_tested, seroprevalence_unadjusted))

# cumulative deaths
deathsdf <- readr::read_tsv("data/raw/cumulative_deaths.tsv") %>%
  dplyr::select(-c("ref", "notes")) %>%
  dplyr::mutate(date_start_survey = lubridate::ymd(date_start_survey), # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format
                date_end_survey = lubridate::ymd(date_end_survey))

# care home deaths
deaths_ch <- readr::read_csv("data/raw/care_home_deaths.csv")

# demography (non-US Census data)
populationdf <- readr::read_tsv("data/raw/population.tsv") %>%
  dplyr::select(-c("reference")) %>%
  dplyr::mutate(age_low = ifelse(age_low == 0 & age_high == 0, 1, age_low),
                age_high = ifelse(age_low == 1 & age_high == 0, 1, age_high))  # liftover "zero" year olds to be 1, 1 as well

#..................................................................................
#---- Preprocess Latin America Data  #-----
#..................................................................................
#............................................................
#---- BRA1 #----
#...........................................................
#......................
# pre-process Brazil
#......................
bra_cumdeaths <- readRDS("data/raw/Brazil_state_age_sex_deaths.rds") %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::rename(gender = sex) %>%
  dplyr::mutate(age = ifelse(age >= 100, 100, age), # liftover some very old ages -- capping at 100
                age = factor(age, levels = c(0:100), labels = c(0:100))) %>% # ugly code to fill in dates
  dplyr::group_by(region, age, gender, .drop = FALSE) %>%
  dplyr::summarise( deaths = sum(count) ) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate(
    age = as.numeric(as.character(age)),
    country = "BRA",
    study_id = "BRA1",
    age_low = age,
    age_high = age + 1, # make 1-based for cuts
    age_breakdown = 1,
    gender_breakdown = 1,
    for_regional_analysis = 1) %>%
  dplyr::ungroup(.) %>%
  dplyr::rename(n_deaths = deaths) %>%
  dplyr::arrange(region)


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
BRA.regions.dat <- process_data4(cum_tp_deaths = bra_cumdeaths,
                                 time_series_totdeaths_df = JHUdf,
                                 time_series_totdeaths_geocode = "BRA",
                                 population = bra_populationdf,
                                 sero_val = sero_valdf,
                                 seroprev = sero_prevdf,
                                 get_descriptive_dat = TRUE,
                                 groupingvar = "region",
                                 study_ids = "BRA1",
                                 origin = lubridate::ymd("2020-01-01"))


#......................
# ages
#......................
BRA.agebands.dat <- process_data4(cum_tp_deaths = bra_cumdeaths,
                                  time_series_totdeaths_df = JHUdf,
                                  time_series_totdeaths_geocode = "BRA",
                                  population = bra_populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  get_descriptive_dat = TRUE,
                                  groupingvar = "ageband",
                                  study_ids = "BRA1",
                                  origin = lubridate::ymd("2020-01-01"),
                                  death_agebreaks = c(0, 4, 9,
                                                      19, 29, 39,
                                                      49, 59, 69,
                                                      79, 999))


#..................................................................................
#............................................................
#---- BRA1 CITIES #----
#...........................................................
#......................
# pre-process Brazil cities data.
#......................
bra_sero<-read.csv("data/raw/bra1_city_sero.csv") %>%
  dplyr::rename(city=City)
bra_city_d<-read.csv("data/raw/bra1_city_deaths.csv")
bra_city<-left_join(bra_city_d,bra_sero,by="city")

#......................
# MANUAL ADJUSTMENTS
#......................
# none needed with linelist

#......................
# save out
#......................
dir.create("data/derived/BRA1/", recursive = T)
saveRDS(BRA.regions.dat, "data/derived/BRA1/BRA1_regions.RDS")
saveRDS(BRA.agebands.dat, "data/derived/BRA1/BRA1_agebands.RDS")


#..................................................................................
#---- Preprocess European Data #----
#..................................................................................
#............................................................
#---- CHE1, Geneva #----
#...........................................................
CHE1TimeSeries <- readr::read_tsv("data/raw/deaths_time_series_subnat.tsv") %>%
  dplyr::filter(study_id == "CHE1")
# sanity check
identical(CHE1TimeSeries$date_start_survey, CHE1TimeSeries$date_end_survey)
CHE1TimeSeries <- CHE1TimeSeries %>%
  dplyr::rename(date = date_end_survey,
                deaths = n_deaths,
                georegion = region) %>%
  dplyr::select(c("date", "georegion", "deaths")) %>%
  dplyr::mutate(date = lubridate::dmy(date)) # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format



#......................
# ages
#......................
CHE1.agebands.dat <- process_data4(cum_tp_deaths = deathsdf,
                                   time_series_totdeaths_df = CHE1TimeSeries,
                                   time_series_totdeaths_geocode = "Geneva", ## use study id in case we get more studies later.
                                   population = populationdf,
                                   sero_val = sero_valdf,
                                   seroprev = sero_prevdf,
                                   get_descriptive_dat = TRUE,
                                   study_ids = "CHE1",
                                   groupingvar = "ageband",
                                   death_agebreaks = c(0, 10, 20, 30,
                                                       40, 50, 60, 70,
                                                       80, 999)) # for pop splits

#......................
# MANUAL ADJUSTMENTS
#......................
# Assume:
# (1) seroprevalence 5-9 is representative of ages 0-10
# (2) seroprevalence 10-19 is representative of ages 10-20
# (3) seroprevalence 20-49 is representative of ages 20-30 and 30-40 and 40-50
# (4) seroprevalence 50-64 is representative of ages 50-60 and 60-70
# (5) seroprevalence 65+ is representative of ages 70-80 and 80+
ageband <- unique(CHE1.agebands.dat$deaths_propMCMC$ageband)
che_adj_seroprev <- tibble::tibble(
  ObsDaymin = unique(CHE1.agebands.dat$seroprevMCMC$ObsDaymin),
  ObsDaymax = unique(CHE1.agebands.dat$seroprevMCMC$ObsDaymax))
che_adj_seroprev <- tidyr::expand_grid(che_adj_seroprev, ageband) %>%
  dplyr::mutate(age_low = as.numeric(stringr::str_split_fixed(ageband, "-[0-9]+", n=2)[,1]),
                age_high= as.numeric(stringr::str_split_fixed(ageband, "[0-9]+-", n=2)[,2])) %>%
  dplyr::arrange(age_low)

# pull out original
che_org_seroprev <- CHE1.agebands.dat$seroprev_group %>%
  dplyr::select(c("ObsDaymin", "ObsDaymax", "age_low", "age_high", "n_positive", "n_tested", "seroprevalence_unadjusted"))

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
  dplyr::select(c("ObsDaymin", "ObsDaymax", "ageband", "n_positive", "n_tested", "SeroPrev")) %>%
  dplyr::arrange(ObsDaymin, ObsDaymax, ageband)

# bring together
CHE1.agebands.dat$seroprevMCMC <- che_adj_seroprev


#......................
# save out
#......................
dir.create("data/derived/CHE1", recursive = T)
saveRDS(CHE1.agebands.dat, "data/derived/CHE1/CHE1_agebands.RDS")



#............................................................
#---- CHE2, Zurich #----
#...........................................................
CHE2TimeSeries <- readr::read_tsv("data/raw/deaths_time_series_subnat.tsv") %>%
  dplyr::filter(study_id == "CHE2")
# sanity check
identical(CHE2TimeSeries$date_start_survey, CHE2TimeSeries$date_end_survey)
CHE2TimeSeries <- CHE2TimeSeries %>%
  dplyr::rename(date = date_end_survey,
                deaths = n_deaths,
                georegion = region) %>%
  dplyr::select(c("date", "georegion", "deaths")) %>%
  dplyr::mutate(date = lubridate::dmy(date)) # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format



#......................
# ages
#......................
CHE2.agebands.dat <- process_data4(cum_tp_deaths = deathsdf,
                                   time_series_totdeaths_df = CHE2TimeSeries,
                                   time_series_totdeaths_geocode = "Zurich",
                                   population = populationdf,
                                   sero_val = sero_valdf,
                                   seroprev = sero_prevdf,
                                   get_descriptive_dat = TRUE,
                                   groupingvar = "ageband",
                                   study_ids = "CHE2",
                                   death_agebreaks = c(0, 10, 20, 30,
                                                       40, 50, 60, 70,
                                                       80, 999)) # for pop splits

#......................
# MANUAL ADJUSTMENTS
#......................
# Assume:
# (1) seroprevalence 18-75 is representative of all ages
ageband <- unique(CHE2.agebands.dat$deaths_propMCMC$ageband)
che2_adj_seroprev <- tibble::tibble(
  ObsDaymin = unique(CHE2.agebands.dat$seroprevMCMC$ObsDaymin),
  ObsDaymax = unique(CHE2.agebands.dat$seroprevMCMC$ObsDaymax),
  n_positive = unique(CHE2.agebands.dat$seroprevMCMC$n_positive),
  n_tested = unique(CHE2.agebands.dat$seroprevMCMC$n_tested),
  SeroPrev = unique(CHE2.agebands.dat$seroprevMCMC$SeroPrev))
che2_adj_seroprev <- tidyr::expand_grid(che2_adj_seroprev, ageband) %>%
  dplyr::mutate(age_low = as.numeric(stringr::str_split_fixed(ageband, "-[0-9]+", n=2)[,1]),
                age_high= as.numeric(stringr::str_split_fixed(ageband, "[0-9]+-", n=2)[,2])) %>%
  dplyr::arrange(age_low)

che2_adj_seroprev <- che2_adj_seroprev %>%
  dplyr::arrange(ObsDaymin, ObsDaymax, ageband) %>%
  dplyr::select(-c("age_low", "age_high"))

# bring together
CHE2.agebands.dat$seroprevMCMC <- che2_adj_seroprev

#......................
# save out
#......................
dir.create("data/derived/CHE2", recursive = T)
#saveRDS(CHE2.region.dat, "data/derived/CHE2/CHE2_region.RDS")
saveRDS(CHE2.agebands.dat, "data/derived/CHE2/CHE2_agebands.RDS")


#............................................................
#---- DNK1 #----
#...........................................................
# #......................
# # regions
# #......................
dnk_sero<-sero_prevdf %>%
  dplyr::filter(study_id=="DNK1") %>%
  dplyr::mutate(for_regional_analysis=ifelse(region=="all",0,1))

DNK.regions.dat <- process_data4(cum_tp_deaths = deathsdf,
                                 time_series_totdeaths_df = JHUdf,
                                 time_series_totdeaths_geocode = "DNK",
                                 population = populationdf,
                                 sero_val = sero_valdf,
                                 seroprev = dnk_sero,
                                 get_descriptive_dat = TRUE,
                                 groupingvar = "region",
                                 study_ids = "DNK1")

#......................
# ages
#......................
DNK.agebands.dat <- process_data4(cum_tp_deaths = deathsdf,
                                  time_series_totdeaths_df = JHUdf,
                                  time_series_totdeaths_geocode = "DNK",
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  get_descriptive_dat = TRUE,
                                  groupingvar = "ageband",
                                  study_ids = "DNK1",
                                  death_agebreaks = c(0, 59, 69, 79, 999),
                                  sero_agebreaks = c(0, 59, 69, 79, 999))
#......................
# ages
#......................
## Use age sero data for age region model
dnk_age_sero<-read.csv("data/raw/DNK1_age.csv") %>%
  dplyr::select(-c("ref", "notes")) %>%
  dplyr::mutate(date_start_survey = lubridate::ymd(date_start_survey), # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format
                date_end_survey = lubridate::ymd(date_end_survey),
                seroprevalence_unadjusted = ifelse(is.na(seroprevalence_unadjusted), n_positive/n_tested, seroprevalence_unadjusted))

DNK.agebands_age_sero.dat <- process_data4(cum_tp_deaths = deathsdf,
                                  time_series_totdeaths_df = JHUdf,
                                  time_series_totdeaths_geocode = "DNK",
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = dnk_age_sero,
                                  get_descriptive_dat = TRUE,
                                  groupingvar = "ageband",
                                  study_ids = "DNK1",
                                  death_agebreaks = c(0, 59, 69, 79, 999),
                                  sero_agebreaks = c(0, 59, 69, 79, 999))


#......................
# MANUAL ADJUSTMENTS
#......................
# Assume DNK blood donors 17-69 years representative of all age groups
dnkagebands <- unique(DNK.agebands.dat$deaths_propMCMC$ageband)
dnk_adj_seroprev <- tibble::tibble(
  ObsDaymin = DNK.agebands.dat$seroprevMCMC$ObsDaymin,
  ObsDaymax = DNK.agebands.dat$seroprevMCMC$ObsDaymax,
  n_positive = NA,
  n_tested = NA)

dnk_adj_seroprev <- lapply(dnkagebands, function(x){
  dnk_adj_seroprev %>%
    dplyr::mutate(ageband = x)}) %>%
  dplyr::bind_rows() %>%
  dplyr::arrange(., ObsDaymin, ObsDaymax, ageband)

DNKnatprev <- DNK.agebands.dat$seroprev_group %>%
  dplyr::select(c("ObsDaymin", "ObsDaymax", "seroprevalence_unadjusted", "range_sero_low", "range_sero_high")) %>%
  dplyr::rename(SeroPrev = seroprevalence_unadjusted,
                SeroLCI = range_sero_low,
                SeroUCI = range_sero_high)

dnk_adj_seroprev <- dplyr::left_join(dnk_adj_seroprev, DNKnatprev, by = c("ObsDaymin", "ObsDaymax"))

# overwrite
DNK.agebands.dat$seroprevMCMC <- dnk_adj_seroprev

#......................
# save out
#......................
dir.create("data/derived/DNK1", recursive = T)
saveRDS(DNK.regions.dat, "data/derived/DNK1/DNK1_regions.RDS")
saveRDS(DNK.agebands.dat, "data/derived/DNK1/DNK1_agebands.RDS")
### save for age region analysis:
saveRDS(DNK.agebands_age_sero.dat, "data/derived/DNK1/DNK1_agebands_age_sero.dat.RDS")



#............................................................
#--- ESP1-2 #----
#...........................................................
sero_prevdfESP <- sero_prevdf %>%
  dplyr::mutate(age_high = ifelse(study_id == "ESP1-2" & age_low == 0 & age_high == 0, 0.99, # to make the cut easier
                                  age_high))

#......................
# regions
#......................
ESP.regions.dat <- process_data4(cum_tp_deaths = deathsdf,
                                 time_series_totdeaths_df = JHUdf,
                                 time_series_totdeaths_geocode = "ESP",
                                 population = populationdf,
                                 sero_val = sero_valdf,
                                 seroprev = sero_prevdfESP,
                                 get_descriptive_dat = TRUE,
                                 groupingvar = "region",
                                 study_ids = "ESP1-2",
                                 filtRegions = c("Andalucía", "Aragón", "Asturias, Principado de",  "Balears, Illes",  "Comunitat Valenciana",
                                                 "Canarias", "Cantabria", "Castilla La Mancha", "Castilla y León", "Cataluña", "Extremadura",
                                                 "Galicia",   "Rioja, La",  "Madrid, Comunidad de",
                                                 "Murcia, Región de", "Navarra, Comunidad Foral de", "País Vasco")) # limit to mainland Spain



#......................
# agebands
#......................
ESP.agebands.dat <- process_data4(cum_tp_deaths = deathsdf,
                                  time_series_totdeaths_df = JHUdf,
                                  time_series_totdeaths_geocode = "ESP",
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdfESP,
                                  get_descriptive_dat = TRUE,
                                  groupingvar = "ageband",
                                  study_ids = "ESP1-2",
                                  death_agebreaks = c(0, 10, 20, 30, 40,
                                                      50, 60, 70, 80, 90, 999),
                                  sero_agebreaks = c(0, 10, 20, 30, 40,
                                                     50, 60, 70, 80, 90, 999))
#......................
# MANUAL ADJUSTMENTS
#......................
# No adjustments to serology as there is alignment of deaths and seroprevalence age groups.

#......................
# save out
#......................
dir.create("data/derived/ESP1-2/", recursive = T)
saveRDS(ESP.agebands.dat, "data/derived/ESP1-2/ESP1-2_agebands.RDS")
saveRDS(ESP.regions.dat, "data/derived/ESP1-2/ESP1-2_regions.RDS")


#............................................................
#---- GBR3 #----
#...........................................................
GBR3TimeSeries <- readr::read_tsv("data/raw/deaths_time_series_subnat.tsv") %>%
  dplyr::filter(study_id == "GBR3") %>%
  dplyr::rename(date = date_end_survey,
                deaths = n_deaths) %>%
  dplyr::mutate(date = lubridate::dmy(date), # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format
                georegion = "ENG") %>%
  dplyr::select(c("date", "georegion", "deaths"))
#......................
# ages
#......................
GBR3.agebands.dat <- process_data4(cum_tp_deaths = deathsdf,
                                   time_series_totdeaths_df = GBR3TimeSeries,
                                   time_series_totdeaths_geocode = "ENG",
                                   population = populationdf,
                                   sero_val = sero_valdf,
                                   seroprev = sero_prevdf,
                                   get_descriptive_dat = TRUE,
                                   groupingvar = "ageband",
                                   study_ids = "GBR3",
                                   death_agebreaks = c(0, 44, 64, 74, 999),
                                   sero_agebreaks = c(0, 44, 64, 74, 999))
#......................
# regions
#......................
GBR3.regions.dat <- process_data4(cum_tp_deaths = deathsdf,
                                  time_series_totdeaths_df = GBR3TimeSeries,
                                  time_series_totdeaths_geocode = "ENG",
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  get_descriptive_dat = TRUE,
                                  groupingvar = "region",
                                  study_ids = "GBR3",
                                  death_agebreaks = c(0, 999), # right now GBR3 popdf not refined at all
                                  filtRegions = c("North East",
                                                  "North West",
                                                  "Yorkshire",
                                                  "East Midlands",
                                                  "West Midlands",
                                                  "East of England",
                                                  "London",
                                                  "South East",
                                                  "South West"))

#......................
# MANUAL ADJUSTMENTS
#......................
# none needed with linelist

#......................
# save out
#......................
dir.create("data/derived/GBR3", recursive = T)
saveRDS(GBR3.agebands.dat, "data/derived/GBR3/GBR3_agebands.RDS")
saveRDS(GBR3.regions.dat, "data/derived/GBR3/GBR3_regions.RDS")


#............................................................
#---- ITA1 #----
#...........................................................
#......................
# ages
#......................
ITA.agebands.dat <- process_data4(cum_tp_deaths = deathsdf,
                                  time_series_totdeaths_df = JHUdf,
                                  time_series_totdeaths_geocode = "ITA",
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  get_descriptive_dat = TRUE,
                                  groupingvar = "ageband",
                                  study_ids = "ITA1",
                                  death_agebreaks = c(0,9,19,29,39,49,59,69,79,89, 999),
                                  sero_agebreaks = c(0, 17,  34,  49,  59,  69, 999))
#......................
# regions
#......................
# regional population data for ITA available
ITA.regions.dat <- process_data4(cum_tp_deaths = deathsdf,
                                  time_series_totdeaths_df = JHUdf,
                                  time_series_totdeaths_geocode = "ITA",
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  get_descriptive_dat = TRUE,
                                  groupingvar = "region",
                                  study_ids = "ITA1")




#......................
# MANUAL ADJUSTMENTS
#......................
# some lack of overlap between serology and deaths data.
# also no seropositives numbers, so must use logit method
# use sero 0-17 for 9-19 yr olds
# use the mean of 18-34 and 35-49 for the 30-39 group (they are v similar anyway)

agebands <- unique(ITA.agebands.dat$deaths_propMCMC$ageband)
ita_adj_seroprev <- tibble::tibble(
  ObsDaymin = unique(ITA.agebands.dat$seroprevMCMC$ObsDaymin),
  ObsDaymax = unique(ITA.agebands.dat$seroprevMCMC$ObsDaymax),
  ageband = agebands,
  n_positive = NA,
  n_tested = NA,
  SeroPrev = NA,
  SeroLCI = NA,
  SeroUCI = NA)

ita_adj_seroprev$SeroLCI[1:2] <- ITA.agebands.dat$seroprev_group$range_sero_low[1]
ita_adj_seroprev$SeroLCI[3] <- ITA.agebands.dat$seroprev_group$range_sero_low[2]
ita_adj_seroprev$SeroLCI[4] <- mean(ITA.agebands.dat$seroprev_group$range_sero_low[2:3])
ita_adj_seroprev$SeroLCI[5] <- ITA.agebands.dat$seroprev_group$range_sero_low[3]
ita_adj_seroprev$SeroLCI[6:7] <- ITA.agebands.dat$seroprev_group$range_sero_low[4:5]
ita_adj_seroprev$SeroLCI[8:10] <- ITA.agebands.dat$seroprev_group$range_sero_low[6]
ita_adj_seroprev$SeroUCI[1:2] <- ITA.agebands.dat$seroprev_group$range_sero_high[1]
ita_adj_seroprev$SeroUCI[3] <- ITA.agebands.dat$seroprev_group$range_sero_high[2]
ita_adj_seroprev$SeroUCI[4] <- mean(ITA.agebands.dat$seroprev_group$range_sero_high[2:3])
ita_adj_seroprev$SeroUCI[5] <- ITA.agebands.dat$seroprev_group$range_sero_high[3]
ita_adj_seroprev$SeroUCI[6:7] <- ITA.agebands.dat$seroprev_group$range_sero_high[4:5]
ita_adj_seroprev$SeroUCI[8:10] <- ITA.agebands.dat$seroprev_group$range_sero_high[6]
ita_adj_seroprev$SeroPrev[1:2] <- ITA.agebands.dat$seroprev_group$seroprevalence_unadjusted[1]
ita_adj_seroprev$SeroPrev[3] <- ITA.agebands.dat$seroprev_group$seroprevalence_unadjusted[2]
ita_adj_seroprev$SeroPrev[4] <- mean(ITA.agebands.dat$seroprev_group$seroprevalence_unadjusted[2:3])
ita_adj_seroprev$SeroPrev[5] <- ITA.agebands.dat$seroprev_group$seroprevalence_unadjusted[3]
ita_adj_seroprev$SeroPrev[6:7] <- ITA.agebands.dat$seroprev_group$seroprevalence_unadjusted[4:5]
ita_adj_seroprev$SeroPrev[8:10] <- ITA.agebands.dat$seroprev_group$seroprevalence_unadjusted[6]


# overwrite
ITA.agebands.dat$seroprevMCMC <- ita_adj_seroprev

#......................
# save out
#......................
dir.create("data/derived/ITA1", recursive = T)
saveRDS(ITA.agebands.dat, "data/derived/ITA1/ITA1_agebands.RDS")
saveRDS(ITA.regions.dat, "data/derived/ITA1/ITA1_regions.RDS")



#............................................................
#---- LUX1 #----
#...........................................................

#......................
# ages
#......................
LUX.agebands.dat <- process_data4(cum_tp_deaths = deathsdf,
                                  time_series_totdeaths_df = JHUdf,
                                  time_series_totdeaths_geocode = "LUX",
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  get_descriptive_dat = TRUE,
                                  groupingvar = "ageband",
                                  study_ids = "LUX1",
                                  death_agebreaks = c(0, 29, 39, 49, 59, 69, 79, 999))


#......................
# MANUAL ADJUSTMENTS
#......................
# Only one ageband - assume all the for all age groups
agebands <- unique(LUX.agebands.dat$deaths_propMCMC$ageband)
lux_adj_seroprev <- tibble::tibble(
  ObsDaymin = unique(LUX.agebands.dat$seroprevMCMC$ObsDaymin),
  ObsDaymax = unique(LUX.agebands.dat$seroprevMCMC$ObsDaymax),
  ageband = agebands,
  n_positive = NA,
  n_tested = NA,
  SeroPrev = NA)
lux_adj_seroprev$n_positive <- LUX.agebands.dat$seroprev_group$n_positive
lux_adj_seroprev$n_tested <- LUX.agebands.dat$seroprev_group$n_tested
lux_adj_seroprev$SeroPrev <- LUX.agebands.dat$seroprev_group$seroprevalence_unadjusted

LUX.agebands.dat$seroprevMCMC <- lux_adj_seroprev

#......................
# save out
#......................
dir.create("data/derived/LUX1", recursive = T)
saveRDS(LUX.agebands.dat, "data/derived/LUX1/LUX1_agebands.RDS")




#............................................................
#---- NLD1 #-----
#...........................................................
#......................
# ages
#......................
NLD.agebands.dat <- process_data4(cum_tp_deaths = deathsdf,
                                  time_series_totdeaths_df = JHUdf,
                                  time_series_totdeaths_geocode = "NLD",
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  get_descriptive_dat = TRUE,
                                  groupingvar = "ageband",
                                  study_ids = "NLD1",
                                  death_agebreaks = c(0, 49, 59, 69, 79, 999),
                                  sero_agebreaks = c(0, 49, 59, 69, 79, 999))


#......................
# MANUAL ADJUSTMENTS
#......................
# Netherlands seroprevalence and deaths not perfectly aligned.
# Assumptions:
#   Let 18-30, 31-40, and 41-50 stand in for 0-49 ageband
#   Let 51-60 stand in for 49-59 ageband
#   Let 60-72 stand in for 59-69, 69-79, 79++ ageband

timepoint1 <- NLD.agebands.dat$seroprevMCMC %>%
  dplyr::filter(ObsDaymax == 106)
agebands <- unique(NLD.agebands.dat$deaths_propMCMC$ageband)
nld_adj_seroprev <- tibble::tibble(
  ObsDaymin = 92, # only timepoint 1
  ObsDaymax = 106,
  ageband = c("79-999"),
  n_positive = 47,
  n_tested = 1742) %>%
  dplyr::mutate(SeroPrev = n_positive/n_tested)

nld_adj_seroprev <- dplyr::bind_rows(timepoint1,
                                     nld_adj_seroprev)

# time point2
timepoint2 <- tibble::tibble(
  ObsDaymin = 131,
  ObsDaymax = 141,
  ageband = nld_adj_seroprev$ageband,
  n_positive = 385,
  n_tested = 7000) %>%
  dplyr::mutate(SeroPrev = n_positive/n_tested)

nld_adj_seroprev <- dplyr::bind_rows(nld_adj_seroprev, timepoint2)

# overwrite
NLD.agebands.dat$seroprevMCMC <- nld_adj_seroprev

#......................
# save out
#......................
dir.create("data/derived/NLD1", recursive = T)
saveRDS(NLD.agebands.dat, "data/derived/NLD1/NLD1_agebands.RDS")


#............................................................
#--- SWE1 #----
#...........................................................
SWEdeathsdf <- deathsdf %>%
  dplyr::filter(study_id == "SWE1" & region == "all" & gender == "both") %>%
  dplyr::mutate(region = "SWE",
                age_low = 0,
                age_high = 999,
                n_deaths = sum(n_deaths),
                for_regional_analysis = 1,
                age_breakdown = 1,
                gender_breakdown = 1) %>%
  dplyr::filter(!duplicated(.))

#......................
# ages
#......................
SWE.agebands.dat <- process_data4(cum_tp_deaths = deathsdf,
                                  time_series_totdeaths_df = JHUdf,
                                  time_series_totdeaths_geocode = "SWE",
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  get_descriptive_dat = TRUE,
                                  groupingvar = "ageband",
                                  study_ids = "SWE1")

#......................
# MANUAL ADJUSTMENTS
#......................
ndays <- length(SWE.agebands.dat$seroprevMCMC$ObsDaymin)
agelen <- length(SWE.agebands.dat$prop_pop$ageband)
# assume blood donors stand in for all ages
swe_adj_seroprev <- tibble::tibble(
  ObsDaymin = SWE.agebands.dat$seroprevMCMC$ObsDaymin,
  ObsDaymax = SWE.agebands.dat$seroprevMCMC$ObsDaymax,
  n_positive = NA,
  n_tested = NA)
swe_adj_seroprev <- lapply(SWE.agebands.dat$prop_pop$ageband, function(x){
  swe_adj_seroprev %>%
    dplyr::mutate(ageband = x)}) %>%
  dplyr::bind_rows() %>%
  dplyr::arrange(ObsDaymin, ObsDaymax, ageband)

SWEnatprev <- SWE.agebands.dat$seroprev_group %>%
  dplyr::select(c("ObsDaymin", "ObsDaymax", "seroprevalence_unadjusted", "range_sero_low", "range_sero_high")) %>%
  dplyr::rename(SeroPrev = seroprevalence_unadjusted,
                SeroLCI = range_sero_low,
                SeroUCI = range_sero_high)


swe_adj_seroprev <- dplyr::left_join(swe_adj_seroprev,
                                     SWEnatprev, by = c("ObsDaymin", "ObsDaymax"))

# overwrite
SWE.agebands.dat$seroprevMCMC <- swe_adj_seroprev

#......................
# save out
#......................
dir.create("data/derived/SWE1", recursive = T)
saveRDS(SWE.agebands.dat, "data/derived/SWE1/SWE1_agebands.RDS")




#..................................................................................
#---- Preprocess Asian Data  #-----
#..................................................................................
#............................................................
#---- CHN1 #----
#...........................................................
CHN1TimeSeries <- readr::read_tsv("data/raw/deaths_time_series_subnat.tsv") %>%
  dplyr::filter(study_id == "CHN1") %>%
  dplyr::rename(date = date_end_survey,
                deaths = n_deaths) %>%
  dplyr::mutate(date = lubridate::dmy(date), # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format
                georegion = "CHN") %>%
  dplyr::select(c("date", "georegion", "deaths"))

#......................
# ages
#......................
CHN.agebands.dat <- process_data4(cum_tp_deaths = deathsdf,
                                  time_series_totdeaths_df = CHN1TimeSeries,
                                  time_series_totdeaths_geocode = "CHN",
                                  population = populationdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  get_descriptive_dat = TRUE,
                                  groupingvar = "ageband",
                                  study_ids = "CHN1",
                                  death_agebreaks = c(0, 9, 19, 29, 39, 49, 59, 69, 79, 999),
                                  sero_agebreaks = c(0, 19, 24, 29, 34, 39, 44, 49, 54, 59, 999))


#......................
# MANUAL ADJUSTMENTS
#......................
agebands <- unique(CHN.agebands.dat$deaths_propMCMC$ageband)
chn_adj_seroprev <- tibble::tibble(
  ObsDaymin = unique(CHN.agebands.dat$seroprevMCMC$ObsDaymin),
  ObsDaymax = unique(CHN.agebands.dat$seroprevMCMC$ObsDaymax),
  ageband = agebands,
  n_positive = NA,
  n_tested = NA,
  SeroPrev = NA)

chn_adj_seroprev$n_positive[1:2] <- round(mean(CHN.agebands.dat$seroprev_group$n_positive))
chn_adj_seroprev$n_positive[3] <- sum(CHN.agebands.dat$seroprev_group$n_positive[2:3])
chn_adj_seroprev$n_positive[4] <- sum(CHN.agebands.dat$seroprev_group$n_positive[4:5])
chn_adj_seroprev$n_positive[5] <- sum(CHN.agebands.dat$seroprev_group$n_positive[6:7])
chn_adj_seroprev$n_positive[6] <- sum(CHN.agebands.dat$seroprev_group$n_positive[8:9])
chn_adj_seroprev$n_positive[7:9] <- CHN.agebands.dat$seroprev_group$n_positive[10]

chn_adj_seroprev$n_tested[1:2] <- round(mean(CHN.agebands.dat$seroprev_group$n_tested))
chn_adj_seroprev$n_tested[3] <- sum(CHN.agebands.dat$seroprev_group$n_tested[2:3])
chn_adj_seroprev$n_tested[4] <- sum(CHN.agebands.dat$seroprev_group$n_tested[4:5])
chn_adj_seroprev$n_tested[5] <- sum(CHN.agebands.dat$seroprev_group$n_tested[6:7])
chn_adj_seroprev$n_tested[6] <- sum(CHN.agebands.dat$seroprev_group$n_tested[8:9])
chn_adj_seroprev$n_tested[7:9] <- CHN.agebands.dat$seroprev_group$n_tested[10]

chn_adj_seroprev$SeroPrev[1:2] <- mean(CHN.agebands.dat$seroprev_group$seroprevalence_unadjusted)
chn_adj_seroprev$SeroPrev[3] <- sum(CHN.agebands.dat$seroprev_group$n_positive[2:3]) / sum(CHN.agebands.dat$seroprev_group$n_tested[2:3])
chn_adj_seroprev$SeroPrev[4] <- sum(CHN.agebands.dat$seroprev_group$n_positive[4:5]) / sum(CHN.agebands.dat$seroprev_group$n_tested[4:5])
chn_adj_seroprev$SeroPrev[5] <- sum(CHN.agebands.dat$seroprev_group$n_positive[6:7]) / sum(CHN.agebands.dat$seroprev_group$n_tested[6:7])
chn_adj_seroprev$SeroPrev[6] <- sum(CHN.agebands.dat$seroprev_group$n_positive[8:9]) / sum(CHN.agebands.dat$seroprev_group$n_tested[8:9])
chn_adj_seroprev$SeroPrev[7:9] <- sum(CHN.agebands.dat$seroprev_group$n_positive[10]) / sum(CHN.agebands.dat$seroprev_group$n_tested[10])


CHN.agebands.dat$seroprevMCMC <- chn_adj_seroprev

#......................
# save out
#......................
dir.create("data/derived/CHN1", recursive = T)
saveRDS(CHN.agebands.dat, "data/derived/CHN1/CHN1_agebands.RDS")


#..................................................................................
#---- Preprocess Africa Data  #-----
#..................................................................................
#............................................................
#---- KEN1 #----
#...........................................................
KENpop <- populationdf %>%
  dplyr::filter(study_id == "KEN1") %>%
  dplyr::mutate(age_high = ifelse(age_low == 0 & age_high == 0, 0.99, # to make the cut easier
                                  age_high))
#......................
# ages
#......................
KEN.agebands.dat <- process_data4(cum_tp_deaths = deathsdf,
                                  time_series_totdeaths_df = JHUdf,
                                  time_series_totdeaths_geocode = "KEN",
                                  population = KENpop,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  get_descriptive_dat = TRUE,
                                  groupingvar = "ageband",
                                  study_ids = "KEN1") # age breaks for demography/population df
#......................
# MANUAL ADJUSTMENTS
#......................
# Assume seroprevalence from 15-24 yr is appropriate for 0-9 and 9-19
# Assume seroprev from 15-24 and 25-34 averaged is appropriate for 19-29
# Assume seroprev from 25-34 and 35-44 averaged is appropriate for 29-39
# Assume seroprev from 35-44 and 45-54 averaged is appropriate for 39-49
# Assume seroprev from 45-54 and 55-64 averaged is appropriate for 49-59
# Assume seroprev from 55-64 is appropriate for 59-69
agebands <- unique(KEN.agebands.dat$deaths_propMCMC$ageband)
ken_adj_seroprev <- tibble::tibble(
  ObsDaymin = unique(KEN.agebands.dat$seroprevMCMC$ObsDaymin),
  ObsDaymax = unique(KEN.agebands.dat$seroprevMCMC$ObsDaymax),
  ageband = agebands,
  n_positive = NA,
  n_tested = NA,
  SeroPrev = NA)

ken_adj_seroprev$n_positive[1:2] <- KEN.agebands.dat$seroprev_group$n_positive[1]
ken_adj_seroprev$n_positive[3] <- round(mean(KEN.agebands.dat$seroprev_group$n_positive[1:2]))
ken_adj_seroprev$n_positive[4] <- round(mean(KEN.agebands.dat$seroprev_group$n_positive[2:3]))
ken_adj_seroprev$n_positive[5] <- round(mean(KEN.agebands.dat$seroprev_group$n_positive[3:4]))
ken_adj_seroprev$n_positive[6] <- round(mean(KEN.agebands.dat$seroprev_group$n_positive[4:5]))
ken_adj_seroprev$n_positive[7] <- KEN.agebands.dat$seroprev_group$n_positive[5]

ken_adj_seroprev$n_tested[1:2] <- KEN.agebands.dat$seroprev_group$n_tested[1]
ken_adj_seroprev$n_tested[3] <- round(mean(KEN.agebands.dat$seroprev_group$n_tested[1:2]))
ken_adj_seroprev$n_tested[4] <- round(mean(KEN.agebands.dat$seroprev_group$n_tested[2:3]))
ken_adj_seroprev$n_tested[5] <- round(mean(KEN.agebands.dat$seroprev_group$n_tested[3:4]))
ken_adj_seroprev$n_tested[6] <- round(mean(KEN.agebands.dat$seroprev_group$n_tested[4:5]))
ken_adj_seroprev$n_tested[7] <-  KEN.agebands.dat$seroprev_group$n_tested[5]
ken_adj_seroprev$SeroPrev <- ken_adj_seroprev$n_positive /  ken_adj_seroprev$n_tested
KEN.agebands.dat$seroprevMCMC <- ken_adj_seroprev

ken_adj_seroprev$SeroPrev<-ken_adj_seroprev$n_positive/ken_adj_seroprev$n_tested
KEN.agebands.dat$seroprevMCMC <- ken_adj_seroprev

#......................
# save out
#......................
dir.create("data/derived/KEN1", recursive = T)
saveRDS(KEN.agebands.dat, "data/derived/KEN1/KEN1_agebands.RDS")


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

# JHU data for new recast df
JHUdf <- readr::read_csv("data/raw/time_series_covid19_deaths_US.csv") %>%
  dplyr::select(-c("UID", "iso2", "code3", "FIPS", "Country_Region", "Lat", "Long_", "Combined_Key", "Population")) %>%
  tidyr::pivot_longer(., cols = -c("iso3", "Admin2", "Province_State"),
                      names_to = "date", values_to = "cumdeaths") %>%
  dplyr::filter(!is.na(Admin2)) %>%
  dplyr::mutate(Admin2sp = sub(" ", "-", Admin2),
                georegion = paste0(Province_State, "_", Admin2sp)) %>%
  dplyr::mutate(date = lubridate::mdy(date)) %>% # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format
  dplyr::group_by(georegion) %>% # group by region for daily deaths
  dplyr::mutate(deaths = cumdeaths - dplyr::lag(cumdeaths),
                deaths = ifelse(is.na(deaths), 0, deaths), # take care of first value
                deaths = ifelse(deaths < 1, 0, deaths)) %>% # take care of cumulative deaths corrections
  dplyr::select(c("date", "georegion", "deaths")) %>%
  dplyr::filter(date <= lubridate::mdy("08-17-2020")) %>%
  dplyr::ungroup(.)

# fill in from origin
georegions <- dplyr::group_keys(dplyr::group_by(JHUdf, georegion)) %>% dplyr::pull(.)
origindf <- tibble::as_tibble(expand.grid(date = seq(lubridate::ymd("2020-01-01"), (min(JHUdf$date)-1), by = "1 day"),
                                          georegion = list(georegions),
                                          deaths = 0)) %>%
  tidyr::unnest(cols = georegion)
# combine now
JHUdf <- dplyr::bind_rows(JHUdf, origindf) %>%
  dplyr::arrange(georegion, date)


#............................................................
#---- LA_CA1 #----
# Los Angeles, CA
#...........................................................
LACA_TSdeathsdf <- JHUdf %>%
  dplyr::filter(georegion %in% "California_Los-Angeles")


LA_CA.agebands.dat <- process_data4(cum_tp_deaths = deathsdf,
                                   time_series_totdeaths_df = LACA_TSdeathsdf,
                                   time_series_totdeaths_geocode = "California_Los-Angeles",
                                   population = populationdf,
                                   sero_val = sero_valdf,
                                   seroprev = sero_prevdf,
                                   groupingvar = "ageband",
                                   study_ids = "LA_CA1")

#......................
# MANUAL ADJUSTMENTS
#......................
# assume blood group donors are representative for all ages
laca_adj_seroprev <- tibble::tibble(
  ObsDaymin = LA_CA.agebands.dat$seroprevMCMC$ObsDaymin,
  ObsDaymax = LA_CA.agebands.dat$seroprevMCMC$ObsDaymax,
  ageband = LA_CA.agebands.dat$prop_pop$ageband,
  n_positive = LA_CA.agebands.dat$seroprevMCMC$n_positive,
  n_tested = LA_CA.agebands.dat$seroprevMCMC$n_tested,
  SeroPrev = LA_CA.agebands.dat$seroprevMCMC$SeroPrev)

# overwrite
LA_CA.agebands.dat$seroprevMCMC <- laca_adj_seroprev

#......................
# save out
#......................
dir.create("data/derived/USA/", recursive = TRUE)
saveRDS(LA_CA.agebands.dat, "data/derived/USA/LA_CA1_agebands.RDS")



#............................................................
#---- NYS_1: New York State #-----
#...........................................................
#......................
# preprocess New York State
#......................
# groupings of seroprevalence in the study:
wr <- c("New York_Westchester|New York_Rockland")
long_island <- c("New York_Nassau|New York_Suffolk")
nyc <- "New York_New-York|New York_Kings|New York_Bronx|New York_Queens|New York_Richmond"
upstate <- c("Albany|Allegany|Broome|Cattaraugus|Cayuga|Chautauqua|Chemung|Chenango|Clinton|Columbia|Cortland|Delaware|Dutchess|Erie|Essex|Franklin|Fulton|Genessee|Greene|Hamilton|Herkimer|Jefferson|Lewis|Livingston|Madison|Monroe|Montgomery|Niagara|Oneida|Onondaga|Ontario|Orange|Orleans|Oswego|Otsego|Putnam|Rensselaer|St. Lawrence|Saratoga|Schenectady|Schoharie|Schuyler|Seneca|Steuben|Sullivan|Tioga|Tompkins|Ulster|Warren|Washington|Wayne|Wyoming|Yates")

NYSJHU <- JHUdf %>%
  dplyr::filter(grepl("New York_", georegion))

wrd <- NYSJHU %>%
  dplyr::filter(grepl(wr, georegion)) %>%
  dplyr::mutate(georegion = "Westchester and Rockland")
long_islandd <- NYSJHU %>%
  dplyr::filter(grepl(long_island, georegion)) %>%
  dplyr::mutate(georegion = "Long Island")
nycd <- NYSJHU %>%
  dplyr::filter(grepl(nyc, georegion)) %>%
  dplyr::mutate(georegion = "New York City")
upstated <- NYSJHU %>%
  dplyr::filter(grepl(upstate, georegion)) %>%
  dplyr::mutate(georegion = "Upstate New York")
NYSJHU <- rbind(wrd,nycd, long_islandd,upstated)
NYSJHU <- NYSJHU %>%
  dplyr::group_by(date, georegion) %>%
  dplyr::summarise(deaths = sum(deaths))

NYSdeathsrgndf <- NYSJHU %>%
  dplyr::group_by(georegion) %>%
  dplyr::summarise(n_deaths = sum(deaths)) %>%
  dplyr::mutate(country = "USA",
                study_id = "NYS1",
                age_low = 0,
                age_high = 999,
                region = georegion,
                gender = "both",
                age_breakdown = 0,
                for_regional_analysis = 1,
                gender_breakdown = 0)



NYSJHU_tseries<- NYSJHU %>%
  dplyr::mutate(georegion="NYS") %>%
  dplyr::group_by(date,georegion) %>%
  dplyr::summarise(deaths = sum(deaths)) %>%
  dplyr::ungroup()

NYSpopdf <- populationdf %>%
  dplyr::filter(grepl("New York_", region))
wrp <- NYSpopdf %>%
  dplyr::filter(grepl(wr,region)) %>%
  dplyr::mutate(region="Westchester and Rockland")
long_islandp<-NYSpopdf %>%
  dplyr::filter(grepl(long_island,region)) %>%
  dplyr::mutate(region="Long Island")
nycp <- NYSpopdf %>%
  dplyr::filter(grepl(nyc,region)) %>%
  dplyr::mutate(region="New York City")
upstatep <- NYSpopdf %>%
  dplyr::filter(grepl(upstate,region)) %>%
  dplyr::mutate(region="Upstate New York")

NYpopdf <- rbind(wrp,long_islandp,nycp,upstatep)
NYpopdf <- NYpopdf %>%
  dplyr::group_by(region, age_low, age_high, gender) %>%
  dplyr::summarise(population = sum(population)) %>%
  dplyr::mutate(
    country = "USA",
    study_id = "NYS1",
    age_breakdown = 1,
    for_regional_analysis = 1,
    gender_breakdown = 1) %>%
  dplyr::ungroup(.)

NYSdeathsdf <- deathsdf %>%
  dplyr::filter(study_id == "NYS1") %>%
  dplyr::mutate(
    country = "USA",
    study_id = "NYS1",
    region = "NYS",
    #age_low = 0,
    #age_high = 999,
    gender = "both",
    age_breakdown = 1,
    gender_breakdown = 1,
    for_regional_analysis = 1)


#......................
# age
#......................
NYS.age.dat <- process_data4(cum_tp_deaths = NYSdeathsdf,
                                  time_series_totdeaths_df = NYSJHU_tseries,
                                  time_series_totdeaths_geocode = "NYS",
                                  population = NYpopdf,
                                  sero_val = sero_valdf,
                                  seroprev = sero_prevdf,
                                  get_descriptive_dat = TRUE,
                                  groupingvar = "ageband",
                                  study_ids = "NYS1",
                                  death_agebreaks = c(0, seq(9, 89, 10), 999))

#......................
# regions
#......................
NYS.region.dat <- process_data4(cum_tp_deaths = NYSdeathsrgndf,
                                time_series_totdeaths_df = NYSJHU_tseries,
                                time_series_totdeaths_geocode = "NYS",
                                population = NYpopdf,
                                sero_val = sero_valdf,
                                seroprev = sero_prevdf,
                                get_descriptive_dat = TRUE,
                                groupingvar = "region",
                                study_ids = "NYS1",
                                death_agebreaks = c(0, 999))


#......................
# MANUAL ADJUSTMENTS
#......................
# Age bands deaths and seroprev exactly non overlapping. Set all ages to be the weighted population average except:
# assume 60+ has seroprevalence of 55+
# weighted pop average from the paper for the whole study area is 0.125
nys_adj_seroprev <- tibble::tibble(
  ObsDaymin = unique(NYS.age.dat$seroprevMCMC$ObsDaymin),
  ObsDaymax = unique(NYS.age.dat$seroprevMCMC$ObsDaymax),
  ageband = unique(NYS.age.dat$deaths_propMCMC$ageband),
  n_positive = NA,
  n_tested = NA,
  SeroPrev = NA)
nys_adj_seroprev <- nys_adj_seroprev %>%
  dplyr::mutate(SeroPrev = ifelse(ageband=="59-69" | ageband=="69-79" | ageband=="79-89" |
                                           ageband=="89-999",NYS.age.dat$seroprevMCMC$SeroPrev[4],0.125),
                n_tested = ifelse(ageband=="59-69" | ageband=="69-79" | ageband=="79-89" |
                                           ageband=="89-999",NYS.age.dat$seroprevMCMC$n_tested[4],
                                         sum(NYS.age.dat$seroprevMCMC$n_tested[1:3])),
                n_positive=round(n_tested*SeroPrev))

# overwrite
NYS.age.dat$seroprevMCMC <- nys_adj_seroprev

#......................
# save out
#......................
saveRDS(NYS.region.dat, "data/derived/USA/NYS1_regions.RDS")
saveRDS(NYS.age.dat, "data/derived/USA/NYS1_agebands.RDS")





#..................................................................................
#---- Care Home Data Processing  #-----
#..................................................................................
dir.create("data/derived/carehomes/", recursive = TRUE)

### ESP1-2
ESP.agebands_noCH.dat <- remove_ch_deaths(ESP.agebands.dat,"ESP1-2")
saveRDS(ESP.agebands_noCH.dat, "data/derived/carehomes/ESP1-2_agebands_noCH.RDS")

### GBR3
GBR3.agebands_noCH.dat <- remove_ch_deaths(GBR3.agebands.dat,"GBR3")
saveRDS(GBR3.agebands_noCH.dat, "data/derived/carehomes/GBR3_agebands_noCH.RDS")

### CHE1
CHE1.agebands_noCH.dat <- remove_ch_deaths(CHE1.agebands.dat,"CHE1")
saveRDS(CHE1.agebands_noCH.dat, "data/derived/carehomes/CHE1_agebands_noCH.RDS")

### CHE2
CHE2.agebands_noCH.dat <- remove_ch_deaths(CHE2.agebands.dat,"CHE2")
saveRDS(CHE2.agebands_noCH.dat, "data/derived/carehomes/CHE2_agebands_noCH.RDS")

### NYC_NY_1
NYC_NY_1.agebands_noCH.dat <- remove_ch_deaths(NYC_NY_1.agebands.dat,"NYC_NY_1")
saveRDS(NYC_NY_1.agebands_noCH.dat, "data/derived/carehomes/NYC_NY_1_agebands_noCH.RDS")

#..................................................................................
#---- Confirmed Deaths Data Processing  #-----
#..................................................................................
# cumulative deaths
deathsdf <- readr::read_tsv("data/raw/cumulative_deaths.tsv") %>%
  dplyr::select(-c("ref", "notes")) %>%
  dplyr::mutate(date_start_survey = lubridate::ymd(date_start_survey), # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format
                date_end_survey = lubridate::ymd(date_end_survey))

# demography (non-US Census data)
populationdf <- readr::read_tsv("data/raw/population.tsv") %>%
  dplyr::select(-c("reference")) %>%
  dplyr::mutate(age_low = ifelse(age_low == 0 & age_high == 0, 1, age_low),
                age_high = ifelse(age_low == 1 & age_high == 0, 1, age_high))  # liftover "zero" year olds to be 1, 1 as well

dir.create("data/derived/confirmeddeaths/", recursive = TRUE)
#......................
# GBR3
#......................
gbr3TimeSeries <- readr::read_tsv("data/raw/deathsconfirmed_time_series.tsv") %>%
  dplyr::filter(study_id == "GBR3") %>%
  dplyr::rename(date = date_end_survey,
                deaths = n_deaths) %>%
  dplyr::mutate(georegion = "GBR3") %>%
  dplyr::select(c("date", "georegion", "deaths")) %>%
  dplyr::mutate(date = lubridate::ymd(date)) # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format


# note, using proportion of deaths from general pop deaths (missing confirmed deaths by age)
GBR3_confirmed_deaths <-  process_data4(cum_tp_deaths = deathsdf,
                                        time_series_totdeaths_df = gbr3TimeSeries,
                                        time_series_totdeaths_geocode = "GBR3",
                                        population = populationdf,
                                        sero_val = sero_valdf,
                                        seroprev = sero_prevdf,
                                        get_descriptive_dat = TRUE,
                                        groupingvar = "ageband",
                                        study_ids = "GBR3",
                                        death_agebreaks = c(0, 44, 64, 74, 999),
                                        sero_agebreaks = c(0, 44, 64, 74, 999))

#......................
# NYC1
#......................
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

nyc_cumdeaths <- readr::read_tsv("data/raw/cumulative_deathsconfirmed.tsv") %>%
  dplyr::select(c("country", "study_id", "age_low", "age_high", "n_deaths", "date_start_survey", "date_end_survey")) %>%
  dplyr::mutate(date_start_survey = lubridate::ymd(date_start_survey), # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format
                date_end_survey = lubridate::ymd(date_end_survey),
                region = "NYC_NY_1",
                gender= "both",
                age_breakdown = 1,
                for_regional_analysis = 0,
                gender_breakdown = 0)


nycTimeSeries <- readr::read_tsv("data/raw/deathsconfirmed_time_series.tsv") %>%
  dplyr::filter(study_id == "NYC_NY_1") %>%
  dplyr::rename(date = date_end_survey,
                deaths = n_deaths) %>%
  dplyr::mutate(georegion = "NYC_NY_1") %>%
  dplyr::select(c("date", "georegion", "deaths")) %>%
  dplyr::mutate(date = lubridate::ymd(date)) # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format


# note, using proportion of deaths from general pop deaths (missing confirmed deaths by age)
NYC_confirmed_deaths <-  process_data4(cum_tp_deaths = nyc_cumdeaths,
                                       time_series_totdeaths_df = nycTimeSeries,
                                       time_series_totdeaths_geocode = "NYC_NY_1",
                                       population = populationdf,
                                       sero_val = sero_valdf,
                                       seroprev = sero_prevdf,
                                       get_descriptive_dat = TRUE,
                                       groupingvar = "ageband",
                                       study_ids = "NYC_NY_1",
                                       death_agebreaks = c(0, 18, 45, 65, 75, 999),
                                       sero_agebreaks = c(0, 18, 45, 65, 75, 999))

# write over from above for seroprev
NYC_confirmed_deaths$seroprevMCMC <- nyc_adj_seroprev


#......................
# save out
#......................
saveRDS(GBR3_confirmed_deaths, "data/derived/confirmeddeaths/GBR3_confirmed_deaths.rds")
saveRDS(NYC_confirmed_deaths, "data/derived/confirmeddeaths/NYC_confirmed_deaths.rds")

