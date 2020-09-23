##................................................................................................
## Purpose:
##
## Notes:
##................................................................................................
library(tidyverse)
library(COVIDCurve)
source("R/my_themes.R")
source("R/covidcurve_helper_functions.R")
# colors
study_cols <- readr::read_csv("data/plot_aesthetics/color_studyid_map.csv")
mycolors <- study_cols$cols
names(mycolors) <- study_cols$study_id


#............................................................
# read in data
#...........................................................
dscdat <- readRDS("results/descriptive_results/descriptive_results_datamap.RDS") %>%
  dplyr::mutate(drop = ifelse(study_id %in% c("DNK1", "GBR3", "ESP1-2", "BRA1", "CHE1", "ITA1")
                              & breakdown == "region", T, # drop regions that aren't basic
                              ifelse(grepl("_nch", study_id), T, F)) # drop carehomes
  ) %>%
  dplyr::filter(!drop)

#......................
# results
#......................
regrets <- list.files("results/Modfits_noserorev/", full.names = T)
regretmap <- tibble::tibble(study_id = toupper(stringr::str_split(basename(regrets), "_age|_rgn", simplify = T)[,1]),
                            lvl =  ifelse(grepl("age", basename(regrets)), "Age-Band", "Region"),
                            sero = "reg",
                            path = regrets)

serorevrets <- list.files("results/ModFits_SeroRev/", full.names = T)
serorevretmap <- tibble::tibble(study_id = toupper(stringr::str_split(basename(serorevrets), "_age|_rgn", simplify = T)[,1]),
                                lvl =  ifelse(grepl("age", basename(serorevrets)), "Age-Band", "Region"),
                                sero = "serorev",
                                path = serorevrets)
# bring together
retmap <- dplyr::bind_rows(regretmap, serorevretmap) %>%
  dplyr::mutate(overallIFRret = purrr::map(path, get_overall_IFRs)) %>%
  tidyr::unnest(cols = "overallIFRret")

#......................
# extra data
#......................
# descriptive data from our study
dscdat <- readRDS("data/derived/descriptive_results_datamap.RDS")

# our world in data
owid_covid <- readr::read_csv("data/raw/owid-covid-data.csv") %>%
  dplyr::filter(iso_code %in% c("BRA", "CHE", "CHN", "DNK", "ESP", "GBR", "ITA", "LUX", "NLD", "SWE", "USA")) %>%
  dplyr::mutate(date = lubridate::ymd(date),
                newcases_per_hospbed = new_cases/(hospital_beds_per_thousand * 1e3),
                newdeaths_per_hospbed = new_deaths/(hospital_beds_per_thousand * 1e3),
                totcases_per_hospbed = total_cases/(hospital_beds_per_thousand * 1e3),
                totdeaths_per_hospbed = total_deaths/(hospital_beds_per_thousand * 1e3),
                newcases_smooth_per_hospbed = new_cases_smoothed/(hospital_beds_per_thousand * 1e3),
                newdeaths_smooth_per_hospbed = new_deaths_smoothed/(hospital_beds_per_thousand * 1e3)) %>%
  dplyr::group_by(iso_code) %>%
  dplyr::summarise(maxnewdeaths_per_hospbed = max(newdeaths_per_hospbed, na.rm = T))

# plot aesthetics
plotmap <- readr::read_csv("data/plot_aesthetics/color_studyid_map.csv") %>%
  dplyr::rename(iso_code = national_georegion)

#............................................................
#---- Age Structure #----
#...........................................................
# get proportion over 65
get_65prop <- function(path) {
  # read data
  datin <- readRDS(path)
  # get total pop
  totpop <- datin$prop_pop %>%
    dplyr::summarise(totpop = sum(popN)) %>%
    dplyr::pull("totpop")
  # prop over 65 by mean of age split
  datin$prop_pop %>%
    dplyr::mutate(age_mid = purrr::map_dbl(ageband, function(x){
      nums <- as.numeric(stringr::str_split_fixed(x, "-", n = 2))
      nums[nums == 999] <- 100
      return(mean(nums))})) %>%
    dplyr::filter(age_mid >= 65) %>%
    dplyr::summarise(pop65plus = sum(popN)) %>%
    dplyr::mutate(prop65 = pop65plus/totpop) %>%
    dplyr::pull("prop65")
}

prop65 <- dscdat %>%
  dplyr::filter(breakdown == "ageband") %>% # basic models are 0-999 so don't know 65+s
  dplyr::select(c("study_id", "relpath")) %>%
  dplyr::mutate(prop65 = purrr::map_dbl(relpath, get_65prop)) %>%
  dplyr::select(c("study_id", "prop65")) %>%
  dplyr::filter(!duplicated(.))

#......................
# plotObj
#......................
prop65rgnIFR_plotObj <- retmap %>%
  dplyr::filter(lvl == "Age-Band") %>%
  dplyr::left_join(., y = prop65, by = "study_id") %>%
  dplyr::left_join(., plotmap, by = "study_id") %>%
  dplyr::mutate(sero = factor(sero, levels = c("reg", "serorev"), labels = c("SeroInf.", "SeroRev"))) %>%
  ggplot() +
  geom_pointrange(aes(x = prop65, y = median, ymin = LCI, ymax = UCI,
                      color =  study_id), alpha = 0.5) +
  facet_wrap(.~sero) +
  scale_color_manual("Study ID", values = unique(plotmap$cols)) +
  ylab("Overall IFR (95% Cred. Int.)") + xlab("Prop. of Population Over 65-Years") +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#............................................................
#---- Heatlh Care Capacity  #----
# get proxy for healthcare capacity through OWID
#...........................................................
# https://github.com/owid/covid-19-data/tree/master/public/data
# https://github.com/owid/covid-19-data/blob/master/public/data/owid-covid-codebook.csv
#......................
#plotObj
#......................
health_cap_rgnIFR_plotObj <- retmap %>%
  dplyr::mutate(drop = ifelse(study_id %in% c("DNK1", "GBR3", "ESP1-2", "BRA1")
                              & lvl == "Region", T, F),
                sero = factor(sero, levels = c("reg", "serorev"), labels = c("SeroInf.", "SeroRev."))) %>%
  dplyr::filter(!drop) %>%
  dplyr::left_join(., plotmap, by = "study_id") %>%
  dplyr::inner_join(., y = owid_covid, by = c("iso_code")) %>%
  ggplot() +
  geom_pointrange(aes(x = maxnewdeaths_per_hospbed, y = median, ymin = LCI, ymax = UCI,
                      color =  study_id), alpha = 0.5) +
  facet_wrap(.~sero) +
  scale_color_manual("Study ID", values = unique(plotmap$cols)) +
  ylab("Overall IFR (95% Cred. Int.)") + xlab("Max. Ratio of Daily Deaths per Hospital Bed") +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#............................................................
# Main Figure
#...........................................................

# bring together
(main_fig <- cowplot::plot_grid(prop65rgnIFR_plotObj, health_cap_rgnIFR_plotObj,
                                 align = "h", ncol = 1, nrow = 2, labels = c("(A)", "(B)")))
# out
jpeg("figures/final_figures/Figure_explain_IFRs.jpg", width = 10, height = 7, units = "in", res = 500)
plot(main_fig)
graphics.off()

