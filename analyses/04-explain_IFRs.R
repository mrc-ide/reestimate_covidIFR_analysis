##................................................................................................
## Purpose:
##
## Notes:
##................................................................................................
library(tidyverse)
library(COVIDCurve)
source("R/my_themes.R")
get_overall_IFRs <- function(path) {
  modout <- readRDS(path)
  out <- COVIDCurve::get_globalIFR_cred_intervals(IFRmodel_inf = modout,
                                                  whichrung = "rung1",
                                                  by_chain = FALSE)
  return(out)
}


#............................................................
# read in data
#...........................................................
#......................
# results
#......................
regrets <- list.files("results/ModFits/", full.names = T)
regretmap <- tibble::tibble(study_id = toupper(stringr::str_split(basename(regrets), "_age|_rgn", simplify = T)[,1]),
                            lvl =  ifelse(grepl("age", basename(regrets)), "Age-Band", "Region"),
                            sero = "reg",
                            path = regrets)

serorevrets <- list.files("results/ModFits_SeroRev//", full.names = T)
serorevretmap <- tibble::tibble(study_id = toupper(stringr::str_split(basename(regrets), "_age|_rgn", simplify = T)[,1]),
                                lvl =  ifelse(grepl("age", basename(regrets)), "Age-Band", "Region"),
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
  dplyr::left_join(., y = prop65, by = "study_id") %>%
  dplyr::left_join(., plotmap, by = "study_id") %>%
  ggplot() +
  geom_pointrange(aes(x = prop65, y = median, ymin = LCI, ymax = UCI,
                      color =  study_id), alpha = 0.5) +
  scale_color_manual("Study ID", values = unique(plotmap$cols)) +
  ylab("Regional IFR (95% Cred. Int.)") + xlab("Prop. of Population Over 65-Years") +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#............................................................
#---- Heatlh Care Capacity  #----
# get proxy for healthcare capacity through OWID
#...........................................................

#......................
#plotObj
#......................
health_cap_rgnIFR_plotObj <- retmap %>%
  dplyr::left_join(., plotmap, by = "study_id") %>%
  dplyr::inner_join(., y = owid_covid, by = c("iso_code")) %>%
  ggplot() +
  geom_pointrange(aes(x = maxnewdeaths_per_hospbed, y = median, ymin = LCI, ymax = UCI,
                      group = region, color =  study_id), alpha = 0.5) +
  scale_color_manual("Study ID", values = unique(retmap_age$cols)) +
  ylab("Regional IFR (95% Cred. Int.)") + xlab("Max. Prop. of Daily Deaths per Hospital Bed") +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#............................................................
# Main Figure
#...........................................................
# get common legend
commlegend <- cowplot::get_legend(health_cap_rgnIFR_plotObj +  guides(color = guide_legend(nrow = 1)) +
                                    theme(legend.box.margin = margin(12, 0, 0, 0),
                                          legend.position = "bottom"))

# remove legends
prop65rgnIFR_plotObj <- prop65rgnIFR_plotObj + theme(legend.position = "none")
health_cap_rgnIFR_plotObj <- health_cap_rgnIFR_plotObj + theme(legend.position = "none")

# bring together
(main_fig <- cowplot::plot_grid(prop65rgnIFR_plotObj, health_cap_rgnIFR_plotObj,
                                 align = "h", nrow = 1, labels = c("(A)", "(B)")))
# out
jpeg("figures/final_figures/Figure4.jpg", width = 10, height = 7, units = "in", res = 500)
plot(main_fig)
graphics.off()




#............................................................
#---- IFR-Death Comparisons #----
#...........................................................
# TODO --> don't have excess or care home deaths? --> but tell Lucy easiest way to link is by study_id
percap_deaths <- dscdat %>%
  dplyr::filter(!duplicated(study_id)) %>%
  dplyr::select(c("study_id", "data")) %>%
  group_by(study_id) %>%
  dplyr::mutate(deaths = purrr::map(data, "deaths_TSMCMC"),
                deaths = purrr::map_dbl(deaths, function(x){sum(x$deaths)}),
                popN = purrr::map(data, "prop_pop"),
                popN = purrr::map_dbl(popN, function(x){sum(x$popN)})) %>%
  tidyr::unnest(cols = c("deaths", "popN")) %>%
  dplyr::summarise(cumdeaths = sum(deaths),
                   popN = sum(popN)) %>%
  dplyr::mutate(capita = cumdeaths/popN) %>%
  dplyr::select(c("study_id", "capita"))


# bring together others deaths
ext_deaths <- percap_deaths %>%
  #dplyr::left_join(., y = excess_deaths, by = "study_id") %>%
  #dplyr::left_join(., y = carehome_deaths, by = "study_id") %>%
  tidyr::pivot_longer(., cols = -c("study_id"), names_to = "param", values_to = "est") %>%
  dplyr::mutate(param = factor(param,
                               levels = c("crude", "capita", "excess"),
                               labels = c("Crude", "Per Capita", "Excess")
  ))

# plot out
(death_type_plot <- ggplot() +
  geom_point(data = ext_deaths, aes(x = study_id, y = est, fill = param, shape = param),
             size = 6.5, alpha = 0.8) +
  geom_pointrange(data = retmap,
                  aes(x = study_id, ymin = LCI, y = median, ymax = UCI), size = 0.75, alpha = 0.5) +
  scale_shape_manual("Death Type", values = c(22, 23, 24, 25)) +
  scale_fill_manual("Death Type", values = c(wesanderson::wes_palette("IsleofDogs2", type = "discrete"))) +
  ylab("Infection Fatality Rate") +
  coord_flip() +
  theme(
    plot.title =  element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(family = "Helvetica", face = "bold", vjust = 0.5, hjust = 0.5, size = 16),
    axis.text.x = element_text(family = "Helvetica", color = "#000000", vjust = 0.5, hjust = 0.5, size = 14),
    axis.text.y = element_text(family = "Helvetica", color = "#000000", face = "bold", vjust = 0.5, hjust = 1, size = 14),
    legend.title = element_text(family = "Helvetica", face = "bold", vjust = 0.5, hjust = 0.5, size = 16),
    legend.text = element_text(family = "Helvetica", vjust = 0.8, hjust = 0.5, size = 14, angle = 0),
    legend.position = "right",
    axis.line.x = element_line(color = "black", size = 1.5),
    axis.line.y = element_line(color = "black", size = 1.5),
    axis.ticks.y = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent"),
    panel.grid = element_blank(),
    panel.border = element_blank()
  ))

# out
jpeg("figures/final_figures/Figure5.jpg", width = 10, height = 7, units = "in", res = 500)
plot(death_type_plot)
graphics.off()

