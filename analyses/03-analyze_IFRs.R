##................................................................................................
## Purpose:
##
## Notes:
##................................................................................................
library(tidyverse)
library(COVIDCurve)
source("R/my_themes.R")
#......................
# quick functions
#......................
get_IFRs <- function(path) {
  modout <- readRDS(path)
  stratachar <- ifelse(grepl("age", basename(path)), "ageband",
                       ifelse(grepl("rgn", basename(path)), "region", NA))
  dictkey <- modout$inputs$IFRmodel$IFRdictkey %>%
    dplyr::rename(param = Strata)
  colnames(dictkey)[colnames(dictkey) == stratachar] <- "strata"
  # get ifrs
  ifrs <- COVIDCurve::get_cred_intervals(IFRmodel_inf = modout, whichrung = paste0("rung", 1),
                                         what = "IFRparams", by_chain = FALSE)
  # out
  dplyr::left_join(dictkey, ifrs) %>%
    dplyr::mutate(param = factor(param, levels = paste0("ma", 1:nrow(dictkey))),
                  strata = forcats::fct_reorder(strata, as.numeric(param)))
}

#......................
# read in data
#......................
rets <- list.files("results/ModFits/", full.names = T)
retmap <- tibble::tibble(study_id = stringr::str_split_fixed(basename(rets), "_", n = 3)[,1],
                         lvl =  ifelse(grepl("age", basename(rets)), "ageband", ifelse(grepl("rgn", basename(rets)), "region", NA)),
                         IFRret = purrr::map(rets, get_IFRs)) %>%
  tidyr::unnest(cols = IFRret) %>%
  dplyr::mutate(study_id = ifelse(study_id == "NYC", "NYC_NY_1", study_id))


#............................................................
# Age Specific Results
#...........................................................
colormap <- readr::read_csv("data/plot_aesthetics/color_studyid_map.csv")

retmap_age <- retmap %>%
  dplyr::filter(lvl == "ageband") %>%
  dplyr::left_join(., colormap, by = "study_id") %>%
  dplyr::mutate(age_mid = purrr::map_dbl(strata, function(x){
    nums <- as.numeric(stringr::str_split_fixed(x, "-", n = 2))
    nums[nums == 999] <- 100
    return(mean(nums))}))
ageIFR_plotObj <- ggplot() +
  geom_pointrange(data = retmap_age, aes(x = age_mid, y = median, ymin = LCI, ymax = UCI, color =  study_id), alpha = 0.5) +
  geom_line(data = retmap_age, aes(x = age_mid, y = median, group = study_id, color =  study_id),
            alpha = 0.4, size = 0.5, show.legend = FALSE) +
  scale_color_manual("Study ID", values = unique(retmap_age$cols)) +
  ylab("Age-Specific IFR (95% Cred. Int.)") + xlab("Mid. Age") +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#............................................................
# Region Specific Results
#...........................................................
#......................
# get proportion over 65
#......................
get_65prop <- function(path) {
  stratachar <- ifelse(grepl("ageband", basename(path)), "ageband",
                       ifelse(grepl("region", basename(path)), "region", NA))
  if (stratachar == "region") {
    datin <- readRDS(path)
    totpopdf <- datin$prop_pop %>%
      dplyr::group_by(region) %>%
      dplyr::summarise(totpop = sum(popN))

    datin$prop_pop %>%
      dplyr::mutate(age_mid = purrr::map_dbl(ageband, function(x){
        nums <- as.numeric(stringr::str_split_fixed(x, "-", n = 2))
        nums[nums == 999] <- 100
        return(mean(nums))})) %>%
      dplyr::filter(age_mid >= 65) %>%
      dplyr::group_by(region) %>%
      dplyr::summarise(pop65plus = sum(popN)) %>%
      dplyr::left_join(., totpopdf) %>%
      dplyr::mutate(prop65 = pop65plus/totpop) %>%
      dplyr::select(c("region", "prop65"))

  } else {
    return(NA)
  }
}

rgnprop65 <- readRDS("data/derived/descriptive_results_datamap.RDS") %>%
  dplyr::select(c("study_id", "breakdown", "relpath")) %>%
  dplyr::mutate(prop65 = purrr::map(relpath, get_65prop)) %>%
  dplyr::filter(breakdown == "region") %>%
  dplyr::select(c("study_id", "prop65")) %>%
  tidyr::unnest(cols = prop65)

#......................
# get proxy for healthcare capacity
#......................
study_id_liftover <- tibble::tibble(study_id = unique(retmap$study_id),
                                    iso_code = stringr::str_extract(unique(retmap$study_id), ".+?(?=[0-9]+)"))
owid_covid <- readr::read_csv("data/raw/owid-covid-data.csv") %>%
  dplyr::filter(iso_code %in% c("BRA", "CHE", "CHN", "DNK", "ESP", "GBR", "ITA", "LUX", "NLD", "SWE", "USA")) %>%
  dplyr::mutate(date= lubridate::ymd(date),
                newcases_per_hospbed = new_cases/(hospital_beds_per_thousand * 1e3),
                newdeaths_per_hospbed = new_deaths/(hospital_beds_per_thousand * 1e3),
                totcases_per_hospbed = total_cases/(hospital_beds_per_thousand * 1e3),
                totdeaths_per_hospbed = total_deaths/(hospital_beds_per_thousand * 1e3),
                newcases_smooth_per_hospbed = new_cases_smoothed/(hospital_beds_per_thousand * 1e3),
                newdeaths_smooth_per_hospbed = new_deaths_smoothed/(hospital_beds_per_thousand * 1e3)) %>%
  dplyr::inner_join(., y = study_id_liftover, by = "iso_code") %>%
  dplyr::group_by(study_id) %>%
  dplyr::summarise(maxnewdeaths_per_hospbed = max(newdeaths_per_hospbed, na.rm = T))

#......................
# make IFR plots for regions by prop age
#......................
retmap_rgn <- retmap %>%
  dplyr::filter(lvl == "region")

prop65rgnIFR_plotObj <- retmap_rgn %>%
  dplyr::rename(region = strata) %>%
  dplyr::left_join(., y = rgnprop65, by = c("study_id", "region")) %>%
  dplyr::left_join(., colormap, by = "study_id") %>%
  ggplot() +
  geom_pointrange(aes(x = prop65, y = median, ymin = LCI, ymax = UCI,
                      group = region, color =  study_id), alpha = 0.5) +
  scale_color_manual("Study ID", values = unique(retmap_age$cols)) +
  ylab("Regional IFR (95% Cred. Int.)") + xlab("Prop. of Population Over 65-Years") +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#......................
# make IFR plots for regions by prop age
#......................
health_cap_rgnIFR_plotObj <- retmap_rgn %>%
  dplyr::rename(region = strata) %>%
  dplyr::left_join(., y = owid_covid, by = c("study_id")) %>%
  dplyr::left_join(., colormap, by = "study_id") %>%
  dplyr::mutate(jitter = purrr::map_dbl(1:nrow(.), function(x){rnorm(1, 1e-2, 5e-3)}),
                maxnewdeaths_per_hospbed_jitt = maxnewdeaths_per_hospbed + jitter) %>%
  ggplot() +
  geom_pointrange(aes(x = maxnewdeaths_per_hospbed_jitt, y = median, ymin = LCI, ymax = UCI,
                      group = region, color =  study_id), alpha = 0.5) +
  scale_color_manual("Study ID", values = unique(retmap_age$cols)) +
  ylab("Regional IFR (95% Cred. Int.)") + xlab("Max. Prop. of Daily Deaths per Hospital Bed") +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#............................................................
# Main Figure
#...........................................................
# get common legend
commlegend <- cowplot::get_legend(ageIFR_plotObj +  guides(color = guide_legend(nrow = 1)) +
                       theme(legend.box.margin = margin(12, 0, 0, 0),
                             legend.position = "bottom"))

# remove legends
ageIFR_plotObj <- ageIFR_plotObj + theme(legend.position = "none")
prop65rgnIFR_plotObj <- prop65rgnIFR_plotObj + theme(legend.position = "none")
health_cap_rgnIFR_plotObj <- health_cap_rgnIFR_plotObj + theme(legend.position = "none")

# make some white space
whitespace <- tibble::tibble(
  xmin = -Inf, ymin = 0,
  xmax = Inf, ymax = 5
) %>%
  ggplot() +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymin), color = "#000000")  +
  theme_bw() +
  theme(
    plot.background = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank()
  )

# bring together
bottow_row <- cowplot::plot_grid(prop65rgnIFR_plotObj, health_cap_rgnIFR_plotObj,
                                align = "h", nrow = 1, labels = c("(B)", "(C)"))
main_fig <- cowplot::plot_grid(ageIFR_plotObj, whitespace, bottow_row, commlegend, labels = c("(A)", ""),
                               nrow = 4, rel_heights = c(1, 0.1, 1, 0.1))
main_fig
