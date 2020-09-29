##................................................................................................
## Purpose: Analyze care home specific IFRs
##
## Notes:
##................................................................................................
library(tidyverse)
library(COVIDCurve)
source("R/my_themes.R")
source("R/covidcurve_helper_functions.R")
source("R/delta_method.R")
# colors
study_cols <- readr::read_csv("data/plot_aesthetics/color_studyid_map.csv")
mycolors <- study_cols$cols
names(mycolors) <- study_cols$study_id
studyidnames <- study_cols %>%
  dplyr::select(c("study_id", "names")) %>%
  dplyr::filter(!is.na(names))

# get seroprevalence posterior for carehome denominator
get_seroprev_est_se <- function(path, date) {
  ret <- get_strata_seroprevs(path)
  oldest <- get_data_dict(path)
  oldest <- oldest$Strata[length(oldest$Strata)]
  # get inferred observed seroprev
  ret %>%
    dplyr::select(c("ObsDay", "sim", paste0("RG_pd_seroprev_", oldest))) %>%
    dplyr::rename(est =  paste0("RG_pd_seroprev_", oldest)) %>%
    dplyr::mutate(obsdate = ObsDay + lubridate::ymd("2020-01-01") - 1) %>%
    dplyr::filter(obsdate == date) %>%
    dplyr::summarise(
      LCI = quantile(est, 0.025),
      median = median(est),
      UCI = quantile(est, 0.975),
      SE = (COVIDCurve:::logit(UCI + .Machine$double.xmin) - COVIDCurve:::logit(LCI + .Machine$double.xmin)) / (2*1.96)
    )
}


#............................................................
# read in data
#...........................................................
# care home orig df
carehomedf <- readr::read_csv("data/raw/care_home_deaths.csv")
# observed data
dscdat <- readRDS("results/descriptive_results/descriptive_results_datamap.RDS")
dsc_agedat <- dscdat %>%
  dplyr::filter(breakdown == "ageband") %>%
  dplyr::filter(!grepl("_nch", study_id)) %>%
  dplyr::filter(study_id %in% carehomedf$study_id) %>%
  dplyr::select(c("study_id", "std_deaths")) %>%
  tidyr::unnest(cols = "std_deaths") %>%
  dplyr::mutate(obsdate = obsday + lubridate::ymd("20200101") - 1) %>%
  dplyr::group_by(study_id, obsdate) %>%
  dplyr::summarise(cumdeaths = sum(cumdeaths),
                   popn = sum(popn)) # sum over agebands



#......................
# Deal with deaths
#......................
deathdat <- dplyr::left_join(carehomedf, dsc_agedat, by = "study_id") %>%
  dplyr::mutate(date = lubridate::dmy(date)) %>%
  dplyr::filter(date == obsdate) %>%
  dplyr::mutate(no_ch_cum_deaths = cumdeaths - care_home_deaths) %>%
  dplyr::select(c("study_id", "no_ch_cum_deaths", "popn"))

#......................
# Deal with seroprevs
#......................
regrets <- list.files("results/Modfits_noserorev/", full.names = T)
regretmap <- tibble::tibble(study_id = toupper(stringr::str_split(basename(regrets), "_age", simplify = T)[,1]),
                            sero = "reg",
                            path = regrets) %>%
  dplyr::filter(study_id %in% carehomedf$study_id) %>%
  dplyr::left_join(., carehomedf, by = "study_id") %>%
  dplyr::mutate(date = lubridate::dmy(date))
regretmap$seroprevs <- purrr::pmap(regretmap[, c("path", "date")], get_seroprev_est_se)
# get bits I need
regretmap_serosumm <- regretmap %>%
  dplyr::group_by(study_id) %>%
  dplyr::mutate(median = purrr::map_dbl(seroprevs, "median"),
                SE = purrr::map_dbl(seroprevs, "SE")) %>%
  dplyr::select(c("study_id", "median", "SE")) %>%
  dplyr::rename(seroprev = median) %>%
  dplyr::ungroup(.)

#......................
# get IFRs
#......................
care_home_crude_IFRs <- dplyr::left_join(deathdat, regretmap_serosumm) %>%
  dplyr::group_by(study_id) %>%
  dplyr::mutate(IFRcalc = no_ch_cum_deaths  / (seroprev * popn + no_ch_cum_deaths),
                IFRbound = purrr::map(seroprev, get_delta_CI_vals, deaths = no_ch_cum_deaths, popN = popn, SE = SE, tol = 1e-4),
                lower_ci = purrr::map_dbl(IFRbound, "lower.ci"),
                upper_ci = purrr::map_dbl(IFRbound, "upper.ci"))
# out
care_home_crude_IFRs %>%
  dplyr::mutate(no_ch_cum_deaths = prettyNum(no_ch_cum_deaths, big.mark=",", scientific=FALSE),
                IFRcalc = round(IFRcalc * 100, 2),
                lower_ci = round(lower_ci * 100 , 2),
                upper_ci = round(upper_ci * 100, 2),
                carehome_col = paste0(IFRcalc, " (", lower_ci, ", ", upper_ci, ")")) %>%
  dplyr::select(c("study_id", "no_ch_cum_deaths", "carehome_col")) %>%
  readr::write_tsv(., path = "tables/final_tables/carehome_crude_data.tsv")



