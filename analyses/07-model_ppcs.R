##................................................................................................
## Purpose: Analyze differences in IFRs
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

#............................................................
# read in results
#...........................................................
# no serorev
mod_NOserorev_paths <- list.files("results/Modfits_noserorev/", full.names = T)
mod_NOserorev_retmap <- tibble::tibble(study_id = toupper(stringr::str_split(basename(mod_NOserorev_paths), "_age|_carehomes", simplify = T)[,1]),
                                       lvl = "NoSeroRev",
                                       paths = mod_NOserorev_paths)

# yes serorev
mod_serorev_paths <- list.files("results/Modfits_serorev/", full.names = T)
mod_serorev_retmap <- tibble::tibble(study_id = toupper(stringr::str_split(basename(mod_serorev_paths), "_age|_carehomes", simplify = T)[,1]),
                                     lvl = "SeroRev",
                                     paths = mod_serorev_paths)
# carehomes
mod_carehome_paths <- list.files("results/Modfits_carehomes//", full.names = T)
mod_carehomes_retmap <- tibble::tibble(study_id = toupper(stringr::str_split(basename(mod_carehome_paths), "_age|_carehomes", simplify = T)[,1]),
                                       lvl = "CareHomes",
                                       paths = mod_carehome_paths)

#.....................
# come together
#......................
datmap <- dplyr::bind_rows(mod_NOserorev_retmap, mod_serorev_retmap, mod_carehomes_retmap)
datmap <- datmap %>%
  dplyr::mutate(modout = purrr::map(paths, readRDS))

#............................................................
# internal functions to get result tables
#...........................................................
# sens/spec
get_study_seros_posts <- function(modout) {
  COVIDCurve::get_cred_intervals(IFRmodel_inf = modout, what = "Serotestparams",
                                 whichrung = "rung1", by_chain = F) %>%
    dplyr::mutate_if(is.numeric, round, 2) %>%
    dplyr::mutate(outcol = paste0(median, " (", LCI, ", ", UCI, ")" )) %>%
    dplyr::filter(param %in% c("sens", "spec")) %>%
    dplyr::select(c("param", "outcol")) %>%
    tidyr::pivot_wider(., id_cols = "param", names_from = "param",
                       values_from = "outcol")
}

# delays
get_study_delays_posts <- function(modout) {
  sero_delays <- COVIDCurve::get_cred_intervals(IFRmodel_inf = modout, what = "Serotestparams",
                                 whichrung = "rung1", by_chain = F) %>%
    dplyr::filter(!param %in% c("sens", "spec"))
  death_delays <- COVIDCurve::get_cred_intervals(IFRmodel_inf = modout, what = "DeathDelayparams",
                                            whichrung = "rung1", by_chain = F)

    dplyr::bind_rows(sero_delays, death_delays) %>%
    dplyr::mutate_if(is.numeric, round, 2) %>%
    dplyr::mutate(outcol = paste0(median, " (", LCI, ", ", UCI, ")" )) %>%
    dplyr::select(c("param", "outcol")) %>%
    tidyr::pivot_wider(., id_cols = "param", names_from = "param",
                       values_from = "outcol")

}


#............................................................
#---- Posterior Sero Characteristics and Noise Effects Table  #----
#...........................................................
sero_char_tab <- datmap %>%
  dplyr::mutate(tbl = purrr::map(modout, get_study_seros_posts))

#............................................................
#---- Posterior Delay Param Table  #----
#...........................................................
delay_char_tab <- datmap %>%
  dplyr::mutate(tbl = purrr::map(modout, get_study_delays_posts))






#............................................................
# internal functions for PPC plots
#...........................................................
# draw posterior seroprev simply
get_seroprev_ppcs <- function(modout) {

  # demog information
  demog <- modout$inputs$IFRmodel$demog %>%
    dplyr::left_join(., modout$inputs$IFRmodel$IFRdictkey) %>%
    dplyr::rename(param = Strata)

  # seroprev points
  seropnts <- COVIDCurve::draw_posterior_sero_curves(IFRmodel_inf = modout,
                                                     dwnsmpl = 1e2,
                                                     by_chain = F)
  # serocurve data
  serocurvedat <- seropnts %>%
    dplyr::select(c("sim", "ObsDay", dplyr::starts_with("RG_pd_"),
                    dplyr::starts_with("crude_pd_"))) %>%
    tidyr::pivot_longer(., cols = -c("sim", "ObsDay"),
                        names_to = "seroprev_strata_lvl", values_to = "seroprev") %>%
    dplyr::mutate(seroprevlvl = ifelse(stringr::str_detect(seroprev_strata_lvl, "RG_"), "RG Corr.", "Crude"),
                  param = stringr::str_extract(seroprev_strata_lvl, "ma[0-9]+")) %>%
    dplyr::left_join(., demog)

  # dwnsmpl to oldest age group to not overwhelm figure
  oldest <- modout$inputs$IFRmodel$IFRdictkey$ageband[length(modout$inputs$IFRmodel$IFRdictkey$ageband)]
  serocurvedat %>%
    dplyr::filter(ageband == oldest)

}


# draw posterior deaths simply
get_death_ppcs <- function(modout) {
  #......................
  # get deaths posterior pred check
  #......................
  postdat <- COVIDCurve::posterior_check_infxns_to_death(IFRmodel_inf = modout,
                                                         dwnsmpl = 1e2,
                                                         by_chain = FALSE)
  postdat_long <- postdat %>%
    dplyr::select(c("sim", "time", dplyr::starts_with("deaths"))) %>%
    tidyr::gather(., key = "Strata", value = "deaths", 3:ncol(.)) %>%
    dplyr::mutate(Strata = gsub("deaths_", "", Strata)) %>%
    dplyr::left_join(., y = modout$inputs$IFRmodel$IFRdictkey) %>%
    dplyr::rename(inf_deaths = deaths)
  #......................
  # get deaths data
  #......................
  # recast deaths
  proplist <- split(modout$inputs$IFRmodel$data$prop_deaths, 1:nrow(modout$inputs$IFRmodel$data$prop_deaths))

  deathrecast <- lapply(proplist,
                        function(x){
                          tibble::tibble(
                            Strata = x$Strata,
                            Deaths = x$PropDeaths * modout$inputs$IFRmodel$data$obs_deaths$Deaths)}) %>%
    dplyr::bind_rows(.) %>%
    dplyr::group_by(Strata) %>%
    dplyr::mutate(ObsDay = 1:nrow(modout$inputs$IFRmodel$data$obs_deaths)) %>%
    dplyr::ungroup(.)

  datclean <-  deathrecast %>%
    dplyr::left_join(., modout$inputs$IFRmodel$IFRdictkey, by = "Strata") %>%
    dplyr::rename(obs_deaths = Deaths)

  #......................
  # combine
  # note, truth is being duplicated many times...
  #......................
  # dwnsmpl to oldest age group to not overwhelm figure
  oldest <- modout$inputs$IFRmodel$IFRdictkey$ageband[length(modout$inputs$IFRmodel$IFRdictkey$ageband)]
  dplyr::left_join(postdat_long, datclean) %>%
    dplyr::filter(ageband == oldest)


}


#............................................................
#---- PPC Sero Plots  #----
#...........................................................
dsc_agedat <- readRDS("results/descriptive_results/descriptive_results_datamap.RDS") %>%
  dplyr::filter(breakdown == "ageband")
# crude data
crudedat <- dsc_agedat %>%
  dplyr::select(c("study_id", "plotdat")) %>%
  tidyr::unnest(cols = "plotdat")

SeroPrevObs <- crudedat %>%
  dplyr::select(c("study_id", "ageband", "obsdaymin", "obsdaymax", "seroprev")) %>%
  dplyr::mutate(obsmidday = (obsdaymin + obsdaymax)/2)

#..................................
# standard model
#...................................
noserorev_ppc_seroPlotObj

datmap %>%
  dplyr::filter(lvl == "NoSeroRev") %>%
  dplyr::mutate(plotdat = purrr::map(modout, get_seroprev_ppcs)) %>%
  tidyr::unnest(cols = plotdat) %>%
  dplyr::left_join(., SeroPrevObs)

