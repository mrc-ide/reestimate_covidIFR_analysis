##................................................................................................
## Purpose: Get Posteriors for Age Based Model
##
## Notes:
##................................................................................................
library(tidyverse)
library(COVIDCurve)
source("R/my_themes.R")
source("R/covidcurve_helper_functions.R")
# colors now based on location
locatkey <- readr::read_csv("data/plot_aesthetics/color_studyid_map.csv")
mycolors <- locatkey$cols
names(mycolors) <- locatkey$location
# order
order <- readr::read_csv("data/plot_aesthetics/study_id_order.csv")

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

#.....................
# come together
#......................
datmap <- dplyr::bind_rows(mod_NOserorev_retmap, mod_serorev_retmap)
datmap <- datmap %>%
  dplyr::mutate(modout = purrr::map(paths, readRDS))

#............................................................
# internal functions to get result tables
#...........................................................
# sens/spec
get_study_seros_posts <- function(modout) {
  COVIDCurve::get_cred_intervals(IFRmodel_inf = modout, what = "Serotestparams",
                                 whichrung = "rung50", by_chain = F) %>%
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
                                                whichrung = "rung50", by_chain = F) %>%
    dplyr::filter(!param %in% c("sens", "spec"))
  death_delays <- COVIDCurve::get_cred_intervals(IFRmodel_inf = modout, what = "DeathDelayparams",
                                                 whichrung = "rung50", by_chain = F)

  dplyr::bind_rows(sero_delays, death_delays) %>%
    dplyr::mutate_if(is.numeric, round, 2) %>%
    dplyr::mutate(outcol = paste0(median, " (", LCI, ", ", UCI, ")" )) %>%
    dplyr::select(c("param", "outcol")) %>%
    tidyr::pivot_wider(., id_cols = "param", names_from = "param",
                       values_from = "outcol")

}


#............................................................
#---- Posterior Sero Characteristics  #----
#...........................................................
#......................
# regional model
#......................
rgnmod <- list.files("results/prior_inputs/", pattern = "_reg_age.csv", full.names = T)
rgnmod <- tibble::tibble(study_id = stringr::str_split_fixed(basename(rgnmod), "_", n = 2)[,1],
                         path = rgnmod) %>%
  dplyr::mutate(lvl = ifelse(grepl("sens", basename(path)), "sens", "spec")) %>%
  dplyr::mutate(ret = purrr::map(path, readr::read_csv)) %>%
  dplyr::mutate(median = purrr::map_dbl(ret, function(x){round(median(x$x) * 100, 2)}),
                LCI = purrr::map_dbl(ret, function(x){round(quantile(x$x, 0.025) * 100, 2)}),
                UCI = purrr::map_dbl(ret, function(x){round(quantile(x$x, 0.975) * 100, 2)}),
                rgn_sens_spec = paste0(median, " (", LCI, ", ", UCI, ")")) %>%
  dplyr::select(c("study_id", "lvl", "rgn_sens_spec")) %>%
  tidyr::pivot_wider(., names_from = "lvl", values_from = "rgn_sens_spec") %>%
  magrittr::set_colnames(c("study_id", "rgn_sens", "rgn_spec"))


#......................
# age based model
#......................
datmap %>%
  dplyr::mutate(sero_spec_sens = purrr::map(paths, get_sens_spec)) %>%
  tidyr::unnest(cols = "sero_spec_sens") %>%
  dplyr::mutate(median = round(median * 100, 2),
                LCI = round(LCI * 100, 2),
                UCI = round(UCI * 100, 2)) %>%
  dplyr::mutate(sens_spec = paste0(median, " (", LCI, ", ", UCI, ")"),
                lvl = paste0(lvl, "_", param)) %>%
  dplyr::select("study_id", "lvl", "sens_spec") %>%
  tidyr::pivot_wider(., names_from = "lvl", values_from = "sens_spec") %>%
  dplyr::select(c("study_id", "NoSeroRev_sens", "NoSeroRev_spec", "SeroRev_sens", "SeroRev_spec")) %>%
  dplyr::left_join(., rgnmod, by = "study_id") %>%
  dplyr::select(c("study_id", "rgn_sens", "rgn_spec", dplyr::everything())) %>%
  dplyr::left_join(locatkey, ., by = "study_id") %>%
  dplyr::left_join(., order, by = "study_id") %>%
  dplyr::arrange(order) %>%
  dplyr::select(-c("order", "study_id", "cols")) %>%
  readr::write_tsv(., path = "tables/final_tables/overall_sens_spec_for_all.tsv")


#...........................................................
#---- Posterior Delay Param Table  #----
#...........................................................
datmap %>%
  dplyr::mutate(tbl = purrr::map(modout, get_study_delays_posts)) %>%
  dplyr::left_join(locatkey, .) %>%
  dplyr::left_join(order, .) %>%
  dplyr::arrange(order) %>%
  dplyr::select(c("location", "lvl", "tbl")) %>%
  tidyr::unnest(cols = "tbl") %>%
  dplyr::mutate(lvl = factor(lvl,
                             levels = c("NoSeroRev", "SeroRev"),
                             labels = c("Without Serorev.", "With Serorev."))) %>%
  readr::write_tsv(., path = "tables/final_tables/onset_delay_params_posterior_tbl.tsv")

