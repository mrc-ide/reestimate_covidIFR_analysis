##................................................................................................
## Purpose:
##
## Notes:
##................................................................................................
library(tidyverse)
library(COVIDCurve)
source("R/my_themes.R")
# quick functions
get_strata_IFRs <- function(path) {
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

get_overall_IFRs <- function(path) {
  modout <- readRDS(path)
  out <- COVIDCurve::get_globalIFR_cred_intervals(IFRmodel_inf = modout,
                                                  whichrung = "rung1",
                                                  by_chain = FALSE)
  return(out)
}




#......................
# read in data
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
  dplyr::mutate(strataIFRret = purrr::map(path, get_strata_IFRs),
                overallIFRret = purrr::map(path, get_overall_IFRs))


#............................................................
#---- Overall IFR Results #----
#...........................................................
tbl2 <- retmap %>%
  dplyr::select(c("study_id", "sero", "overallIFRret")) %>%
  tidyr::unnest(cols = overallIFRret) %>%
  dplyr::select(-c(dplyr::starts_with("Geweke"))) %>%
  dplyr::mutate(sero = factor(sero, levels = c("reg", "serorev"), labels = c("SeroInf", "SeroRev."))) %>%
  dplyr::arrange(study_id) %>%
  magrittr::set_colnames(
    c("Study ID", "Seroreversion", "Min.", "LCI", "Med.", "Mean", "UCI", "Max", "ESS"))

dir.create("tables/final_tables", recursive = T)
readr::write_csv(tbl2, path = "tables/final_tables/overall_IFRs.csv")

#............................................................
#---- Age Specific Results #----
#...........................................................
colormap <- readr::read_csv("data/plot_aesthetics/color_studyid_map.csv")
retmap_age <- retmap %>%
  tidyr::unnest(cols = "strataIFRret") %>%
  dplyr::filter(lvl == "ageband") %>%
  dplyr::left_join(., colormap, by = "study_id") %>%
  dplyr::mutate(age_mid = purrr::map_dbl(strata, function(x){
    nums <- as.numeric(stringr::str_split_fixed(x, "-", n = 2))
    nums[nums == 999] <- 100
    return(mean(nums))})) %>%
  dplyr::mutate(sero = factor(sero, levels = c("reg", "serorev"), labels = c("SeroInf", "SeroRev")))

ageIFR_plotObj <- ggplot() +
  geom_pointrange(data = retmap_age, aes(x = age_mid, y = median, ymin = LCI, ymax = UCI, color =  study_id), alpha = 0.5) +
  #geom_line(data = retmap_age, aes(x = age_mid, y = median, group = study_id, color =  study_id),
  #          alpha = 0.4, size = 0.5, show.legend = FALSE) +
  facet_wrap(.~sero) +
  scale_color_manual("Study ID", values = unique(retmap_age$cols)) +
  ylab("Age-Specific IFR (95% Cred. Int.)") + xlab("Mid. Age") +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


