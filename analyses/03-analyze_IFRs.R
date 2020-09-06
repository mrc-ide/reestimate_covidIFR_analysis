##................................................................................................
## Purpose:
##
## Notes:
##................................................................................................
library(tidyverse)
library(COVIDCurve)
source("R/my_themes.R")
source("R/covidcurve_helper_functions.R")

#......................
# read in data
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
tbl2 %>%
  dplyr::mutate_if(is.numeric, round, 2) %>%
  readr::write_csv(., path = "tables/final_tables/overall_IFRs.csv")

#............................................................
#---- Age Specific Results #----
#...........................................................
colormap <- readr::read_csv("data/plot_aesthetics/color_studyid_map.csv")
retmap_age <- retmap %>%
  tidyr::unnest(cols = "strataIFRret") %>%
  dplyr::filter(lvl == "Age-Band") %>%
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


