##................................................................................................
## Purpose: Analyze care home specific IFRs
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
#---- Read in Fitted Care Home Data #----
#...........................................................
chrets <- list.files("results/Modfits_carehomes/", full.names = T)
chretmap <- tibble::tibble(study_id = toupper(stringr::str_split(basename(chrets), "_carehomes", simplify = T)[,1]),
                           path = chrets) %>%
  dplyr::mutate(overallIFRret = purrr::map(path, get_overall_IFRs,
                                           whichstandard = "arpop")) %>%
  tidyr::unnest(cols = "overallIFRret")

carehome_column <- chretmap %>%
  dplyr::mutate(median = round(median * 100, 2),
                LCI = round(LCI * 100, 2),
                UCI = round(UCI * 100, 2)) %>%
  dplyr::mutate(mod_carehome_ifr_col = paste0(median, " (", LCI, ", ", UCI, ")")) %>%
  dplyr::select(c("study_id", "mod_carehome_ifr_col"))


# bring together
readr::write_tsv(carehome_column, path = "tables/final_tables/carehome_excluded_ifr_data.tsv")



