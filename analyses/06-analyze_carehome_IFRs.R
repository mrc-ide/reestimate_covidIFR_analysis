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
                                           whichstandard = "pop")) %>%
  tidyr::unnest(cols = "overallIFRret")

carehome_column <- chretmap %>%
  dplyr::mutate(median = round(median * 100, 2),
                LCI = round(LCI * 100, 2),
                UCI = round(UCI * 100, 2)) %>%
  dplyr::mutate(mod_carehome_ifr_col = paste0(median, " (", LCI, ", ", UCI, ")")) %>%
  dplyr::select(c("study_id", "mod_carehome_ifr_col"))


# bring together
readr::write_tsv(carehome_column, path = "tables/final_tables/carehome_excluded_ifr_data.tsv")




#............................................................
#---- Fig of Age Structure versus IFR #----
#...........................................................
# read in observed data
dscdat <- readRDS("results/descriptive_results/descriptive_results_datamap.RDS")
dsc_agedat <- dscdat %>%
  dplyr::filter(breakdown == "ageband") %>%
  dplyr::filter(!grepl("_nch", study_id))

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

prop65 <- dsc_agedat %>%
  dplyr::select(c("study_id", "relpath")) %>%
  dplyr::mutate(prop65 = purrr::map_dbl(relpath, get_65prop)) %>%
  dplyr::select(c("study_id", "prop65")) %>%
  dplyr::filter(!duplicated(.))

#......................
# plotObj
#......................
prop65_modelIFR_plotObj <- chretmap %>%
  dplyr::mutate(median = round(median * 100, 2),
                LCI = round(LCI * 100, 2),
                UCI = round(UCI * 100, 2)) %>%
  dplyr::left_join(., locatkey, by = "study_id") %>%
  dplyr::left_join(., y = prop65, by = "study_id") %>%
  ggplot() +
  geom_pointrange(aes(x = prop65, y = median, ymin = LCI, ymax = UCI,
                      color =  location), alpha = 0.95, size = 1.2) +
  scale_color_manual("Location", values = mycolors) +
  ylab("Care-Home Excluded IFR (95% CrI)") + xlab("Prop. of Population Over 65-Years") +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  theme(plot.margin = unit(c(0.05, 0.05, 0.05, 1),"cm"))


#............................................................
#---- Fig of Attack Rate in Pop over 65#----
#...........................................................
# get seroprevalence over 65
seroprevdat <- dsc_agedat %>%
  dplyr::select(c("study_id", "location", "seroprev_adjdat")) %>%
  tidyr::unnest(cols = "seroprev_adjdat") %>%
  dplyr::group_by(study_id) %>%
  dplyr::filter(seromidpt == max(seromidpt)) %>% # subset to latest serostudy
  dplyr::ungroup()  %>%
  dplyr::mutate(age_mid = purrr::map_dbl(ageband, function(x){
    nums <- as.numeric(stringr::str_split_fixed(x, "-", n = 2))
    nums[nums == 999] <- 100
    return(mean(nums))}))

# subset to studies w/ counts
seroprevdat_65subn <- seroprevdat %>%
  dplyr::filter(!is.na(n_positive) & !is.na(n_tested)) %>%
  dplyr::filter(age_mid >= 65) %>%
  dplyr::group_by(study_id) %>%
  dplyr::summarise(n_65pos = sum(n_positive),
                   n_65test = sum(n_tested)) %>%
  dplyr::mutate(crude_seroprev_obj = purrr::map2(n_65pos, n_65test, .f = function(x,n){ binom.test(x,n) }),
                crude_seroprev_CI = purrr::map(crude_seroprev_obj, "conf.int"),
                seroprev_65_LCI = purrr::map_dbl(crude_seroprev_CI, function(x){x[[1]]}),
                seroprev_65_UCI = purrr::map_dbl(crude_seroprev_CI, function(x){x[[2]]}),
                seroprev_65 = purrr::map_dbl(crude_seroprev_obj, "estimate"))



# subset to studies w/ CIs
seroprevdat_65CIs <- seroprevdat %>%
  dplyr::filter(is.na(n_positive) & is.na(n_tested)) %>%
  dplyr::filter(age_mid >= 65) %>%
  dplyr::group_by(study_id) %>%
  dplyr::summarise(seroprev_65 = mean(seroprev),
                   seroprev_65_LCI  = mean(serolci),
                   seroprev_65_UCI  = mean(serouci))

# bring together
seroprevdat_65 <- dplyr::bind_rows(seroprevdat_65subn, seroprevdat_65CIs)


#......................
# plotObj
#......................
attackrate65_modelIFR_plotObj <- chretmap %>%
  dplyr::left_join(., locatkey, by = "study_id") %>%
  dplyr::left_join(., y = seroprevdat_65, by = "study_id") %>%
  ggplot() +
  geom_pointrange(aes(x = seroprev_65, y = median, ymin = LCI, ymax = UCI,
                      color =  location), alpha = 0.95, size = 1.2) +
  geom_errorbar(aes(x = seroprev_65, y = median, xmin = seroprev_65_LCI, xmax = seroprev_65_UCI,
                    color =  location), alpha = 0.95, size = 1.2) +
  scale_color_manual("Location", values = mycolors) +
  ylab("Overall IFR (95% CrI)") + xlab("Over 65-Years Observed Seroprevalence (95% CI)") +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  theme(plot.margin = unit(c(0.05, 0.05, 0.05, 1),"cm"))


#......................
# bring together
#......................
FigA <- prop65_modelIFR_plotObj + theme(legend.position = "none")
FigB <- attackrate65_modelIFR_plotObj + theme(legend.position = "none")
mainFig <- cowplot::plot_grid(FigA, FigB,
                              align = "h", ncol = 2, nrow = 1,
                              labels = c("(A)", "(B)"))
legend <- cowplot::get_legend(attackrate65_modelIFR_plotObj +
                                guides(color = guide_legend(nrow = 2)) +
                                theme(legend.position = "bottom"))
mainFig <- cowplot::plot_grid(mainFig, legend,
                              nrow = 2, ncol = 1, rel_heights = c(0.9, 0.1))

jpeg("figures/final_figures/Carehome_vs_Over65_years_plots.jpg",
     width = 11, height = 8, units = "in", res = 500)
plot(mainFig)
graphics.off()


