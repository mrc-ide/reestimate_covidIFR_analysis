##................................................................................................
## Purpose: Analyze differences in IFRs
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
regrets <- list.files("results/Modfits_noserorev/", full.names = T)
regretmap <- tibble::tibble(study_id = toupper(stringr::str_split(basename(regrets),
                                                                  "_age", simplify = T)[,1]),
                            sero = "reg",
                            path = regrets)

serorevrets <- list.files("results/ModFits_SeroRev/", full.names = T)
serorevretmap <- tibble::tibble(study_id = toupper(stringr::str_split(basename(serorevrets),
                                                                      "_age", simplify = T)[,1]),
                                sero = "serorev",
                                path = serorevrets)
# bring together
retmap <- dplyr::bind_rows(regretmap, serorevretmap) %>%
  dplyr::mutate(overallIFRret = purrr::map(path,
                                           get_overall_IFRs,
                                           whichstandard = "pop")) %>%
  tidyr::unnest(cols = "overallIFRret")





#............................................................
#---- Table of Overall IFR Standardized by Pop  #----
#...........................................................
#.........................................
# read in observed data
#.........................................
dscdat <- readRDS("results/descriptive_results/descriptive_results_datamap.RDS")

dsc_agedat <- dscdat %>%
  dplyr::filter(breakdown == "ageband") %>%
  dplyr::filter(!grepl("_nch", study_id))

#......................
# get colum death data
#......................
death_col <- dsc_agedat %>%
  dplyr::select(c("study_id", "location", "plotdat")) %>%
  tidyr::unnest(cols = "plotdat") %>%
  dplyr::group_by(study_id) %>%
  dplyr::filter(seromidpt == max(seromidpt)) %>% # subset to latest serostudy
  dplyr::filter(obsday == seromidpt) %>% # subset day of observation to sero midpoint
  dplyr::ungroup() %>%
  dplyr::group_by(study_id) %>% # group across age bands
  dplyr::summarise(totdeaths = sum(cumdeaths))

#......................
# get simple seropred
#......................
simp_seroprevdat <- dsc_agedat %>%
  dplyr::select(c("study_id", "location", "seroprev_adjdat")) %>%
  tidyr::unnest(cols = "seroprev_adjdat") %>%
  dplyr::group_by(study_id) %>%
  dplyr::filter(seromidpt == max(seromidpt)) %>% # subset to latest serostudy
  dplyr::mutate(sero_start = obsdaymin + lubridate::ymd("2020-01-01") - 1,
                sero_end = obsdaymax + lubridate::ymd("2020-01-01") - 1) %>%
  group_by(study_id, sero_start, sero_end) %>%
  dplyr::summarise(n_totpos = sum(n_positive),
                   n_tottest = sum(n_tested)) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate(seroprev = n_totpos/n_tottest)
# manual liftover for studies that only reported CIs
ci_simp_seroprevdat <- dsc_agedat %>%
  dplyr::select(c("study_id", "location", "seroprev_adjdat")) %>%
  tidyr::unnest(cols = "seroprev_adjdat") %>%
  dplyr::group_by(study_id) %>%
  dplyr::filter(seromidpt == max(seromidpt)) %>% # subset to latest serostudy
  dplyr::filter(is.na(n_positive) & is.na(n_tested)) %>%
  dplyr::summarise(seroprev_ci = mean(seroprev))

# fix and take to simple column
simp_seroprevdat <- simp_seroprevdat %>%
  dplyr::left_join(., ci_simp_seroprevdat, by = "study_id") %>%
  dplyr::mutate(seroprev = ifelse(is.na(seroprev), seroprev_ci, seroprev))

seroprev_column <- simp_seroprevdat %>%
  dplyr::mutate(seroprev = round(seroprev * 100, 2),
                serocol = paste0(seroprev, "%", " (", sero_start, " - ", sero_end, ")"),
                serocol = gsub("2020-", "", serocol),
                serocol = str_replace_all(serocol,
                                          c("01-" = "Jan. ",
                                            "02-" = "Feb. ",
                                            "03-" = "Mar. ",
                                            "04-" = "Apr. ",
                                            "05-" = "May ",
                                            "06-" = "Jun. ",
                                            "07-" = "Jul. ",
                                            "08-" = "Aug. "))) %>%
  dplyr::select(c("study_id", "serocol"))

#......................
# crude IFRs
#......................
source("R/monte_carlo_cis.R")
# calculate CIs for binomial
crude_IFRs_binomial <- dsc_agedat %>%
  dplyr::select(c("study_id", "plotdat")) %>%
  dplyr::filter(!study_id %in% c("ITA1", "SWE1", "DNK1")) %>%
  tidyr::unnest(cols = plotdat) %>%
  dplyr::group_by(study_id) %>%
  dplyr::filter(seromidpt == max(seromidpt)) %>% # latest serostudy
  dplyr::filter(obsday == seromidpt) %>% # latest serostudy
  dplyr::summarise(cumdeaths = sum(cumdeaths), # sum over agebands
                   popn = sum(popn),
                   n_positive = sum(n_positive),
                   n_tested = sum(n_tested)) %>%
  dplyr::select(c("study_id", "cumdeaths", "popn", "n_positive", "n_tested")) %>%
  dplyr::group_by(study_id) %>%
  dplyr::mutate(seroprev = n_positive/n_tested,
                ifr_range = purrr::map(cumdeaths, get_binomial_monte_carlo_cis, popN = popn,
                                       npos = n_positive, ntest = n_tested, iters = 1e5),
                crudeIFR = cumdeaths/((seroprev * popn) + cumdeaths),
                lower_ci = purrr::map_dbl(ifr_range, quantile, 0.025),
                upper_ci = purrr::map_dbl(ifr_range, quantile, 0.975)) %>%
  dplyr::select(c("study_id", "seroprev", "crudeIFR", "lower_ci", "upper_ci")) %>%
  dplyr::ungroup(.)

# calculate CIs for logit
crude_IFRs_logit <- dsc_agedat %>%
  dplyr::select(c("study_id", "plotdat")) %>%
  dplyr::filter(study_id %in% c("ITA1", "SWE1", "DNK1")) %>%
  tidyr::unnest(cols = plotdat) %>%
  dplyr::group_by(study_id) %>%
  dplyr::filter(seromidpt == max(seromidpt)) %>% # latest serostudy
  dplyr::filter(obsday == seromidpt) %>% # latest serostudy
  dplyr::summarise(cumdeaths = sum(cumdeaths), # sum over agebands
                   popn = sum(popn),
                   seroprev = mean(seroprev),
                   serolci = mean(serolci),
                   serouci = mean(serouci)) %>%
  dplyr::select(c("study_id", "cumdeaths", "popn", "seroprev",  "serolci", "serouci")) %>%
  dplyr::group_by(study_id) %>%
  dplyr::mutate(SE = (COVIDCurve:::logit(serouci) - COVIDCurve:::logit(serolci))/(1.96 * 2))  %>%
  dplyr::mutate(ifr_range = purrr::map(cumdeaths, get_normal_monte_carlo_cis, popN = popn,
                                       mu = seroprev, sigma = SE, iters = 1e5),
                crudeIFR = cumdeaths/((seroprev * popn) + cumdeaths),
                lower_ci = purrr::map_dbl(ifr_range, quantile, 0.025),
                upper_ci = purrr::map_dbl(ifr_range, quantile, 0.975)) %>%
  dplyr::select(c("study_id", "seroprev", "crudeIFR", "lower_ci", "upper_ci"))

# out
crudeIFRs_CI <- dplyr::bind_rows(crude_IFRs_binomial, crude_IFRs_logit)
crude_IFR_column <- crudeIFRs_CI %>%
  dplyr::mutate(crudeIFR = round(crudeIFR * 100, 2),
                lower_ci = round(lower_ci * 100, 2),
                upper_ci = round(upper_ci * 100, 2)) %>%
  dplyr::mutate(crude_ifr_col = paste0(crudeIFR, " (", lower_ci, ", ", upper_ci, ")")) %>%
  dplyr::select(c("study_id", "crude_ifr_col"))


#......................
# Sens Spec
#......................
no_serorev_sensspec <- retmap %>%
  dplyr::select(c("study_id", "sero", "path")) %>%
  dplyr::mutate(sero_spec_sens = purrr::map(path, get_sens_spec)) %>%
  tidyr::unnest(cols = "sero_spec_sens") %>%
  dplyr::filter(sero == "reg") %>%
  dplyr::mutate(median = round(median * 100, 2),
                LCI = round(LCI * 100, 2),
                UCI = round(UCI * 100, 2)) %>%
  dplyr::mutate(mod_no_serorev_sero_spec_sens_col = paste0(median, " (", LCI, ", ", UCI, ")")) %>%
  dplyr::select("study_id", "param", "mod_no_serorev_sero_spec_sens_col") %>%
  tidyr::pivot_wider(., names_from = "param", values_from = "mod_no_serorev_sero_spec_sens_col") %>%
  dplyr::select(c("study_id", "sens", "spec"))

#......................
# Modeled IFRs
#......................
no_serorev_ifrs <- retmap %>%
  dplyr::filter(sero == "reg") %>%
  dplyr::mutate(median = round(median * 100, 2),
                LCI = round(LCI * 100, 2),
                UCI = round(UCI * 100, 2)) %>%
  dplyr::mutate(mod_no_serorev_ifr_col = paste0(median, " (", LCI, ", ", UCI, ")")) %>%
  dplyr::select(c("study_id", "mod_no_serorev_ifr_col"))

serorev_ifrs <- retmap %>%
  dplyr::filter(sero == "serorev") %>%
  dplyr::mutate(median = round(median * 100, 2),
                LCI = round(LCI * 100, 2),
                UCI = round(UCI * 100, 2)) %>%
  dplyr::mutate(mod_serorev_ifr_col = paste0(median, " (", LCI, ", ", UCI, ")")) %>%
  dplyr::select(c("study_id", "mod_serorev_ifr_col"))


#......................
# bring together
#......................
dir.create("tables/final_tables/", recursive = TRUE)
dplyr::left_join(death_col, seroprev_column, by = "study_id") %>%
  dplyr::left_join(., crude_IFR_column, by = "study_id") %>%
  dplyr::left_join(., no_serorev_sensspec, by = "study_id") %>%
  dplyr::left_join(., no_serorev_ifrs, by = "study_id") %>%
  dplyr::left_join(., serorev_ifrs, by = "study_id") %>%
  dplyr::left_join(., order, by = "study_id") %>%
  dplyr::arrange(order) %>%
  dplyr::select(-c("order")) %>%
  dplyr::mutate(totdeaths = prettyNum(totdeaths, big.mark=",", scientific=FALSE)) %>%
  readr::write_tsv(., path = "tables/final_tables/overall_ifr_data_popstandardized.tsv")


#............................................................
#---- Supp Table of Overall IFR Standardized by AR and Pop  #----
#...........................................................
# run new standardization
alt_stand_retmap <- dplyr::bind_rows(regretmap, serorevretmap) %>%
  dplyr::mutate(overallIFRret = purrr::map(path,
                                           get_overall_IFRs,
                                           whichstandard = "arpop")) %>%
  tidyr::unnest(cols = "overallIFRret")

#......................
# Modeled IFRs
#......................
alt_no_serorev_ifrs <- alt_stand_retmap %>%
  dplyr::filter(sero == "reg") %>%
  dplyr::mutate(median = round(median * 100, 2),
                LCI = round(LCI * 100, 2),
                UCI = round(UCI * 100, 2)) %>%
  dplyr::mutate(alt_mod_no_serorev_ifr_col = paste0(median, " (", LCI, ", ", UCI, ")")) %>%
  dplyr::select(c("study_id", "alt_mod_no_serorev_ifr_col"))

alt_serorev_ifrs <- alt_stand_retmap %>%
  dplyr::filter(sero == "serorev") %>%
  dplyr::mutate(median = round(median * 100, 2),
                LCI = round(LCI * 100, 2),
                UCI = round(UCI * 100, 2)) %>%
  dplyr::mutate(alt_mod_serorev_ifr_col = paste0(median, " (", LCI, ", ", UCI, ")")) %>%
  dplyr::select(c("study_id", "alt_mod_serorev_ifr_col"))

# alternative overall IFR out
dplyr::left_join(crude_IFR_column, no_serorev_ifrs, by = "study_id") %>%
  dplyr::left_join(., serorev_ifrs, by = "study_id") %>%
  dplyr::left_join(., alt_no_serorev_ifrs) %>%
  dplyr::left_join(., alt_serorev_ifrs, by = "study_id") %>%
  dplyr::left_join(., order, by = "study_id") %>%
  dplyr::arrange(order) %>%
  dplyr::select(-c("order")) %>%
  readr::write_tsv(., path = "tables/final_tables/overall_ifr_data_attackrate_pop_standardized.tsv")




#............................................................
#---- Fig of IFR-Death Comparisons #----
#...........................................................

#......................
# crude IFRs
#......................
crude_IFRs <- dsc_agedat %>%
  dplyr::select(c("study_id", "location", "plotdat")) %>%
  tidyr::unnest(cols = "plotdat") %>%
  dplyr::group_by(study_id, location) %>%
  dplyr::filter(seromidpt == max(seromidpt)) %>% # latest serostudy
  dplyr::filter(obsday == seromidpt) %>%  # sero obs day
  dplyr::summarise(cumdeaths = sum(cumdeaths),
                   popn = sum(popn)) %>% # sum across agebands
  dplyr::left_join(., simp_seroprevdat, by = "study_id") %>%
  dplyr::group_by(study_id, location) %>%
  dplyr::mutate(crude = cumdeaths  / (seroprev * popn + cumdeaths),
                crude = crude * 100) %>%
  dplyr::ungroup() %>%
  dplyr::select(c("study_id", "location", "crude"))


#......................
# per capita deaths
#......................
percap_deaths <- dsc_agedat %>%
  dplyr::select(c("study_id", "location", "data")) %>%
  group_by(study_id, location) %>%
  dplyr::mutate(deaths = purrr::map(data, "deaths_TSMCMC"),
                deaths = purrr::map_dbl(deaths, function(x){sum(x$deaths)}),
                popN = purrr::map(data, "prop_pop"),
                popN = purrr::map_dbl(popN, function(x){sum(x$popN)})) %>%
  tidyr::unnest(cols = c("deaths", "popN")) %>%
  dplyr::summarise(cumdeaths = sum(deaths),
                   popN = sum(popN)) %>%
  dplyr::mutate(capita = cumdeaths/popN) %>%
  dplyr::select(c("study_id", "location", "capita"))

#......................
# bring together others deaths
#......................
extra_deaths <- percap_deaths %>%
  dplyr::mutate(capita = capita * 100) %>%
  dplyr::left_join(., crude_IFRs, by = c("study_id", "location")) %>%
  #dplyr::left_join(., y = excess_deaths, by = "study_id") %>%
  #dplyr::left_join(., y = carehome_deaths, by = "study_id") %>%
  tidyr::pivot_longer(., cols = -c("study_id", "location"), names_to = "param", values_to = "est") %>%
  dplyr::mutate(param = factor(param,
                               levels = c("crude", "capita", "excess", "carehome"),
                               labels = c("Crude", "Per Capita", "Excess", "Care-Home")
  ))
#......................
# modelled IFRs
#......................
studylocats <- crude_IFRs %>%
  dplyr::select(c("study_id", "location"))
locatlvls <- rev(sort(unique(crude_IFRs$location)))
studylocats$location <- factor(studylocats$location, levels = locatlvls)
extra_deaths$location <- factor(extra_deaths$location, levels = locatlvls)

# modelled
modelled_IFRs <- retmap %>%
  dplyr::left_join(studylocats, ., by = "study_id") %>%
  dplyr::mutate(median = median * 100,
                LCI = LCI * 100,
                UCI = UCI * 100,
                upperbound_inf = ifelse(UCI >= 3 & sero == "reg", 3, NA),
                upperbound_rev = ifelse(UCI >= 3 & sero == "serorev", 3, NA),
                UCI = ifelse(UCI >= 3, 3, UCI),
                sero = factor(sero, levels = c("reg", "serorev"),
                              labels = c("Without Serorev.", "With Serorev.")))

#......................
# plot out
#......................
death_type_plot <- ggplot() +
  geom_point(data = extra_deaths, aes(x = location, y = est, fill = param, shape = param),
             size = 6.5, alpha = 0.8) +
  geom_pointrange(data = modelled_IFRs,
                  aes(x = location, ymin = LCI, y = median, ymax = UCI, group = sero, color = sero),
                  size = 0.75, alpha = 0.5, position = position_dodge(width = 0.5)) +
  geom_point(data = modelled_IFRs,
             aes(x = location, y = upperbound_inf),
             color = "#4285F4", size = 3, alpha = 0.9, shape = 3,
             position = position_nudge(x = -0.125)) +
  geom_point(data = modelled_IFRs,
             aes(x = location, y = upperbound_rev),
             color = "#EA4335", size = 3, alpha = 0.9, shape = 3,
             position = position_nudge(x = 0.125)) +
  scale_color_manual("", values = c("#4285F4", "#EA4335")) +
  scale_shape_manual("", values = c(23, 22, 24, 25)) +
  scale_fill_manual("", values = c(wesanderson::wes_palette("IsleofDogs2", type = "discrete"))) +
  ylab("Infection Fatality Ratio (%)") +
  coord_flip() +
  theme(
    plot.title =  element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(family = "Helvetica", face = "bold", vjust = 0.5, hjust = 0.5, size = 16),
    axis.text.x = element_text(family = "Helvetica", color = "#000000", vjust = 0.5, hjust = 0.5, size = 14),
    axis.text.y = element_text(family = "Helvetica", color = "#000000", face = "bold", vjust = 0.5, hjust = 1, size = 14),
    legend.title = element_blank(),
    legend.text = element_text(family = "Helvetica", vjust = 0.8, hjust = 0.5, size = 14, angle = 0),
    legend.position = "right",
    legend.key = element_blank(),
    axis.line.x = element_line(color = "black", size = 1.5),
    axis.line.y = element_line(color = "black", size = 1.5),
    axis.ticks.y = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent"),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )
# out
jpeg("figures/final_figures/Figure_IFR_ranges.jpg", width = 10, height = 7, units = "in", res = 500)
plot(death_type_plot)
graphics.off()




#............................................................
#---- Fig of Age Structure versus IFR #----
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

prop65 <- dsc_agedat %>%
  dplyr::select(c("study_id", "relpath")) %>%
  dplyr::mutate(prop65 = purrr::map_dbl(relpath, get_65prop)) %>%
  dplyr::select(c("study_id", "prop65")) %>%
  dplyr::filter(!duplicated(.))

#......................
# plotObj
#......................
prop65_modelIFR_plotObj <- retmap %>%
  dplyr::left_join(., locatkey, by = "study_id") %>%
  dplyr::left_join(., y = prop65, by = "study_id") %>%
  dplyr::filter(sero == "reg") %>% # only no serorev models
  ggplot() +
  geom_pointrange(aes(x = prop65, y = median, ymin = LCI, ymax = UCI,
                      color =  location), alpha = 0.95, size = 1.2) +
  scale_color_manual("Location", values = mycolors) +
  ylab("Overall IFR (95% CrI)") + xlab("Prop. of Population Over 65-Years") +
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
attackrate65_modelIFR_plotObj <- retmap %>%
  dplyr::left_join(., locatkey, by = "study_id") %>%
  dplyr::left_join(., y = seroprevdat_65, by = "study_id") %>%
  dplyr::filter(sero == "reg") %>% # only no serorev models
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

jpeg("figures/final_figures/Over65_years_plots.jpg",
     width = 11, height = 8, units = "in", res = 500)
plot(mainFig)
graphics.off()


