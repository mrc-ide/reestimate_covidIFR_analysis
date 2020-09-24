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
# read in data
#...........................................................
dscdat <- readRDS("results/descriptive_results/descriptive_results_datamap.RDS")

dsc_agedat <- dscdat %>%
  dplyr::filter(breakdown == "ageband") %>%
  dplyr::filter(!grepl("_nch", study_id))

#......................
# results
#......................
regrets <- list.files("results/Modfits_noserorev/", full.names = T)
regretmap <- tibble::tibble(study_id = toupper(stringr::str_split(basename(regrets), "_age", simplify = T)[,1]),
                            sero = "reg",
                            path = regrets)

serorevrets <- list.files("results/ModFits_SeroRev/", full.names = T)
serorevretmap <- tibble::tibble(study_id = toupper(stringr::str_split(basename(serorevrets), "_age", simplify = T)[,1]),
                                sero = "serorev",
                                path = serorevrets)
# bring together
retmap <- dplyr::bind_rows(regretmap, serorevretmap) %>%
  dplyr::mutate(overallIFRret = purrr::map(path, get_overall_IFRs)) %>%
  tidyr::unnest(cols = "overallIFRret")



#............................................................
#---- Table of Overall IFR  #----
#...........................................................
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
  group_by(study_id, seromidpt) %>%
  dplyr::summarise(n_totpos = sum(n_positive),
                   n_tottest = sum(n_tested)) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate(sero_midday = seromidpt + lubridate::ymd("2020-01-01"),
                seroprev = n_totpos/n_tottest)
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
                serocol = paste0(seroprev, "%", " (", sero_midday, ")")) %>%
  dplyr::select(c("study_id", "serocol"))

#......................
# crude IFRs
#......................
# get standard errors
# Delta method needs standard error of seroprev SE(p)
# where SE(p) is the standard error of the binomial proportion for all studies except ITA
# for ITA, SWE, DNK we use SE(p) as the logit transformed SE gleaned from the provided CIs

SE_subn <- dsc_agedat %>%
  dplyr::select(c("study_id", "location", "seroprev_adjdat")) %>%
  tidyr::unnest(cols = "seroprev_adjdat") %>%
  dplyr::group_by(study_id) %>%
  dplyr::filter(seromidpt == max(seromidpt)) %>% # latest serostudy
  dplyr::summarise(n_positive = sum(n_positive),
                   n_tested = sum(n_tested)) %>%
  dplyr::filter(!is.na(n_positive) & !is.na(n_tested)) %>%
  dplyr::mutate(seroprev = n_positive/n_tested,
                binom_se = sqrt(seroprev * (1-seroprev))) %>%
  dplyr::select(c("study_id", "binom_se"))

SE_cis <- dsc_agedat %>%
  dplyr::select(c("study_id", "location", "seroprev_adjdat")) %>%
  tidyr::unnest(cols = "seroprev_adjdat") %>%
  dplyr::group_by(study_id) %>%
  dplyr::filter(seromidpt == max(seromidpt)) %>%
  dplyr::filter(study_id %in% c("ITA1", "DNK1", "SWE1")) %>%
  dplyr::summarise(serolci = mean(serolci),
                   serouci = mean(serouci)) %>%
  dplyr::mutate(binom_se = (COVIDCurve:::logit(serouci) - COVIDCurve:::logit(serolci))/(1.96 * 2))  %>%
  dplyr::select(c("study_id", "binom_se"))

deltaSE <- dplyr::bind_rows(SE_subn, SE_cis)


#......................
# subset regions to parts we need
# and perform calculation
#......................
delta_IFR <- dsc_agedat %>%
  dplyr::select(c("study_id", "location", "plotdat")) %>%
  tidyr::unnest(cols = "plotdat") %>%
  dplyr::group_by(study_id) %>%
  dplyr::filter(seromidpt == max(seromidpt)) %>% # latest serostudy
  dplyr::filter(obsday == seromidpt) %>%  # sero obs day
  dplyr::summarise(cumdeaths = sum(cumdeaths),
                   popn = sum(popn)) %>% # sum across agebands
  dplyr::left_join(., simp_seroprevdat, by = "study_id") %>%
  dplyr::select(c("study_id", "cumdeaths", "popn", "seroprev")) %>%
  dplyr::left_join(., deltaSE, by = "study_id") %>%
  dplyr::group_by(study_id) %>%
  dplyr::mutate(IFRcalc = cumdeaths  / (seroprev * popn + cumdeaths),
                IFRbound = purrr::map(seroprev, get_delta_CI_vals, deaths = cumdeaths, popN = popn, SE = binom_se, tol = 1e-4),
                lower_ci = purrr::map_dbl(IFRbound, "lower.ci"),
                upper_ci = purrr::map_dbl(IFRbound, "upper.ci"))


delta_IFR_column <- delta_IFR %>%
  dplyr::mutate(IFRcalc = round(IFRcalc * 100, 2),
                lower_ci = round(lower_ci * 100, 2),
                upper_ci = round(upper_ci * 100, 2)) %>%
  dplyr::mutate(delta_ifr_col = paste0(IFRcalc, " (", lower_ci, ", ", upper_ci, ")")) %>%
  dplyr::select(c("study_id", "delta_ifr_col"))


#......................
# Modeled IFRs
#......................
no_serorev <- retmap %>%
  dplyr::filter(sero == "reg") %>%
  dplyr::mutate(median = round(median * 100, 2),
                LCI = round(LCI * 100, 2),
                UCI = round(UCI * 100, 2)) %>%
  dplyr::mutate(mod_no_serorev_ifr_col = paste0(median, " (", LCI, ", ", UCI, ")")) %>%
  dplyr::select(c("study_id", "mod_no_serorev_ifr_col"))

serorev <- retmap %>%
  dplyr::filter(sero == "serorev") %>%
  dplyr::mutate(median = round(median * 100, 2),
                LCI = round(LCI * 100, 2),
                UCI = round(UCI * 100, 2)) %>%
  dplyr::mutate(mod_serorev_ifr_col = paste0(median, " (", LCI, ", ", UCI, ")")) %>%
  dplyr::select(c("study_id", "mod_serorev_ifr_col"))


#......................
# bring together
#......................
dplyr::left_join(death_col, seroprev_column, by = "study_id") %>%
  dplyr::left_join(., delta_IFR_column, by = "study_id") %>%
  dplyr::left_join(., no_serorev, by = "study_id") %>%
  dplyr::left_join(., serorev, by = "study_id") %>%
  dplyr::mutate(totdeaths = prettyNum(totdeaths, big.mark=",", scientific=FALSE)) %>%
  readr::write_tsv(., path = "tables/final_tables/overall_ifr_data.tsv")




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
                UCI = ifelse(UCI >= 3, 3, UCI))

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
  scale_color_manual("Sero. Type", values = c("#4285F4", "#EA4335")) +
  scale_shape_manual("Death Type", values = c(23, 22, 24, 25)) +
  scale_fill_manual("Death Type", values = c(wesanderson::wes_palette("IsleofDogs2", type = "discrete"))) +
  ylab("Infection Fatality Ratio (%)") +
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
  dplyr::left_join(., y = prop65, by = "study_id") %>%
  dplyr::filter(sero == "reg") %>% # only no serorev models
  ggplot() +
  geom_pointrange(aes(x = prop65, y = median, ymin = LCI, ymax = UCI,
                      color =  study_id), alpha = 0.95) +
  scale_color_manual("Study ID", values = mycolors) +
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
  dplyr::left_join(., y = seroprevdat_65, by = "study_id") %>%
  dplyr::filter(sero == "reg") %>% # only no serorev models
  ggplot() +
  geom_pointrange(aes(x = seroprev_65, y = median, ymin = LCI, ymax = UCI,
                      color =  study_id), alpha = 0.95) +
  geom_errorbar(aes(x = seroprev_65, y = median, xmin = seroprev_65_LCI, xmax = seroprev_65_UCI,
                    color =  study_id), alpha = 0.95) +
  scale_color_manual("Study ID", values = mycolors) +
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
                                theme(legend.position = "bottom"))
mainFig <- cowplot::plot_grid(mainFig, legend,
                              nrow = 2, ncol = 1, rel_heights = c(0.9, 0.1))

jpeg("figures/final_figures/Over65_years_plots.jpg",
     width = 11, height = 8, units = "in", res = 500)
plot(mainFig)
graphics.off()



#............................................................
#---- Fig of Heath Care Capacity beds  #----
# get proxy for healthcare capacity through OWID
#...........................................................
# our world in data for national hospital beds
# https://github.com/owid/covid-19-data/tree/master/public/data
# https://github.com/owid/covid-19-data/blob/master/public/data/owid-covid-codebook.csv
nat_hosp_beds <- readr::read_csv("data/raw/owid-covid-data.csv") %>%
  dplyr::filter(iso_code %in% c("BRA", "DNK", "ESP", "ITA", "KEN", "LUX", "NLD", "SWE")) %>% # national studies
  dplyr::select(c("iso_code", "hospital_beds_per_thousand")) %>%
  dplyr::mutate(hospital_beds = hospital_beds_per_thousand * 1e3) %>%
  dplyr::filter(!duplicated(.))

# get data for deaths as well
# replicated from above
simpdeaths <- dsc_agedat %>%
  dplyr::select(c("study_id", "location", "plotdat")) %>%
  tidyr::unnest(cols = "plotdat") %>%
  dplyr::group_by(study_id) %>%
  dplyr::filter(seromidpt == max(seromidpt)) %>% # subset to latest serostudy
  dplyr::filter(obsday == seromidpt) %>% # subset day of observation to sero midpoint
  dplyr::ungroup() %>%
  dplyr::group_by(study_id) %>% # group across age bands
  dplyr::summarise(totdeaths = sum(cumdeaths))

hlth_care_plotdat <- retmap %>%
  dplyr::mutate(iso_code = stringr::str_extract_all(study_id, "[A-Z]+", simplify = T)[,1]) %>%
  dplyr::left_join(., nat_hosp_beds, by = "iso_code") %>%
  dplyr::filter(sero == "reg") %>% # only no serorev models
  dplyr::left_join(., simpdeaths, by = "study_id") %>%
  dplyr::mutate(death_bed_ratio = totdeaths/hospital_beds)

#......................
#plotObj
#......................
health_cap_modIFR_plotObj <- hlth_care_plotdat %>%
  ggplot() +
  geom_pointrange(aes(x = death_bed_ratio, y = median,
                      ymin = LCI, ymax = UCI,
                      color = study_id)) +
  scale_color_manual("Study ID", values = mycolors) +
  xlab("Ratio of Cumulative Deaths to Hospital Beds") + ylab("IFR (95% CrI)") +
  xyaxis_plot_theme +
  theme(legend.position = "bottom")

# out
jpeg("figures/final_figures/hosp_bed_versus_ifr.jpg",
     width = 10, height = 7, units = "in", res = 500)
plot(health_cap_modIFR_plotObj)
graphics.off()




