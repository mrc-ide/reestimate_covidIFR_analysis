##................................................................................................
## Purpose:
##
## Notes:
##................................................................................................
library(tidyverse)
library(COVIDCurve)
source("R/my_themes.R")
source("R/covidcurve_helper_functions.R")

# colors
study_cols <- readr::read_csv("data/plot_aesthetics/color_studyid_map.csv")
mycolors <- study_cols$cols
names(mycolors) <- study_cols$study_id
studyidnames <- study_cols %>%
  dplyr::select(c("study_id", "names")) %>%
  dplyr::filter(!is.na(names))

#......................
# read in fitted data
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


#......................
# read in crude data
#......................
dscdat <- readRDS("results/descriptive_results/descriptive_results_datamap.RDS") %>%
  dplyr::mutate(drop = ifelse(study_id %in% c("DNK1", "GBR3", "ESP1-2", "BRA1", "CHE1", "ITA1")
                              & breakdown == "region", T, # drop regions that aren't basic
                              ifelse(grepl("_nch", study_id), T, F)) # drop carehomes
  ) %>%
  dplyr::filter(!drop) %>%
  dplyr::mutate(prop_pop = purrr::map(data, "prop_pop")) %>%
  dplyr::select(c("study_id", "prop_pop", "plotdat"))

# proportion of population
prop_pop_dscdat <- dscdat %>%
  dplyr::select(c("study_id", "prop_pop")) %>%
  tidyr::unnest(cols = "prop_pop") %>%
  dplyr::select(-c("region"))
# data we need
plotdat_dscdat <- dscdat %>%
  dplyr::select(c("study_id", "plotdat")) %>%
  tidyr::unnest(cols = "plotdat") %>%
  dplyr::select(-c("region"))
# get descriptive IFR
dsc_overall_ifr <- dplyr::left_join(plotdat_dscdat, prop_pop_dscdat, by = c("study_id", "ageband")) %>%
  dplyr::mutate(infxns = popn * seroprev,
                crudeIFR =  cumdeaths/(infxns + cumdeaths),
                crudeIFR = ifelse(infxns == 0, NA, crudeIFR)) %>% # if infxns are 0, can't calculate
  dplyr::group_by(study_id, ageband) %>%
  dplyr::filter(obsday == seromidpt) %>% # want serosurvey midpoint
  dplyr::filter(seromidpt == max(seromidpt)) %>%  # want latest date
  dplyr::ungroup(.) %>%
  dplyr::mutate(pop_prop = ifelse(study_id %in% c(c("BRA4", "BRA5", "LA_CA1", "SF_CA1")),
                                  1, pop_prop)) %>%
  dplyr::group_by(study_id) %>%
  dplyr::summarise(crudeIFR = sum(crudeIFR * pop_prop)) %>%
  dplyr::select(c("study_id", "crudeIFR")) %>%
  dplyr::mutate(crudeIFR = crudeIFR * 100)

# descriptive age-based IFR
dscdat_age <- dscdat %>%
  dplyr::select(c("study_id", "plotdat")) %>%
  tidyr::unnest(cols = "plotdat") %>%
  dplyr::mutate(infxns = popn * seroprev,
                crudeIFR =  cumdeaths/(infxns + cumdeaths),
                crudeIFR = ifelse(infxns == 0, NA, crudeIFR)) %>% # if infxns are 0, can't calculate
  dplyr::group_by(study_id, ageband) %>%
  dplyr::filter(obsday == seromidpt) %>% # want serosurvey midpoint
  dplyr::filter(seromidpt == max(seromidpt)) %>%  # want latest date
  dplyr::ungroup(.)  %>%
  dplyr::select(c("study_id", "age_mid", "crudeIFR")) %>%
  dplyr::mutate(crudeIFR = crudeIFR * 100)


#............................................................
#---- Overall IFR Results #----
#...........................................................
# fitted ifrs
tbl2 <- retmap %>%
  dplyr::select(c("study_id", "sero", "overallIFRret")) %>%
  tidyr::unnest(cols = overallIFRret) %>%
  dplyr::select(-c(dplyr::starts_with("Geweke"))) %>%
  dplyr::mutate(sero = factor(sero, levels = c("reg", "serorev"), labels = c("SeroInf.", "SeroRev."))) %>%
  dplyr::arrange(study_id)

pretty_tbl2 <- tbl2 %>%
  dplyr::mutate(median = median * 100,
                LCI = LCI * 100,
                UCI = UCI * 100) %>%
  dplyr::mutate_if(is.numeric, round, 2) %>%
  dplyr::mutate(IFR = paste0(median, " ", "(", LCI, ", ", UCI, ")")) %>%
  dplyr::select(c("study_id", "sero", "IFR")) %>%
  tidyr::pivot_wider(., names_from = "sero", values_from = "IFR")


# combine crude w/ fitted
pretty_tbl2 <- pretty_tbl2 %>%
  dplyr::left_join(., dsc_overall_ifr, by = "study_id") %>%
  dplyr::left_join(., studyidnames, by = "study_id") %>%
  dplyr::select(c("names", "crudeIFR", "SeroInf.", "SeroRev."))



dir.create("tables/final_tables", recursive = T)
pretty_tbl2 %>%
  dplyr::mutate_if(is.numeric, round, 2) %>%
  readr::write_csv(., path = "tables/final_tables/overall_IFRs.csv")

#............................................................
#---- Age Specific Results #----
#...........................................................
retmap_age <- retmap %>%
  tidyr::unnest(cols = "strataIFRret") %>%
  dplyr::filter(lvl == "Age-Band") %>%
  dplyr::left_join(., studyidnames, by = "study_id") %>%
  dplyr::mutate(age_mid = purrr::map_dbl(strata, function(x){
    nums <- as.numeric(stringr::str_split_fixed(x, "-", n = 2))
    nums[nums == 999] <- 100
    return(mean(nums))})) %>%
  dplyr::select(c("study_id", "age_mid", "median", "LCI", "UCI", "sero")) %>%
  dplyr::mutate(median = median * 100,
                LCI = LCI * 100,
                UCI = UCI * 100)


# overlay crude
dscdat_age <- dscdat_age %>%
  dplyr::mutate(median = crudeIFR,
                LCI = crudeIFR,
                UCI = crudeIFR,
                sero = "crude") %>%
  dplyr::select(-c("crudeIFR")) %>%
  dplyr::filter(study_id %in% unique(retmap_age$study_id))
# bring together
age_IFR_plotdat <- dplyr::bind_rows(retmap_age, dscdat_age) %>%
  dplyr::mutate(sero = factor(sero, levels = c("crude", "reg", "serorev"), labels = c("Crude", "SeroInf.", "SeroRev."))) %>%
  dplyr::inner_join(studyidnames, ., by = "study_id")

# plot
ageIFR_plotObj <- ggplot() +
  geom_pointrange(data = age_IFR_plotdat,
                  aes(x = age_mid, y = median, ymin = LCI, ymax = UCI, color =  study_id),
                  alpha = 0.5, shape = 16, size = 0.9) +
  facet_wrap(.~sero) +
  scale_color_manual("Study ID", values = mycolors) +
  ylab("Age-Specific IFR (95% Cred. Int.)") + xlab("Mid. Age") +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

jpeg("figures/final_figures/Figure_age_IFR.jpg", width = 10, height = 7, units = "in", res = 600)
plot(ageIFR_plotObj)
graphics.off()



#............................................................
#---- IFR-Death Comparisons #----
#...........................................................
percap_deaths <- dscdat %>%
  dplyr::filter(!duplicated(study_id)) %>%
  dplyr::select(c("study_id", "data")) %>%
  group_by(study_id) %>%
  dplyr::mutate(deaths = purrr::map(data, "deaths_TSMCMC"),
                deaths = purrr::map_dbl(deaths, function(x){sum(x$deaths)}),
                popN = purrr::map(data, "prop_pop"),
                popN = purrr::map_dbl(popN, function(x){sum(x$popN)})) %>%
  tidyr::unnest(cols = c("deaths", "popN")) %>%
  dplyr::summarise(cumdeaths = sum(deaths),
                   popN = sum(popN)) %>%
  dplyr::mutate(capita = cumdeaths/popN) %>%
  dplyr::select(c("study_id", "capita")) %>%
  dplyr::left_join(studyidnames, ., by = "study_id") %>%
  dplyr::select(-c("study_id"))

crude_deaths <- pretty_tbl2 %>%
  dplyr::select(c("names", "crudeIFR")) %>%
  dplyr::rename(crude = crudeIFR)

# bring together others deaths
extra_deaths <- percap_deaths %>%
  dplyr::mutate(capita = capita * 100) %>%
  dplyr::left_join(., crude_deaths, by = "names") %>%
  #dplyr::left_join(., y = excess_deaths, by = "study_id") %>%
  #dplyr::left_join(., y = carehome_deaths, by = "study_id") %>%
  tidyr::pivot_longer(., cols = -c("names"), names_to = "param", values_to = "est") %>%
  dplyr::mutate(param = factor(param,
                               levels = c("crude", "capita", "excess", "carehome"),
                               labels = c("Crude", "Per Capita", "Excess", "Care-Home")
  ))


# modelled
seroinf <- tbl2 %>%
  dplyr::left_join(studyidnames, ., by = "study_id") %>%
  dplyr::mutate(median = median * 100,
                LCI = LCI * 100,
                UCI = UCI * 100,
                upperbound_inf = ifelse(UCI >= 3 & sero == "SeroInf.", 3, NA),
                upperbound_rev = ifelse(UCI >= 3 & sero == "SeroRev.", 3, NA),
                UCI = ifelse(UCI >= 3, 3, UCI))

# plot out
death_type_plot <- ggplot() +
    geom_point(data = extra_deaths, aes(x = names, y = est, fill = param, shape = param),
               size = 6.5, alpha = 0.8) +
    geom_pointrange(data = seroinf,
                    aes(x = names, ymin = LCI, y = median, ymax = UCI, group = sero, color = sero),
                    size = 0.75, alpha = 0.5, position = position_dodge(width = 0.5)) +
    geom_point(data = seroinf,
               aes(x = names, y = upperbound_inf),
               color = "#4285F4", size = 3, alpha = 0.9, shape = 3,
               position = position_nudge(x = -0.125)) +
    geom_point(data = seroinf,
               aes(x = names, y = upperbound_rev),
               color = "#EA4335", size = 3, alpha = 0.9, shape = 3,
               position = position_nudge(x = 0.125)) +
    scale_color_manual("Sero. Type", values = c("#4285F4", "#EA4335")) +
    scale_shape_manual("Death Type", values = c(23, 22, 24, 25)) +
    scale_fill_manual("Death Type", values = c(wesanderson::wes_palette("IsleofDogs2", type = "discrete"))) +
    ylab("Infection Fatality Rate (%)") +
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
