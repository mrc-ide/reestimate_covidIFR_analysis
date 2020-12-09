##................................................................................................
## Purpose: Get Age Based Model PPCs
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
mod_NOserorev_retmap <- tibble::tibble(study_id = toupper(stringr::str_split(basename(mod_NOserorev_paths), "_age", simplify = T)[,1]),
                                       lvl = "NoSeroRev",
                                       paths = mod_NOserorev_paths)

# yes serorev
mod_serorev_paths <- list.files("results/Modfits_serorev/", full.names = T)
mod_serorev_retmap <- tibble::tibble(study_id = toupper(stringr::str_split(basename(mod_serorev_paths), "_age", simplify = T)[,1]),
                                     lvl = "SeroRev",
                                     paths = mod_serorev_paths)

#.....................
# come together
#......................
datmap <- dplyr::bind_rows(mod_NOserorev_retmap, mod_serorev_retmap)
datmap <- datmap %>%
  dplyr::mutate(modout = purrr::map(paths, readRDS))


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
                                                     whichrung = "rung50",
                                                     by_chain = F)
  # serocurve data
  serocurvedat <- seropnts %>%
    dplyr::select(c("sim", "ObsDay", dplyr::starts_with("RG_pd_"),
                    dplyr::starts_with("crude_pd_"))) %>%
    tidyr::pivot_longer(., cols = -c("sim", "ObsDay"),
                        names_to = "seroprev_strata_lvl", values_to = "seroprev") %>%
    dplyr::mutate(seroprevlvl = ifelse(stringr::str_detect(seroprev_strata_lvl, "RG_"), "RGCorr", "Crude"),
                  param = stringr::str_extract(seroprev_strata_lvl, "ma[0-9]+")) %>%
    dplyr::select(-c("seroprev_strata_lvl")) %>%
    tidyr::pivot_wider(., names_from = "seroprevlvl", values_from = "seroprev") %>%
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
                                                         whichrung = "rung50",
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
  postdat_long <- postdat_long %>%
    dplyr::filter(ageband == oldest)
  datclean <- datclean %>%
    dplyr::filter(ageband == oldest)
  # out
  out <- list(postdat_long = postdat_long,
              datclean = datclean)
  return(out)

}


#............................................................
#---- PPC Sero Plots  #----
#...........................................................
dsc_agedat <- readRDS("results/descriptive_results/descriptive_results_datamap.RDS") %>%
  dplyr::filter(breakdown == "ageband") %>%
  dplyr::filter(!grepl("_nch", study_id))
# crude data
crudedat <- dsc_agedat %>%
  dplyr::select(c("study_id", "location", "plotdat")) %>%
  tidyr::unnest(cols = "plotdat")

SeroPrevObs <- crudedat %>%
  dplyr::mutate(obsmidday = (obsdaymin + obsdaymax)/2) %>%
  dplyr::rename(crude_obs_seroprev = seroprev) %>%
  dplyr::group_by(study_id) %>%
  dplyr::mutate(age_high= as.numeric(stringr::str_extract(ageband, "[0-9]+?(?=])"))) %>%
  dplyr::filter(age_high == max(age_high)) %>%
  dplyr::filter(obsday == seromidpt) %>%
  dplyr::ungroup(.)
#......................
# get uncertainty
#......................
SeroPrevObs_binom <- SeroPrevObs %>%
  dplyr::filter(!is.na(n_positive)) %>%
  dplyr::filter(!is.na(n_tested)) %>%
  dplyr::mutate(crude_seroprev_obj = purrr::map2(n_positive, n_tested, .f = function(x,n){ binom.test(x,n) }),
                crude_seroprev_CI = purrr::map(crude_seroprev_obj, "conf.int"),
                crude_obs_seroprevLCI = purrr::map_dbl(crude_seroprev_CI, function(x){x[[1]]}),
                crude_obs_seroprevUCI = purrr::map_dbl(crude_seroprev_CI, function(x){x[[2]]}),
                crude_obs_seroprev = purrr::map_dbl(crude_seroprev_obj, "estimate"))
# observed 95% CI
SeroPrevObs_logit <- SeroPrevObs %>%
  dplyr::filter(study_id %in% c("ITA1", "SWE1", "DNK1")) %>%
  dplyr::mutate(crude_seroprev_obj = NA,
                crude_seroprev_CI = NA,
                crude_obs_seroprevLCI = serolci,
                crude_obs_seroprevUCI = serouci)


SeroPrevObs <- dplyr::bind_rows(SeroPrevObs_binom, SeroPrevObs_logit)%>%
  dplyr::select(c("study_id", "location", "ageband", "obsdaymin", "obsmidday", "obsdaymax",
                  "crude_obs_seroprev", "crude_obs_seroprevLCI", "crude_obs_seroprevUCI"))

#..................................
# standard model
#...................................
noserorev_ppc_seroPlotObj <- datmap %>% # NB with these functions have already subsetted to oldest age group
  dplyr::filter(lvl == "NoSeroRev") %>%
  dplyr::mutate(plotdat = purrr::map(modout, get_seroprev_ppcs)) %>%
  dplyr::select(-c("modout")) %>% # for memory
  tidyr::unnest(cols = plotdat) %>%
  dplyr::select(c("study_id", "sim", "ObsDay",
                  "RGCorr", "Crude")) %>%
  dplyr::left_join(., locatkey) %>%
  ggplot() +
  geom_line(aes(x = ObsDay, y = Crude, group = sim, color = "Crude"), alpha = 0.5) +
  geom_line(aes(x = ObsDay, y = RGCorr, group = sim, color = "RGCorr"), alpha = 0.5) +
  geom_rect(data = SeroPrevObs,
            aes(xmin = obsdaymin, xmax = obsdaymax, ymin = -Inf, ymax = Inf),
            fill = "#d9d9d9", alpha = 0.4) +
  geom_pointrange(data = SeroPrevObs,
             aes(x = obsmidday, y = crude_obs_seroprev,
                 ymin = crude_obs_seroprevLCI, max = crude_obs_seroprevUCI),
             color = "#000000", size = 0.75, alpha = 0.6) +
  facet_wrap(.~location, scales = "free_y") +
  scale_color_manual("Seroprev. \n Adjustment",
                     values = c("Crude" = "#FFD301", "RGCorr" = "#246BCF"),
                     labels = c("Inferred 'Truth'", "Inferred 'Observed' - \n Rogan-Gladen Corrected")) +
  ggtitle("Seroprevalence in Oldest Age-Group Modelled without Seroreversion") +
  xyaxis_plot_theme +
  ylab("Seropev.") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.90, hjust= 1, face = "bold"),
        legend.position = "bottom",
        legend.key = element_rect(fill = "#252525"))

# out
jpeg("figures/final_figures/ppc_seroprev_NOSeroRev.jpg",
     width = 11, height = 8, units = "in", res = 600)
plot(noserorev_ppc_seroPlotObj)
graphics.off()




#..................................
# SeroRev model
#...................................
serorev_ppc_seroPlotObj <- datmap %>% # NB with these functions have already subsetted to oldest age group
  dplyr::filter(lvl == "SeroRev") %>%
  dplyr::mutate(plotdat = purrr::map(modout, get_seroprev_ppcs)) %>%
  dplyr::select(-c("modout")) %>% # for memory
  tidyr::unnest(cols = plotdat) %>%
  dplyr::select(c("study_id", "sim", "ObsDay",
                  "RGCorr", "Crude")) %>%
  dplyr::left_join(., locatkey) %>%
  ggplot() +
  geom_line(aes(x = ObsDay, y = Crude, group = sim, color = "Crude"), alpha = 0.5) +
  geom_line(aes(x = ObsDay, y = RGCorr, group = sim, color = "RGCorr"), alpha = 0.5) +
  geom_rect(data = SeroPrevObs,
            aes(xmin = obsdaymin, xmax = obsdaymax, ymin = -Inf, ymax = Inf),
            fill = "#d9d9d9", alpha = 0.4) +
  geom_pointrange(data = SeroPrevObs,
                  aes(x = obsmidday, y = crude_obs_seroprev,
                      ymin = crude_obs_seroprevLCI, max = crude_obs_seroprevUCI),
             color = "#000000", size = 0.75, alpha = 0.6) +
  facet_wrap(.~location, scales = "free_y") +
  scale_color_manual("Seroprev. \n Adjustment",
                     values = c("Crude" = "#FFD301", "RGCorr" = "#246BCF"),
                     labels = c("Inferred 'Truth'", "Inferred 'Observed' - \n Rogan-Gladen Corrected")) +
  xyaxis_plot_theme +
  ggtitle("Seroprevalence in Oldest Age-Group Modelled with Seroreversion") +
  ylab("Seropev.") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.90, hjust= 1, face = "bold"),
        legend.position = "bottom",
        legend.key = element_rect(fill = "#252525"))


# out
jpeg("figures/final_figures/ppc_seroprev_SeroRev.jpg",
     width = 11, height = 8, units = "in", res = 600)
plot(serorev_ppc_seroPlotObj)
graphics.off()



#............................................................
#---- PPC Death Data  #----
#...........................................................
#..................................
# without seroreversion model
#...................................
death_datmap <- datmap %>% # NB with these functions have already subsetted to oldest age group
  dplyr::filter(lvl == "NoSeroRev") %>%
  dplyr::mutate(plotdat = purrr::map(modout, get_death_ppcs)) %>%
  dplyr::select(-c("modout"))
ObsDeathDat <- death_datmap %>%
  dplyr::mutate(datclean = purrr::map(plotdat, "datclean")) %>%
  dplyr::select(c("study_id", "datclean")) %>%
  tidyr::unnest(cols = "datclean") %>%
  dplyr::left_join(., locatkey)

InfDeathDat <- death_datmap %>%
  dplyr::mutate(postdat_long = purrr::map(plotdat, "postdat_long")) %>%
  dplyr::select(c("study_id", "postdat_long")) %>%
  tidyr::unnest(cols = "postdat_long") %>%
  dplyr::left_join(., locatkey)


noserorev_ppc_DeathPlotObj <- ggplot() +
  geom_line(data = InfDeathDat, aes(x= time, y = inf_deaths, group = sim),
            size = 1.2, color = "#bdbdbd") +
  geom_line(data = ObsDeathDat, aes(x = ObsDay, y = obs_deaths),
            color = "#3182bd") +
  facet_wrap(.~location, scales = "free_y") +
  theme_bw() +
  xlab("Time") + ylab("Deaths") +
  ggtitle("Deaths in Oldest Age-Group Modelled without Seroreversion") +
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, vjust = 0.5))

# out
jpeg("figures/final_figures/ppc_Deaths_NOSeroRev.jpg",
     width = 8, height = 11, units = "in", res = 600)
plot(noserorev_ppc_DeathPlotObj)
graphics.off()



#..................................
# seroreversion model
#...................................
serorev_death_datmap <- datmap %>% # NB with these functions have already subsetted to oldest age group
  dplyr::filter(lvl == "SeroRev") %>%
  dplyr::mutate(plotdat = purrr::map(modout, get_death_ppcs)) %>%
  dplyr::select(-c("modout"))
ObsDeathDat <- death_datmap %>%
  dplyr::mutate(datclean = purrr::map(plotdat, "datclean")) %>%
  dplyr::select(c("study_id", "datclean")) %>%
  tidyr::unnest(cols = "datclean") %>%
  dplyr::left_join(., locatkey)

serorev_DeathDat <- serorev_death_datmap %>%
  dplyr::mutate(postdat_long = purrr::map(plotdat, "postdat_long")) %>%
  dplyr::select(c("study_id", "postdat_long")) %>%
  tidyr::unnest(cols = "postdat_long") %>%
  dplyr::left_join(., locatkey)


serorev_ppc_DeathPlotObj <- ggplot() +
  geom_line(data = serorev_DeathDat, aes(x= time, y = inf_deaths, group = sim),
            size = 1.2, color = "#bdbdbd") +
  geom_line(data = ObsDeathDat, aes(x = ObsDay, y = obs_deaths),
            color = "#3182bd") +
  facet_wrap(.~location, scales = "free_y") +
  theme_bw() +
  xlab("Time") + ylab("Deaths") +
  ggtitle("Deaths in Oldest Age-Group Modelled with Seroreversion") +
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, vjust = 0.5))

# out
jpeg("figures/final_figures/ppc_Deaths_SeroRev.jpg",
     width = 8, height = 11, units = "in", res = 600)
plot(serorev_ppc_DeathPlotObj)
graphics.off()
