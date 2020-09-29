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
  dplyr::left_join(., order, by = "study_id") %>%
  dplyr::arrange(order) %>%
  dplyr::select(-c("order")) %>%
  readr::write_tsv(., path = "tables/final_tables/overall_sens_spec_for_all.tsv")


#...........................................................
#---- Posterior Delay Param Table  #----
#...........................................................
datmap %>%
  dplyr::mutate(tbl = purrr::map(modout, get_study_delays_posts)) %>%
  dplyr::left_join(order, .) %>%
  dplyr::arrange(order) %>%
  dplyr::select(c("study_id", "lvl", "tbl")) %>%
  tidyr::unnest(cols = "tbl") %>%
  dplyr::mutate(lvl = factor(lvl,
                             levels = c("NoSeroRev", "SeroRev"),
                             labels = c("Without Serorev.", "With Serorev."))) %>%
  readr::write_tsv(., path = "tables/final_tables/onset_delay_params_posterior_tbl.tsv")

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
  dplyr::mutate(age_high = as.numeric(stringr::str_split_fixed(ageband, "-", n=2)[,2])) %>%
  dplyr::filter(age_high == max(age_high)) %>%
  dplyr::filter(obsday == seromidpt) %>%
  dplyr::ungroup(.) %>%
  dplyr::select(c("study_id", "location", "ageband", "obsdaymin", "obsmidday", "obsdaymax", "crude_obs_seroprev"))

# get location simple
locat <- dsc_agedat %>%
  dplyr::select(c("study_id", "location")) %>%
  dplyr::filter(!duplicated(.))

#..................................
# standard model
#...................................
noserorev_ppc_seroPlotObj <- datmap %>% # NB with these functions have already subsetted to oldest age group
  dplyr::filter(lvl == "NoSeroRev") %>%
  dplyr::mutate(plotdat = purrr::map(modout, get_seroprev_ppcs)) %>%
  dplyr::select(-c("modout")) %>% # for memory
  tidyr::unnest(cols = plotdat) %>%
  dplyr::rename(model_seroprev = seroprev) %>%
  dplyr::select(c("study_id", "sim", "ObsDay",
                  "seroprevlvl", "model_seroprev")) %>%
  dplyr::left_join(., locat) %>%
  ggplot() +
    geom_line(aes(x = ObsDay, y = model_seroprev, color = seroprevlvl), alpha = 0.5) +
    geom_rect(data = SeroPrevObs,
              aes(xmin = obsdaymin, xmax = obsdaymax, ymin = -Inf, ymax = Inf),
              fill = "#d9d9d9", alpha = 0.4) +
    geom_point(data = SeroPrevObs,
               aes(x = obsmidday, y = crude_obs_seroprev),
               color = "#000000", size = 1.2, alpha = 0.6) +
    facet_wrap(.~location, scales = "free_y") +
    scale_color_manual("Seroprev. \n Adjustment", values = c("#FFD301", "#246BCF"),
                       labels = c("Inferred 'Truth'", "Inferred 'Observed' - \n Rogan-Gladen Corrected")) +
    xyaxis_plot_theme +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.90, hjust= 1, face = "bold"),
          legend.position = "bottom",
          legend.key = element_rect(fill = "#252525"))

# out
jpeg("figures/final_figures/ppc_seroprev_NOSeroRev.jpg",
     width = 8, height = 11, units = "in", res = 600)
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
  dplyr::rename(model_seroprev = seroprev) %>%
  dplyr::select(c("study_id", "sim", "ObsDay",
                  "seroprevlvl", "model_seroprev")) %>%
  dplyr::left_join(., locat) %>%
  ggplot() +
  geom_line(aes(x = ObsDay, y = model_seroprev, color = seroprevlvl), alpha = 0.5) +
  geom_rect(data = SeroPrevObs,
            aes(xmin = obsdaymin, xmax = obsdaymax, ymin = -Inf, ymax = Inf),
            fill = "#d9d9d9", alpha = 0.4) +
  geom_point(data = SeroPrevObs,
             aes(x = obsmidday, y = crude_obs_seroprev),
             color = "#000000", size = 1.2, alpha = 0.6) +
  facet_wrap(.~location, scales = "free_y") +
  scale_color_manual("Seroprev. \n Adjustment", values = c("#FFD301", "#246BCF"),
                     labels = c("Inferred 'Truth'", "Inferred 'Observed' - \n Rogan-Gladen Corrected")) +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.90, hjust= 1, face = "bold"),
        legend.position = "bottom",
        legend.key = element_rect(fill = "#252525"))

# out
jpeg("figures/final_figures/ppc_seroprev_SeroRev.jpg",
     width = 8, height = 11, units = "in", res = 600)
plot(serorev_ppc_seroPlotObj)
graphics.off()



#............................................................
#---- PPC Death Data  #----
#...........................................................
#..................................
# standard model
#...................................

death_datmap <- datmap %>% # NB with these functions have already subsetted to oldest age group
  dplyr::filter(lvl == "NoSeroRev") %>%
  dplyr::mutate(plotdat = purrr::map(modout, get_death_ppcs)) %>%
  dplyr::select(-c("modout"))
ObsDeathDat <- death_datmap %>%
  dplyr::mutate(datclean = purrr::map(plotdat, "datclean")) %>%
  dplyr::select(c("study_id", "datclean")) %>%
  tidyr::unnest(cols = "datclean") %>%
  dplyr::left_join(locat, .)

InfDeathDat <- death_datmap %>%
  dplyr::mutate(postdat_long = purrr::map(plotdat, "postdat_long")) %>%
  dplyr::select(c("study_id", "postdat_long")) %>%
  tidyr::unnest(cols = "postdat_long") %>%
  dplyr::left_join(locat, .)


noserorev_ppc_DeathPlotObj <- ggplot() +
  geom_line(data = InfDeathDat, aes(x= time, y = inf_deaths, group = sim),
            size = 1.2, color = "#bdbdbd") +
  geom_line(data = ObsDeathDat, aes(x = ObsDay, y = obs_deaths),
            color = "#3182bd") +
  facet_wrap(.~location, scales = "free_y") +
  theme_bw() +
  xlab("Time") + ylab("Deaths") +
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5),
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
  dplyr::left_join(locat, .)

serorev_InfDeathDat <- serorev_death_datmap %>%
  dplyr::mutate(postdat_long = purrr::map(plotdat, "postdat_long")) %>%
  dplyr::select(c("study_id", "postdat_long")) %>%
  tidyr::unnest(cols = "postdat_long") %>%
  dplyr::left_join(locat, .)


serorev_ppc_DeathPlotObj <- ggplot() +
  geom_line(data = serorev_InfDeathDat, aes(x= time, y = inf_deaths, group = sim),
            size = 1.2, color = "#bdbdbd") +
  geom_line(data = ObsDeathDat, aes(x = ObsDay, y = obs_deaths),
            color = "#3182bd") +
  facet_wrap(.~location, scales = "free_y") +
  theme_bw() +
  xlab("Time") + ylab("Deaths") +
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, vjust = 0.5))

# out
jpeg("figures/final_figures/ppc_Deaths_SeroRev.jpg",
     width = 8, height = 11, units = "in", res = 600)
plot(serorev_ppc_DeathPlotObj)
graphics.off()
