##................................................................................................
## Purpose: Collate Simulations
##
## Notes:
##................................................................................................
library(tidyverse)
library(COVIDCurve)
source("R/my_themes.R")

#............................................................
# lightweight function to make plot
#...........................................................
get_sim_IFR_curve <- function(curve, modout) {
  # pull out  true curve
  truecurve <- tibble::tibble(time = 1:length(curve),
                              truecurve = curve)
  # pull out inferred infection curve
  infxncurve <- COVIDCurve::draw_posterior_infxn_cubic_splines(IFRmodel_inf = modout,
                                                               whichrung = "rung50",
                                                               dwnsmpl = 1e2,
                                                               by_chain = TRUE,
                                                               by_strata = FALSE)$curvedata %>%
    dplyr::left_join(., truecurve, by = "time") # note this will be serially duplicated
  return(infxncurve)
}


get_sim_IFR_est <- function(modout) {
  COVIDCurve::get_cred_intervals(IFRmodel_inf = modout,
                                 by_chain = FALSE,
                                 whichrung = "rung50",
                                 what = "IFRparams") %>%
    # subset to oldest age group for space
    dplyr::filter(param == "ma5")
}


#............................................................
# read in data
#...........................................................
sim_noserorev_parammap <- readRDS("data/param_map/SimCurves_noserorev/simfit_param_map.RDS")
sim_serorev_parammap <- readRDS("data/param_map/SimCurves_noserorev/simfit_param_map.RDS")
# no serorev
noserorev_paths <- list.files("results/SimCurves_noserorev/", full.names = TRUE)
noserorev_paths <- tibble::tibble(sim = gsub("_NoSeroRev.RDS", "", basename(noserorev_paths)),
                                  path = noserorev_paths) %>%
  dplyr::mutate(modout = purrr::map(path, readRDS)) %>%
  dplyr::left_join(sim_noserorev_parammap, ., by = "sim") %>%
  dplyr::mutate(serorev_accounted = FALSE)

# serorev
serorev_paths <- list.files("results/SimCurves_serorev/", full.names = TRUE)
serorev_paths <- tibble::tibble(sim = gsub("_SeroRev.RDS", "", basename(serorev_paths)),
                                path = serorev_paths) %>%
  dplyr::mutate(modout = purrr::map(path, readRDS)) %>%
  dplyr::left_join(sim_serorev_parammap, ., by = "sim") %>%
  dplyr::mutate(serorev_accounted = TRUE)

#.....................
# come together
#......................
datmap <- dplyr::bind_rows(noserorev_paths, serorev_paths)
datmap$infxncurve <- purrr::pmap(datmap[, c("curve", "modout")], get_sim_IFR_curve)
datmap$ifrests <- purrr::pmap(datmap[, c("modout")], get_sim_IFR_est)

#.....................
# plot out
#......................
collated_infxn_curve_NOserorev_plotObj <- datmap %>%
  dplyr::filter(!serorev_accounted) %>% # no seroreversion first
  dplyr::select(-c("sim")) %>%
  dplyr::mutate(simtype = purrr::map_chr(nm, function(x){ switch(x,
                                                                 "expgrowth" = {"Exp. Growth"},
                                                                 "intervene" = {"Outbreak Control"},
                                                                 "secondwave" = {"Second Wave"})})) %>%
  tidyr::unnest(cols = "infxncurve") %>%
  dplyr::select(c("simtype", "sens", "spec", "sim", "chain", "time", "totinfxns", "truecurve")) %>%
  dplyr::mutate(serotestchar = paste0("Sens: ", sens, " ; Spec: ", spec)) %>%
  ggplot() +
  geom_line(aes(x = time, y = totinfxns, group = sim), color = "#9ecae1", alpha = 0.8, size = 0.8) +
  geom_line(aes(x = time, y = truecurve), color = "#969696", size = 1) +
  facet_grid(simtype ~ serotestchar, scales = "free") +
  scale_y_continuous(labels = scales::comma) +
  xlab("Time") + ylab("Num. Infxns") +
  xyaxis_plot_theme +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold"),
        plot.margin = unit(c(0.05, 0.05, 0.05, 1),"cm"))
# ifr ests
collated_Noserorev_plotObj <-  datmap %>%
  dplyr::filter(!serorev_accounted) %>% # no seroreversion first
  dplyr::select(-c("sim")) %>%
  dplyr::mutate(simtype = purrr::map_chr(nm, function(x){ switch(x,
                                                                 "expgrowth" = {"Exp. Growth"},
                                                                 "intervene" = {"Outbreak Control"},
                                                                 "secondwave" = {"Second Wave"})})) %>%
  tidyr::unnest(cols = "ifrests") %>%
  dplyr::select(c("simtype", "sens", "spec", "median", "LCI", "UCI")) %>%
  dplyr::mutate(serotestchar = paste0("Sens: ", sens, " ; Spec: ", spec),
                param = "Oldest Age Group",
                median = round(median * 100, 2),
                LCI = round(LCI * 100, 2),
                UCI = round(UCI * 100, 2)) %>%
  ggplot() +
  geom_pointrange(aes(x = param, y = median, ymin = LCI, ymax = UCI),
                  color = "#3182bd", size = 1.25) +
  geom_hline(yintercept = 20.0, color = "#969696", size = 1.5, linetype = "dashed") +
  facet_grid(simtype ~ serotestchar, scales = "free") +
  ylab("IFR (95% CrI)") +
  xyaxis_plot_theme +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold"),
        plot.margin = unit(c(0.05, 0.05, 0.05, 1),"cm"))

# come togther
together_NOserorev_plotObj <- cowplot::plot_grid(collated_infxn_curve_NOserorev_plotObj,
                                                 collated_Noserorev_plotObj,
                                                 nrow = 1, labels = c("(A)", "(B)"),
                                                 rel_widths = c(0.6, 0.4))
# out
jpeg("figures/final_figures/collated_simulation_noserorev.jpg",
     width = 11, height = 8, units = "in", res = 600)
plot(together_NOserorev_plotObj)
graphics.off()


# with serorev
collated_infxn_curve_serorev_plotObj <- datmap %>%
  dplyr::filter(serorev_accounted) %>% # seroreversion SECOND
  dplyr::select(-c("sim")) %>%
  dplyr::mutate(simtype = purrr::map_chr(nm, function(x){ switch(x,
                                                                 "expgrowth" = {"Exp. Growth"},
                                                                 "intervene" = {"Outbreak Control"},
                                                                 "secondwave" = {"Second Wave"})})) %>%
  tidyr::unnest(cols = "infxncurve") %>%
  dplyr::select(c("simtype", "sens", "spec", "sim", "chain", "time", "totinfxns", "truecurve")) %>%
  dplyr::mutate(serotestchar = paste0("Sens: ", sens, " ; Spec: ", spec)) %>%
  ggplot() +
  geom_line(aes(x = time, y = totinfxns, group = sim), color = "#9ecae1", alpha = 0.8, size = 0.8) +
  geom_line(aes(x = time, y = truecurve), color = "#969696", size = 1) +
  facet_grid(simtype ~ serotestchar, scales = "free") +
  scale_y_continuous(labels = scales::comma) +
  xlab("Time") + ylab("Num. Infxns") +
  xyaxis_plot_theme +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold"),
        plot.margin = unit(c(0.05, 0.05, 0.05, 1),"cm"))
# ifr ests
collated_serorev_plotObj <-  datmap %>%
  dplyr::filter(serorev_accounted) %>% # seroreversion SECOND
  dplyr::select(-c("sim")) %>%
  dplyr::mutate(simtype = purrr::map_chr(nm, function(x){ switch(x,
                                                                 "expgrowth" = {"Exp. Growth"},
                                                                 "intervene" = {"Outbreak Control"},
                                                                 "secondwave" = {"Second Wave"})})) %>%
  tidyr::unnest(cols = "ifrests") %>%
  dplyr::select(c("simtype", "sens", "spec", "median", "LCI", "UCI")) %>%
  dplyr::mutate(serotestchar = paste0("Sens: ", sens, " ; Spec: ", spec),
                param = "Oldest Age Group",
                median = round(median * 100, 2),
                LCI = round(LCI * 100, 2),
                UCI = round(UCI * 100, 2)) %>%
  ggplot() +
  geom_pointrange(aes(x = param, y = median, ymin = LCI, ymax = UCI),
                  color = "#3182bd", size = 1.25) +
  geom_hline(yintercept = 20.0, color = "#969696", size = 1.5, linetype = "dashed") +
  facet_grid(simtype ~ serotestchar, scales = "free") +
  ylab("IFR (95% CrI)") +
  xyaxis_plot_theme +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold"),
        plot.margin = unit(c(0.05, 0.05, 0.05, 1),"cm"))

# come togther
together_serorev_plotObj <- cowplot::plot_grid(collated_infxn_curve_serorev_plotObj,
                                               collated_serorev_plotObj,
                                               nrow = 1, labels = c("(A)", "(B)"),
                                               rel_widths = c(0.6, 0.5))


# out
jpeg("figures/final_figures/collated_simulation_withserorev.jpg",
     width = 11, height = 8, units = "in", res = 600)
plot(together_serorev_plotObj)
graphics.off()



