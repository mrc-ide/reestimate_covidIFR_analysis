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
                                                               dwnsmpl = 1e2,
                                                               by_chain = TRUE,
                                                               by_strata = FALSE)$curvedata %>%
    dplyr::left_join(., truecurve, by = "time") # note this will be serially duplicated
  return(infxncurve)
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
  geom_line(aes(x = time, y = totinfxns, group = sim, color = chain), alpha = 0.8, size = 0.8) +
  geom_line(aes(x = time, y = truecurve), color = "#bdbdbd",
            linetype = "dashed", size = 1) +
  facet_grid(simtype ~ serotestchar, scales = "free") +
  scale_color_viridis_d("Chain") +
  labs(caption = "Viridis is inferred, Grey Line is Simulated Infection Curve") +
  xlab("Time") + ylab("Num. Infxns") +
  xyaxis_plot_theme +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 11, face = "bold"))

# with serorev
collated_infxn_curve_serorev_plotObj <- datmap %>%
  dplyr::filter(serorev_accounted) %>% # no seroreversion first
  dplyr::select(-c("sim")) %>%
  dplyr::mutate(simtype = purrr::map_chr(nm, function(x){ switch(x,
                                                                 "expgrowth" = {"Exp. Growth"},
                                                                 "intervene" = {"Outbreak Control"},
                                                                 "secondwave" = {"Second Wave"})})) %>%
  tidyr::unnest(cols = "infxncurve") %>%
  dplyr::select(c("simtype", "sens", "spec", "sim", "chain", "time", "totinfxns", "truecurve")) %>%
  dplyr::mutate(serotestchar = paste0("Sens: ", sens, " ; Spec: ", spec)) %>%
  ggplot() +
  geom_line(aes(x = time, y = totinfxns, group = sim, color = chain), alpha = 0.8, size = 0.8) +
  geom_line(aes(x = time, y = truecurve), color = "#bdbdbd",
            linetype = "dashed", size = 1) +
  facet_grid(simtype ~ serotestchar, scales = "free") +
  scale_color_viridis_d("Chain") +
  labs(caption = "Viridis is inferred, Grey Line is Simulated Infection Curve") +
  xlab("Time") + ylab("Num. Infxns") +
  xyaxis_plot_theme +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 11, face = "bold"))

# out
jpeg("figures/final_figures/collated_simulation_noserorev.jpg",
     width = 10, height = 7, units = "in", res = 600)
plot(collated_infxn_curve_NOserorev_plotObj)
graphics.off()

jpeg("figures/final_figures/collated_simulation_withserorev.jpg",
     width = 10, height = 7, units = "in", res = 600)
plot(collated_infxn_curve_serorev_plotObj)
graphics.off()



