####################################################################################
## Purpose: Plot for Seroday concept
##
## Notes:
####################################################################################
set.seed(48)
library(COVIDCurve)
library(tidyverse)
source("R/covidcurve_helper_functions.R")
source("R/extra_plotting_functions.R")
source("R/crude_plot_summ.R")
source("R/my_themes.R")

#............................................................
# read results in
#...........................................................
dat_map <- tibble::tibble(lvl = c("oneday", "twoday"),
                          mod = c("data/param_map/SeroDays_Concept/OneDay_mod_rung50_burn10000_smpl10000.RDS",
                                  "data/param_map/SeroDays_Concept/TwoDays_mod_rung50_burn10000_smpl10000.RDS")) %>%
  dplyr::mutate(mod = purrr::map(mod, readRDS)) %>%
  tidyr::unnest(cols = mod)

fits <- tibble::tibble(lvl = c("oneday", "twoday"),
                       fit = c("results/SeroDays_Concept/OneDay_mod_rung50_burn10000_smpl10000.RDS",
                               "results/SeroDays_Concept/TwoDays_mod_rung50_burn10000_smpl10000.RDS")) %>%
  dplyr::mutate(fit = purrr::map(fit, readRDS))

# bring together
param_map <- dplyr::left_join(dat_map, fits, by = "lvl")


#............................................................
# check that simulations worked and converged
#...........................................................
param_map$fit[[1]]$mcmcout$output %>%
  dplyr::filter(rung == "rung1") %>%
  dplyr::filter(stage == "sampling") %>%
  dplyr::select(c("chain", "iteration", "loglikelihood", "logprior")) %>%
  tidyr::gather(., key = "like", value = "val", 3:4) %>%
  ggplot() +
  geom_line(aes(x = iteration, y = val, color = chain), size = 0.25, alpha = 0.8) +
  scale_color_viridis_d() +
  ylab("") + xlab("Iteration") +
  facet_wrap(.~like, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.title = element_blank())

param_map$fit[[2]]$mcmcout$output %>%
  dplyr::filter(rung == "rung1") %>%
  dplyr::filter(stage == "sampling") %>%
  dplyr::select(c("chain", "iteration", "loglikelihood", "logprior")) %>%
  tidyr::gather(., key = "like", value = "val", 3:4) %>%
  ggplot() +
  geom_line(aes(x = iteration, y = val, color = chain), size = 0.25, alpha = 0.8) +
  scale_color_viridis_d() +
  ylab("") + xlab("Iteration") +
  facet_wrap(.~like, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.title = element_blank())


drjacoby::plot_mc_acceptance(param_map$fit[[1]]$mcmcout)
drjacoby::plot_mc_acceptance(param_map$fit[[2]]$mcmcout)
quick_sero_diagnostics(param_map$fit[[1]])
quick_sero_diagnostics(param_map$fit[[2]])
