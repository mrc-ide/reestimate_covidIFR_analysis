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
                          mod = c("data/param_map/SeroDays_Concept/OneDay_mod_rung50_burn10000_smpl20000.RDS",
                                  "data/param_map/SeroDays_Concept/TwoDays_mod_rung50_burn10000_smpl20000.RDS")) %>%
  dplyr::mutate(mod = purrr::map(mod, readRDS)) %>%
  tidyr::unnest(cols = mod)

fits <- tibble::tibble(lvl = c("oneday", "twoday"),
                       fit = c("results/SeroDays_Concept/OneDay_mod_rung50_burn10000_smpl20000.RDS",
                               "results/SeroDays_Concept/TwoDays_mod_rung50_burn10000_smpl20000.RDS")) %>%
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
COVIDCurve::get_gelman_rubin_diagnostic(param_map$fit[[1]])
COVIDCurve::get_gelman_rubin_diagnostic(param_map$fit[[2]])


#............................................................
# plot out
#...........................................................
#......................
# get ifrs
#......................
param_map$ifrs <- purrr::map(param_map$fit, COVIDCurve::get_cred_intervals,
                             whichrung = "rung1",
                             what = "IFRparams", by_chain = FALSE)
#......................
# get incidence curve
#......................
param_map$infxncurve <- purrr::map(param_map$fit, COVIDCurve::draw_posterior_infxn_cubic_splines,
                                                             dwnsmpl = 1e2,
                                                             by_chain = FALSE,
                                                             by_strata = TRUE)

#......................
# get crude ifrs from sim
#......................
fatalitydata <- tibble::tibble(param = c("ma1", "ma2", "ma3"),
                               IFR = c(1e-3, 0.05, 0.1))

#......................
# one seroday
#......................
oneday_ifrs <- ggplot() +
  geom_pointrange(data = param_map$ifrs[param_map$name == "OneDay_mod"][[1]],
                  aes(x = param, ymin = LCI, ymax = UCI, y = median),
                  color = "#969696", size = 1.2) +
  geom_point(data = fatalitydata,
             aes(x = param, y = IFR),
             color = "#000000", shape = 2,
             size = 3, alpha = 0.75, show.legend = F) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.90, hjust= 1, face = "bold"),
        legend.position = "right") +
  xlab("") + ylab("Median (95% CrIs)") +
  ylim(0, 0.125) +
  xyaxis_plot_theme +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 1),"cm"))
#......................
# two seroday
#......................
twoday_ifrs <- ggplot() +
  geom_pointrange(data = param_map$ifrs[param_map$name == "TwoDays_mod"][[1]],
                  aes(x = param, ymin = LCI, ymax = UCI, y = median),
                  color = "#969696", size = 1.2) +
  geom_point(data = fatalitydata,
             aes(x = param, y = IFR),
             color = "#000000", shape = 2,
             size = 3, alpha = 0.75, show.legend = F) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.90, hjust= 1, face = "bold"),
        legend.position = "right") +
  xlab("") + ylab("Median (95% CrIs)") +
  ylim(0, 0.125) +
  xyaxis_plot_theme +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 1),"cm"))


#......................
# come togehter
#......................
mainFig <- cowplot::plot_grid(oneday_ifrs, twoday_ifrs, labels = c("(A)", "(B)", ncol = 1))
jpeg("figures/final_figures/seroday_comparinson.jpg",
     width = 8, height = 6, units = "in", res = 500)
plot(mainFig)
graphics.off()


