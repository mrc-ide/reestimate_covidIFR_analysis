####################################################################################
## Purpose: Plot for Figure 1 Showing Delays and Inference Framework
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
#----- Simulation Results #-----
#...........................................................
dat_map <- tibble::tibble(lvl = c("reg", "serorev"),
                          mod = c("data/param_map/Fig_ConceptualFits/reg_mod_rung50_burn10000_smpl10000.RDS",
                                  "data/param_map/Fig_ConceptualFits/serorev_mod_rung50_burn10000_smpl10000.RDS")) %>%
  dplyr::mutate(mod = purrr::map(mod, readRDS)) %>%
  tidyr::unnest(cols = mod)

fits <- tibble::tibble(lvl = c("reg", "serorev"),
                       fit = c("results/Fig_ConceptualFits/reg_mod_rung50_burn10000_smpl10000.RDS",
                               "results/Fig_ConceptualFits/serorev_mod_rung50_burn10000_smpl10000.RDS")) %>%
  dplyr::mutate(fit = purrr::map(fit, readRDS))

# bring together
param_map <- dplyr::left_join(dat_map, fits, by = "lvl")


#............................................................
# check that simulations worked and converged
#...........................................................
param_map$fit[[1]]$mcmcout$output %>%
  dplyr::filter(rung == "rung50") %>%
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
  dplyr::filter(rung == "rung50") %>%
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
COVIDCurve::get_cred_intervals(param_map$fit[[1]], whichrung = "rung50", what = "IFRparams", by_chain = F)
COVIDCurve::get_cred_intervals(param_map$fit[[2]], whichrung = "rung50", what = "IFRparams", by_chain = F)

COVIDCurve::get_cred_intervals(param_map$fit[[1]], whichrung = "rung50", what = "Serotestparams", by_chain = F)
COVIDCurve::get_cred_intervals(param_map$fit[[2]], whichrung = "rung50", what = "Serotestparams", by_chain = F)

#......................
# quick plot
#......................
fatalitydata <- tibble::tibble(param = c("ma1", "ma2", "ma3"),
                               IFR = c(1e-3, 0.05, 0.1))

#......................
# get ifrs quick plot to confirm
#......................
ifrs_reg <- COVIDCurve::get_cred_intervals(IFRmodel_inf = param_map$fit[[1]],
                                       whichrung = "rung50",
                                       what = "IFRparams", by_chain = FALSE) %>%
  dplyr::left_join(., fatalitydata, by = "param")

ifrs_reg %>%
  ggplot() +
  geom_pointrange(aes(x = param, ymin = LCI, ymax = UCI, y = median),
                  color = "#969696", size = 1.2) +
  geom_hline(aes(yintercept = IFR),
             color = "#3182bd", size = 3, alpha = 0.75, show.legend = F) +
  facet_wrap(~param, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right") +
  xlab("") + ylab("Median (95% CrIs)")


ifrs_serorev <- COVIDCurve::get_cred_intervals(IFRmodel_inf = param_map$fit[[2]],
                                           whichrung = "rung50",
                                           what = "IFRparams", by_chain = FALSE) %>%
  dplyr::left_join(., fatalitydata, by = "param")

ifrs_serorev %>%
  ggplot() +
  geom_pointrange(aes(x = param, ymin = LCI, ymax = UCI, y = median),
                  color = "#969696", size = 1.2) +
  geom_hline(aes(yintercept = IFR),
             color = "#3182bd", size = 3, alpha = 0.75, show.legend = F) +
  facet_wrap(~param, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right") +
  xlab("") + ylab("Median (95% CrIs)")

#......................
# gelman
#......................

COVIDCurve::get_gelman_rubin_diagnostic(param_map$fit[[1]])
COVIDCurve::get_gelman_rubin_diagnostic(param_map$fit[[2]])





#...........................................................
# get IFR over time
#...........................................................
# IFR over time
get_data_IFR_longit <- function(simdat, modelobj, fit, sens, spec, dwnsmpl = 1e2) {

  #......................
  # cumulativde data through time
  #......................
  cumdat <-  simdat$StrataAgg_Seroprev %>%
    dplyr::left_join(simdat$StrataAgg_TimeSeries_Death, ., by = c("ObsDay", "Strata")) %>%
    dplyr::left_join(., modelobj$demog, by = "Strata") %>%
    dplyr::group_by(Strata) %>%
    dplyr::mutate(
      cumdeaths = cumsum(Deaths),
      RGIFR = cumdeaths/(popN * rogan_gladen(obs_prev = ObsPrev, sens = sens, spec = spec)),
      CrudeIFR = cumdeaths/(popN * ObsPrev)
    ) %>%
    dplyr::ungroup(.)

  #......................
  # get posterior IFRs
  #......................
  mcmcout.nodes <- fit$mcmcout$output
  mcmcout.nodes <- mcmcout.nodes %>%
    dplyr::mutate(logposterior = loglikelihood + logprior)
  probs <- COVIDCurve:::convert_post_probs(mcmcout.nodes$logposterior)
  # downsample
  dwnsmpl_rows <- sample(1:nrow(mcmcout.nodes), size = dwnsmpl,
                         prob = probs)
  dwnsmpl_rows <- sort(dwnsmpl_rows)
  mcmcout.nodes <- mcmcout.nodes[dwnsmpl_rows, ]

  # inferred IFR
  infIFR <- mcmcout.nodes %>%
    dplyr::select(c("iteration", "ma3"))

  #......................
  # bring together
  #......................
  ret <- list(cumdat = cumdat, infIFR = infIFR)
  return(ret)
}

# run
param_map$plotdat <- purrr::pmap(param_map[, c("simdat", "modelobj", "fit")],
                                 get_data_IFR_longit,
                                 dwnsmpl = 1e2, sens = 0.85, spec = 0.95) # sens and spec from sim


#............................................................
#----- Simulation Figure Crude IFR vs. Modelled #-----
#...........................................................
#............................
# Plot No SeroReversion
#............................
# truth from simulation
fatalitydata_intercept <- 0.1

cumplotdat <- param_map$plotdat[[1]]$cumdat %>%
  dplyr::filter(ObsDay >= 50 & Strata == "ma3") %>%
  dplyr::select(c("ObsDay", "RGIFR", "CrudeIFR")) %>%
  tidyr::pivot_longer(., cols = c("RGIFR", "CrudeIFR"), names_to = "IFRlvl", values_to = "IFR") %>%
  dplyr::mutate(IFRlvl = factor(IFRlvl, levels = c("RGIFR", "CrudeIFR"), labels = c("Test-Adj.", "Crude")))

infplotdat <- param_map$plotdat[[1]]$infIFR


no_serorev_infIFR_plotObj <- ggplot() +
  geom_hline(data = infplotdat, aes(yintercept = ma3), color = "#d9d9d9", alpha = 0.8, size = 0.9) +
  geom_hline(yintercept = fatalitydata_intercept, color = "#252525", size = 1.2,
             linetype = "dashed") +
  geom_line(data = cumplotdat,
            aes(x = ObsDay, y = IFR, color = IFRlvl), size = 1.1) +
  scale_color_manual("IFR Calc.", values = c("#F5390D", "#EDAC2C")) +
  xlab("Time (days)") +
  theme(plot.title = element_blank(),
        axis.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 16),
        axis.text = element_text(family = "Helvetica", hjust = 0.5, size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(family = "Helvetica", hjust = 0.5, vjust = 0.5, size = 15),
        legend.key=element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "#000000", size = 1))


#................................
# Plot w. SeroReversion
#................................
cumplotdat <- param_map$plotdat[[2]]$cumdat %>%
  dplyr::filter(ObsDay >= 50 & Strata == "ma3") %>%
  dplyr::select(c("ObsDay", "RGIFR", "CrudeIFR")) %>%
  tidyr::pivot_longer(., cols = c("RGIFR", "CrudeIFR"), names_to = "IFRlvl", values_to = "IFR") %>%
  dplyr::mutate(IFRlvl = factor(IFRlvl, levels = c("RGIFR", "CrudeIFR"), labels = c("Test Adj.", "Crude")))

infplotdat <- param_map$plotdat[[2]]$infIFR

serorev_infIFR_plotObj <- ggplot() +
  geom_hline(data = infplotdat, aes(yintercept = ma3), color = "#d9d9d9", alpha = 0.8, size = 0.9) +
  geom_hline(yintercept = fatalitydata_intercept, color = "#252525", size = 1.2,
             linetype = "dashed") +
  geom_line(data = cumplotdat,
            aes(x = ObsDay, y = IFR, color = IFRlvl), size = 1.1) +
  scale_color_manual("IFR Calc.", values = c("#F5390D", "#EDAC2C")) +
  xlab("Time (days)") +
  theme(plot.title = element_blank(),
        axis.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 16),
        axis.text = element_text(family = "Helvetica", hjust = 0.5, size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(family = "Helvetica", hjust = 0.5, vjust = 0.5, size = 15),
        legend.key=element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "#000000", size = 1))





#............................................................
#----- Figure of Crude/RG vs Modelled IFR #-----
#...........................................................
popN <- param_map$modelobj[[1]]$demog %>%
  dplyr::filter(Strata == "ma3") %>%
  dplyr::pull("popN")

#......................
# tidy up and combine
#......................
cumdeaths <- param_map$simdat[[1]]$StrataAgg_TimeSeries_Death %>%
  dplyr::filter(Strata == "ma3") %>%
  dplyr::mutate(cumDeaths = cumsum(Deaths)/popN) %>%
  dplyr::rename(time = ObsDay)


serodf_norev <- param_map$simdat[[1]]$StrataAgg_Seroprev %>%
  dplyr::filter(Strata == "ma3") %>%
  dplyr::select(c("ObsDay", "TruePrev", "ObsPrev")) %>%
  dplyr::rename(time = ObsDay,
                regTruePrev = TruePrev,
                regObsPrev = ObsPrev)

serodf_rev <- param_map$simdat[[2]]$StrataAgg_Seroprev %>%
  dplyr::filter(Strata == "ma3") %>%
  dplyr::select(c("ObsDay", "ObsPrev")) %>%
  dplyr::rename(time = ObsDay,
                revObsPrev = ObsPrev)


# combine
datdf <- dplyr::left_join(cumdeaths, serodf_norev, by = "time") %>%
  dplyr::left_join(., serodf_rev, by = "time")

# long
plotdatdf <- datdf %>%
  dplyr::select(-c("Strata", "Deaths")) %>%
  tidyr::pivot_longer(., cols = -c("time"), names_to = "datlevel", values_to = "prop") %>%
  dplyr::mutate(datlevel = factor(datlevel,
                                  levels = c("regTruePrev", "cumDeaths",
                                              "regObsPrev", "revObsPrev"),
                                  labels = c("Infected",
                                             "Deaths",
                                             "Obs. Seroprev.",
                                             "Obs. Serorev."
                                  )))

#......................
# labels and arrows
#......................
arrows <- tibble::tibble(
  lvl =  c("mod", "serocon", "sens", "spec", "serorev"),
  x =    c(102,    172,       250,    10,     290),
  xend = c(136,    187,       250,    10,     290),
  y =    c(0.02,   0.5,       0.66,  0,      0.675),
  yend = c(0.02,   0.5,       0.58,   0.05,   0.445)
)


labels <- tibble::tibble(
  lvl =    c("mod",       "serocon",    "sens",    "spec",  "serorev"),
  label =  c("I-D Delay", "I-S Delay",  "Sens.",   "Spec.",  "I-R Delay"),
  x =      c(168,          149,          227,       12,       257),
  y =      c(0.02,         0.5,         0.595,    0.09,      0.50),
)



#......................
# plot
#......................
delay_plotObj <- plotdatdf %>%
  ggplot() +
  geom_line(aes(x = time, y = prop, color = datlevel), size = 2) +
  geom_segment(data = arrows, aes(x = x, xend = xend, y = y, yend = yend),
               size = 1.8, color = "#000000",
               arrow = arrow(length = unit(0.025, "npc"))) +
  geom_text(data = labels, aes(x = x, y = y, label = label),
            size = 5,
            family = "Helvetica",
            fontface = "bold",
            vjust = 0.5,
            hjust = 0.5) +
  scale_color_manual(values = c("#224B8B", "#06CDF4", "#3BABFD", "#B2DDFF")) +
  xlab("") + ylab("Proportion") +
  theme(plot.title = element_blank(),
        axis.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 16),
        axis.text = element_text(family = "Helvetica", hjust = 0.5, size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(family = "Helvetica", hjust = 0.5, vjust = 0.5, size = 15),
        legend.key=element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "#000000", size = 1))

# view
delay_plotObj


#..........................
# Inset of Infxn Curve
#..........................
# original curve that we simulated from w/ poisson distribution
infxnsdf <- tibble::tibble(time = 1:length(dat_map$infxns[[1]]),
                           infxns = dat_map$infxns[[1]])

infxninset_plotObj <- ggplot() +
  geom_line(data = infxnsdf, aes(x = time, y = infxns), color = "#252525", size = 1.1) +
  xlab("Time (days)") + ylab("Daily Infections") +
  theme(axis.title = element_text(family = "Helvetica", hjust = 0.5, size = 7.5),
        axis.text = element_text(family = "Helvetica", hjust = 0.5, size = 6),
        legend.position = "right",
        legend.title = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "#000000", size = 0.25),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "#000000", size = 1))


#......................
# put inset in
#......................
(toprow <- cowplot::ggdraw() +
    cowplot::draw_plot(delay_plotObj, x = 0, y = 0, width = 1, height = 1, scale = 1) +
    cowplot::draw_plot(infxninset_plotObj, x = 0.1, y= 0.65, width = 0.2, height = 0.3))


#......................
# bottom row
#......................
no_serorev_infIFR_plotObj <- no_serorev_infIFR_plotObj +
  theme(legend.position = "none",
        plot.margin = unit(c(0, 0, 0.1, 1),"cm"))
serorev_infIFR_plotObj <- serorev_infIFR_plotObj +
  theme(plot.margin = unit(c(0, 0, 0.1, 1),"cm"))

bottomrow <- cowplot::plot_grid(no_serorev_infIFR_plotObj, serorev_infIFR_plotObj,
                                align = "h", labels = c("(B)", "(C)"), nrow = 1, ncol = 2,
                                rel_widths = c(0.7, 1))
#......................
# bring together
#......................
(mainfig <- cowplot::plot_grid(toprow, bottomrow, labels = c("(A)", ""),
                               nrow = 2, rel_heights = c(0.8, 1)))

dir.create("figures/final_figures/", recursive = TRUE)
jpeg("figures/final_figures/Fig_crude_RG_versus_modelled_diagram.jpg",
     height = 8, width = 8, units = "in", res = 800)
plot(mainfig)
graphics.off()
