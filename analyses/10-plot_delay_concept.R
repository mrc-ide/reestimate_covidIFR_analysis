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
    dplyr::mutate(
      cumdeaths = cumsum(Deaths),
      RGIFR = cumdeaths/(popN * rogan_gladen(obs_prev = ObsPrev, sens = sens, spec = spec)),
      CrudeIFR = cumdeaths/(popN * ObsPrev)
    ) %>%
    dplyr::mutate(RGIFR = ifelse(is.nan(RGIFR), 0, RGIFR)) # if denom is 0, just set to 0

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
    dplyr::select(c("iteration", "ma1"))

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
  dplyr::select(c("ObsDay", "RGIFR", "CrudeIFR")) %>%
  dplyr::filter(ObsDay >= 50) %>%
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
  dplyr::select(c("ObsDay", "RGIFR", "CrudeIFR")) %>%
  dplyr::filter(ObsDay >= 50) %>%
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
  dplyr::pull("popN")

#......................
# tidy up and combine
#......................
cuminxns <-  tibble::tibble(time = 1:length(param_map$infxns[[1]]),
                            infxns = param_map$infxns[[1]])%>%
  dplyr::mutate(cumincidence = (cumsum(infxns) * 1/3 )/popN) %>% # rho infections were split evenly across
  dplyr::select(-c("infxns"))


cumdeaths <- param_map$simdat[[1]]$StrataAgg_TimeSeries_Death %>%
  dplyr::mutate(cumDeaths = cumsum(Deaths)/popN) %>%
  dplyr::rename(time = ObsDay)


serodf_norev <- param_map$simdat[[1]]$StrataAgg_Seroprev %>%
  dplyr::select(c("ObsDay", "TruePrev", "ObsPrev")) %>%
  dplyr::rename(time = ObsDay,
                regTruePrev = TruePrev,
                regObsPrev = ObsPrev)

serodf_rev <- param_map$simdat[[2]]$StrataAgg_Seroprev %>%
  dplyr::select(c("ObsDay", "ObsPrev")) %>%
  dplyr::rename(time = ObsDay,
                revObsPrev = ObsPrev)


# combine
datdf <- dplyr::left_join(cumdeaths, cuminxns, by = "time") %>%
  dplyr::left_join(., serodf_norev, by = "time") %>%
  dplyr::left_join(., serodf_rev, by = "time")

# long
plotdatdf <- datdf %>%
  dplyr::select(-c("Deaths")) %>%
  tidyr::pivot_longer(., cols = -c("time"), names_to = "datlevel", values_to = "prop") %>%
  dplyr::mutate(datlevel = factor(datlevel,
                                  levels = c("cumincidence", "cumDeaths",
                                             "regTruePrev", "regObsPrev", "revObsPrev"),
                                  labels = c("Infected",
                                             "Deaths",
                                             "True Seroprev.", "Obs. Seroprev.",
                                             "Obs. Serorev."
                                  )))

#......................
# labels and arrows
#......................
arrows <- tibble::tibble(
  lvl =  c("mod", "serocon", "sens", "spec", "serorev"),
  x =    c(87,    139,       250,    10,     290),
  xend = c(160,    162,        250,    10,        290),
  y =    c(0.02,    0.4,       0.715,   0,     0.725),
  yend = c(0.02,    0.4,       0.615,  0.05,   0.29)
)


labels <- tibble::tibble(
  lvl =    c("mod",       "serocon",    "sens",    "spec",  "serorev"),
  label =  c("O-D Delay", "O-S Delay",  "Sens.",   "Spec.",  "O-R Delay"),
  x =      c(170,          110,          230,       12,       260),
  y =      c(0.05,          0.4,         0.65,    0.07,      0.50),
)


#......................
# plot
#......................
delay_plotObj <- plotdatdf %>%
  dplyr::filter(datlevel != "True Seroprev.") %>% # true prev was decided to be excluded from fig
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

delay_plotObj

#..........................
# Inset of Infxn Curve
#..........................
true_infxns <- tibble::tibble(time = 1:length(param_map$infxns[[1]]),
                                 infxns = param_map$infxns[[1]])

infxninset_plotObj <- ggplot() +
  geom_line(data = true_infxns, aes(x = time, y = infxns), color = "#252525", size = 1.1, linetype = "dashed") +
  xlab("Time (days)") + ylab("Daily INfections") +
  theme(axis.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 10),
        axis.text = element_text(family = "Helvetica", hjust = 0.5, size = 8),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(family = "Helvetica", hjust = 0.5, vjust = 0.5, size = 10),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "#000000", size = 0.25),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "#000000", size = 1))


#......................
# bottom row
#......................
no_serorev_infIFR_plotObj <- no_serorev_infIFR_plotObj +
  theme(legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 1),"cm"))
serorev_infIFR_plotObj <- serorev_infIFR_plotObj +
  theme(plot.margin = unit(c(0, 0, 0, 1),"cm"))

bottomrow <- cowplot::plot_grid(no_serorev_infIFR_plotObj, serorev_infIFR_plotObj,
                                align = "h", labels = c("(B)", "(C)"), nrow = 1, ncol = 2,
                                rel_widths = c(0.7, 1))
#......................
# bring together
#......................
(mainfig <- cowplot::plot_grid(delay_plotObj, bottomrow, labels = c("(A)", ""),
                               nrow = 2, rel_heights = c(0.6, 1)))

dir.create("figures/final_figures/", recursive = TRUE)
jpeg("figures/final_figures/Fig_crude_RG_versus_modelled_diagram.jpg",
     height = 9, width = 8, units = "in", res = 800)
plot(mainfig)
graphics.off()