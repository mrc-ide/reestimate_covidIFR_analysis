####################################################################################
## Purpose: Plot for Figure 1 Showing Delays and Inference Framework
##
## Notes:
####################################################################################
set.seed(48)
library(COVIDCurve)
library(tidyverse)
source("R/covidcurve_helper_functions.R")
source("R/crude_plot_summ.R")
source("R/my_themes.R")


#............................................................
# read results in
#...........................................................
dat_map <- tibble::tibble(lvl = c("reg", "serorev"),
                            mod = c("data/param_map/Fig1_ConceptualFits/reg_mod_rung50_burn10000_smpl10000.RDS",
                                    "data/param_map/Fig1_ConceptualFits/serorev_mod_rung50_burn10000_smpl10000.RDS")) %>%
  dplyr::mutate(mod = purrr::map(mod, readRDS)) %>%
  tidyr::unnest(cols = mod)

fits <- tibble::tibble(lvl = c("reg", "serorev"),
                       fit = c("results/Fig1_ConceptualFits/reg_mod_rung50_burn10000_smpl10000.RDS",
                               "results/Fig1_ConceptualFits/serorev_mod_rung50_burn10000_smpl10000.RDS")) %>%
  dplyr::mutate(fit = purrr::map(fit, readRDS))

# bring together
param_map <- dplyr::left_join(dat_map, fits, by = "lvl")


#...........................................................
# get IFR over time
#...........................................................
# IFR over time
get_data_IFR_longit <- function(simdat, modelobj, fit, sens, spec, dwnsmpl = 1e2) {

  # need popN for crude IFR
  demog <- modelobj$demog

  #......................
  # cumulativde data through time
  #......................
  cumdat <-  simdat$StrataAgg_Seroprev %>%
    dplyr::left_join(simdat$StrataAgg_TimeSeries_Death, ., by = c("ObsDay", "Strata")) %>%
    dplyr::left_join(., demog, by = "Strata") %>%
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
  # Log-Sum-Exp trick
  convert_post_probs <- function(logpost) {
    exp(logpost - (log(sum(exp(logpost - max(logpost)))) + max(logpost)))
  }
  probs <- convert_post_probs(mcmcout.nodes$logposterior)
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
# Plot No SeroReversion
#...........................................................
# truth from simulation
fatalitydata_intercept <- 0.1

cumplotdat <- param_map$plotdat[[1]]$cumdat %>%
  dplyr::select(c("ObsDay", "RGIFR", "CrudeIFR")) %>%
  dplyr::filter(ObsDay >= 50) %>%
  tidyr::pivot_longer(., cols = c("RGIFR", "CrudeIFR"), names_to = "IFRlvl", values_to = "IFR") %>%
  dplyr::mutate(IFRlvl = factor(IFRlvl, levels = c("RGIFR", "CrudeIFR"), labels = c("Rogan-Gladen", "Crude")))

infplotdat <- param_map$plotdat[[1]]$infIFR


no_serorev_infIFR_plotObj <- ggplot() +
  geom_hline(data = infplotdat, aes(yintercept = ma1), color = "#d9d9d9", alpha = 0.8, size = 0.9) +
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
# Plot w. SeroReversion
#...........................................................

cumplotdat <- param_map$plotdat[[2]]$cumdat %>%
  dplyr::select(c("ObsDay", "RGIFR", "CrudeIFR")) %>%
  dplyr::filter(ObsDay >= 50) %>%
  tidyr::pivot_longer(., cols = c("RGIFR", "CrudeIFR"), names_to = "IFRlvl", values_to = "IFR") %>%
  dplyr::mutate(IFRlvl = factor(IFRlvl, levels = c("RGIFR", "CrudeIFR"), labels = c("Rogan-Gladen", "Crude")))

infplotdat <- param_map$plotdat[[2]]$infIFR

serorev_infIFR_plotObj <- ggplot() +
  geom_hline(data = infplotdat, aes(yintercept = ma1), color = "#d9d9d9", alpha = 0.8, size = 0.9) +
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
#----- Conceptual Diagram #-----
#...........................................................
popN <- param_map$modelobj[[1]]$demog$popN # only one strata

#......................
# tidy up and combine
#......................
cuminxns <-  param_map$infxns[[1]] %>%
  dplyr::mutate(cumincidence = cumsum(infxns)/popN) %>%
  dplyr::select(-c("infxns"))


cumdeaths <- param_map$simdat[[1]]$Agg_TimeSeries_Death %>%
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
  x =    c(112,    147,       250,    10,     290),
  xend = c(160,    170.5,     250,    10,     290),
  y =    c(0.1,    0.5,       0.71,   0,      0.72),
  yend = c(0.1,    0.5,       0.625,  0.05,   0.62)
)


labels <- tibble::tibble(
  lvl =    c("mod",       "serocon",    "sens",    "spec",  "serorev"),
  label =  c("O-D Delay", "O-S Delay",  "Sens.",   "Spec.",  "O-R Delay"),
  x =      c(185,          195,          265,       12,       220),
  y =      c(0.1,          0.5,         0.6675,    0.07,      0.55),
)


#......................
# plot
#......................
delay_plotObj <- plotdatdf %>%
  dplyr::filter(datlevel != "True Seroprev.") %>% # discuss wheter true goes into figure or not
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

#............................................................
#----- Inset with Posteriors #-----
#...........................................................
true_incidence <- param_map$infxns[[1]] %>%
  dplyr::filter(infxns != -1) %>%  # remove missing
  dplyr::mutate(incidence = infxns/popN)

infxninset_plotObj <- ggplot() +
  geom_line(data = true_incidence, aes(x = time, y = incidence), color = "#252525", size = 1.1, linetype = "dashed") +
  xlab("Time (days)") + ylab("Daily Incidence") +
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
# Top row figure
#......................
(toprow <- cowplot::ggdraw() +
    cowplot::draw_plot(delay_plotObj, x = 0, y = 0, width = 1, height = 1, scale = 1) +
    cowplot::draw_plot(infxninset_plotObj, x = 0.085, y= 0.65, width = 0.3, height = 0.3))

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
(mainfig <- cowplot::plot_grid(toprow, bottomrow, labels = c("(A)", ""),
                              nrow = 2, rel_heights = c(1, 0.6)))

dir.create("figures/final_figures/", recursive = TRUE)
jpeg("figures/final_figures/Fig1.jpg",
     height = 9, width = 8, units = "in", res = 500)
plot(mainfig)
graphics.off()
