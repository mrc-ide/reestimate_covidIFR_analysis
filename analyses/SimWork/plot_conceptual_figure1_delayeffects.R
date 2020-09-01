####################################################################################
## Purpose: Plot for Figure 1 Showing Delays and Inference Framework
##
## Notes:
####################################################################################
set.seed(48)
library(COVIDCurve)
library(tidyverse)
source("R/covidcurve_helper_functions.R")
source("R/my_themes.R")


#............................................................
# read results in
#...........................................................

#...........................................................
# Plot Results
#...........................................................
# IFR over time
cumdat <-  dat$AggSeroPrev %>%
  dplyr::rename(ObsDay = event_obs_day) %>%
  dplyr::left_join(dat$AggDeath, ., by = c("ObsDay", "Strata"))

cumdat <- cumdat %>%
  dplyr::mutate(
    cumdeaths = cumsum(Deaths),
    RGIFR = cumdeaths/TrueSeroCount,
    CrudeIFR = cumdeaths/(popN * ObsPrev)
  )


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
dwnsmpl_rows <- sample(1:nrow(mcmcout.nodes), size = 1e2,
                       prob = probs)
dwnsmpl_rows <- sort(dwnsmpl_rows)
mcmcout.nodes <- mcmcout.nodes[dwnsmpl_rows, ]

# inferred IFR
infIFR <- mcmcout.nodes %>%
  dplyr::select(c("iteration", "ma1"))

# truth
fatalitydata_intercept <- fatalitydata %>%
  dplyr::pull(c("IFR"))

#......................
# come together for plot right
#......................
infIFR_plotObj <- cumdat %>%
  dplyr::select(c("ObsDay", "RGIFR", "CrudeIFR")) %>%
  dplyr::filter(ObsDay >= 50) %>%
  tidyr::gather(., key = "IFRlvl", value = "IFR", 2:ncol(.)) %>%
  dplyr::mutate(IFRlvl = factor(IFRlvl, levels = c("RGIFR", "CrudeIFR"), labels = c("Rogan-Gladen \n Adj.", "Crude"))) %>%
  ggplot() +
  geom_hline(data = infIFR, aes(yintercept = ma1), color = "#d9d9d9", alpha = 0.8, size = 0.9) +
  geom_hline(yintercept = fatalitydata_intercept, color = "#252525", size = 1.2,
             linetype = "dashed") +
  geom_line(aes(x = ObsDay, y = IFR, color = IFRlvl), size = 1.1) +
  scale_color_manual("IFR Calc.", values = c("#F5390D", "#EDAC2C")) +
  xlab("Time") +
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
#......................
# tidy up and combine
#......................
cuminxns <- infxns %>%
  dplyr::mutate(cumincidence = cumsum(infxns)/popN) %>%
  dplyr::select(-c("infxns"))


cumdeaths <- dat$AggDeath %>%
  dplyr::mutate(cumDeaths = cumsum(Deaths)/popN) %>%
  dplyr::select(-c("Deaths", "Strata")) %>%
  dplyr::rename(time = ObsDay)


serodf <- dat$AggSeroPrev %>%
  dplyr::select(c("event_obs_day", "TruePrev", "ObsPrev")) %>%
  dplyr::rename(time = event_obs_day)
# combine
datdf <- dplyr::left_join(cumdeaths, cuminxns, by = "time") %>%
  dplyr::left_join(., serodf, by = "time")


# long
plotdatdf <- datdf %>%
  tidyr::gather(., key = "datlevel", value = "prop", 2:ncol(.)) %>%
  dplyr::mutate(datlevel = factor(datlevel,
                                  levels = c("cumincidence", "TruePrev", "ObsPrev", "cumDeaths"),
                                  labels = c("Cum. Prop. \n of Infected", "True Seroprev.", "Obs. Seroprev.", "Cum. Prop. \n of Deaths")))

#......................
# labels and arrows
#......................
arrows <- tibble::tibble(
  lvl =  c("mod", "serocon", "sens", "spec"),
  x =    c(112,    147,       250,    10),
  xend = c(160,    170.5,       250,    10),
  y =    c(0.1,    0.5,       0.71,    0),
  yend = c(0.1,    0.5,       0.625,   0.05)
)


labels <- tibble::tibble(
  lvl =    c("mod",       "serocon",    "sens",    "spec"),
  label =  c("O-D Delay", "O-S Delay",  "Sens.",   "Spec."),
  x =      c(185,          195,          265,       12),
  y =      c(0.1,          0.5,         0.6675,    0.07),
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
            vjust = c(0.5, 0.5, 0.5, 0.5),
            hjust = c(0.5, 0.5, 0.5, 0.5)
  ) +
  scale_color_manual(values = c("#224B8B", "#3BABFD", "#06CDF4")) +
  xlab("Time (days)") + ylab("Proportion") +
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
true_incidence <- infxns %>%
  dplyr::filter(infxns != -1) %>%  # remove missing
  dplyr::mutate(incidence = infxns/popN)

post_infxn_draws <- COVIDCurve::draw_posterior_infxn_cubic_splines(IFRmodel_inf = fit,
                                                                   whichrung = "rung1",
                                                                   by_chain = FALSE,
                                                                   by_strata = FALSE,
                                                                   dwnsmpl = 100)$curvedata %>%
  dplyr::mutate(inf_incidence = totinfxns/popN)

infxninset_plotObj <- ggplot() +
  geom_line(data = post_infxn_draws, aes(x = time, y = inf_incidence, group = sim), color = "#d9d9d9", size = 0.9, alpha = 0.8) +
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
# out
#......................
(toprow <- cowplot::ggdraw() +
    cowplot::draw_plot(delay_plotObj, x = 0, y = 0, width = 1, height = 1, scale = 1) +
    cowplot::draw_plot(infxninset_plotObj, x = 0.085, y= 0.65, width = 0.3, height = 0.3))

(fig1 <- cowplot::plot_grid(toprow, infIFR_plotObj, ncol = 1, align = "h",
                            labels = c("(A)", "(B)"), rel_heights = c(1, 0.5)))

dir.create("figures/SimCurves/", recursive = T)
jpeg("figures/SimCurves/conceptual_fig_1.jpg", width = 11, height = 8, units = "in", res = 500)
fig1
graphics.off()
saveRDS(fig1, "figures/SimCurves/conceptual_fig_1.RDS")


