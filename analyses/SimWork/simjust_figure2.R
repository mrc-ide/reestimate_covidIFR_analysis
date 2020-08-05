source("R/my_themes.R")
library(tidyverse)
library(drjacoby)


#............................................................
# data wrangle
#...........................................................
modoutpath <- list.files("results/Aug3/SimCurves/", pattern = ".RDS", full.names = T)
fit_map <- readRDS("data/param_map/SimCurves/small_param_map.RDS") %>%
  dplyr::select(c("sim", "curve", "sens", "spec", "fatalitydata", "demog"))

nick <- fit_map$curve[[1]]


sim <- readRDS("results/Aug3/SimCurves/sim1.RDS")
sim$inputs$IFRmodel$data

#......................................................................
# Look at IFR through time
#......................................................................
dwnsmpl <- 1e3
#............................................................
# infection curve
#...........................................................
infxncurve <- COVIDCurve::draw_posterior_infxn_cubic_splines(IFRmodel_inf = modout,
                                                             dwnsmpl = dwnsmpl,
                                                             by_chain = FALSE,
                                                             by_strata = FALSE)

plotleft <- infxncurve$plotObj +
  geom_line(data = infxns, aes(x = time, y = infxns), color = "#3BABFD", size = 1.1) +
  labs(caption = "") +
  theme(plot.title = element_blank(),
        axis.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 12),
        axis.text.x =  element_text(family = "Helvetica", angle = 45, hjust = 0.5, vjust = 0.5, size = 11),
        axis.text.y =  element_text(family = "Helvetica", hjust = 0.5, vjust = 0.5, size = 11),
        axis.ticks = element_line(color = "#000000"),
        plot.caption = element_text(family = "Helvetica", hjust = 0.5, vjust = 0.5, size = 10), # so spacing equal
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())

#............................................................
# IFR over time
#...........................................................
seroprevdat <-  dat$AggSeroPrev %>%
  dplyr::rename(ObsDay = event_obs_day)

cumdat <- dplyr::left_join(x = dat$AggDeath, y = seroprevdat, by = c("ObsDay", "Strata"))
cumdat <- cumdat %>%
  dplyr::filter(Strata == "ma6") %>%
  dplyr::mutate(
    cumdeaths = cumsum(Deaths),
    RGIFR = cumdeaths/TrueSeroCount,
    CrudeIFR = cumdeaths/(popN * ObsPrev)
  )


#......................
# get posterior IFRs
#......................
mcmcout.nodes <- modout$mcmcout$output
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
  dplyr::select(c("iteration", "ma6"))

# truth
fatalitydata_intercept <- fatalitydata %>%
  dplyr::filter(Strata == "ma6") %>%
  dplyr::pull(c("IFR"))

#......................
# come together for plot right
#......................
plotright <- cumdat %>%
  dplyr::select(c("ObsDay", "RGIFR", "CrudeIFR")) %>%
  tidyr::gather(., key = "IFRlvl", value = "IFR", 2:ncol(.)) %>%
  dplyr::mutate(IFRlvl = factor(IFRlvl, levels = c("RGIFR", "CrudeIFR"), labels = c("Rogan-Gladen \n Adj.", "Crude"))) %>%
  ggplot() +
  geom_line(aes(x = ObsDay, y = IFR, color = IFRlvl), size = 1.1) +
  geom_hline(data = infIFR, aes(yintercept = ma6), color = "#d9d9d9", alpha = 0.5) +
  geom_hline(yintercept = fatalitydata_intercept, color = "#3BABFD", size = 1.2,
             linetype = "dashed", alpha = 0.9) +
  scale_color_manual("IFR Calc.", values = c("#F5390D", "#EDAC2C")) +
  labs(caption = c("Grey lines are posterior draws of the IFR; Blue dashed line is the true IFR")) +
  xlab("Time") +
  xyaxis_plot_theme +
  theme(plot.caption = element_text(family = "Helvetica", hjust = 0.5, vjust = 0.5, size = 10),
        axis.text.x =  element_text(family = "Helvetica", angle = 45, hjust = 0.5, vjust = 0.5, size = 11))


#............................................................
# main plot
#...........................................................
main_plotObj <- cowplot::plot_grid(plotleft, plotright, ncol = 2, nrow = 1)
jpeg("~/Desktop/temp.jpg", width = 8, height = 6, units = "in", res = 300)
main_plotObj
graphics.off()


