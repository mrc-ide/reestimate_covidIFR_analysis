## ....................................................................................................
## Purpose: Plot for Figure Showing Delays and Inference Framework
##
## Notes:
## ....................................................................................................
set.seed(48)
library(COVIDCurve)
library(tidyverse)
source("R/covidcurve_helper_functions.R")
#............................................................
#---- Read in and Simulate Data #----
#...........................................................
infxn_shapes <- readr::read_csv("data/simdat/infxn_curve_shapes.csv")
interveneflat <- infxn_shapes$intervene
# note need more infxns for sensitivity to be apparent on conceptual diagrams
interveneflat <- interveneflat * 1.5
interveneflat <- c(interveneflat, round(seq(from = interveneflat[200],
                                            to = 10, length.out = 100)))

# read in fitted rate of seroreversion parameter
weibullparams <- readRDS("results/prior_inputs/weibull_params.RDS")
weibullparams$wscale <- weibullparams$wscale - 13.3 # account for delay in onset of symptoms to seroconversion

# make up fatality data
fatalitydata <- tibble::tibble(Strata = c("ma1"),
                               IFR = 0.1,
                               Rho = 1)
demog <- tibble::tibble(Strata = "ma1",
                        popN = 3e6)

# run COVIDCurve sims for no seroreversion and seroreversion
# using whole population data -- just for smooth curve/conceputal
dat <- COVIDCurve::Agesim_infxn_2_death(
  fatalitydata = fatalitydata,
  demog = demog,
  m_od = 19.8,
  s_od = 0.85,
  curr_day = 300,
  infections = interveneflat,
  simulate_seroreversion = FALSE,
  smplfrac = 1,
  sens = 0.85,
  spec = 0.95,
  sero_delay_rate = 18.3,
  return_linelist = FALSE)

# with seroreversion
serorevdat <- COVIDCurve::Agesim_infxn_2_death(
  fatalitydata = fatalitydata,
  demog = demog,
  m_od = 19.8,
  s_od = 0.85,
  curr_day = 300,
  infections = interveneflat,
  simulate_seroreversion = TRUE,
  sero_rev_shape = weibullparams$wshape,
  sero_rev_scale = weibullparams$wscale,
  smplfrac = 1,
  sens = 0.85,
  spec = 0.95,
  sero_delay_rate = 18.3,
  return_linelist = FALSE)

#............................................................
#----- Conceptual Diagram #-----
#...........................................................
# read-in in Panel A
panelA <- readRDS("figures/final_figures/weibull_survplot.RDS")  +
  theme(plot.margin = unit(c(0.05, 0.05, 0.05, 1),"cm"))

#............................................................
# Work for Panel B
#...........................................................
#......................
# tidy up and combine
#......................
cuminxns <-  tibble::tibble(time = 1:length(interveneflat),
                            infxns = interveneflat)%>%
  dplyr::mutate(cumprevalence = cumsum(infxns)/demog$popN)


cumdeaths <- dat$Agg_TimeSeries_Death %>%
  dplyr::mutate(cumDeaths = cumsum(Deaths)/demog$popN) %>%
  dplyr::rename(time = ObsDay)


serodf_norev <- dat$StrataAgg_Seroprev %>%
  dplyr::select(c("ObsDay", "TruePrev", "ObsPrev")) %>%
  dplyr::rename(time = ObsDay,
                regTruePrev = TruePrev,
                regObsPrev = ObsPrev)

serodf_rev <- serorevdat$StrataAgg_Seroprev %>%
  dplyr::select(c("ObsDay", "ObsPrev")) %>%
  dplyr::rename(time = ObsDay,
                revObsPrev = ObsPrev)

# combine
datdf <- dplyr::left_join(cumdeaths, cuminxns, by = "time") %>%
  dplyr::left_join(., serodf_norev, by = "time") %>%
  dplyr::left_join(., serodf_rev, by = "time")

# long
plotdatdf <- datdf %>%
  dplyr::select(-c("Deaths", "infxns")) %>%
  tidyr::pivot_longer(., cols = -c("time"), names_to = "datlevel", values_to = "prop") %>%
  dplyr::mutate(datlevel = factor(datlevel,
                                  levels = c("cumprevalence", "cumDeaths",
                                             "regTruePrev", "regObsPrev", "revObsPrev"),
                                  labels = c("Infected",
                                             "Deaths",
                                             "True Seroprev.", "Seroprev. without Serorev.",
                                             "Seroprev. with Serorev."
                                  )))

#......................
# labels and arrows
#......................
arrows <- tibble::tibble(
  lvl =  c("mod", "serocon", "sens", "spec", "serorev"),
  x =    c(87,    139,       250,    10,     290),
  xend = c(136,    160,        250,    10,        290),
  y =    c(0.02,    0.4,       0.715,   0,     0.725),
  yend = c(0.02,    0.4,       0.615,  0.05,   0.265)
)


labels <- tibble::tibble(
  lvl =    c("mod",       "serocon",    "sens",    "spec",  "serorev"),
  label =  c("I-D Delay", "I-S Delay",  "Sens.",   "Spec.",  "I-R Delay"),
  x =      c(172,          102,          225,       12,       262),
  y =      c(0.02,          0.4,         0.65,    0.09,      0.50),
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
  xlab("Time (Days)") + ylab("Proportion") +
  theme(plot.title = element_blank(),
        axis.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 12),
        axis.text = element_text(family = "Helvetica", hjust = 0.5, size = 11),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(family = "Helvetica", hjust = 0.5, vjust = 0.5, size = 12),
        legend.key=element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "#000000", size = 1))  +
  theme(plot.margin = unit(c(0.05, 0.05, 0.05, 1),"cm"))

delay_plotObj


#..........................
# Inset of Infxn Curve
#..........................

infxninset_plotObj <- ggplot() +
  geom_line(data = datdf, aes(x = time, y = infxns), color = "#252525", size = 1.1, linetype = "dashed") +
  xlab("Time (days)") + ylab("Daily Infections") +
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
# Panel B row figure
#......................
delay_plotObj_nolegend <- delay_plotObj +
  theme(legend.position = "none")
legend <- cowplot::get_legend(delay_plotObj)
(panelB <- cowplot::ggdraw() +
    cowplot::draw_plot(delay_plotObj_nolegend, x = 0, y = 0, width = 1, height = 1, scale = 1) +
    cowplot::draw_plot(infxninset_plotObj, x = 0.2, y= 0.65, width = 0.3, height = 0.3))


#......................
# bring together
#......................
mainfig <- cowplot::plot_grid(panelA, panelB, labels = c("(A)", "(B)"),
                               nrow = 1, rel_widths = c(0.6, 1))
(mainfig <- cowplot::plot_grid(mainfig, legend, ncol = 1, rel_heights = c(1, 0.1)))


dir.create("figures/final_figures/", recursive = TRUE)
jpeg("figures/final_figures/Fig_concept_diagram.jpg",
     height = 8, width = 11, units = "in", res = 800)
plot(mainfig)
graphics.off()







