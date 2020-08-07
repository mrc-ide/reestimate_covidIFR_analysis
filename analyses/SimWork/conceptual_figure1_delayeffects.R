####################################################################################
## Purpose: Simple plot so show the delay of onsets effects
##
## Notes:
####################################################################################
set.seed(48)
library(COVIDCurve)
library(tidyverse)
source("R/simple_seir_model.R")

#......................
# run simple SEIR
#......................
# make infxns from exponential growth followed by intervetions
nsims <- 100
popN <- 3e6
infxns <- lapply(1:nsims, function(x){
  run_simple_seir(N = popN,
                  E0 = 50,
                  R0 = 0,
                  betas = c(0.33, 0.13, 0.12, 0.11),
                  beta_changes = c(1, 130, 140, 150),
                  sigma = 0.2,
                  gamma = 0.2,
                  time = 300)
})
infxns <- infxns %>%
  dplyr::bind_rows(.) %>%
  dplyr::group_by(step) %>%
  dplyr::summarise(
    infxns = mean(I)
  ) %>%
  dplyr::mutate_if(is.numeric, round, 0) %>%
  dplyr::rename(time = step)


# make up fatality data
fatalitydata <- tibble::tibble(Strata = "ma1",
                               IFR = 0.2,
                               Rho = 1,
                               Ne = 1)
demog <- tibble::tibble(Strata = "ma1",
                        popN = popN)

#..................
# run COVIDCurve sim
#..................
dat <- COVIDCurve::Aggsim_infxn_2_death(
  fatalitydata = fatalitydata,
  demog = demog,
  m_od = 14.26,
  s_od = 0.19,
  curr_day = 300,
  infections = infxns$infxns,
  simulate_seroprevalence = TRUE,
  sens = 0.86,
  spec = 0.95,
  sero_delay_rate = 13.3)

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
# plot & inset
#......................
delay_plotObj <- plotdatdf %>%
  dplyr::filter(datlevel != "True Prev.") %>% # discuss wheter true goes into figure or not
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



infxninet_plotObj <- infxns %>%
  dplyr::filter(infxns != -1) %>%  # remove missing
  dplyr::mutate(incidence = infxns/popN) %>%
  ggplot() +
  geom_line(aes(x = time, y = incidence), color = "#252525", size = 1.1, linetype = "dashed") +
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
jpeg("figures/conceptual_fig_1.jpg", width = 11, height = 8, units = "in", res = 500)
cowplot::ggdraw() +
  cowplot::draw_plot(delay_plotObj, x = 0, y = 0, width = 1, height = 1, scale = 1) +
  cowplot::draw_plot(infxninet_plotObj, x = 0.085, y= 0.75, width = 0.3, height = 0.2)
graphics.off()


