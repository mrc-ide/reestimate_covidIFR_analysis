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
                  time = 200)
})
infxns <- infxns %>%
  dplyr::bind_rows(.) %>%
  dplyr::group_by(step) %>%
  dplyr::summarise(
    infxns = mean(I)
  ) %>%
  dplyr::mutate_if(is.numeric, round, 0) %>%
  dplyr::rename(time = step)

#......................
# Trim Infection Curve
#......................
infxns <- infxns %>%
  dplyr::filter(time > 50) %>% # trim beginning
  dplyr::mutate(infxns = ifelse(time > 170, 0, infxns), # "stop" time to show seroprev and deaths continue
                time = 1:nrow(.)) # reset time



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
  m_od = 18.8,
  s_od = 0.45,
  curr_day = 150,
  level = "Time-Series",
  infections = infxns$infxns,
  simulate_seroprevalence = TRUE,
  sens = 0.86,
  spec = 0.95,
  sero_delay_rate = 10)

#......................
# tidy up and combine
#......................
infxns <- infxns %>%
  dplyr::mutate(infxns = ifelse(infxns == 0 & time > 50, -1, infxns))

cuminxns <- infxns %>%
  dplyr::filter(infxns != -1) %>%
  dplyr::mutate(cumincidence = cumsum(infxns)/popN) %>%
  dplyr::select(-c("infxns"))


cumdeaths <- dat$AggDat %>%
  dplyr::mutate(cumDeaths = cumsum(Deaths)/popN) %>%
  dplyr::select(-c("Deaths", "Strata")) %>%
  dplyr::rename(time = ObsDay)


serodf <- dat$seroprev %>%
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
                                  labels = c("Cum. Infxns", "True Prev.", "Obs. Prev", "Cum. Deaths"))
  )


#......................
# plot & inset
#......................
delay_plotObj <- ggplot() +
  geom_vline(xintercept = 120, color = "#a50f15", linetype = "dashed", size = 0.5, alpha = 0.8) +
  geom_line(data = plotdatdf, aes(x = time, y = prop, color = datlevel, size = datlevel)) +
  # geom_segment(data = reshape(DATA, v.names="VALUE", idvar = "NAME", timevar = "YEAR", direction = "wide"),
  #              aes(x=VALUE.2011, xend=VALUE.2016, y=NAME, yend=NAME), size = 2,
  #              arrow = arrow(length = unit(0.5, "cm")))
  scale_size_manual(values = c(1.3, 1.3, 0.9, 0.9)) +
  scale_color_manual(values = c("#252525", "#0868ac", "#74a9cf", "#969696")) +
  xlab("Time (days)") + ylab("Prop. of Population Affected") +
  theme(plot.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 14),
        axis.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 12),
        axis.text = element_text(family = "Helvetica", hjust = 0.5, size = 11),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(family = "Helvetica", hjust = 0.5, vjust = 0.5, size = 10),
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
jpeg("figures/temp_fig_1.jpg", width = 11, height = 8, units = "in", res = 500)
cowplot::ggdraw() +
  cowplot::draw_plot(delay_plotObj, x = 0, y = 0, width = 1, height = 1, scale = 1) +
  cowplot::draw_plot(infxninet_plotObj, x = 0.085, y= 0.75, width = 0.3, height = 0.2)
graphics.off()





# sanity

