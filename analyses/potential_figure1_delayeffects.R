####################################################################################
## Purpose: Simple plot so show the delay of onsets effects
##
## Notes:
####################################################################################
set.seed(1234)
library(COVIDCurve)
library(tidyverse)

# make infxns from exp growth
infxns <- data.frame(time = 1:100)
offset <- 50
expgrowth <- function(r0, k, t){r0*exp(k*t)}
infxns$infxns <- sapply((infxns$time + offset), expgrowth, r0 = 1.05, k = 0.095)
infxns$infxns[76:100] <- 0


cumsum(infxns$infxns)

# make up fatality data
popN <- 3e6
fatalitydata <- tibble::tibble(Strata = "ma1",
                               IFR = 0.2,
                               Rho = 1,
                               Ne = 1)
demog <- tibble::tibble(Strata = "ma1",
                        popN = popN)

#..................
# run sim
#..................
dat <- COVIDCurve::Aggsim_infxn_2_death(
  fatalitydata = fatalitydata,
  demog = demog,
  m_od = 18.8,
  s_od = 0.45,
  curr_day = 100,
  level = "Time-Series",
  infections = infxns$infxns,
  simulate_seroprevalence = TRUE,
  sens = 0.86,
  spec = 0.99,
  sero_delay_rate = 10)

#......................
# tidy up and combine
#......................
cuminxns <- infxns %>%
  dplyr::mutate(cumincidence = cumsum(infxns)/popN) %>%
  dplyr::select(-c("infxns")) %>%
  dplyr::filter(time <= 75)


cumdeaths <- dat$AggDat %>%
  dplyr::mutate(cumDeaths = cumsum(Deaths)/popN) %>%
  dplyr::select(-c("Deaths", "Strata")) %>%
  dplyr::rename(time = ObsDay)


serodf <- dat$seroprev %>%
  dplyr::select(c("event_obs_day", "TruePrev", "ObsPrev")) %>%
  dplyr::rename(time = event_obs_day) %>%
  dplyr::filter(time <= 75)

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



jpeg("figures/temp_fig_1.jpg", width = 11, height = 8, units = "in", res = 500)
ggplot() +
  geom_vline(xintercept = 75, color = "#a50f15", linetype = "dashed", size = 0.75, alpha = 0.8) +
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

graphics.off()





# sanity

