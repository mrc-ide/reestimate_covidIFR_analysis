####################################################################################
## Purpose: Plot Spain Data
##
## Author: Nick Brazeau
##
## Date: 27 May, 2020
####################################################################################
library(drjacoby)
library(tidyverse)

#............................................................
# read in data
#...........................................................
r_mcmc_out.ageband <- readRDS("data/derived/ESP/ESP_mcmc_agebands_fit.rds")

#............................................................
# plots
#...........................................................
plot_par(r_mcmc_out.ageband$mcmcout, "ma1")
plot_par(r_mcmc_out.ageband$mcmcout, "ma2")
plot_par(r_mcmc_out.ageband$mcmcout, "ma3")
plot_par(r_mcmc_out.ageband$mcmcout, "ma4")
plot_par(r_mcmc_out.ageband$mcmcout, "ma5")
plot_par(r_mcmc_out.ageband$mcmcout, "ma6")
plot_par(r_mcmc_out.ageband$mcmcout, "ma7")
plot_par(r_mcmc_out.ageband$mcmcout, "ma8")
plot_par(r_mcmc_out.ageband$mcmcout, "ma9")
plot_par(r_mcmc_out.ageband$mcmcout, "ma10")
plot_par(r_mcmc_out.ageband$mcmcout, "y1")
plot_par(r_mcmc_out.ageband$mcmcout, "y2")
plot_par(r_mcmc_out.ageband$mcmcout, "y3")
plot_par(r_mcmc_out.ageband$mcmcout, "y4")
plot_par(r_mcmc_out.ageband$mcmcout, "y5")
plot_par(r_mcmc_out.ageband$mcmcout, "y6")
plot_par(r_mcmc_out.ageband$mcmcout, "y7")
plot_par(r_mcmc_out.ageband$mcmcout, "y8")
plot_par(r_mcmc_out.ageband$mcmcout, "y9")
plot_par(r_mcmc_out.ageband$mcmcout, "y10")
plot_par(r_mcmc_out.ageband$mcmcout, "x1")
plot_par(r_mcmc_out.ageband$mcmcout, "x2")
plot_par(r_mcmc_out.ageband$mcmcout, "x3")
plot_par(r_mcmc_out.ageband$mcmcout, "x4")
plot_par(r_mcmc_out.ageband$mcmcout, "x5")
plot_par(r_mcmc_out.ageband$mcmcout, "x6")
plot_par(r_mcmc_out.ageband$mcmcout, "x7")
plot_par(r_mcmc_out.ageband$mcmcout, "x8")
plot_par(r_mcmc_out.ageband$mcmcout, "x9")
plot_par(r_mcmc_out.ageband$mcmcout, "x10")

plot_par(r_mcmc_out.ageband$mcmcout, "sens")
plot_par(r_mcmc_out.ageband$mcmcout, "spec")
plot_par(r_mcmc_out.ageband$mcmcout, "sero_day")


# cred intervals
(ifr <- COVIDCurve::get_cred_intervals(IFRmodel_inf = r_mcmc_out.ageband,
                                       whichrung = paste0("rung", 1),
                                       what = "IFRparams", by_chain = F))

# infxn curve
r_mcmc_out.ageband.infxncurve <- COVIDCurve::draw_posterior_infxn_points_cubic_splines(IFRmodel_inf = r_mcmc_out.ageband,
                                                                                       CIquant = 0.95,
                                                                                       by_chain = FALSE)
postdat <- COVIDCurve::posterior_check_infxns_to_death(IFRmodel_inf = r_mcmc_out.ageband,
                                                       CIquant = 0.95,
                                                       by_chain = FALSE)
#............................................................
# make summary plots for curves and IFR
#...........................................................
# plot out

liftover_table <- data.frame(ageband = c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-999"),
                             param = c(paste0("ma", 1:9), "ma10"))
jpeg("~/Desktop/SPAIN_posterior_curve_draws.jpg", width = 11, height = 8, units = "in", res = 500)

infxnfatalitydataplot <- ifr %>%
  dplyr::left_join(liftover_table, ., by = "param") %>%
  dplyr::mutate(ageband = factor(ageband, levels =  c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-999")))

plot1 <- ggplot() +
  geom_pointrange(data = infxnfatalitydataplot, aes(x = ageband, ymin = LCI, ymax = UCI, y = median, color = ageband)) +
  scale_color_viridis_d() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust= 0.5, face = "bold"))

plot2 <- r_mcmc_out.ageband.infxncurve$plotObj

cowplot::plot_grid(plot1, plot2, ncol = 1, nrow = 2)

graphics.off()




#......................
# get deaths posterior pred check
#......................
liftover_table <- data.frame(ageband = c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-999"),
                             param = c(paste0("ma", 1:9), "ma10")) %>%
  dplyr::mutate(param = paste0("deaths_", param))
datclean <-  r_mcmc_out.ageband$inputs$IFRmodel$data$obs_deaths %>%
  dplyr::filter(Deaths != -1)

postdat_long <- postdat %>%
  dplyr::select(c("sim", "time", dplyr::starts_with("deaths"))) %>%
  tidyr::gather(., key = "param", value = "deaths", 3:ncol(.)) %>%
  dplyr::left_join(., y = liftover_table, by = "param") %>%
  dplyr::filter(sim %in% c(1:1e3)) # downsample sims


jpeg("~/Desktop/post_dat_spain.jpg", width = 11, height = 8, units = "in", res = 500)
ggplot() +
  geom_line(data = postdat_long, aes(x= time, y = deaths, group = ageband, color = ageband), size = 1.2) +
  scale_color_viridis_d() +
  geom_line(data = datclean, aes(x=ObsDay, y = Deaths, group = ageband), color = "#bdbdbd") +
  facet_wrap(.~ageband) +
  theme_bw() +
  ggtitle("Posterior Predictive Check", subtitle = "Grey Lines are ECDC Data, Viridis Lines are Draws from Posterior") +
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, vjust = 0.5))
graphics.off()

