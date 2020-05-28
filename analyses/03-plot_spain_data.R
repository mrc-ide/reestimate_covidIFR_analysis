####################################################################################
## Purpose: Plot Spain Data
##
## Author: Nick Brazeau
##
## Date: 27 May, 2020
####################################################################################
library(tidyverse)
source("R/summarize_mcmc_infxns.R")

#............................................................
# read in data
#...........................................................
r_mcmc_out.ageband <- readRDS("data/derived/ESP/ESP_mcmc_agebands_25rungs_5chains_2GTI.rds")
mod1 <- readRDS("data/derived/ESP/ESP_modinf_agebands.rds")

#............................................................
# plots
#...........................................................
plot_par(r_mcmc_out.ageband, "r1")
plot_par(r_mcmc_out.ageband, "r2")
plot_par(r_mcmc_out.ageband, "r3")
plot_par(r_mcmc_out.ageband, "r4")
plot_par(r_mcmc_out.ageband, "r5")
plot_par(r_mcmc_out.ageband, "r6")
plot_par(r_mcmc_out.ageband, "r7")
plot_par(r_mcmc_out.ageband, "r8")
plot_par(r_mcmc_out.ageband, "r9")
plot_par(r_mcmc_out.ageband, "ma10")
plot_par(r_mcmc_out.ageband, "y1")
plot_par(r_mcmc_out.ageband, "y2")
plot_par(r_mcmc_out.ageband, "y3")
plot_par(r_mcmc_out.ageband, "y4")
plot_par(r_mcmc_out.ageband, "y5")
plot_par(r_mcmc_out.ageband, "y6")
plot_par(r_mcmc_out.ageband, "y7")
plot_par(r_mcmc_out.ageband, "y8")
plot_par(r_mcmc_out.ageband, "y9")
plot_par(r_mcmc_out.ageband, "y10")

plot_par(r_mcmc_out.ageband, "sens")
plot_par(r_mcmc_out.ageband, "spec")
plot_par(r_mcmc_out.ageband, "sero_date")

plot_mc_acceptance(r_mcmc_out.ageband)
plot_rung_loglike(r_mcmc_out.ageband, y_axis_type = 1)
plot_rung_loglike(r_mcmc_out.ageband, y_axis_type = 2)
plot_rung_loglike(r_mcmc_out.ageband, y_axis_type = 3)


# infxn curve
r_mcmc_out.ageband.infxncurve <- get_infxn_curve(rmcmcout = r_mcmc_out.ageband,
                                                 modinf = mod1,
                                                 CIquant = 0.9)
r_mcmc_out.ageband.infxncurve$plotObj

plotObj <- r_mcmc_out.ageband.infxncurve$plotdata %>%
  dplyr::filter(time < 125) %>%
  ggplot() +
  geom_line(mapping = aes(time, infxns, group = sim), alpha = 0.25,
            lwd = 0.5, color = "#deebf7") +
  geom_vline(xintercept = mod1$knots, color = "#cb181d", lwd = 0.25, linetype = "dashed", alpha = 0.8) +
  xlab("Time") + ylab("Num. Infxns")  +
  labs(title = "Posterior Draws of the Infection Curve", subtitle = "Red lines are position of knots") +
  theme_minimal() +
  theme(
    plot.title = element_text(family = "Helvetica", face = "bold", vjust = 0.5,  hjust = 0.5, size = 18),
    plot.subtitle = element_text(family = "Helvetica", face = "bold", vjust = 0.5,  hjust = 0.5, size = 12),
    axis.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, vjust = 0.5, size = 16),
    axis.text.x = element_text(family = "Helvetica", angle = 45, hjust = 0.5, vjust = 0.5, size = 15),
    axis.text.y = element_text(family = "Helvetica", hjust = 0.5, vjust = 0.5, size = 15),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(color = "#000000", size = 1.2),
    legend.position = "none")

# summary plots
r_mcmc_out.ageband.summary <- get_param_summaries(rmcmcout = r_mcmc_out.ageband,
                                                  modinf = mod1)

#......................
# make summary plots for IFR
#......................
liftover_table <- data.frame(agebands = c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-999"),
                             param = c(paste0("r", 1:9), "ma10"))
agebandsPlotObj <- r_mcmc_out.ageband.summary$IFRparams %>%
  dplyr::left_join(., liftover_table, by = "param") %>%
  ggplot() +
  geom_pointrange(aes(x = agebands, ymin = LCI, ymax = UCI,
                      y = median, color = chain),
                  alpha = 0.5, size = 1.5) +
  scale_color_viridis_d("Chain") +
  xlab("Age Bands") + ylab("IFR") +
  theme(axis.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 12),
        axis.text.x = element_text(family = "Helvetica", hjust = 1, size = 11, angle = 45),
        axis.text.y = element_text(family = "Helvetica", hjust = 0.5, size = 11),
        legend.position = "bottom",
        legend.title = element_text(family = "Helvetica", face = "bold", vjust = 0.85, size = 12),
        legend.text = element_text(family = "Helvetica", hjust = 0.5, vjust = 0.5, size = 10),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        panel.grid.minor.y = element_line(size = 0.5, linetype = "dashed", color = "#bdbdbd"),
        panel.border = element_blank(),
        axis.line = element_line(color = "#000000", size = 1))


jpeg("~/Desktop/preliminary_esp_result.jpg", units = "in", height = 8, width = 11, res = 500)
cowplot::plot_grid(agebandsPlotObj, plotObj,
                   labels = c("(A)", "(B)"), nrow=2, ncol=1)
graphics.off()



summ <- r_mcmc_out.ageband.summary$IFRparams %>%
  dplyr::left_join(., liftover_table, by = "param")
summ %>%
  dplyr::mutate_if(is.numeric, round, 2) %>%
readr::write_csv(., "~/Desktop/preliminary_esp_dat_results.csv")

