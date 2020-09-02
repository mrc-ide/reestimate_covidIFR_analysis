#........................................................................................
## Purpose: Use Drake to fit IFR Models
##
## Notes:
#........................................................................................
library(tidyverse)
library(COVIDCurve)
library(drjacoby)
source("R/my_themes.R")

#............................................................
# cpp liftover functions
#...........................................................
source("analyses/Rgn_Age_Mod/src/draw_rgn_age_mod_seroprev_posterior.R")


#...................................................................................
# Read In Rgn Fit Map
#.................................................................................
rgn_fit_map <- tibble::tibble(study_id = "ESP1-2",
                              rgnmod = "analyses/Rgn_Age_Mod/temp_ESP/mcmcout_newrgnl_mod.RDS",
                              agemod = "results/Sep1/ModFits/ESP1-2_age_rung50_burn10000_smpl10000.RDS") %>%
  dplyr::mutate(rgnmod = purrr::map(rgnmod, readRDS),
                agemod = purrr::map(agemod, readRDS))

#............................................................
# Get Seroprevalences for Regional Models
#...........................................................
rgn_fit_map$seroprev <- purrr::map(rgn_fit_map$rgnmod, RgnAgeMod_draw_posterior_sero_curves,
                                   whichrung = "rung1", dwnsmpl = 1e2, by_chain = FALSE)

#............................................................
# Compare Specificities
#...........................................................
plot_par(rgnmod, "spec")

plot_par(agemod$mcmcout, "spec")

pA <- gbr_curvedat %>%
  ggplot() +
  geom_line(aes(x = time, y = totinfxns, group = sim, color = lvl), alpha = 0.8) +
  scale_color_manual("Marginal Model \n Inference", values = c("#4285F4", "#EA4335")) +
  xlab("") + ylab("Tot. Infxns") +
  xyaxis_plot_theme



#............................................................
# Read in and get global IFR for single unit models (i.e. agebands summed up)
#...........................................................
ab_mod_fits <- readRDS("")
# use what we were doing in individual reports
