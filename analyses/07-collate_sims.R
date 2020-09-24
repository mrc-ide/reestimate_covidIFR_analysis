##................................................................................................
## Purpose: Collate Simulations
##
## Notes:
##................................................................................................
library(tidyverse)
library(COVIDCurve)
source("R/my_themes.R")

#............................................................
# lightweight function to make plot
#...........................................................
get_sim_IFR_curve <- function(sim_param_map, modout) {
  trueseroprev <- sim_param_map %>%
    dplyr::pull("simdat")
  trueseroprev <- trueseroprev[[1]]$StrataAgg_Seroprev


  infxncurve <- COVIDCurve::draw_posterior_infxn_cubic_splines(IFRmodel_inf = modout,
                                                               dwnsmpl = 1e2,
                                                               by_chain = TRUE,
                                                               by_strata = FALSE)$curvedata
  infxncurve %>%
    dplyr::select(c("chain", "sim", "time", "totinfxns")) %>%
    ggplot() +
    geom_line(aes(x = time, y = totinfxns, group = sim, color = chain), alpha = 0.8, size = 0.8) +
    geom_line(data = truecurve, aes(x = time, y = infxns), color = "#000000", linetype = "dashed", size = 1.1) +
    scale_color_viridis_d("Chain") +
    labs(caption = "Viridis is inferred, Grey Line is Simulated Infection Curve") +
    xlab("Time") + ylab("Num. Infxns") +
    xyaxis_plot_theme
}


#............................................................
# read in data
#...........................................................
sim_noserorev_parammap <- readRDS("data/param_map/SimCurves_noserorev/simfit_param_map.RDS")
sim_serorev_parammap <- readRDS("data/param_map/SimCurves_noserorev/simfit_param_map.RDS")
# no serorev
noserorev_paths <- list.files("results/SimCurves_noserorev/", full.names = TRUE)
noserorev_paths <- tibble::tibble(sim = gsub("_NoSeroRev.RDS", "", basename(noserorev_paths)),
                                  path = noserorev_paths) %>%
  dplyr::mutate(modout = purrr::map(path, readRDS)) %>%
  dplyr::left_join(sim_noserorev_parammap, ., by = "sim") %>%
  dplyr::mutate(serorev_accounted = FALSE)

# serorev
serorev_paths <- list.files("results/SimCurves_serorev/", full.names = TRUE)
serorev_paths <- tibble::tibble(sim = gsub("_SeroRev.RDS", "", basename(serorev_paths)),
                                path = serorev_paths) %>%
  dplyr::mutate(modout = purrr::map(path, readRDS)) %>%
  dplyr::left_join(sim_serorev_parammap, ., by = "sim") %>%
  dplyr::mutate(serorev_accounted = TRUE)

#.....................
# come together
#......................
datmap <- dplyr::bind_rows()
