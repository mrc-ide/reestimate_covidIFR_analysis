#........................................................................................
## Purpose: Use Drake to fit IFR Models
##
## Notes:
#........................................................................................
library(tidyverse)
library(COVIDCurve)
library(drjacoby)
source("R/my_themes.R")
source("R/covidcurve_helper_functions.R")

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
# Look at Rnes for RgnMod (param we are basically trying to infer)
#...........................................................

# TODO fix this / make proper function if we keep
liftover <- tibble::tibble(Rne = paste0("Rne", 1:length(rgn_fit_map$rgnmod[[1]]$rgnnames)),
                           region = rgn_fit_map$rgnmod[[1]]$rgnnames)

rgn_fit_map$rgnmod[[1]]$mcmcout$output %>%
  dplyr::filter(rung == "rung1" & stage == "sampling") %>%
  dplyr::select(c("iteration", dplyr::starts_with("Rne"))) %>%
  tidyr::pivot_longer(., cols = -c("iteration"),
                      names_to = "Rne", values_to = "est") %>%
  dplyr::group_by(Rne) %>%
  dplyr::summarise(
    min = min(est),
    LCI = quantile(est, 0.025),
    median = median(est),
    mean = mean(est),
    UCI = quantile(est, 0.975),
    max = max(est),
    ESS = coda::effectiveSize(coda::as.mcmc(est)),
    GewekeZ = coda::geweke.diag(coda::as.mcmc(est))[[1]],
    GewekeP = dnorm(GewekeZ)
  ) %>%
  dplyr::ungroup(.) %>%
  dplyr::left_join(., liftover, by = "Rne") %>%
  dplyr::arrange(median) %>%
  knitr::kable(.)


#............................................................
# Get Seroprevalences for Regional Models
#...........................................................
rgn_fit_map$seroprev <- purrr::map(rgn_fit_map$rgnmod, RgnAgeMod_draw_posterior_sero_curves,
                                   whichrung = "rung1", dwnsmpl = 1e2, by_chain = FALSE)

#............................................................
# Compare Specificities
#...........................................................
drjacoby::plot_par(rgn_fit_map$rgnmod[[1]]$mcmcout, "spec")
drjacoby::plot_par(rgn_fit_map$agemod[[1]]$mcmcout, "spec")

# get region spec
spec_rgn <- rgn_fit_map$rgnmod[[1]]$mcmcout$output %>%
  dplyr::filter(rung == "rung1" & stage == "sampling") %>%
  dplyr::select(c("iteration", "spec")) %>%
  dplyr::mutate(lvl = "rgn")

# get age spec
spec_age <- rgn_fit_map$agemod[[1]]$mcmcout$output %>%
  dplyr::filter(rung == "rung1" & stage == "sampling") %>%
  dplyr::select(c("iteration", "spec")) %>%
  dplyr::mutate(lvl = "age")

# bring together
specdat <- dplyr::bind_rows(spec_rgn, spec_age)

specPlotObj <- specdat %>%
  ggplot() +
  geom_histogram(aes(x = spec, y = ..density..,
                     group = lvl, fill = lvl),
                 alpha = 0.5, position = "identity") +
  scale_fill_manual("Model Level", values = c("#4285F4", "#EA4335")) +
  xlab("Specificity") + ylab("Density") +
  xyaxis_plot_theme
specPlotObj

#......................
# for regional models extract relevenat seroprevs
#......................
rgn_mod_seroprevs <- rgn_fit_map %>%
  dplyr::select(c("study_id", "seroprev")) %>%
  tidyr::unnest(cols = c("seroprev")) %>%
  tidyr::pivot_longer(., cols = -c("study_id", "ObsDay", "sim"),
                      names_to = "sero_type", values_to = "sero_est") %>%
  dplyr::mutate(strata = stringr::str_extract(sero_type, "R[0-9]+"),
                sero_type = stringr::str_split_fixed(sero_type, "_R[0-9]+", n = 2)[,1]) %>%
  dplyr::group_by(study_id, sero_type, strata, ObsDay) %>%
  dplyr::summarise(
    min = min(sero_est),
    LCI = quantile(sero_est, 0.025),
    median = median(sero_est),
    mean = mean(sero_est),
    UCI = quantile(sero_est, 0.975),
    max = max(sero_est),
    ESS = coda::effectiveSize(coda::as.mcmc(sero_est)),
    GewekeZ = coda::geweke.diag(coda::as.mcmc(sero_est))[[1]],
    GewekeP = dnorm(GewekeZ)
  ) %>%
  dplyr::ungroup(.)


#............................................................
# Read in and get global IFR for single unit models (i.e. agebands summed up)
#...........................................................
ab_mod_paths <- list.files("results/Modfits_noserorev/", pattern = ".RDS", full.names = TRUE)
ab_mod_fits <- tibble::tibble(
  study_id = toupper(stringr::str_split(basename(ab_mod_paths), "_age|_rgn", simplify = T)[,1]),
  path = ab_mod_paths) %>%
  dplyr::mutate(seroprev = purrr::map(path, get_seroprevs))

# get demog for weights
# TODO

age_mod_seroprevs <- ab_mod_fits %>%
  dplyr::select(c("study_id", "seroprev")) %>%
  tidyr::unnest(cols = c("seroprev")) %>%
  tidyr::pivot_longer(., cols = -c("study_id", "ObsDay", "sim"),
                      names_to = "sero_type", values_to = "sero_est") %>%
  dplyr::mutate(sero_type = stringr::str_split_fixed(sero_type, "_ma[0-9]+", n = 2)[,1]) %>%
  dplyr::group_by(study_id, sero_type, ObsDay) %>%
  dplyr::summarise(
    min = min(sero_est),
    LCI = quantile(sero_est, 0.025),
    median = median(sero_est),
    mean = mean(sero_est),
    UCI = quantile(sero_est, 0.975),
    max = max(sero_est),
    ESS = coda::effectiveSize(coda::as.mcmc(sero_est)),
    GewekeZ = coda::geweke.diag(coda::as.mcmc(sero_est))[[1]],
    GewekeP = dnorm(GewekeZ)
  )

#............................................................
# Read in and get crude IFR
#...........................................................
# TODO look at how we are determining cum deaths
# descriptive data from our study
dscdat <- readRDS("data/derived/descriptive_results_datamap.RDS")
plotdat <- dscdat %>%
  dplyr::select(c("study_id", "plotdat")) %>%
  tidyr::unnest(cols = "plotdat")


# TODO fix this / make proper function if we keep
liftover <- tibble::tibble(strata = paste0("R", 1:length(rgn_fit_map$rgnmod[[1]]$rgnnames)),
                           region = rgn_fit_map$rgnmod[[1]]$rgnnames)

rgn_mod_seroprevs <- dplyr::left_join(rgn_mod_seroprevs, liftover, by = "strata")

# TODO
midpt <- max(plotdat$seromidpt)
rgn_mod_seroprevs_mid <- rgn_mod_seroprevs %>%
  dplyr::filter(ObsDay == midpt)


plotdat_mid <- plotdat %>%
  dplyr::filter(obsday == midpt)
#......................
# plot
#......................
rgn_mod_seroprevs_mid %>%
  dplyr::left_join(., plotdat_mid) %>%
  dplyr::filter(sero_type == "RG_pd_seroprev") %>%
  # TODO fix this -- this should not be happening -- why is plot dat duplicating
  dplyr::select(c("region", "std_cum_deaths", "median", "LCI", "UCI")) %>%
  dplyr::filter(!duplicated(.)) %>%
  ggplot() +
  geom_pointrange(aes(x = std_cum_deaths, y = median, ymin = LCI, ymax = UCI)) +
  ggrepel::geom_text_repel(aes(x = std_cum_deaths, y = median, label = region)) +
  xlim(c(0, 2000)) + ylim(c(0, 0.1)) +
  coord_flip() +
  ylab("Inferred True Seroprevalence") + xlab("Deaths per Million") +
  xyaxis_plot_theme
