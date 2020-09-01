library(tidyverse)
library(drjacoby)
source("R/my_themes.R")
source("R/extra_plotting_functions.R")

#............................................................
# data wrangle
#...........................................................
modoutpath <- list.files("results/SimCurves/", pattern = ".RDS", full.names = T)
modoutpath <- modoutpath[!grepl("simfit_param_map.RDS", modoutpath)]
modoutpath <- tibble::tibble(sim = sub(".RDS", "", basename(modoutpath)),
                             path = modoutpath)
fit_map <- readRDS("results/SimCurves/simfit_param_map.RDS") %>%
  dplyr::select(c("sim", "nm", "simdat", "curve", "sens", "spec", "mod", "sero_rate", "fatalitydata", "demog")) %>%
  dplyr::left_join(., modoutpath, by = "sim")

fit_map$modout <- purrr::map(fit_map$path, readRDS)


#......................................................................
# Look at IFR through time
#......................................................................
analyze_fits <- function(sim, simdat, curve, modout, fatalitydata,
                         spec, sens, mod, sero_rate, dwnsmpl = 1e2) {

  # get title for throughout
  title <- cowplot::ggdraw() + cowplot::draw_label(
    paste0("Simulation run: ", sim),
    fontface = 'bold', hjust = 0.5, vjust = 0.5) +
    theme(plot.margin = margin(0, 5, 0, 0))

  #......................
  # like and mc accept
  #......................
  likeplot <- modout$mcmcout$output %>%
    dplyr::filter(rung == "rung1") %>%
    dplyr::filter(stage == "sampling") %>%
    dplyr::select(c("chain", "iteration", "loglikelihood", "logprior")) %>%
    tidyr::gather(., key = "like", value = "val", 3:4) %>%
    ggplot() +
    geom_line(aes(x = iteration, y = val, color = chain), size = 0.25, alpha = 0.8) +
    scale_color_viridis_d() +
    ylab("") + xlab("Iteration") +
    facet_wrap(.~like, scales = "free_y") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          legend.title = element_blank())
  mcaccplot <- drjacoby::plot_mc_acceptance(modout$mcmcout)

  diagplot <- cowplot::plot_grid(title, likeplot, mcaccplot, ncol = 1, rel_heights = c(0.1, 1, 0.4))

  #......................
  # sero mixing
  #......................
  maxma <- modout$inputs$IFRmodel$maxMa
  maxmachain <- drjacoby::plot_par(modout$mcmcout, maxma, display = FALSE)
  spechain <- drjacoby::plot_par(modout$mcmcout, "spec", display = FALSE)
  senschain <- drjacoby::plot_par(modout$mcmcout, "sens", display = FALSE)
  modchain <- drjacoby::plot_par(modout$mcmcout, "mod", display = FALSE)
  sodchain <- drjacoby::plot_par(modout$mcmcout, "sod", display = FALSE)
  seroratechain <- drjacoby::plot_par(modout$mcmcout, "sero_rate", display = FALSE)

  maxmachain <- maxmachain[[1]][["trace"]] + theme(legend.position = "none")
  spechain <- spechain[[1]][["trace"]] + geom_hline(yintercept = spec, linetype = "dashed", size = 1.25) +
    labs(caption = paste0("Prior was Beta(",
                         paste(modout$inputs$IFRmodel$paramdf[modout$inputs$IFRmodel$paramdf$name == "spec", c("dsc1", "dsc2")], collapse = ","),
                         ")")) +
    theme(legend.position = "none")
  senschain <- senschain[[1]][["trace"]] + geom_hline(yintercept = sens, linetype = "dashed", size = 1.25) +
    labs(caption = paste0("Prior was Beta(",
                          paste(modout$inputs$IFRmodel$paramdf[modout$inputs$IFRmodel$paramdf$name == "sens", c("dsc1", "dsc2")], collapse = ","),
                          ")")) +
    theme(legend.position = "none")
  modchain <- modchain[[1]][["trace"]] + geom_hline(yintercept = mod, linetype = "dashed", size = 1.25) +
    labs(caption = paste0("Prior was Norm+(",
                          paste(modout$inputs$IFRmodel$paramdf[modout$inputs$IFRmodel$paramdf$name == "mod", c("dsc1", "dsc2")], collapse = ","),
                          ")")) +
    theme(legend.position = "none")
  sodchain <- sodchain[[1]][["trace"]] +
    labs(caption = paste0("Prior was Beta(",
                          paste(modout$inputs$IFRmodel$paramdf[modout$inputs$IFRmodel$paramdf$name == "sod", c("dsc1", "dsc2")], collapse = ","),
                          ")")) +
    theme(legend.position = "none")
  seroratechain <- seroratechain[[1]][["trace"]] + geom_hline(yintercept = sero_rate, linetype = "dashed", size = 1.25) +
    labs(caption = paste0("Prior was Norm+(",
                          paste(modout$inputs$IFRmodel$paramdf[modout$inputs$IFRmodel$paramdf$name == "sero_rate", c("dsc1", "dsc2")], collapse = ","),
                          ")")) +
    theme(legend.position = "none")
  # get legend
  legend_bt <- cowplot::get_legend(maxmachain + theme(legend.position = "bottom",
                                                      legend.title = element_blank()))
  # out
  topp <- cowplot::plot_grid(spechain, senschain, modchain, sodchain, maxmachain, seroratechain,
                             ncol = 2, nrow = 3)
  sero_plot_dets <- cowplot::plot_grid(title, topp, legend_bt, nrow = 3, rel_heights = c(0.1, 1, 0.1))

  #.......................
  # infection curve
  #.......................
  infxncurve <- COVIDCurve::draw_posterior_infxn_cubic_splines(IFRmodel_inf = modout,
                                                               dwnsmpl = dwnsmpl,
                                                               by_chain = TRUE,
                                                               by_strata = TRUE)$curvedata
  pA <- infxncurve %>%
    dplyr::select(c("chain", "sim", "time", dplyr::starts_with("infxns_"))) %>%
    tidyr::pivot_longer(., cols =  dplyr::starts_with("infxns_"),
                        names_to = "strata", values_to = "infxns") %>%
    dplyr::mutate(strata = sub("infxns_", "", strata)) %>%
    ggplot() +
    geom_line(aes(x = time, y = infxns, group = sim, color = chain), alpha = 0.8, size = 0.8) +
    geom_line(data = curve, aes(x = time, y = infxns), color = "#6baed6", linetype = "dashed", size = 1.1) +
    scale_color_viridis_d("Chain") +
    xlab("Time") + ylab("Num. Infxns") +
    xyaxis_plot_theme

  pB <- infxncurve %>%
    dplyr::select(c("chain", "sim", "time", "spec", dplyr::starts_with("infxns_"))) %>%
    tidyr::pivot_longer(., cols =  dplyr::starts_with("infxns_"),
                        names_to = "strata", values_to = "infxns") %>%
    dplyr::mutate(strata = sub("infxns_", "", strata)) %>%
    ggplot() +
    geom_line(aes(x = time, y = infxns, group = sim, color = spec), alpha = 0.8, size = 0.8) +
    geom_line(data = curve, aes(x = time, y = infxns), color = "#6baed6", linetype = "dashed", size = 1.1) +
    scale_color_viridis_c("Spec.") +
    xlab("Time") + ylab("Num. Infxns") +
    labs(caption = paste("True Specificity is", spec)) +
    xyaxis_plot_theme

  infxnplot <- cowplot::plot_grid(pA, pB, nrow = 1, align = "h")
  infxnplot <- cowplot::plot_grid(title, infxnplot, nrow = 2, rel_heights = c(0.1, 1))


  #..................
  # out
  #..................
  out <- list(diagplot = diagplot,
              sero_plot_dets = sero_plot_dets,
              infxnplot = infxnplot)
  return(out)
}


#......................
# run analysis function of fits
#......................
fit_map$retplots <- purrr::pmap(fit_map[, c("sim", "simdat", "curve", "modout", "fatalitydata",
                                            "spec", "sens", "mod", "sero_rate")],
                                analyze_fits, dwnsmpl = 1e2)

# combine plots
fit_map_sm <- fit_map %>%
  dplyr::mutate(diagplot = purrr::map(retplots, "diagplot"),
                sero_plot_dets = purrr::map(retplots, "sero_plot_dets"),
                infxnplot = purrr::map(retplots, "infxnplot")) %>%
  dplyr::select(c("sim", "nm", "spec", "sens", "mod", "sero_rate", "diagplot", "sero_plot_dets", "infxnplot"))


fit_map_sm_diagplots <- fit_map_sm %>%
  dplyr::select(-c("sero_plot_dets", "infxnplot")) %>%
  dplyr::rename(plotObjs = diagplot)
fit_map_sm_sero_plot_dets <- fit_map_sm %>%
  dplyr::select(-c("infxnplot", "diagplot")) %>%
  dplyr::rename(plotObjs = sero_plot_dets)
fit_map_sm_infxnplot <- fit_map_sm %>%
  dplyr::select(-c("sero_plot_dets", "diagplot")) %>%
  dplyr::rename(plotObjs = infxnplot)

fit_map_sm <- dplyr::bind_rows(fit_map_sm_diagplots, fit_map_sm_sero_plot_dets, fit_map_sm_infxnplot) %>%
  dplyr::arrange(sim, spec, sens, mod, sero_rate)

# save out
dir.create("figures/SimCurves/PowerAnalysis", recursive = T)
pdf("figures/SimCurves/PowerAnalysis/Fits_comparisons.pdf", width = 8, height = 11)
fit_map_sm$plotObjs
graphics.off()
