#' @title Quick Diagnostic Plot for MCMC Framework
#' @details Plot the MC acceptance as well as the chains that are likely to mix the slowest

library(ggplot2)
quick_mc_diagnostics <- function(modout) {
  # get mcaccplot
  mcaccplot <- drjacoby::plot_mc_acceptance(modout$mcmcout)
  mcacclogplot1 <- drjacoby::plot_rung_loglike(modout$mcmcout, x_axis_type = 2, y_axis_type = 2)
  mcacclogplot2 <- drjacoby::plot_rung_loglike(modout$mcmcout, x_axis_type = 2, y_axis_type = 3)
  # find max ma (will mix slowest)
  maxma <- modout$inputs$IFRmodel$maxMa
  maxmachain <- drjacoby::plot_par(modout$mcmcout, maxma, display = FALSE)
  spechain <- drjacoby::plot_par(modout$mcmcout, "spec", display = FALSE)
  serodaychain <- drjacoby::plot_par(modout$mcmcout, "sero_day", display = FALSE)
  seroratechain <- drjacoby::plot_par(modout$mcmcout, "sero_rate", display = FALSE)

  maxmachain <- maxmachain[[1]][["trace"]] + theme(legend.position = "none")
  spechain <- spechain[[1]][["trace"]] + theme(legend.position = "none")
  serodaychain <- serodaychain[[1]][["trace"]] + theme(legend.position = "none")
  seroratechain <- seroratechain[[1]][["trace"]] + theme(legend.position = "bottom")
  # out
  lftside <- cowplot::plot_grid(mcaccplot, mcacclogplot1, mcacclogplot2,
                                nrow = 3)
  rightside <- cowplot::plot_grid(maxmachain, spechain, serodaychain, seroratechain,
                                  nrow = 4, rel_heights = c(1,1,1,1.5))
  cowplot::plot_grid(lftside, rightside, ncol = 2)

}
