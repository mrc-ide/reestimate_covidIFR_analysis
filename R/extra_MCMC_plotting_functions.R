#' @title Quick Diagnostic Plot for MCMC Framework
#' @details Plot the MC acceptance as well as the chains that are likely to mix the slowest

library(ggplot2)
quick_mc_diagnostics <- function(modout) {
  # get mcaccplot
  mcaccplot <- drjacoby::plot_mc_acceptance(modout$mcmcout)
  # find max ma (will mix slowest)
  maxma <- modout$inputs$IFRmodel$maxMa
  maxmachain <- drjacoby::plot_par(modout$mcmcout, maxma, display = FALSE)
  spechain <- drjacoby::plot_par(modout$mcmcout, "spec", display = FALSE)
  serodaychain <- drjacoby::plot_par(modout$mcmcout, "sero_day", display = FALSE)

  maxmachain <- maxmachain[[1]][["trace"]] + theme(legend.position = "none")
  spechain <- spechain[[1]][["trace"]] + theme(legend.position = "none")
  serodaychain <- serodaychain[[1]][["trace"]] + theme(legend.position = "bottom")

  # out
  cowplot::plot_grid(mcaccplot, maxmachain,
                     spechain, serodaychain,
                     nrow = 2, rel_heights = c(1,1,1,3))

}
