source(paste0(here::here(), "/R/assertions_v5.R"))

#' @title Make a nice DT table from a dataframe
#' @import dplyr, DT
pretty_DT_tab <- function(df, pageLength = 10) {
  df %>%
    dplyr::mutate_if(is.numeric, round, 2) %>%
    DT::datatable(., extensions='Buttons',
                  options = list(
                    searching = T,
                    pageLength = pageLength,
                    dom = 'Bfrtip',
                    buttons = c('csv')))
}

#' @title Quick Jpeg
jpgsnapshot <- function(outpath, plot, type = "wide",width_wide=11,height_wide=8) {
  assert_in(type, c("long", "wide"))
  if (type == "long") {
    jpeg(outpath, width = 8, height = 11, units = "in", res = 500)
    plot(plot)
    graphics.off()
  } else if (type == "wide") {
    jpeg(filename = outpath, width = width_wide, height = height_wide, units = "in", res = 500)
    plot(plot)
    graphics.off()
  }
}


#' @title Quick Diagnostic Plot for MCMC Framework
#' @details Plot the serological chains that are highly correlated

library(ggplot2)
quick_sero_diagnostics <- function(modout) {

  # find max ma (will mix slowest)
  maxma <- modout$inputs$IFRmodel$maxMa
  maxmachain <- drjacoby::plot_par(modout$mcmcout, maxma, display = FALSE)
  spechain <- drjacoby::plot_par(modout$mcmcout, "spec", display = FALSE)
  senschain <- drjacoby::plot_par(modout$mcmcout, "sens", display = FALSE)
  modchain <- drjacoby::plot_par(modout$mcmcout, "mod", display = FALSE)
  sodchain <- drjacoby::plot_par(modout$mcmcout, "sod", display = FALSE)
  seroratechain <- drjacoby::plot_par(modout$mcmcout, "sero_con_rate", display = FALSE)

  maxmachain <- maxmachain[[1]][["trace"]] + theme(legend.position = "none")
  spechain <- spechain[[1]][["trace"]] + theme(legend.position = "none")
  senschain <- senschain[[1]][["trace"]] + theme(legend.position = "none")
  modchain <- modchain[[1]][["trace"]] + theme(legend.position = "none")
  sodchain <- sodchain[[1]][["trace"]] + theme(legend.position = "none")
  seroratechain <- seroratechain[[1]][["trace"]] + theme(legend.position = "none")
  # get legend
  legend_bt <- cowplot::get_legend(maxmachain + theme(legend.position = "bottom",
                                                      legend.title = element_blank()))

  if (modout$inputs$account_seroreversion) {
    revratechain <- drjacoby::plot_par(modout$mcmcout, "sero_rev_rate", display = FALSE)
    revratechain <- revratechain[[1]][["trace"]] + theme(legend.position = "none")

    topp <- cowplot::plot_grid(spechain, senschain, modchain, sodchain, maxmachain, seroratechain,
                               revratechain,
                               ncol = 2, nrow = 4)
  } else {
    # out
    topp <- cowplot::plot_grid(spechain, senschain, modchain, sodchain, maxmachain, seroratechain,
                               ncol = 2, nrow = 3)

  }
  cowplot::plot_grid(topp, legend_bt, nrow = 2, rel_heights = c(1, 0.1))

}
