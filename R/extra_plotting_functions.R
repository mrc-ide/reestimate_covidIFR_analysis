source(paste0(here::here(), "/R/assertions_v5.R"))

#' @title Specific Function for Summarizing Simulated Runs
summarize_plot_simulated_runs <- function(modoutpath, fit_map) {
  # read in
  modout <- readRDS(modoutpath)
  fitdat <- fit_map %>%
    dplyr::filter(sim == sub(".RDS", "", basename(modoutpath)))
  # extract pieces
  spec <- fitdat$spec
  sens <- fitdat$sens

  infxns <- fitdat %>%
    dplyr::select(c("curve")) %>%
    tidyr::unnest(cols = "curve")

  demog <- fitdat %>%
    dplyr::select(c("demog")) %>%
    tidyr::unnest(cols = "demog")
  fatalitydata <- fitdat %>%
    dplyr::select(c("fatalitydata")) %>%
    tidyr::unnest(cols = "fatalitydata")
  fitdat <- dplyr::left_join(fatalitydata, demog, by = "Strata")

  #............................................................
  # IFR plot
  #...........................................................
  ifrs <- COVIDCurve::get_cred_intervals(IFRmodel_inf = modout, whichrung = paste0("rung", 1),
                                         what = "IFRparams", by_chain = FALSE) %>%
    dplyr::rename(Strata = param)
  infxncurve <- COVIDCurve::draw_posterior_infxn_cubic_splines(IFRmodel_inf = modout,
                                                               dwnsmpl = 1e3,
                                                               by_chain = FALSE)
  # make ifr and incidence plots
  plot1 <- ggplot() +
    geom_pointrange(data = ifrs, aes(x = Strata, ymin = LCI, ymax = UCI, y = median),
                    color = "#969696", size = 1.2) +
    geom_point(data = fitdat, aes(x = Strata, y = IFR), color = "#000000", size = 2) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.90, hjust= 1, face = "bold"),
          legend.position = "right") +
    xlab("") + ylab("Median (95% CIs)") +
    labs(caption = "Grey and black points represent model fits and true IFR estimates, respectively")
  plot2 <- infxncurve$plotObj +
    geom_line(data = infxns, aes(x = time, y = infxns))
  # main plot
  mainIFR_plotObj <- cowplot::plot_grid(plot1, plot2, ncol = 1, nrow = 2)

  #............................................................
  # Specifitity
  #...........................................................
  # specPlot <- drjacoby::plot_par(modout$mcmcout, "spec", display = F)
  # top <- specPlot[[1]][[1]]
  # left <- specPlot[[1]][[2]] + geom_vline(xintercept = spec, size = 2, color = "#F61843")
  # right <- specPlot[[1]][[3]]
  # bottom <- cowplot::plot_grid(left, right, nrow = 1)
  # specPlot <- cowplot::plot_grid(top, bottom, ncol = 1)
  specPlot <- drjacoby::plot_par(modout$mcmcout, "spec", display = F)
  specPlot <- specPlot[[1]][[2]] + geom_vline(xintercept = spec, size = 2, color = "#F61843")
  #............................................................
  # Sensitivity
  #...........................................................
  # sensPlot <- drjacoby::plot_par(modout$mcmcout, "sens", display = F)
  # top <- sensPlot[[1]][[1]]
  # left <- sensPlot[[1]][[2]] + geom_vline(xintercept = sens, size = 2, color = "#F61843")
  # right <- sensPlot[[1]][[3]]
  # bottom <- cowplot::plot_grid(left, right, nrow = 1)
  # sensPlot <- cowplot::plot_grid(top, bottom, ncol = 1)
  sensPlot <- drjacoby::plot_par(modout$mcmcout, "sens", display = F)
  sensPlot <- sensPlot[[1]][[2]] + geom_vline(xintercept = sens, size = 2, color = "#F61843")


  #............................................................
  # seroday and mc accept
  #...........................................................
  seroday <- drjacoby::plot_par(modout$mcmcout, "sero_day", display = F)
  seroday <- seroday[[1]][[1]]
  mcacc <- drjacoby::plot_mc_acceptance(modout$mcmcout)

  #............................................................
  # bring together
  #...........................................................
  whitespace <- tibble::tibble(
    xmin = -Inf, ymin = 0,
    xmax = Inf, ymax = 5
  ) %>%
    ggplot() +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymin), color = "#000000")  +
    theme_bw() +
    theme(
      plot.background = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank()
    )

  cowplot::plot_grid(mainIFR_plotObj, specPlot, sensPlot, seroday, mcacc,
                     whitespace, ncol = 1)
  out <- cowplot::plot_grid(mainIFR_plotObj, specPlot, sensPlot, seroday, mcacc, ncol = 1,
                            rel_heights = c(0.75, 0.1, 0.1, 0.1, 0.25))
  return(list(out))
}





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
    revshapechain <- drjacoby::plot_par(modout$mcmcout, "sero_rev_shape", display = FALSE)
    revscalechain <- drjacoby::plot_par(modout$mcmcout, "sero_rev_scale", display = FALSE)
    revshapechain <- revshapechain[[1]][["trace"]] + theme(legend.position = "none")
    revscalechain <- revscalechain[[1]][["trace"]] + theme(legend.position = "none")

    topp <- cowplot::plot_grid(spechain, senschain, modchain, sodchain, maxmachain, seroratechain,
                               revshapechain, revscalechain,
                               ncol = 2, nrow = 4)
  } else {
    # out
    topp <- cowplot::plot_grid(spechain, senschain, modchain, sodchain, maxmachain, seroratechain,
                               ncol = 2, nrow = 3)

  }
  cowplot::plot_grid(topp, legend_bt, nrow = 2, rel_heights = c(1, 0.1))

}
