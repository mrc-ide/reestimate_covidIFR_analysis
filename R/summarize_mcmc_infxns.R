library(magrittr)
get_param_summaries <- function(rmcmcout, modinf) {
  # get param sets
  infxnparams <- modinf$Infxnparams
  ifrparams <- modinf$IFRparams
  seroparams <- modinf$Seroparams
  params <- list(infxnparams, ifrparams, seroparams)
  names(params) <- c("Infxnparams", "IFRparams", "Seroparams")
  # call cred intervals
  credintervals <- function(rmcmcout, params) {
    rmcmcout$output %>%
      dplyr::select(c("chain", "iteration", params)) %>%
      tidyr::gather(., key = "param", value = "est", 3:ncol(.)) %>%
      dplyr::group_by(chain, param) %>%
      dplyr::summarise(
        min = min(est),
        LCI = quantile(est, 0.025),
        median = median(est),
        mean = mean(est),
        UCI = quantile(est, 0.975),
        max = max(est),
        ESS = coda::effectiveSize(coda::as.mcmc(est))
      )
  }

  # out
  lapply(params, credintervals, rmcmcout = rmcmcout)

}


#' @title Draw Posterior Observations from Infection Curve
#' @param rmcmcout DrJacoby output; MCMC from Dr. Jacoby
#' @param modinf R6 Class; modinf from COVIDCuvre
#' @param CIquant numeric; CI quantile to draw from posterior log-likelihood greater than or equal to
get_infxn_curve <- function(rmcmcout, modinf, CIquant) {
  #......................
  # get gradients for slopes
  #......................
  liftover_infxn_curve <- function(mcmcout, knots){
    mcmcout$output$nodegrad1 <- (rmcmcout$output$y2 - rmcmcout$output$y1)/(knots[2] - knots[1])
    mcmcout$output$nodegrad2 <- (rmcmcout$output$y3 - rmcmcout$output$y2)/(knots[3] - knots[2])
    return(mcmcout)
  }
  r_mcmc_out.infxncurve <- liftover_infxn_curve(rmcmcout, knots = modinf$knots)

  #......................
  # make infection curve
  #......................
  make_infxn_curve <- function(y1, nodegrad1, nodegrad2, knots){
    curr_day <- knots[length(knots)] - knots[1] + 1
    ret <- rep(NA, times = curr_day)
    ret[1] <- y1
    for (i in 2:curr_day) {
      lvl <- cut(i, breaks = knots, labels = paste0("k", 1:(length(knots)-1)))
      switch(as.character(lvl),
             "k1" = {
               ret[i] <- nodegrad1 + ret[i-1]
             },
             "k2" = {
               ret[i] <- nodegrad2 + ret[i-1]
             }
      )

    }
    ret <- data.frame(time = 1:curr_day, infxns = exp(ret))
    return(ret)
  }
  knots <- modinf$knots
  r_mcmc_out.infxncurve$output$knots <- lapply(1:nrow(r_mcmc_out.infxncurve$output), function(x) return(unlist(knots)))
  r_mcmc_out.infxncurve$output$infxncurves <- purrr::pmap(r_mcmc_out.infxncurve$output[,c("y1", "nodegrad1", "nodegrad2", "knots")],
                                                          make_infxn_curve)

  #......................
  # sample by CI limit
  #......................
  upperci <- quantile(r_mcmc_out.infxncurve$output$loglikelihood, probs = CIquant)
  plotdat <- r_mcmc_out.infxncurve$output %>%
    dplyr::filter(stage == "sampling") %>%
    dplyr::filter(loglikelihood >= upperci) %>%
    dplyr::select(c("chain", "infxncurves")) %>%
    dplyr::group_by(chain) %>%
    dplyr::mutate(sim = 1:dplyr::n()) %>%
    dplyr::ungroup(chain) %>%
    tidyr::unnest(cols = "infxncurves")

  # plot
  plotObj <- ggplot() +
    geom_line(data = plotdat, mapping = aes(time, infxns, group = sim), alpha = 0.25,
              lwd = 0.5, color = "#d9d9d9") +
    geom_vline(xintercept = knots, color = "#cb181d", lwd = 0.25, linetype = "dashed", alpha = 0.5) +
    xlab("Time") + ylab("Num. Infxns")  +
    labs(title = "Posterior Draws of the Infection Curve") +
    facet_wrap(. ~ chain) +
    theme_minimal() +
    theme(
      plot.title = element_text(family = "Helvetica", face = "bold", vjust = 0.5,  hjust = 0.5, size = 18),
      plot.subtitle = element_text(family = "Helvetica", face = "bold", vjust = 0.5,  hjust = 0.5, size = 18),
      axis.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, vjust = 0.5, size = 16),
      axis.text.x = element_text(family = "Helvetica", angle = 45, hjust = 0.5, vjust = 0.5, size = 15),
      axis.text.y = element_text(family = "Helvetica", hjust = 0.5, vjust = 0.5, size = 15),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.line = element_line(color = "#000000", size = 1.2),
      legend.position = "none")

  #......................
  # out
  #......................
  ret <- list(
    plotdata = plotdat,
    plotObj = plotObj
  )
  return(ret)
}
