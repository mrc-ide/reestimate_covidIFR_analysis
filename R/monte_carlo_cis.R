
#' @title Binomial Monte Carlo Draws
get_binomial_monte_carlo_cis <- function(deaths, popN, npos, ntest, iters) {
  # IFR calc out
  deaths/((rbinom(n = iters, size = ntest, prob = npos/ntest)/ntest) * popN + deaths)
}

get_normal_monte_carlo_cis <- function(deaths, popN, mu, sigma, iters) {
  # internal functions
  logit <- function(x, tol=1e-4){
    return( log(((x+tol)/(1-x+tol))) )
  }
  expit <- function(x, tol=1e-4){
    return( 1/(1+exp(-x + tol)) )
  }
  # IFR calc out
  deaths/(expit(rnorm(n = iters, mean = logit(mu), sd = sigma)) * popN + deaths)
}

