
#' @title Binomial Monte Carlo Draws
get_binomial_monte_carlo_cis <- function(deaths, popN, npos, ntest, iters) {
  # IFR calc out
  ret <- deaths/((rbinom(n = iters, size = ntest, prob = npos/ntest)/ntest) * popN + deaths)
  # catch when denominator goes to zero -- saying that there are no infections but deaths, so inf
  if(any(is.nan(ret))) {
    ret[is.nan(ret)] <- Inf
  }
  return(ret)
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
  ret <- deaths/(expit(rnorm(n = iters, mean = logit(mu), sd = sigma)) * popN + deaths)
  # catch when denominator goes to zero -- saying that there are no infections but deaths, so inf
  if(any(is.nan(ret))) {
    ret[is.nan(ret)] <- Inf
  }
  return(ret)
}

