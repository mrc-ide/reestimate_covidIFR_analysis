# Notes (adapted from Kylie Ainslie)
# the calculation of the lower and upper bounds above is based on the Delta method, where
# p_adj is the adjusted prevalence
# se_adj is adjusted standard error
# Briefly, the math for the delta method is as follows
# we say IFR = g(p) = D/(Np + D), where D = deaths, p = seroprevalence, N = population size
# by the Delta method we can calculate the Var(g(p)) = [g'(p)]^2 * Var(p), where g'(p) is the first derivative
# Thus,
# Var(g(p)) = [-(D * N) / (Np + D)^2]^2 * Var(p) = [D^2 N^2]/[(Np + D)^4] * Var(p)
# We can then calculate the 95% confidence intervals as
# IFR +- 1.96 * (D * N) / (Np + D)^2 * SE(p)

#' @return critical value for CI

get_delta_CI_val <- function(deaths, seroprev, popN, SE) {
  1.96 * ((deaths * popN) / (seroprev * popN)^2) * SE
}



#----------------------------------------------------------------------------------------------------
# proportion CIs
#----------------------------------------------------------------------------------------------------
#' @param pt numeric; point estimate
#' @param crit numeric; critical value
#' @param tol numeric; tolerance for logit transformation

getCI_from_logit_transfrom <- function(pt, crit, tol, alpha){
  # basic transforms
  logit <- function(x, tol=1e-4){
    return( log(((x+tol)/(1-x+tol))) )
  }

  expit <- function(x, tol=1e-4){
    return( 1/(1+exp(-x + tol)) )
  }

  # calculations
  p.logit <- logit(pt, tol=tol)

  LL <- expit( p.logit - crit, tol = tol)
  UL <- expit( p.logit + crit, tol = tol)

  ret <- c(LL, UL)
  names(ret) <- c("lower.ci", "upper.ci")
  return(ret)
}



