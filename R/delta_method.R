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

get_delta_CI_vals <- function(deaths, seroprev, popN, SE, tol) {
  # basic transforms
  logit <- function(x, tol=1e-4){
    return( log(((x+tol)/(1-x+tol))) )
  }

  expit <- function(x, tol=1e-4){
    return( 1/(1+exp(-x + tol)) )
  }
  # calculate critical value from delta method
  crit_value <- 1.96 * ((deaths * popN) / (seroprev * popN)^2) * SE

  # IFR
  IFRcalc <- deaths  / (seroprev * popN + deaths)
  # calculations in transformed space to account for binomial vs. normal
  LL <- expit( logit(IFRcalc, tol=tol) - crit_value, tol = tol)
  UL <- expit( logit(IFRcalc, tol=tol) + crit_value, tol = tol)

  # out
  ret <- c(LL, UL)
  names(ret) <- c("lower.ci", "upper.ci")
  return(ret)
}


