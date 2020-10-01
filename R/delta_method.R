# Notes (adapted from Kylie Ainslie)
# the calculation of the lower and upper bounds above is based on the Delta method, where
# p_adj is the adjusted prevalence
# se_adj is adjusted standard error
# Briefly, the math for the delta method is as follows
# we say IFR = g(p) = D/(Np + D), where D = deaths, p = seroprevalence, N = population size
# by the Delta method we can calculate the Var(g(p)) = [g'(p)]^2 * Var(p), where g'(p) is the first derivative
# Thus,
# Var(g(p)) = [-(DN) / (Np + D)^2]^2 * Var(p) = ()
# We can then calculate the 95% confidence intervals as
# IFR +- 1.96 * (DN) / (Np + D)^2 * SE(p)

get_delta_CI_vals <- function(deaths, seroprev, popN, SE, tol) {

  # calculate critical value from delta method
  crit_value <- 1.96 * ((popN * deaths) / (seroprev * popN + deaths)^2) * SE

  # IFR and cis
  IFRcalc <- deaths  / (seroprev * popN + deaths)
  LL <- IFRcalc - crit_value
  UL <- IFRcalc + crit_value

  # out
  ret <- c(LL, UL)
  names(ret) <- c("lower.ci", "upper.ci")
  return(ret)
}

