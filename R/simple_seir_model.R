library(odin)
#' no export
seir_step <- function(S, E, I, R, N, Beta, mu_EI, mu_IR, mu_L, delta.t){
  SE_events <- rbinom(n = 1, size = S, prob = 1 - exp((-mu_L -Beta*(I/N))*delta.t))
  SE_deaths <- rbinom(n = 1, size = SE_events, prob = 1 - exp(-mu_L/(mu_L +Beta*(I/N))))
  infected <- SE_events - SE_deaths
  EI_events <- rbinom(n = 1, size = E, prob = 1 - exp(-mu_L -mu_EI*delta.t))
  EI_deaths <- rbinom(n = 1, size = EI_events, prob = 1 - exp(-mu_L/(mu_L + mu_EI)))
  infectious <- EI_events - EI_deaths
  IR_events <- rbinom(n = 1, size = I, prob = 1 - exp(-mu_L -mu_IR*delta.t))
  IR_deaths <- rbinom(n = 1, size = IR_events, prob = 1 - exp(-mu_L/(mu_L + mu_IR)))
  recovered <- IR_events - IR_deaths
  recovered_deaths <- rbinom(n = 1, size = R, prob = 1 - exp(-mu_L))

  S <- S - infected + (EI_deaths + IR_deaths + recovered_deaths)
  E <- E + infected - infectious - EI_deaths
  I <- I + infectious - recovered - IR_deaths
  R <- R + recovered - recovered_deaths
  N <- S + E + I + R
  return(c("N" = N, "S" = S, "E" = E, "I" = I, "R" = R))
}



#' @title Simple SEIR Model
#' @param N numeric; population size
#' @param E0 numeric; Number of individuals infected at time 0 (seeding)
#' @param R0 numeric; Number of individuals immune at time 0
#' @param betas numeric vector; beta values corresponding to beta change time points
#' @param beta_change numeric vector; times that correspond to beta changes in the \link{betas} vector
#' @param mu_EI numeric; rate of transition from E-I compartment: (1/duration of incubation/latency period)
#' @param mu_IR numeric; rate of transition from I-R compartment: (1/duration of disease period)
#' @param mu_L numeric; birth/death rate: (1/life expectancy)
#' @param time numeric; total time of epidemic observed, up to but not including this "day"

run_simple_seir <-  function(N, E0, R0, betas, mu_EI, mu_IR, mu_L, time, ...){
  source("R/assertions_v5.R")
  assert_numeric(N)
  assert_numeric(E0)
  assert_numeric(R0)
  assert_numeric(betas)
  assert_numeric(beta_changes)
  assert_numeric(mu_EI)
  assert_numeric(mu_IR)
  assert_numeric(mu_L)
  assert_numeric(time)
  assert_greq(time, max(beta_changes))

  #......................
  # setup
  #......................
  S0 <- N - E0 - R0
  beta_changes_steps <- c((beta_changes - dplyr::lag(beta_changes))[2:length(beta_changes)],
                          (time - beta_changes[length(beta_changes)]))
  betas <- rep(betas, times = beta_changes_steps)
  time.steps <- 1:(time-1)


  # init
  ret <- matrix(NA, nrow = length(time.steps) , ncol = 5)
  ret[1,] <- seir_step(S = S0,
                       E = E0,
                       I = 0,
                       R = R0,
                       N = N,
                       Beta = betas[1],
                       mu_EI = mu_EI,
                       mu_IR = mu_IR,
                       mu_L = mu_L,
                       delta.t = 0)

  for (t in 2:length(time.steps)) {
    ret[t,] <-  seir_step(S = ret[t-1, 2],
                          E = ret[t-1, 3],
                          I = ret[t-1, 4],
                          R = ret[t-1, 5],
                          N = ret[t-1, 1],
                          Beta = betas[t],
                          mu_EI = mu_EI,
                          mu_IR = mu_IR,
                          mu_L = mu_L,
                          delta.t = time.steps[t])
  }
  # out
  ret <- cbind.data.frame(time, ret)
  colnames(ret) <- c("time", "N", "S", "E", "I", "R")
  return(ret)

}
