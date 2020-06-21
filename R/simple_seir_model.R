simple_seir <- odin::odin({
  ## deal with time
  ## Core equations for transitions between compartments:
  update(S) <- S - n_SE
  update(E) <- E + n_SE - n_EI
  update(I) <- I + n_EI - n_IR
  update(R) <- R + n_IR

  ## Individual probabilities of transition:
  p_SE <- 1 - exp(-beta * (I / N))
  p_EI <- 1 - exp(-sigma)
  p_IR <- 1 - exp(-gamma)

  ## Draws from binomial distributions for numbers changing between
  ## compartments:
  n_SE <- rbinom(S, p_SE)
  n_EI <- rbinom(E, p_EI)
  n_IR <- rbinom(I, p_IR)

  ## Total population size
  N <- S + E + I + R

  ## Initial states:
  initial(S) <- S_ini
  initial(E) <- E_ini
  initial(I) <- I_ini
  initial(R) <- R_ini

  ## User defined parameters - default in parentheses:
  S_ini <- user(1000)
  E_ini <- user(1)
  I_ini <- user(1)
  R_ini <- user(0)
  beta <- user(0.2)
  gamma <- user(0.1)
  sigma <- user(0.1)
})


#' @title Simple SEIR Model
#' @param N numeric; population size
#' @param E0 numeric; Number of individuals infected at time 0 (seeding)
#' @param R0 numeric; Number of individuals immune at time 0
#' @param betas numeric vector; beta values corresponding to beta change time points
#' @param beta_change numeric vector; times that correspond to beta changes in the \link{betas} vector
#' @param sigma numeric; rate of transition from E-I compartment: (1/duration of incubation/latency period)
#' @param gamma numeric; rate of transition from I-R compartment: (1/duration of disease period)
#' @param time numeric; total time of epidemic observed, up to but not including this "day"

run_simple_seir <-  function(N, E0, R0, betas, beta_changes, sigma, gamma, time, ...){
  source("R/assertions_v5.R")
  assert_numeric(N)
  assert_numeric(E0)
  assert_numeric(R0)
  assert_numeric(betas)
  assert_numeric(beta_changes)
  assert_numeric(sigma)
  assert_numeric(gamma)
  assert_numeric(time)
  assert_greq(time, max(beta_changes))

  #......................
  # setup
  #......................
  if (length(beta_changes) > 1) {
    beta_changes_steps <- c(1, (beta_changes - dplyr::lag(beta_changes))[2:length(beta_changes)],
                            (time - beta_changes[length(beta_changes)]))
  } else {
    beta_changes_steps <- c(1, time-1)
  }

  beta_changes <- c(beta_changes, time)
  # init empty
  ret <- matrix(NA, nrow = time, ncol = 5)
  # loop through
  mod.prev <- simple_seir(
    beta = betas[1],
    sigma = sigma,
    gamma = gamma,
    S_ini = N,
    E_ini = E0,
    I_ini = 0,
    R_ini = 0)
  mod.prev <- mod.prev$run(1)[1,]

  for (i in 1:length(betas)) {
    mod.new <- simple_seir(
      beta = betas[i],
      sigma = sigma,
      gamma = gamma,
      S_ini = mod.prev["S"],
      E_ini = mod.prev["E"],
      I_ini = mod.prev["I"],
      R_ini = mod.prev["R"])
    mod.new <- mod.new$run(1:(beta_changes_steps[i+1]+1)) # new beta, "new start"
    # store
    ret[beta_changes[i]:beta_changes[i+1], ] <- mod.new
    # overwrite prev
    mod.prev <- mod.new[nrow(mod.new),]
  }
  # clean up "time"
  ret <- as.data.frame(ret)
  colnames(ret) <- colnames(mod.new)
  ret$step <- 1:nrow(ret)
  return(ret)
}

