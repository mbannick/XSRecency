
#' Generates raw data of number positive, negative, and screened
#' based on study parameters
#' that can be used later to simulate infection times.
#'
#' @param n_sims Number of simulations
#' @param n Number of observations per time point
#' @param prevalence Constant prevalence
#' @param times Times at which enrollment happens
generate.raw.data <- function(n_sims, n, prevalence, times=c(0)){

  # number of times
  n_times <- length(times)

  # number of positive subjects
  n_p <- replicate(n_times, rbinom(n=n_sims, size=n, prob=prevalence)) %>% t

  # number of negative subjects
  n_n <- n - n_p

  return(list(
    n=matrix(n, nrow=n_times, ncol=n_sims),
    n_p=n_p,
    n_n=n_n,
    times=matrix(rep(times, n_sims), ncol=n_sims, byrow=FALSE)
  ))
}

#' Simulate recency indicators based on the data, true phi function,
#' and the infection incidence function.
#'
#' @param sim_data Outputs from the `generate.raw.data` function
#' @param infection.function Function that simulates the infection time
#' @param phi.func Test positive function for recency assay
#' @param baseline_incidence Baseline incidence value (time 0)
#' @param prevalence Constant prevalence value
#' @param rho A parameter to the `infection.function` for how quickly incidence changes
simulate.recent <- function(sim_data, infection.function,
                            phi.func, baseline_incidence, prevalence, rho){

  # Get dimensions
  dims <- dim(sim_data$n_p)
  n_times <- dims[1]
  n_sims <- dims[2]

  # get the uniform for the cumulative distribution
  # function for each individual

  times <- as.vector(sim_data$times)
  n_p <- as.vector(sim_data$n_p)

  cdfs <- lapply(n_p, runif)

  # infection time
  t_infect <- mapply(infection.function, e=cdfs, t=times, p=prevalence,
                     lambda_0=baseline_incidence, rho=rho, SIMPLIFY=F)

  # infection duration to pass to the phi hat function
  infect_duration <- mapply(function(t, u) t - u, t=times, u=t_infect, SIMPLIFY=F)

  # probability of recent infection
  recent_probabilities <- lapply(infect_duration, phi.func)

  # indicators
  indicators <- lapply(recent_probabilities,
                       function(x) rbinom(n=length(x), size=1, prob=x))

  # number of recents
  n_r <- lapply(indicators, sum) %>% unlist

  # if we have multiple times, convert it back into a matrix
  if(n_times > 1)  n_r <- matrix(n_r, nrow=n_times, ncol=n_sims, byrow=FALSE)
  if(n_times > 1)  n_p <- sim_data$n_p
  if(n_times > 1)  n_n <- sim_data$n_n
  if(n_times > 1)  times <- sim_data$times
  if(n_times > 1)  n <- sim_data$n

  # if we have only one time (cross-sectional), keep everything
  # in vector form
  if(n_times == 1) n_n <- as.vector(sim_data$n_n)
  if(n_times == 1) n   <- as.vector(sim_data$n)

  aspect_list <- list(
    n=n,
    n_p=n_p,
    n_n=n_n,
    n_r=n_r,
    times=times
  )
  return(aspect_list)
}

#' Create simulations of trials based on a past infection time function
#' and a prevalence and baseline incidence (time 0).
#'
#' @export
#' @param n_sims Number of simulations
#' @param n Number of subjects screened
#' @param infection.function Function that simulates the infection time
#' @param phi.func Test positive function for recency assay
#' @param baseline_incidence Baseline incidence value (time 0)
#' @param prevalence Constant prevalence value
#' @param rho A parameter to the `infection.function` for how quickly incidence changes
#' @param times A vector of times at which to enroll (e.g. c(0, 1, 2, 3))
#' @return A list of vectors for sample size screened across simulations:
#' - n: number screened
#' - n_p: number of positives
#' - n_n: number of negatives
#' - and number of recents (`n_r`)
generate.data <- function(n_sims, n, infection.function,
                          phi.func, baseline_incidence,
                          prevalence, rho, times=c(0)){
  data <- generate.raw.data(n_sims, n, prevalence, times=times)
  data <- simulate.recent(data, infection.function, phi.func,
                          baseline_incidence, prevalence, rho)
  return(data)
}

#' Simulate one example data set with past infection times.
#'
#' @export
#' @param n Number of subjects screened
#' @param infection.function Optional function that simulates the infection time
#' @param phi.func Optional duration-specific test-recent function
#' @param baseline_incidence Baseline incidence value (time 0)
#' @param prevalence Prevalence value (constant over time)
#' @param rho A parameter to the infection fun
#' @return A data frame with subject ID
#' - n: number screened
#' - n_p: number of positives
#' - n_n: number of negatives
#' - and number of recents (`n_r`)

# INFECTION FUNCTIONS

#' Constant incidence function
c.incidence <- function(t, lambda_0, rho=NA) lambda_0

#' Linearly decreasing incidence function
l.incidence <- function(t, lambda_0, rho=1) lambda_0 - rho * t

# Exponentially decreasing incidence function
e.incidence <- function(t, lambda_0, rho=1) lambda_0 * exp(-rho * t)

#' Function to produce infection times
#' from constant incidence.
#' Note that you can work with e or 1 - e since e is Uniform(0, 1)
c.infections <- function(e, t, p, lambda_0, rho=NA){
  infections <- t - p*e / ((1 - p) * lambda_0)
  return(infections)
}

#' Function to produce infection times
#' from linear incidence.
l.infections <- function(e, t, p, lambda_0, rho){
  incidence <- lambda_0
  numerator <- incidence**2 + 2 * rho * p * e / (1 - p)
  numerator <- sqrt(numerator) - incidence

  infections <- t - numerator / rho
  return(infections)
}

#' Function to produce infection times
#' from exponential incidence.
e.infections <- function(e, t, p, lambda_0, rho){
  incidence <- lambda_0
  infections <- t - (1/rho) * log(rho*p*e/((1-p)*incidence) + 1)
  return(infections)
}

# INFECTION FUNCTION FOR ARBITRARY INCIDENCE FUNCTION

#' Function to simulate infection times from an arbitrary incidence
#' function and with constant prevalence.
#'
#' @export
#' @param p Constant prevalence value
#' @param inc.function A function that gives incidence for its single argument time
#' @param nsims Number of infection times to return
#' @param dt The precision for numerical integration
sim.infection.times <- function(p, inc.function, nsims=1000, dt=0.01){
  # Need to vectorize a scalar function
  if(length(inc.function(c(1, 2))) == 1) inc.function <- Vectorize(inc.function)

  # Simulate e's and the solve for T
  e <- runif(nsims)
  lh <- e * p/(1-p)
  ds <- seq(0, 100, by=dt)
  vals <- inc.function(ds)
  int <- cumsum(vals*dt)
  indexes <- sapply(lh, function(x) max(which(int < x)))
  times <- ds[indexes]
  return(times)
}
