
#' Generates raw data based on study parameters
#' that can be used later to simulate infection times
#'
#' @export
#' @param n Number of observations per time point
#' @param prevalence Constant prevalence
generate.raw.data <- function(n_sims, n, prevalence){

  # number of positive subjects
  n_p <- rbinom(n=n_sims, size=n, prob=prevalence)

  # number of negative subjects
  n_n <- n - n_p

  return(list(
    n=rep(n, n_sims),
    n_p=n_p,
    n_n=n_n
  ))
}

#' Simulate recency indicators based on the data, true phi function,
#' and the infection incidence function.
#'
#' @export
#' @param data Outputs from the `generate.raw.data` function
#' @param infection.function Function that simulates the infection time
#' @param model A `glm` model that was fit using data generated from
#'              the same true.phi.function, but then fit with cubic polynomial
#'
simulate.recent <- function(sim_data, infection.function, phi.func, baseline_incidence, prevalence, rho){

  # get the uniform for the cumulative distribution
  # function for each individual
  cdfs <- lapply(sim_data$n_p, runif)

  # infection time
  t_infect <- mapply(infection.function, e=cdfs, t=0, p=prevalence,
                     lambda_0=baseline_incidence, rho=rho)


  # infection duration to pass to the phi hat function
  infect_duration <- mapply(function(t, u) t - u, t=0, u=t_infect)

  # probability of recent infection
  recent_probabilities <- lapply(infect_duration, phi.func)

  # indicators
  indicators <- lapply(recent_probabilities,
                       function(x) rbinom(n=length(x), size=1, prob=x))

  # number of recents
  n_r <- lapply(indicators, sum) %>% unlist

  new.data <- list(
    n_r=n_r
  )

  all.data <- append(data, new.data)
  return(all.data)
}

# INCIDENCE FUNCTIONS

# (1) constant incidence
c.incidence <- function(t, lambda_0, rho=NA) lambda_0
# (2) linear incidence
l.incidence <- function(t, lambda_0, rho=1) lambda_0 - rho * t
# (3) exponential incidence
e.incidence <- function(t, lambda_0, rho=1) lambda_0 * exp(-rho * t)

# infection time functions
# (1) constant incidence infection time
c.infections <- function(e, t, p, lambda_0, rho=NA){
  infections <- t - p*e / ((1 - p) * c.incidence(t=t, lambda_0=lambda_0))
  return(infections)
}

#' Compute the average incidence
#' over an interval.
avg.incidence <- function(inc_function, ts, ...){
  incidence <- lapply(ts, inc_function, ...) %>% unlist
  return(mean(incidence))
}

# INFECTION FUNCTIONS

#' Constant incidence function
c.incidence <- function(t, lambda_0, rho=NA) lambda_0

#' Linearly decreasing incidence function
l.incidence <- function(t, lambda_0, rho=1) lambda_0 - rho * t

# Exponentially decreasing incidence function
e.incidence <- function(t, lambda_0, rho=1) lambda_0 * exp(-rho * t)

#' Function to produce infection times
#' from constant incidence.
c.infections <- function(e, t, p, lambda_0, rho=NA){
  infections <- t - p*e / ((1 - p) * c.incidence(t=t, lambda_0=lambda_0))
  return(infections)
}

#' Function to produce infection times
#' from linear incidence.
l.infections <- function(e, t, p, lambda_0, rho){
  incidence <- l.incidence(t=t, lambda_0=lambda_0, rho=rho)
  numerator <- incidence**2 + 2 * rho * p * e / (1 - p)
  numerator <- sqrt(numerator) - incidence

  infections <- t - numerator / rho
  return(infections)
}

#' Function to produce infection times
#' from exponential incidence.
e.infections <- function(e, t, p, lambda_0, rho){
  incidence <- e.incidence(t=t, rho=rho, lambda_0=lambda_0)
  infections <- t - (1/rho) * log(rho*p*e/((1-p)*incidence) + 1)
  return(infections)
}
