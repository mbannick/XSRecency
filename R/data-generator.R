library(magrittr)

#' @title Simulates cross-sectional data based on epidemiological parameters.
#'
#' @description
#' These are the available incidence functions. \code{t} is relative to time 0.
#' So, to get incidence 2 years ago, \code{t = -2}.
#'
#' @details
#' Constant incidence: \deqn{\lambda(t) = \lambda_0 - \rho}
#' Linear incidence: \deqn{\lambda(t) = \lambda_0 - \rho t}
#' Exponential incidence: \deqn{\lambda(t) = \lambda_0 \exp(-\rho t)}
#' Constant-linear incidence: \deqn{\lambda(t) = \lambda_0 - \rho (t + T^*) I(-t > T^*)}
#'
#' @param n_sims Number of simulations
#' @param n Number of observations per time point
#' @param prevalence Constant prevalence
#' @param baseline_incidence Baseline incidence
#' @param incidence_type Incidence type: constant, linear, exponential, or linear-constant (constant up to \code{bigT} argument, linear after)
#' @param p Constant prevalence
#' @param rho Parameter for linear or exponentially decreasing incidence
#' @param t Time t baseline (defaults to 0)
#'
#' @examples
#' phi.func <- getTestRecentFunc(window=200, shadow=191)
#' sim <- simCrossSect(phi.func=phi.func,
#'                     incidence_type="linear", prevalence=0.29, baseline_incidence=0.032, rho=0.07)
#' sim(10)
#'
#' @export
simCrossSect <- function(phi.func, incidence_type,
                         prevalence, baseline_incidence, rho=NA, bigT=NA){

  sim.screening <- function(n, return_n=TRUE){

    # number of positive subjects
    pos <- rbinom(n=n, size=1, prob=prevalence) %>% sort
    npos <- sum(pos)

    # generate infection time
    infection.func <- simInfections(
      type=incidence_type,
      p=prevalence,
      lambda_0=baseline_incidence,
      rho=rho, bigT=bigT)
    inf.time <- infection.func(npos)

    # recency indicator
    rec <- rbinom(n=npos, size=1, prob=phi.func(-inf.time))

    df <- data.frame(
      di=pos,
      ui=c(rep(NA, n-npos), -inf.time),
      ri=c(rep(NA, n-npos), rec)
    )
    if(!return_n){
      df <- df[df$di == 1,]
    }

    return(df)
  }
  return(sim.screening)
}

#' Generates prior testing data for positive subjects
#' based on a prior testing distribution and prior testing probability
#'
#' @param ptest.dist A prior test result distribution function, which must have
#'   be a function of `u` which is an
#'   infection duration
#'   though it can ignore the `u` argument inside the function.
#'   An example is `rnorm(1, mean=u)`.
#' @param ptest.prob Probability of prior test result being available.
#'   Can instead be a function of `u`, the infection duration.
#' @examples
#' set.seed(101)
#' phi.func <- getTestRecentFunc(window=200, shadow=191)
#' sim <- simCrossSect(phi.func=phi.func,
#'                     incidence_type="linear", prevalence=0.29, baseline_incidence=0.032, rho=0.07)
#' df <- sim(100)
#' sim.pt <- simPriorTests(ptest.dist=function(u) runif(1, 0, 4),
#'                         ptest.prob=function(u) 0.5)
#' sim.pt(df[df$di == 1,])
#'
#' @export
simPriorTests <- function(ptest.dist, ptest.prob, ptest.dist2=NULL){
  sim.pt <- function(df){

    # Timing of test
    df$ti <- sapply(df$ui, ptest.dist)

    # Availability of test
    df$qi <- sapply(df$ui, function(u) rbinom(1, 1, ptest.prob(u)))

    # Apply second testing mechanism
    if(!is.null(ptest.dist2)){
      df$ti2    <- sapply(df$ui, ptest.dist2)
      df$ti     <- ifelse((((df$ti < df$ti2) | (df$ti2 < 0)) & df$qi == 1),
                          df$ti, df$ti2)
      df$qi     <- as.integer(df$ti > 0)
    }

    # Test result
    df$deltai <- as.integer(df$ui >= df$ti)

    df[df$qi == 0, "ti"]     <- NA
    df[df$qi == 0, "deltai"] <- NA

    return(df)
  }
  return(sim.pt)
}

# Modifies the prior testing data for recall bias.
#
# @param t_noise Function to add noise to a single prior testing time.
#   Must be a function of t, the prior testing time.
# @param d_misrep Probability of given a positive prior test, individual
#   reports negative prior test.
# @param p_misrep Probability of not reporting a prior test result,
#   among those with positive prior test.
# @param t_range Range of times in which to keep prior testing times.
#
# @export
# @examples
# infection.func <- simInfections(t=0, type="linear", p=0.29, lambda=0.032, rho=0.07)
# phi.func <- getTestRecentFunc(window=200, shadow=191)
# sim <- simCrossSect(prevalence=0.29, infection.func=infection.func, phi.func=phi.func)
# df <- sim(100)
# sim.pt <- simPriorTests(ptest.dist=function(u) runif(1, 0, 4),
#                          ptest.prob=function(u) 0.5)
# pt.df <- sim.pt(df[df$di == 1,])
# t_noise <- function(t) max(0, t + rnorm(n=1, sd=0.5))
# modifypt <- modifyPriorTests(t_noise=t_noise, d_misrep=1, p_misrep=0, t_range=c(1, 2))
# pt.df2 <- modifypt(pt.df)
#
modifyPriorTests <- function(t_noise=function(t) t, t_range=NULL,
                              d_misrep=0, p_misrep=0){
  modifypt <- function(df){

    if(d_misrep > 0 & p_misrep > 0){
      stop("Cannot have both d_misrep and p_misrep at the same time.")
    }
    # Apply noise to the testing times if desired,
    # but not changing ptest_delta
    if(!is.null(t_noise)){
      df$ti <- sapply(df$ti, function(t) sapply(t, t_noise))
    }
    # Recall bias for those with prior tests
    if(p_misrep > 0){
      p_bin <- rbinom(n=nrow(df), size=1, prob=1-p_misrep)
      df$deltai <- df$deltai * p_bin
    }
    # Recall bias for those with positive prior tests
    if(d_misrep > 0){
      pos_idx <- which(df$deltai == 1)
      if(sum(pos_idx) > 0){
        d_bin <- rbinom(n=length(pos_idx), size=1, prob=1-d_misrep)
        df[pos_idx, "qi"] <- df[pos_idx, "qi"] * d_bin
      }
    }
    # Get rid of test results outside of certain range
    if(!is.null(t_range)){
      rem_idx <- which((df$ti < min(t_range)) | (df$ti > max(t_range)))
      df[rem_idx, "qi"] <- 0
    }

    df[(df$qi == 0), "ti"] <- NA
    df[(df$qi == 0), "deltai"] <- NA

    return(df)
  }
  return(modifypt)
}

# INFECTION FUNCTIONS

# Constant incidence function
#
# @param t Time, float or vector
# @param lambda_0 Incidence constant value
# @return Incidence at time t
# @examples
# incidence.con(0, lambda_0=0.05)
# incidence.con(c(-1, 0, 0), lambda_0=0.05)
incidence.con <- Vectorize(function(t, lambda_0, rho=NA) lambda_0)

# Linearly decreasing incidence function
#
# \deqn{
#   \lambda(t) = \lambda_0 - \rho t
# }
#
# @param t Time, float or vector
# @param lambda_0 Incidence at time 0
# @param rho Linear decrease in incidence
# @return Incidence at time t
# @examples
# incidence.lin(0, lambda_0=0.05, rho=1e-3)
# incidence.lin(c(-1, 0, 1), lambda_0=0.05, rho=1e-3)
incidence.lin <- function(t, lambda_0, rho=1) lambda_0 - rho * t

# Exponentially decreasing incidence function
#
# \deqn{
#   \lambda(t) = \lambda_0 \exp(-\rho t)
# }
#
# @param t Time, float or vector
# @param lambda_0 Incidence at time 0
# @param rho Exponential decrease in incidence
# @return Incidence at time t
# @examples
# incidence.exp(0, lambda_0=0.05, rho=0.07)
# incidence.exp(c(-1, 0, 1), lambda_0=0.05, rho=0.07)
incidence.exp <- function(t, lambda_0, rho=1) lambda_0 * exp(-rho * t)

# Infection times function based on constant incidence.
#
# @param e A number between 0 and 1 (randomly generated), float or vector
# @param t Time, float
# @param p Constant prevalence
# @param lambda_0 Baseline incidence
# @return Infection time
# @examples
# e <- runif(10)
# infections.con(e, t=0, p=0.2, lambda_0=0.05)
# @export
infections.con <- function(e, t, p, lambda_0, rho=NA){
  # Note that you can work with e or 1 - e since e is Uniform(0, 1)
  infections <- t - p*e / ((1 - p) * lambda_0)
  return(infections)
}

# Infection times function based on linearly decreasing incidence.
#
# @param e A number between 0 and 1 (randomly generated), float or vector
# @param t Time, float
# @param p Constant prevalence
# @param lambda_0 Baseline incidence
# @param rho Linearly decreasing incidence parameter
# @return Infection time
# @examples
# e <- runif(10)
# infections.lin(e, t=0, p=0.2, lambda_0=0.05, rho=1e-3)
# @export
infections.lin <- function(e, t, p, lambda_0, rho){
  incidence <- lambda_0
  numerator <- incidence**2 + 2 * rho * p * e / (1 - p)
  numerator <- sqrt(numerator) - incidence

  infections <- t - numerator / rho
  return(infections)
}

# Infection times function based on constant incidence until bigT
# then linearly increasing after. Compared to the other infection functions
# we've switched to 1-e, which doesn't matter because e is uniform.
infections.lincon <- function(e, t, p, lambda_0, rho, bigT){
  incidence <- lambda_0
  estar <- 1-incidence * (1-p) / p * bigT
  efunc <- function(es){
    if(es < estar){
      numerator <- incidence**2 + 2 * rho * (p * (1-es)/(1-p) - incidence * bigT)
      numerator <- sqrt(numerator) - incidence
      infections <- (t - bigT) - numerator / rho
    } else {
      infections <- t - p*(1-es) / ((1 - p) * incidence)
    }
  }
  infections <- sapply(e, efunc)
  return(infections)
}

# Infection times function based on exponential decreasing incidence.
#
# @param e A number between 0 and 1 (randomly generated), float or vector
# @param t Time, float
# @param p Constant prevalence
# @param lambda_0 Baseline incidence
# @param rho Linearly decreasing incidence parameter
# @return Infection time
# @examples
# e <- runif(10)
# infections.exp(e, t=0, p=0.2, lambda_0=0.05, rho=0.07)
#
# @export
infections.exp <- function(e, t, p, lambda_0, rho){
  incidence <- lambda_0
  infections <- t - (1/rho) * log(rho*p*e/((1-p)*incidence) + 1)
  return(infections)
}

# Simulate infection times based on epidemiological parameters.
#
# @param n_sims Number of infection times to simulate
# @param type Incidence type: constant, linear, exponential
# @param p Constant prevalence
# @param lambda_0 Baseline incidence (at time t)
# @param rho Parameter for linear or exponentially decreasing incidence
# @param t Time t baseline (defaults to 0)
#
# @export
#
# @examples
# infs <- simInfections(type="constant", p=0.28, lambda_0=0.03)
# infs(10)
simInfections <- function(type, p, lambda_0, rho=NA, t=0, bigT=NA){

  inf.func <- function(n_sims){
    e <- runif(n_sims)

    if(type == "constant"){
      infs <- infections.con(e=e, p=p, lambda_0=lambda_0, t=t)
    } else {
      if(is.na(rho)) stop("Need to input rho.")
      if(type == "linear"){
        infs <- infections.lin(e=e, p=p, lambda_0=lambda_0, rho=rho, t=t)
      } else if(type == "exponential"){
        infs <- infections.exp(e=e, p=p, lambda_0=lambda_0, rho=rho, t=t)
      } else if(type == "linear-constant"){
        if(!is.na(bigT)) stop("Need to input bigT.")
        infs <- infections.lincon(e=e, p=p, lambda_0=lambda_0, rho=rho, t=t, bigT=bigT)
      } else {
        stop("Unrecognized incidence type. Must be one of constant, linear, or exponential.")
      }
    }
    return(infs)
  }
  return(inf.func)
}

# INFECTION FUNCTION FOR ARBITRARY INCIDENCE FUNCTION

# Function to simulate infection times from an arbitrary incidence
# function and with constant prevalence.
#
# @param p Constant prevalence value
# @param inc.function A function of time that returns incidence
# @param nsims Number of infection times to simulate
# @param dt The precision for numerical integration. May need to increase.
# @return A vector of infection durations
# @examples
# set.seed(1)
# simInfectionsCustom(p=0.2, nsims=15,
#                     inc.function=function(t) 0.03 - 0.001 * t)
# simInfectionsCustom(p=0.2, nsims=15,
#                     inc.function=function(t) 0.03 - 0.001 * t)
simInfectionsCustom <- function(p, inc.function, nsims=1000, dt=0.001){
  # TODO: Change this to infection duration in the generate.data function
  # Need to vectorize a scalar function
  if(length(inc.function(c(1, 2))) == 1) inc.function <- Vectorize(inc.function)

  # Simulate e's and the solve for T
  e <- runif(nsims)
  lh <- e * p/(1-p)
  ds <- seq(0, 100, by=dt)
  vals <- inc.function(ds) # TODO: Shift lambda for longitudinal
  int <- cumsum(vals*dt)
  if(min(lh) < min(int)) stop("Pick a larger dt for the integration")
  indexes <- sapply(lh, function(x) max(which(int < x)))
  times <- ds[indexes]
  return(times)
}
