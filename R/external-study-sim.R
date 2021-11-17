# -------------------------------------------------------------------------
# EXTERNAL STUDY DATA SIMULATION TO ESTIMATE MU, OMEGA, BETA
# -------------------------------------------------------------------------

library(geepack)
library(data.table)
.datatable.aware=TRUE

# number of long infecteds for beta study
N_LONG_INFECT=1500
YEAR2DAY=365.25

#' This simulates the false recency rate using
#' tail probabilities from the phi parameters
#' by sampling from a binomial with a sample size
#' given in the study.
#'
#' @export
simulate.beta <- function(phi.func, minT, maxT, recent=NULL){
  N <- N_LONG_INFECT

  if(is.null(recent)){
    inf_times <- runif(n=N, min=minT, max=maxT)
    recent <- rbinom(n=N, size=1, p=phi.func(inf_times))
  } else {
    # If we already have recency indicators, then don't
    # apply the phi function.

    # However, up-sample or down-sample
    # them so that we have the amount that we need.
    recent <- sample(recent, size=N, replace=TRUE)
  }

  beta <- sum(recent) / N
  beta_var <- beta * (1 - beta) / N
  return(list(est=beta, var=beta_var))
}

#' This simulates many false recency rates using
#' tail probabilities from the phi parameters
#' by sampling from a binomial with a sample size
#' given in the study.
#'
#' @export
simulate.nbeta <- function(nsims, phi.func, minT, maxT, studies=NULL){
  if(is.null(studies)){
    result <- replicate(nsims, simulate.beta(phi.func=phi.func,
                                             minT=minT, maxT=maxT))
  } else {
    recents <- list()
    i <- 1
    for(study in studies){
      recents[[i]] <- study[study$durations >= minT & study$durations <= maxT, "recent"]
      i <- i + 1
    }
    result <- sapply(recents, function(x) simulate.beta(phi.func=phi.func,
                                                        minT=minT,
                                                        maxT=maxT,
                                                        recent=x))
  }

}

#' Helper function for simulating many studies
#' that mimic the data from Duong et al. 2015.
#' Fits a GEE model with piece-wise linear
#' parameterization and Poisson marginal mean / variance.
#'
#' @param df Data frame with sample number column (samp),
#'   days column, and num.samples column.
#' @param knot Knot location for piecewise-linear function
#' @return List of start durations, numbers of samples, and model coefficients
fit.model <- function(df, knot=5){
  df <- data.table(df)
  # Indicator variable for sample 5+
  df[, samp.5 := samp >= knot]

  # Fit GEE model for the gap times
  mod <- geese(gap ~ samp + samp.5 + samp.5*samp, id=id.key, data=df,
             family=poisson(link="log"), corstr="exchangeable")
  # Get only the first days of sampling
  first.day <- df[samp == 1]

  return(list(days=first.day$days,
              num.samples=first.day$num.samples,
              rate=mod$beta))
}

#' Function to simulate one study
#' that mimics the data from Duong et al. 2015
#'
#' @param days Duration on date of first sample
#' @param num.samples Distribution of number of samples
#' @param coefs Coefficients from the piecewise linear model
#' @param knot Knot location for piecewise linear model
#' @param phi.func Optional recency test-positive function
#' @return List of data frames with an id, time, and recency indicator
simulate.study <- function(days, num.samples, coefs, knot=5,
                           phi.func=NULL){

  # Sample the first day
  day1 <- sample(days, size=length(days), replace=TRUE)

  # Sample the number of samples
  nums <- sample(num.samples, size=length(days), replace=TRUE)

  # Get a vector of the sample numbers, starting from the second one
  samp.num <- lapply(nums, function(x) 2:x)

  # Get a sample design matrix
  get.samp.dmat <- function(samp) return(cbind(1, samp, samp>=knot, (samp>=knot)*samp))
  dmats <- lapply(samp.num, get.samp.dmat)

  # Compute the predicted rates
  rates <- lapply(dmats, function(x) exp(x %*% coefs)[, 1])

  # Get the gap-times for each
  gaps <- lapply(rates, function(x) rpois(n=length(x), lambda=x))

  # Add on the start day
  days.gap <- mapply(function(x, y) c(x, y), x=day1, y=gaps)

  # Calculate sampling times
  days.tot <- lapply(days.gap, cumsum) %>% unlist
  time.tot <- days.tot/YEAR2DAY

  # Get ids and put into data frame
  ids <- 1:length(days)
  ids <- rep(ids, nums)

  df <- data.frame(id=ids,
                   durations=time.tot)

  # If there is a phi function, apply it to the durations
  if(!is.null(phi.func)){
    recent <- rbinom(n=length(time.tot), size=1, prob=phi.func(time.tot))
    df$recent <- recent
  }
  return(df)
}

#' Function to simulate many studies
#' that mimic the data from Duong et al. 2015
#'
#' @export
#' @param nsims Number of study simulations
#' @param phi.func Optional recency test-positive function
#' @return List of data frames with an id, time, and recency indicator
simulate.studies <- function(nsims, phi.func=NULL, ext_df=NULL){
  if(!is.null(ext_df)){
    df <- ext_df
  } else {
    df <- XSRecency:::duong
  }
  mod <- fit.model(df)
  sims <- replicate(nsims, simulate.study(mod$days, mod$num.samples, mod$rate,
                                          phi.func=phi.func),
                    simplify=FALSE)
  return(sims)
}

#' Fits a cubic model to indicators and durations
#' and returns the model object.
#'
#' @export
fit.cubic <- function(recent, durations, id){
  durations2 <- durations**2
  durations3 <- durations**3

  # fit a logistic regression with cubic polynomials
  model <- geese(recent ~ durations + durations2 + durations3,
                 family=binomial(link="logit"), id=id, corstr="exchangeable")
  return(model)
}

#' Create phi probability predictions
#' based on a model object and d durations
#'
#' @export
phi.hat <- function(d, model){
  dmat <- cbind(
    rep(1, length(d)),
    d, d**2, d**3
  )
  lin <- dmat %*% matrix(model$beta)
  return(exp(lin) / (1 + exp(lin)))
}

#' Estimate omega_T from a model object
#' Can instead estimate mu hat if you pass in a sufficiently
#' long follow_T parameter (like 9)
#'
#' @export
integrate.phi <- function(model, minT=0, maxT=12){

  # generate a sequence of duration times, starting at 0 and
  # up to max t, just to get the phi_hat
  dt <- 0.001
  ts <- seq(minT, maxT, by=dt)
  ts <- cbind(rep(1, length(ts)),
              ts,
              ts**2,
              ts**3)

  # get the linear predictor
  lin_predictor <- ts %*% matrix(model$beta)

  # expit the linear predictor
  phi <- exp(lin_predictor) / (1 + exp(lin_predictor))

  # integrate the linear predictor over the times
  estimate <- sum(phi * dt)

  # get the variance of omega_t
  # get the gradient of the linear predictor and then multiply
  # by the time diff, so it's the sum of the gradients
  grad <- exp(lin_predictor) / (1 + exp(lin_predictor))**2
  grad <- t(ts) %*% grad
  grad <- grad * dt

  vc <- model$vbeta

  # the variance of omega_t
  variance <- c(t(grad) %*% vc %*% grad)

  return(list(est=estimate, var=variance))
}

#' Use an assay df with recent (yes/no) and duration columns
#' to calculate the mean window period and MDRI of the assay
#'
#' @export
#' @param study Data frame with recent and durations variables
#' @param bigT The T^* time
#' @param tau The maximum time
#' @returns List of properties and their variances
assay.properties.est <- function(study, bigT, tau, last_point=TRUE){
  model <- fit.cubic(recent=study$recent,
                     durations=study$durations,
                     id=study$id)

  # get mu and omega
  if(last_point) tau <- max(study$durations)

  mu_sim <- integrate.phi(model, minT=0, maxT=tau)
  omega_sim <- integrate.phi(model, minT=0, maxT=bigT)

  result <- list(mu_est=mu_sim$est,
                 mu_var=mu_sim$var,
                 omega_est=omega_sim$est,
                 omega_var=omega_sim$var,
                 beta_est=NA,
                 beta_var=NA)
  return(result)
}

#' Function that simulates the mean window period, mean duration of recent
#' infection, and false recent rate (FRR) based on external study data.
#'
#' The external study data for mean window and MDRI is available here from
#' Duong et al. 2015:
#' https://doi.org/10.1371/journal.pone.0114947.s001
#'
#' @export
#' @param n_sims Number of simulations
#' @param phi.func A test-recent positive function of t
#' @param bigT The time cutoff value designating true recent versus false recent
#' @param tau The maximum duration of infection where a subject
#'   could have a false-positive for recent infection
#' @param ext_FRR Whether or not to get FRR from the Duong et al. 2015
#'   study or to calculate it from a simpler, separate study, passed in \code{ext_df}
#' @param ext_FRR A dataset with column "duration" and column for binary "recent" indicator
#' @param max_FRR The maximum duration allowed
#' @param last_point Integrate the mean window period to the last observed
#'   duration in the dataset, rather than tau
#' @return A list of estimated mean window period \eqn{\mu} and its variance,
#'   estimated MDRI \eqn{\Omega_{T^*}} and its variance, and
#'   estimated FRR \eqn{\beta_{T^*}} and its variance.
#' @examples
#' set.seed(0)
#' assay.properties.sim(n_sims=1, phi.func=function(t) 1 - pgamma(t, 1, 1.5),
#'                      bigT=2, tau=12)
assay.properties.nsim <- function(n_sims, phi.func, bigT, tau,
                                  ext_FRR=FALSE, ext_df=NULL,
                                  max_FRR=NULL, last_point=FALSE){
  studies <- simulate.studies(n_sims, phi.func, ext_df=ext_df)
  result <- sapply(studies, function(x) assay.properties.est(study=x,
                                                             bigT=bigT,
                                                             tau=tau,
                                                             last_point=last_point))
  if(ext_FRR){
    frr_studies <- studies
  } else {
    frr_studies <- NULL
  }
  if(!is.null(max_FRR)){
    b.maxT <- max_FRR
  } else {
    b.maxT <- tau
  }
  beta_sim <- simulate.nbeta(nsims=n_sims, phi.func=phi.func,
                             minT=bigT, maxT=b.maxT, studies=frr_studies)
  result[c("beta_est", "beta_var"), ] <- beta_sim

  result <- t(result)
  result <- data.table(result)
  columns <- colnames(result)
  cols <- lapply(columns, function(x) unlist(result[[x]]))
  df <- do.call(cbind, cols) %>% data.table
  names(df) <- columns
  return(df)
}
