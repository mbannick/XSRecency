# -------------------------------------------------------------------------
# EXTERNAL STUDY DATA SIMULATION TO ESTIMATE MU, OMEGA, BETA
# -------------------------------------------------------------------------

library(geepack)
library(data.table)
library(MASS)
.datatable.aware=TRUE

# number of long infecteds for beta study
N_LONG_INFECT=1500
YEAR2DAY=365.25

expit <- function(x) exp(x) / (1 + exp(x))

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

get.cubic.ts <- function(minT, maxT, dt){
  # generate a sequence of duration times, starting at 0 and
  # up to max t, just to get the phi_hat
  ts <- seq(minT, maxT, by=dt)
  ts <- cbind(rep(1, length(ts)),
              ts,
              ts**2,
              ts**3)
  return(ts)
}

#' Create a phi matrix of d x d based on draws from the phi.hat function
#'
#' @param d A vector
#' @param n_approx Number of "draws" for phi.hat
#' @param phi.hat A function that creates draws from a model object using d and any additional args ...
#' @param ... Additional arguments to pass to phi.hat
matrix.phi <- function(model, n_approx, minT, maxT, dt=0.001){

  ts <- get.cubic.ts(minT=minT, maxT=maxT, dt=dt)
  # get the linear predictor
  lin_predictor <- ts %*% matrix(model$beta)

  phi <- expit(lin_predictor)

  lin_var <- ts %*% model$vbeta %*% t(ts)
  delta_g <- exp(lin_predictor) / (1 + exp(lin_predictor))**2
  jacob_g <- diag(c(delta_g))
  phi_var <- jacob_g %*% lin_var %*% jacob_g

  return(list(point=phi, covar=phi_var))
}

#' Estimate omega_T from a model object
#' Can instead estimate mu hat if you pass in a sufficiently
#' long follow_T parameter (like 9)
#'
#' @export
integrate.phi <- function(model, minT=0, maxT=12, dt=0.001){
  ts <- get.cubic.ts(minT=minT, maxT=maxT, dt=dt)

  # get the linear predictor
  lin_predictor <- ts %*% matrix(model$beta)

  # expit the linear predictor
  phi <- expit(lin_predictor)

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
assay.properties.est <- function(study, bigT, tau, last_point=TRUE, dt=0.01,
                                 ptest_times=NULL){
  model <- fit.cubic(recent=study$recent,
                     durations=study$durations,
                     id=study$id)

  # get mu and omega
  if(last_point) tau <- max(study$durations)

  mu_sim <- integrate.phi(model, minT=0, maxT=tau, dt=dt)
  omega_sim <- integrate.phi(model, minT=0, maxT=bigT, dt=dt)

  # If we have prior test results, we need to calculate the
  # following quantities that allow us to calculate the variance
  # of the enhanced adjusted estimator that incorporates prior test results.
  if(!is.null(ptest_times)){
    if(any(is.na(ptest_times))) stop("Not implemented for q < 1.0, or any misspecification scenario.")

    # Convert testing times to indices for integration using dt
    time_indices <- round(-ptest_times/dt)

    # Create covariance matrix between the phi estimates at different times
    rho_mat <- matrix.phi(model, minT=0, maxT=tau, dt=dt)

    # Define function to integrate over rectangles within the rho matrix
    # starting at (1,1).
    integrate_cov <- function(i, j) sum(rho_mat$covar[1:j, 1:i]) * dt^2

    # Calculate the expectation of the integrated covariance matrix
    # for individual i's
    r_Tii <- mapply(integrate_cov, i=time_indices, j=time_indices) %>% mean

    # Create grid over which to integrate for different i < j pairs
    covar_grid <- expand.grid(i=time_indices, j=time_indices)
    covar_grid <- covar_grid[covar_grid$i < covar_grid$j, ]

    # Calculate the expectation of the integrated covariance matrix
    # for i < j pairs.
    r_Tij <- mapply(integrate_cov, i=covar_grid$i, j=covar_grid$j) %>% mean

    # Define function for integrating across
    integrate_point <- function(i) sum(rho_mat$point[1:i]) * dt

    # Calculate the omega_Ti for each prior test time,
    # then get its expectation and variance
    # (this expectation and variance are with respect to the variability
    # in the trial data T_i's, not the external study data.)
    omega_Ti <- sapply(time_indices, integrate_point)
    omega_Ti_est <- mean(omega_Ti)
    omega_Ti_var <- var(omega_Ti)

  } else {
    r_Tii <- NULL
    r_Tij <- NULL
    omega_Ti_est <- NULL
    omega_Ti_var <- NULL
  }

  result <- list(mu_est=mu_sim$est,
                 mu_var=mu_sim$var,
                 omega_est=omega_sim$est,
                 omega_var=omega_sim$var,
                 omega_Ti_est=omega_Ti_est,
                 omega_Ti_var=omega_Ti_var,
                 r_Tii=r_Tii,
                 r_Tij=r_Tij,
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
#' @param ext_df A dataset with column "duration" and column for binary "recent" indicator
#' @param max_FRR The maximum duration allowed
#' @param last_point Integrate the mean window period to the last observed
#'   duration in the dataset, rather than tau
#' @return A list of estimated mean window period \eqn{\mu} and its variance,
#'   estimated MDRI \eqn{\Omega_{T^*}} and its variance, and
#'   estimated FRR \eqn{\beta_{T^*}} and its variance.
#' @examples
#' set.seed(0)
#' assay.properties.nsim(n_sims=1, phi.func=function(t) 1 - pgamma(t, 1, 1.5),
#'                      bigT=2, tau=12)
assay.properties.nsim <- function(n_sims, phi.func, bigT, tau,
                                  ext_FRR=FALSE, ext_df=NULL,
                                  max_FRR=NULL, last_point=FALSE,
                                  ptest_times=NULL){
  studies <- simulate.studies(n_sims, phi.func, ext_df=ext_df)
  if(is.null(ptest_times)){
    result <- sapply(studies, function(x) assay.properties.est(study=x,
                                                               bigT=bigT,
                                                               tau=tau,
                                                               last_point=last_point,
                                                               ptest_times=ptest_times))
  } else {
    mfunc <- function(s, t) assay.properties.est(
      study=s, bigT=bigT, tau=tau, last_point=last_point, ptest_times=t
    )
    result <- mapply(FUN=mfunc, s=studies, t=ptest_times)
  }

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
