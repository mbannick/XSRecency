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

# This simulates the false recency rate using
# tail probabilities from the phi parameters
# by sampling from a binomial with a sample size
# given in the study.
simulate.beta <- function(phi.func, minT, maxT, recent=NULL){
  N <- N_LONG_INFECT

  if(is.null(recent)){
    inf_times <- runif(n=N, min=minT, max=maxT)
    recent <- rbinom(n=N, size=1, prob=phi.func(inf_times))
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

# This simulates many false recency rates using
# tail probabilities from the phi parameters
# by sampling from a binomial with a sample size
# given in the study.
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

# Helper function for simulating many studies
# that mimic the data from Duong et al. 2015.
# Fits a GEE model with piece-wise linear
# parameterization and Poisson marginal mean / variance.
# NOTE: This is not the same thing as fitting the GEE model for the phi function.
# This is simulating the duration of infections in an external dataset for estimating recency assay parameters.
#
# @param df Data frame with sample number column (samp),
#   days column, and num.samples column.
# @param knot Knot location for piecewise-linear function
# @importFrom geepack geese
# @return List of start durations, numbers of samples, and model coefficients
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

# Function to simulate one study
# that mimics the data from Duong et al. 2015
#
# @param days Duration on date of first sample
# @param num.samples Distribution of number of samples
# @param coefs Coefficients from the piecewise linear model
# @param knot Knot location for piecewise linear model
# @param phi.func Optional recency test-positive function
# @return List of data frames with an id, time, and recency indicator
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

#' Function to simulate studies
#' that mimic the infection duration data from Duong et al. 2015 (or another external dataset)
#'
#' @export
#' @param nsims Number of study simulations
#' @param phi.func Optional recency test-positive function
#' @param ext_df A dataset with column "duration" and column for binary "recent" indicator
#' @return List of data frames with an id, time, and recency indicator
#'
#' @examples
#' set.seed(0)
#' study <- simExternal(nsims=1, phi.func=function(t) 1 - pgamma(t, 1, 1.5))
#' head(study)
simExternal <- function(nsims, phi.func=NULL, ext_df=NULL){
  if(!is.null(ext_df)){
    df <- ext_df
  } else {
    df <- XSRecency:::duong
  }
  mod <- fit.model(df)
  sims <- replicate(nsims, simulate.study(mod$days, mod$num.samples, mod$rate,
                                          phi.func=phi.func),
                    simplify=FALSE)

  if(nsims == 1){
    sims <- sims[[1]]
  }
  return(sims)
}

# Fits a cubic model to indicators and durations
# and returns the model object.
fit.cubic <- function(recent, durations, id){
  durations2 <- durations**2
  durations3 <- durations**3

  # fit a logistic regression with cubic polynomials
  model <- geese(recent ~ durations + durations2 + durations3,
                 family=binomial(link="logit"), id=id, corstr="exchangeable")
  return(model)
}

get.ts <- function(minT, maxT, dt){
  index <- seq(minT/dt, maxT/dt, by=1)
  ts <- index*dt
  return(data.table(index=index,
                    ts=ts))
}

get.cubic.ts <- function(ts){
  # generate a sequence of duration times, starting at 0 and
  # up to max t, just to get the phi_hat
  ts.cube <- cbind(rep(1, length(ts)),
              ts,
              ts**2,
              ts**3)
  return(ts.cube)
}

# Estimate omega_T from a model object
# Can instead estimate mu hat if you pass in a sufficiently
# long follow_T parameter (like 9)
#
integrate.phi <- function(model, ts, dt=0.01){
  ts <- get.cubic.ts(ts)

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

ELEM <- list(
  omega_TA     = 0,
  omega_TA_var = 0,
  omega_TAstar = 0,
  den_omega    = 0,
  mu_TA        = 0,
  var_TA       = 0,
  r_TA         = 0,
  r_TAprime    = 0,
  r_TAstar     = 0,
  mu_TB        = 0,
  var_TB       = 0,
  p_A          = 0,
  p_B          = 0
)

# Function that simulates data with a phi function based on external infection duration study data from Duong et al. 2015.
# Data available here from Duong et al. 2015: https://doi.org/10.1371/journal.pone.0114947.s001.
#
# @param n_sims Number of simulations
# @param phi.func A test-recent positive function of t
# @param bigT The time cutoff value designating true recent versus false recent
# @param tau The maximum duration of infection where a subject could have a false-positive for recent infection
# @param ext_FRR Whether or not to get FRR from the Duong et al. 2015
#               study or to calculate it from a simpler, separate study, passed in ext_df
# @param ext_df A dataset with column "duration" and column for binary "recent" indicator
# @param max_FRR The maximum duration allowed
# @param last_point Integrate the mean window period to the last observed duration in the dataset, rather than tau
#
# @import geepack
#
# @returns A list of studies, and FRR estimates
# @examples
# set.seed(0)
# assay.nsim.pt(n_sims=1, phi.func=function(t) 1 - pgamma(t, 1, 1.5),
#               bigT=2, tau=12)
assay.nsim.pt <- function(n_sims, phi.func, tau, bigT,
                          ext_FRR=FALSE, ext_df=NULL,
                          max_FRR=NULL){

  studies <- simExternal(n_sims, phi.func, ext_df=ext_df)

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

  return(list(studies=studies, beta_sim=beta_sim))
}

