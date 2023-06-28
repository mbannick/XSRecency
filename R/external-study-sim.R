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

#' This simulates many false recency rates using
#' tail probabilities from the phi parameters
#' by sampling from a binomial with a sample size
#' given in the study.
#'
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
#' @importFrom geepack geese
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
simulate_studies <- function(nsims, phi.func=NULL, ext_df=NULL){
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

#' Estimate omega_T from a model object
#' Can instead estimate mu hat if you pass in a sufficiently
#' long follow_T parameter (like 9)
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

#' @examples
#' pt.properties.est(ptest_times=runif(100, 1, 4),
#'                   ptest_delta=rbinom(100, 1, 0.5),
#'                   ptest_avail=rbinom(100, 1, 0.5))
pt.properties.est <- function(ptest_times,
                              ptest_delta,
                              ptest_avail,
                              dt=1/365.25){

  elem <- copy(ELEM)

  # Map testing times to indices for integration using dt
  closest <- data.table(
    ts_orig=-ptest_times,
    has_test=ptest_avail
  )
  closest[, id := .I]
  closest[, index := NA]

  get.closest.index <- function(t){
    if(is.na(t)){
      idx <- NA
    } else {
      idx <- ts_index[which.min(abs(t - ts_index$ts)), index]
    }
  }
  closest_index <- sapply(closest$ts_orig, get.closest.index)
  closest[, index := closest_index]
  closest <- merge(closest, ts_index, by="index", all.x=T)

  closest[, Ai := as.numeric((ts <= bigT_index[, ts]) & (has_test))]
  closest[, Bi := as.numeric((ts > bigT_index[, ts]) & (has_test))]

  # Calculate fraction of those with recent/non-recent prior tests
  p_A <- mean(closest$Ai)
  p_B <- mean(closest$Bi)

  elem[["p_A"]] <- p_A
  elem[["p_B"]] <- p_B

  if(p_A > 0){

    # Create point estimates and covariance matrix between
    # the phi estimates at different times
    mat_index <- ts_index[ts <= max(closest[Ai ==1, ts], bigT_index[, ts]),]
    cphi <- matrix.phi(model, mat_index, dt=dt)
    cphi_point <- cphi$cphi
    cphi_covar <- cphi$csum

    # Get ids for those who have a recent prior test, and their
    # corresponding time index
    idmap <- closest[Ai == 1, .(id, index, ts, Mi)]
    # Get unique combinations of the ids
    # for the covariance calculation
    idmap_covar <- combn(idmap$id, 2) %>% t %>% data.table
    setnames(idmap_covar, c("idX", "idY"))
    idmap_covar <- merge(idmap_covar, idmap[, .(id, index)], by.x="idX", by.y="id")
    setnames(idmap_covar, "index", "indexX")
    idmap_covar <- merge(idmap_covar, idmap[, .(id, index)], by.x="idY", by.y="id")
    setnames(idmap_covar, "index", "indexY")

    # Calculate omega_T conditional expectation and variance
    omega_ta     <- merge(idmap, cphi_point, by="index", all.x=TRUE)
    elem[["omega_TA"]] <- mean(omega_ta$phi)
    elem[["omega_TA_var"]] <- var(omega_ta$phi)

    # Calculate omega_T * T conditional expectation
    omega_ta <- omega_ta[, tastar := ts * phi]
    omega_ta <- omega_ta[, taneg := ts - phi]
    elem[["omega_TAstar"]] <- mean(omega_ta$tastar)
    elem[["den_omega"]] <- sum(omega_ta$taneg)

    # Calculate r_TiTi conditional expectation and variance
    idmap2 <- idmap
    idmap2[, indexX := index]
    idmap2[, indexY := index]

    r_Tii  <- merge(idmap2, cphi_covar, by=c("indexX", "indexY"), all.x=T)
    elem[["r_TA"]]   <- mean(r_Tii$rho)
    elem[["var_TA"]] <- var(r_Tii$rho)

    # Calculate r_TiTj conditional expectation
    r_Tij <- merge(idmap_covar, cphi_covar, by=c("indexX", "indexY"), all.x=T)
    elem[["r_TAprime"]] <- mean(r_Tij$rho)

    # Calculate r_TiTstar conditional expectation
    idmap3 <- idmap
    idmap3[, indexX := index]

    idmap3[, indexY := bigT_index[, index]]

    r_Tistar <- merge(idmap3, cphi_covar, by=c("indexX", "indexY"), all.X=T)
    elem[["r_TAstar"]] <- mean(r_Tistar$rho)

    # Calculate the expected time of prior tests
    ta <- closest[Ai == 1, ts_orig]
    elem[["mu_TA"]] <- mean(ta)
    elem[["var_TA"]] <- var(ta)
  }

  if(p_B > 0){

    # Calculate the expected time of prior tests
    tb <- closest[Bi == 1, ts_orig]
    elem[["mu_TB"]] <- mean(tb)
    elem[["var_TB"]] <- var(tb)

  }

  # Return list of elements that were estimated
  return(elem)
}

