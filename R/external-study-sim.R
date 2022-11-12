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
matrix.phi <- function(model, minT, maxT, dt=0.01){

  # Get the set of t's that we want, in design matrix form
  ts <- get.cubic.ts(minT=minT, maxT=maxT, dt=dt)
  # Get the linear predictor
  lin_predictor <- ts %*% matrix(model$beta)
  # Exponentiate the linear predictor
  phi <- expit(lin_predictor)

  # Calculate cumulative phi and convert to a data frame
  cphi <- cumsum(phi * dt)
  cphi_d <- data.table(
    index=1:length(cphi),
    phi=cphi
  )

  # We have by delta method Var(T \beta) = T %*% Var(\beta) %*% T^t
  # Also by delta method, Var(expit(T \beta)) = J %*% T %*% Var(\beta) %*% T^t J
  # where J is the Jacobian of the expit() transformation.

  # For efficiency's sake, we can instead first calculate J %*% T
  # since (J %*% T)^t = T^t %*% J^t = T^t %*% J
  # and then calculate the Cholesky decomposition of Var(\beta) so that we
  # can write Var(\beta) = U^t U where U is the upper triangular
  # matrix from the decomposition (what chol function returns).

  # Then we have, with R = J %*% T %*% U^t,
  # the final variance-covariance matrix is R %*% R^t.

  # 1. Var(\beta) and Cholesky decomposition
  U <- chol(model$vbeta)

  # 2. Jacobial of transformation
  delta_g <- exp(lin_predictor) / (1 + exp(lin_predictor))**2
  J <- diag(c(delta_g))

  # 3. Calculate R matrix
  R <- J %*% ts %*% t(U)

  # 4. Calculate final variance
  phi_var <- R %*% t(R)

  # This is the old way of calculating it -- which is extremely expensive
  # and might (probably will) crash your R.

  # lin_var <- ts %*% model$vbeta %*% t(ts)
  # delta_g <- exp(lin_predictor) / (1 + exp(lin_predictor))**2
  # jacob_g <- diag(c(delta_g))
  # phi_var2 <- jacob_g %*% lin_var %*% jacob_g

  # Calculate the 2-d integral of the variance-covariance matrix
  # at all grid points by using cumsum on the rows and columns.
  csum <- t(apply(apply(phi_var * dt^2, 2, cumsum), 1, cumsum))

  # Reshape the matrix into a data frame so that it's easier to work with
  # and merge onto the id grid.
  csum <- data.table(csum)
  csum <- csum[, indexX := .I]
  csum_d <- melt(csum, id.vars="indexX")
  csum_d <- csum_d[, variable := as.numeric(gsub("V", "", variable))]
  setnames(csum_d, c("variable", "value"), c("indexY", "rho"))

  return(list(cphi=cphi_d, csum=csum_d))
}

#' Estimate omega_T from a model object
#' Can instead estimate mu hat if you pass in a sufficiently
#' long follow_T parameter (like 9)
#'
#' @export
integrate.phi <- function(model, minT=0, maxT=12, dt=0.01){
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
assay.properties.est <- function(study, bigT, tau, last_point=TRUE, dt=1/365.25,
                                 ptest_times=NULL,
                                 ptest_delta=NULL,
                                 ptest_avail=NULL,
                                 ri=ri){
  start.time <- Sys.time()
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
  if(!all(is.null(ptest_times))){

    # Convert testing times to indices for integration using dt
    time_indices <- round(-ptest_times/dt)
    Ai <- as.numeric(((-ptest_times) <= bigT) & (ptest_avail))
    Bi <- as.numeric(((-ptest_times) > bigT) & (ptest_avail))

    # Calculate fraction of those with recent/non-recent prior tests
    p_A <- mean(Ai)
    p_B <- mean(Bi)

    if(p_A == 0){

      omega_TA     <- 0
      omega_TA_var <- 0
      omega_TAstar <- 0

      r_TA         <- 0
      var_TA       <- 0
      r_TAprime    <- 0
      r_TAstar     <- 0

      nATO         <- 0

    } else {

      # Create point estimates and covariance matrix between
      # the phi estimates at different times
      cat("Creating rho matrix\n")
      cphi <- matrix.phi(model, minT=0, maxT=tau, dt=dt)
      cphi_point <- cphi$cphi
      cphi_covar <- cphi$csum

      # Get ids for those who have a test, and their
      # corresponding time index
      has_test <- which(as.logical(Ai))
      idmap <- data.table(
        id=has_test,
        index=time_indices[has_test]
      )
      # Get unique combinations of the ids
      # for the covariance calculation
      idmap_covar <- combn(idmap$id, 2) %>% t %>% data.table
      setnames(idmap_covar, c("idX", "idY"))
      idmap_covar <- merge(idmap_covar, idmap, by.x="idX", by.y="id")
      setnames(idmap_covar, "index", "indexX")
      idmap_covar <- merge(idmap_covar, idmap, by.x="idY", by.y="id")
      setnames(idmap_covar, "index", "indexY")

      # Calculate omega_T conditional expectation and variance
      omega_ta     <- merge(idmap, cphi_point, by="index", all.x=TRUE)
      omega_TA     <- mean(omega_ta$phi)
      omega_TA_var <- var(omega_ta$phi)

      # Calculate omega_T * T conditional expectation
      omega_ta <- omega_ta[, tastar := index * dt * phi]
      omega_TAstar <- mean(omega_ta$tastar)

      # Calculate r_TiTi conditional expectation and variance
      idmap2 <- idmap
      idmap2 <- idmap2[, indexX := index]
      idmap2 <- idmap2[, indexY := index]
      cat("Merge for TiTi\n")
      r_Tii  <- merge(idmap2, cphi_covar, by=c("indexX", "indexY"), all.x=T)
      r_TA   <- mean(r_Tii$rho)
      var_TA <- var(r_Tii$rho)

      # Calculate r_TiTj conditional expectation
      cat("Merge for TiTj\n")
      r_Tij <- merge(idmap_covar, cphi_covar, by=c("indexX", "indexY"), all.x=T)
      r_TAprime <- mean(r_Tij$rho)

      # Calculate r_TiTstar conditional expectation
      idmap3 <- idmap
      idmap3 <- idmap3[, indexX := index]
      idmap3 <- idmap3[, indexY := round(bigT/dt)]
      cat("Merge for TiTstar\n")
      r_Tistar <- merge(idmap3, cphi_covar, by=c("indexX", "indexY"), all.X=T)
      r_TAstar <- mean(r_Tistar$rho)

      # To compare variance terms in debugging
      nATO <- sum(omega_ta$phi)
    }

    if(p_B == 0){

      mu_TB  <- 0
      var_TB <- 0
      nB     <- 0
      nBT    <- 0

    } else {

      # Calculate the expected time of prior tests
      tb <- -ptest_times[as.logical(Bi)]
      mu_TB <- mean(tb)
      var_TB <- var(tb)

      # To compare variance terms for debugging
      nB <- sum(Bi)
      nBT <- sum(-ptest_times[as.logical(Bi)])

    }

  } else {
    r_TA <- NULL
    r_TAprime <- NULL
    r_TAstar <- NULL
    omega_TA <- NULL
    omega_TAstar <- NULL
    mu_TB <- NULL
    var_TB <- NULL
    p_A <- NULL
    p_B <- NULL
    nB <- NULL
    nBT <- NULL
    nATO <- NULL
  }

  result <- list(
    mu_est=mu_sim$est,
    mu_var=mu_sim$var,
    omega_est=omega_sim$est,
    omega_var=omega_sim$var,
    beta_est=NA,
    beta_var=NA,
    r_TA=r_TA,
    r_TAprime=r_TAprime,
    r_TAstar=r_TAstar,
    omega_TA=omega_TA,
    omega_TA_var=omega_TA_var,
    omega_TAstar=omega_TAstar,
    mu_TB=mu_TB,
    var_TB=var_TB,
    p_A=p_A,
    p_B=p_B,
    nB=nB,
    nBT=nBT,
    nATO=nATO
  )
  end.time <- Sys.time()
  print(end.time - start.time)
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
                                  ptest_times=NULL, ptest_delta=NULL,
                                  ptest_avail=NULL, ri=NULL){

  studies <- simulate.studies(n_sims, phi.func, ext_df=ext_df)
  if(is.null(ptest_times)){
    result <- sapply(studies, function(x) assay.properties.est(
      study=x,
      bigT=bigT,
      tau=tau,
      last_point=last_point,
      ptest_times=ptest_times,
      ptest_delta=ptest_delta,
      ptest_avail=ptest_avail,
      ri=ri
      ))
  } else {
    mfunc <- function(s, t, d, a, r) assay.properties.est(
      study=s, bigT=bigT, tau=tau, last_point=last_point,
      ptest_times=t, ptest_delta=d, ptest_avail=a, ri=r
    )
    result <- mapply(FUN=mfunc, s=studies,
                     t=ptest_times, d=ptest_delta, a=ptest_avail, r=ri)
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
