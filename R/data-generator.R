library(magrittr)

#' Generates raw data of number positive, negative, and screened
#' based on study parameters that can be used later to simulate infection times.
#'
#' @param n_sims Number of simulations
#' @param n Number of observations per time point
#' @param prevalence Constant prevalence
#' @param e.func Function to simulate infection time from a 0-1 value
#'
#' @examples
#' e.func <- function(e) infections.lin(e, t=0, p=0.29, lambda=0.032, rho=0.07)
#' params <- get.gamma.params(window=200/365.25, shadow=191/365.25)
#' phi.func <- function(t) 1-pgamma(t, shape = params[1], rate = params[2])
#' sim <- sim.screening.generator(prevalence=0.29, e.func=e.func, phi.func=phi.func)
#' sim(100)
#'
#' @export
sim.screening.generator <- function(prevalence, e.func, phi.func){
  sim.screening <- function(n, return_n=TRUE){

    # number of positive subjects
    pos <- rbinom(n=n, size=1, prob=prevalence) %>% sort
    npos <- sum(pos)

    # generate infection time
    e <- runif(n=npos)
    inf.time <- sapply(e, e.func)

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
#' e.func <- function(e) infections.lin(e, t=0, p=0.29, lambda=0.032, rho=0.07)
#' params <- get.gamma.params(window=200/365.25, shadow=191/365.25)
#' phi.func <- function(t) 1-pgamma(t, shape = params[1], rate = params[2])
#' sim <- sim.screening.generator(prevalence=0.29, e.func=e.func, phi.func=phi.func)
#' df <- sim(100)
#' sim.pt <- sim.pt.generator(ptest.dist=function(u) runif(1, 0, 4),
#'                          ptest.prob=function(u) 0.5)
#' sim.pt(df[df$di == 1,])
#'
#' @export
sim.pt.generator <- function(ptest.dist, ptest.prob, ptest.dist2=NULL){
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

#' Modifies the prior testing data for recall bias.
#'
#' @param t_noise Function to add noise to a single prior testing time.
#'   Must be a function of t, the prior testing time.
#' @param d_misrep Probability of given a positive prior test, individual
#'   reports negative prior test.
#' @param p_misrep Probability of not reporting a prior test result,
#'   among those with positive prior test.
#' @param t_range Range of times in which to keep prior testing times.
#' @examples
#' e.func <- function(e) infections.lin(e, t=0, p=0.29, lambda=0.032, rho=0.07)
#' params <- get.gamma.params(window=200/365.25, shadow=191/365.25)
#' phi.func <- function(t) 1-pgamma(t, shape = params[1], rate = params[2])
#' sim <- sim.screening.generator(prevalence=0.29, e.func=e.func, phi.func=phi.func)
#' df <- sim(100)
#' sim.pt <- sim.pt.generator(ptest.dist=function(u) runif(1, 0, 4),
#'                          ptest.prob=function(u) 0.5)
#' pt.df <- sim.pt(df[df$di == 1,])
#' t_noise <- function(t) max(0, t + rnorm(n=1, sd=0.5))
#' modify.pt <- modify.pt.generator(t_noise=t_noise, d_misrep=1, p_misrep=0, t_range=c(1, 2))
#' pt.df2 <- modify.pt(pt.df)
#'
#' @export
modify.pt.generator <- function(t_noise=function(t) t, t_range=NULL,
                              d_misrep=0, p_misrep=0, t_total_exclude=NULL){
  modify.pt <- function(df){

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
  return(modify.pt)
}

#' Simulate recency indicators based on the data, true phi function,
#' and the infection incidence function.
#'
#' Also possibly simulate prior test results.
#'
#' @param sim_data Outputs from the `generate.raw.data` function
#' @param infection.function Function that simulates the infection time
#' @param phi.func Test positive function for recency assay
#' @param baseline_incidence Baseline incidence value (time 0)
#' @param prevalence Constant prevalence value
#' @param rho A parameter to the `infection.function` for how quickly incidence changes
#' @param incidence.function An incidence function in place of `infection.function`.
#'   Cannot pass both arguments.
#' @param ptest.dist A prior test result distribution function, which must have
#'   be a function of a scalar `n` and vector `u` which is a vector
#'   of infection durations (which will be of length n),
#'   though it can ignore the `u` argument inside the function.
#'   An example is `rnorm(n, mean=u)`.
#' @param ptest.prob Probability of prior test result being available.
#'   Can instead be a function of `u`, the infection duration.
#' @param t_noise Function to add noise to a single prior testing time.
#'   Must be a function of t, the prior testing time.
#' @param d_misrep Probability of given a positive prior test, individual
#'   reports negative prior test.
#' @param q_misrep Probability of not reporting a prior test result.
#' @param ... Additional arguments to the `simulate.infection.times` function if you're
#'   passing an `incidence.function` instead of an `infection.function`.
#' @importFrom data.table data.table
#' @import magrittr
#'
#' @examples
#' dat <- generate.raw.data(n_sims=10, n=100, prevalence=0.5, times=c(0, 0))
#' sim <- simulate.recent(sim_data=dat, infection.function=infections.con,
#'                        baseline_incidence=0.05, prevalence=0.5, rho=NA,
#'                        phi.func=function(t) 1 - pgamma(t, 1, 2),
#'                        times=c(0, 1), summarize=TRUE,
#'                        ptest.dist=function(n, u) runif(n, 0.4, 1.0),
#'                        ptest.prob=1.0, bigT=2,
#'                        t_range=c(0.4, 1.0))
#' sim <- simulate.recent(sim_data=dat, infection.function=infections.con,
#'                        baseline_incidence=0.05, prevalence=0.5, rho=NA,
#'                        phi.func=function(t) 1 - pgamma(t, 1, 2),
#'                        times=c(0, 1), summarize=TRUE,
#'                        ptest.dist=function(n) runif(n, 0.0, 5.0),
#'                        ptest.prob=0.5, bigT=2)
#' dat <- generate.raw.data(n_sims=10, n=100, prevalence=0.5, times=c(0))
#' sim <- simulate.recent(sim_data=dat, infection.function=infections.con,
#'                        baseline_incidence=0.05, prevalence=0.5, rho=NA,
#'                        phi.func=function(t) 1 - pgamma(t, 1, 2),
#'                        summarize=TRUE,
#'                        ptest.dist=function(n) runif(n, 0.0, 5.0),
#'                        ptest.prob=0.5, bigT=2)
simulate.recent <- function(sim_data, infection.function=NULL,
                            phi.func, baseline_incidence, prevalence, rho, summarize=TRUE,
                            incidence.function=NULL,
                            ptest.dist=NULL, ptest.prob=1.0,
                            t_range=NULL, t_noise=NULL,
                            d_misrep=0.0, q_misrep=0.0, p_misrep=0.0,
                            bigT=NULL,
                            ptest.dist2=NULL,
                            exclude_pt_bigT=FALSE){

  if(is.null(incidence.function) & is.null(infection.function)){
    stop("Pass either an incidence or infection function.")
  } else if(!is.null(incidence.function) & !is.null(infection.function)){
    stop("Pass either an incidence or infection function, but not both.")
  }

  # Get dimensions
  dims <- dim(sim_data$n_p)
  n_times <- dims[1]
  n_sims <- dims[2]

  # Convert the matrices to vectors
  # The data are ordered by simulation first, then times
  times <- as.vector(sim_data$times)
  n_p <- as.vector(sim_data$n_p)
  n_n <- as.vector(sim_data$n_n)

  # infection time
  if(!is.null(infection.function)){
    # get the uniform for the cumulative distribution
    # function for each individual
    cdfs <- lapply(n_p, runif)
    t_infect <- mapply(infection.function, e=cdfs, t=times, p=prevalence,
                       lambda_0=baseline_incidence, rho=rho, SIMPLIFY=F)
  } else {
    t_infect <- sim.infection.times(p=prevalence,
                                    inc.function=incidence.function,
                                    ...)
  }

  # infection duration to pass to the phi hat function
  infect_duration <- mapply(function(t, u) t - u, t=times, u=t_infect, SIMPLIFY=F)

  # if there is a recent test positive function provided,
  # get the indicator of recent infection
  if(!is.null(phi.func)){
    # probability of recent infection
    recent_probabilities <- lapply(infect_duration, phi.func)
    # indicators
    indicators <- lapply(recent_probabilities,
                         function(x) rbinom(n=length(x), size=1, prob=x))
  }

  if(!is.null(ptest.dist)){

    # Get the prior testing times for everyone
    test_times <- mapply(function(t, n, u) t - ptest.dist(n, u),
                         t=times, n=n_p, u=infect_duration)

    # See whether or not the test was positive
    ptest_delta <- mapply(function(it, pt) as.integer(it < pt), pt=test_times, it=t_infect)

    # Apply noise to the testing times if desired, but not changing ptest_delta
    if(!is.null(t_noise)){
      test_times <- lapply(test_times, function(t) sapply(t, t_noise))
    }

    # Simulate whether or not those tests are available (were actually taken)
    # Since can pass either a probability or a function,
    # convert the probabilities to a constant function of infection duration.
    if(!is.function(ptest.prob)){
      testprob <- function(u) ptest.prob
    } else {
      testprob <- ptest.prob
    }
    available <- mapply(function(n, u) rbinom(n=n, size=1, prob=testprob(u)),
                        n=n_p, u=infect_duration)

    # Apply a potential second testing mechanism to the first
    if(!is.null(ptest.dist2)){
      # Generate from second testing mechanism
      test_times2 <- mapply(function(t, n, u) t - ptest.dist2(n, u),
                            t=times, n=n_p, u=infect_duration)
      compare <- function(t, t1, t2, a1){
        t <- ifelse((((t1 > t2) | (t2 > 0)) & a1 == 1), t1, t2)
        return(t)
      }
      # Replace the test times and availability by new criteria of availability
      test_times <- mapply(compare,
                           t=times,
                           t1=test_times, t2=test_times2,
                           a1=available)
      available <- lapply(test_times, function(t) as.integer(t <= 0))
    }

    # Apply misreporting of available tests
    if(!is.null(q_misrep)){
      w <- lapply(n_p, function(n) rbinom(n, size=1, prob=1-q_misrep))
      available <- mapply(function(a, w) a * w, a=available, w=w)
    }

    # Get rid of prior testing times that are outside of the range
    if(!is.null(t_range)){
      rem <- function(x){
        x[(x > -t_range[1]) | (x < -t_range[2])] <- NA
        return(x)
      }
      test_times <- lapply(test_times, rem)
    }

    # Generate vector with prior time or NA if not available
    ptest_times <- mapply(function(t, a) ifelse(a, t, NA), t=test_times, a=available)
    available <- mapply(function(t, a) ifelse(is.na(t), 0, a), t=test_times, a=available)

    # Apply misreporting of positive tests
    if(!is.null(d_misrep)){
      z <- lapply(n_p, function(n) rbinom(n, size=1, prob=1-d_misrep))
      ptest_delta <- mapply(function(d, z) d * z, d=ptest_delta, z=z)
    }
    if(!is.null(p_misrep)){

      vs <- lapply(n_p, function(n) rbinom(n, size=1, prob=1-p_misrep))
      replace.func <- function(v, delta){
        v2 <- v
        v2[delta == 0] <- 1
        return(v2)
      }
      vs <- mapply(replace.func, v=vs, delta=ptest_delta)
      available <- mapply(function(a, v) a * v, a=available, v=vs)
    }

    # Apply exclusion for tests that were taken longer ago than bigT
    if(exclude_pt_bigT){
      exclude <- lapply(ptest_times, function(x) x < -bigT)
      available <- mapply(function(a, e) ifelse((!e & !is.na(e)), a, 0), a=available, e=exclude)
    }

    # DO THIS AGAIN
    # Generate vector with prior time or NA if not available
    ptest_times <- mapply(function(t, a) ifelse(a, t, NA), t=ptest_times, a=available)

    # Modify ptest_delta to have NA's when not available
    ptest_delta <- mapply(function(it, pt, d) ifelse(is.na(pt), NA, d),
                          pt=ptest_times,
                          it=t_infect,
                          d=ptest_delta)

    Ai <- list()
    Bi <- list()
    Mi <- list()
    BiTi <- list()

    for(k in 1:n_sims){
      Aik <- rep(0, n_p[k])
      Bik <- rep(0, n_p[k])
      BiTik <- rep(0, n_p[k])

      Cik <- -ptest_times[[k]] <= bigT
      Qik <- as.logical(available[[k]])

      Dik <- ptest_delta[[k]]

      Aik[Cik & Qik] <- 1
      Bik[!Cik & Qik] <- 1

      # This is the old recency indicator
      Rik <- indicators[[k]]

      # This is the new recency indicator
      Mik <- Rik * (1 - Bik * Dik) + (1 - Rik) * Aik * (1 - Dik)
      Mik[!Qik] <- Rik[!Qik]

      BiTik[Qik] <- -ptest_times[[k]][Qik] * Bik[Qik]

      Ai[[k]] <- Aik
      Bi[[k]] <- Bik
      Mi[[k]] <- Mik
      BiTi[[k]] <- BiTik
    }
  }

  if(summarize == TRUE){
    # number of recents
    n_r <- sapply(indicators, sum)
    n_d <- sapply(available, sum)
    q_eff <- n_d / n_p

    # Additional info from prior test results
    if(!is.null(ptest.dist)){
      # This is N^*_{rec}
      n_r_pt <- sapply(Mi, sum)
      # This is \sum 1 - I(T_i > T^*)
      num_beta <- sapply(Bi, function(x) sum(!x))
      # This is I(T_i > T^*) * T_i
      den_beta <- sapply(BiTi, function(x) sum(x))

    } else {
      n_r_pt <- NULL
      num_beta <- NULL
      den_beta <- NULL
      q_eff <- NULL
    }

    if(n_times > 1){
      # if we have multiple times, convert it back into a matrix
      n_r       <- matrix(n_r, nrow=n_times, ncol=n_sims, byrow=FALSE)
      n_r_pt    <- matrix(n_r_pt, nrow=n_times, ncol=n_sims, byrow=FALSE)
      num_beta  <- matrix(num_beta, nrow=n_times, ncol=n_sims, byrow=FALSE)
      den_beta  <- matrix(den_beta, nrow=n_times, ncol=n_sims, byrow=FALSE)
      q_eff     <- matrix(q_eff, nrow=n_times, ncol=n_sims, byrow=FALSE)
      n_p       <- sim_data$n_p
      n_n       <- sim_data$n_n
      times     <- sim_data$times
      n         <- sim_data$n
    } else {
      # if we have only one time (cross-sectional), keep everything
      # in vector form
      n_p   <- as.vector(sim_data$n_p)
      n_n   <- as.vector(sim_data$n_n)
      n     <- as.vector(sim_data$n)
      times <- as.vector(sim_data$times)
    }

    aspect_list <- list(
      n=n,
      n_p=n_p,
      n_n=n_n,
      n_r=n_r,
      times=times,
      n_r_pt=n_r_pt,
      num_beta=num_beta,
      den_beta=den_beta,
      q_eff=q_eff,
      ptest_times=ptest_times,
      ptest_delta=ptest_delta,
      ptest_avail=available,
      ri=Mi
    )
    return(aspect_list)
  } else {

    # Get the simulation number
    sim <- rep(1:n_sims, each=n_times)

    # Get the enrollment times for negative and positive
    time_neg <- rep(times, n_n)
    time_pos <- rep(times, n_p)
    time <- c(time_neg, time_pos)

    # Get simulation number for each negative and positive
    sim_neg <- rep(sim, n_n)
    sim_pos <- rep(sim, n_p)
    sim <- c(sim_neg, sim_pos)

    # Get prevalence positive indicator
    prev_pos <- rep(1, length(sim_pos))
    prev_neg <- rep(0, length(sim_neg))
    prev <- c(prev_neg, prev_pos)

    # Infection time and prior testing time / result
    na_neg <- rep(NA, length(time_neg))
    infect <- c(na_neg, unlist(t_infect))

    if(!is.null(ptest.dist)){
      priorT <- c(na_neg, unlist(ptest_times))
      priorD <- c(na_neg, as.integer(unlist(ptest_delta)))
    } else {
      priorT <- NULL
      priorD <- NULL
    }
    df <- data.table(
      sim=sim,
      time=time,
      pos=prev,
      itime=infect,
      priorT=priorT,
      priorD=priorD
    )

    if(!is.null(phi.func)){

      # Phi function evaluated at the infection times *-1
      probs_pos <- unlist(recent_probabilities)
      probs_neg <- rep(NA, length(time_neg))
      probs <- c(probs_neg, probs_pos)

      # Whether or not they tested positive
      r_pos <- unlist(indicators)
      r_neg <- rep(NA, length(time_neg))
      r <- c(r_neg, r_pos)

      df$probs <- probs
      df$rpos <- r
    }
    setorder(df, sim, time, pos, itime)
    return(df)
  }
}

#' Create simulations of trials based on a past infection time function
#' and a prevalence and baseline incidence (time 0).
#'
#' @export
#' @param n_sims Number of simulations
#' @param n Number of subjects screened
#' @param infection.function Function that simulates the infection time
#' @param phi.func Optional test positive function for recency assay.
#'   If you pass \code{summarize=TRUE}, then you must provide a test positive function.
#' @param baseline_incidence Baseline incidence value (time 0)
#' @param prevalence Constant prevalence value
#' @param rho A parameter to the \code{infection.function} for how quickly incidence changes
#' @param times A vector of times at which to enroll (e.g. c(0, 1, 2, 3))
#' @param summarize Whether to return a summary of the dataset or unit-record data
#' @return If \code{summarize=TRUE}, a list of vectors (or matrices if
#'   length(times > 1) for summary screening quantities
#'   across simulations (and potentially times if length(times) > 1):
#'   \itemize{
#'     \item n: number screened
#'     \item times: times screened
#'     \item n_p: number of positives
#'     \item n_n: number of negatives
#'     \item n_r: number of recents
#'   }
#'   If \code{summarize=FALSE}, then returns a list of data frames by simulation
#'   that have unit-record data with infection times with column names:
#'   \itemize{
#'     \item time: time screened
#'     \item pos: indicator for positive
#'     \item itime: infection time (reference time 0)
#'     \item probs: recency test positive probability based on infection time
#'     \item rpos: indicator for recency test positive, if a phi function is passed
#'   }
#' @examples
#' set.seed(1)
#' generate.data(n_sims=3, n=2, infection.function=infections.con,
#'               baseline_incidence=0.05, prevalence=0.5, rho=NA,
#'               phi.func=function(t) 1 - pgamma(t, 1, 2),
#'               times=c(0, 1), summarize=FALSE)
#' generate.data(n_sims=3, n=100, infection.function=infections.con,
#'               baseline_incidence=0.05, prevalence=0.3, rho=NA,
#'               phi.func=function(t) 1 - pgamma(t, 1, 2),
#'               times=c(0, 1), summarize=TRUE)
#' generate.data(n_sims=3, n=100, infection.function=infections.con,
#'               baseline_incidence=0.05, prevalence=0.3, rho=NA,
#'               phi.func=function(t) 1 - pgamma(t, 1, 2),
#'               times=c(0), summarize=TRUE)
generate.data <- function(n_sims, n, infection.function,
                          baseline_incidence,
                          prevalence, rho, phi.func=NULL, times=c(0), summarize=TRUE, ...){

  # TODO: include some sanity checks
  if(summarize & is.null(phi.func)) stop("Need a phi function
                                       to summarize recents.")

  data <- generate.raw.data(n_sims, n, prevalence, times=times)
  data <- simulate.recent(data, infection.function, phi.func,
                          baseline_incidence, prevalence, rho, summarize=summarize, ...)
  return(data)
}

# INFECTION FUNCTIONS

#' Constant incidence function
#'
#' @param t Time, float or vector
#' @param lambda_0 Incidence constant value
#' @return Incidence at time t
#' @examples
#' incidence.con(0, lambda_0=0.05)
#' incidence.con(c(-1, 0, 0), lambda_0=0.05)
incidence.con <- Vectorize(function(t, lambda_0, rho=NA) lambda_0)

#' Linearly decreasing incidence function
#'
#' \deqn{
#'   \lambda(t) = \lambda_0 - \rho t
#' }
#'
#' @param t Time, float or vector
#' @param lambda_0 Incidence at time 0
#' @param rho Linear decrease in incidence
#' @return Incidence at time t
#' @examples
#' incidence.lin(0, lambda_0=0.05, rho=1e-3)
#' incidence.lin(c(-1, 0, 1), lambda_0=0.05, rho=1e-3)
incidence.lin <- function(t, lambda_0, rho=1) lambda_0 - rho * t

#' Exponentially decreasing incidence function
#'
#' \deqn{
#'   \lambda(t) = \lambda_0 \exp(-\rho t)
#' }
#'
#' @param t Time, float or vector
#' @param lambda_0 Incidence at time 0
#' @param rho Exponential decrease in incidence
#' @return Incidence at time t
#' @examples
#' incidence.exp(0, lambda_0=0.05, rho=0.07)
#' incidence.exp(c(-1, 0, 1), lambda_0=0.05, rho=0.07)
incidence.exp <- function(t, lambda_0, rho=1) lambda_0 * exp(-rho * t)

#' Infection times function based on constant incidence.
#'
#' @param e A number between 0 and 1 (randomly generated), float or vector
#' @param t Time, float
#' @param p Constant prevalence
#' @param lambda_0 Baseline incidence
#' @return Infection time
#' @examples
#' e <- runif(10)
#' infections.con(e, t=0, p=0.2, lambda_0=0.05)
#' @export
infections.con <- function(e, t, p, lambda_0, rho=NA){
  # Note that you can work with e or 1 - e since e is Uniform(0, 1)
  infections <- t - p*e / ((1 - p) * lambda_0)
  return(infections)
}

#' Infection times function based on linearly decreasing incidence.
#'
#' @param e A number between 0 and 1 (randomly generated), float or vector
#' @param t Time, float
#' @param p Constant prevalence
#' @param lambda_0 Baseline incidence
#' @param rho Linearly decreasing incidence parameter
#' @return Infection time
#' @examples
#' e <- runif(10)
#' infections.lin(e, t=0, p=0.2, lambda_0=0.05, rho=1e-3)
#' @export
infections.lin <- function(e, t, p, lambda_0, rho){
  incidence <- lambda_0
  numerator <- incidence**2 + 2 * rho * p * e / (1 - p)
  numerator <- sqrt(numerator) - incidence

  infections <- t - numerator / rho
  return(infections)
}

#' Infection times function based on constant incidence until bigT
#' then linearly increasing after. Compared to the other infection functions
#' we've switched to 1-e, which doesn't matter because e is uniform.
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

#' Infection times function based on exponential decreasing incidence.
#'
#' @param e A number between 0 and 1 (randomly generated), float or vector
#' @param t Time, float
#' @param p Constant prevalence
#' @param lambda_0 Baseline incidence
#' @param rho Linearly decreasing incidence parameter
#' @return Infection time
#' @examples
#' e <- runif(10)
#' infections.exp(e, t=0, p=0.2, lambda_0=0.05, rho=0.07)
#'
#' @export
infections.exp <- function(e, t, p, lambda_0, rho){
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
#' @param inc.function A function of time that returns incidence
#' @param nsims Number of infection times to simulate
#' @param dt The precision for numerical integration. May need to increase.
#' @return A vector of infection durations
#' @examples
#' set.seed(1)
#' sim.infection.times(p=0.2, nsims=15,
#'                     inc.function=function(t) incidence.con(t, 0.5))
#' sim.infection.times(p=0.2, nsims=15,
#'                     inc.function=function(t) incidence.con(t, 0.05))
sim.infection.times <- function(p, inc.function, nsims=1000, dt=0.001){
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
