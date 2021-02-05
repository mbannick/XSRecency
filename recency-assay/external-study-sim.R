# -------------------------------------------------------------------------
# EXTERNAL STUDY DATA SIMULATION TO ESTIMATE MU, OMEGA, BETA
# -------------------------------------------------------------------------
library(tidyr)

#' This simulates the false recency rate using
#' tail probabilities from the phi parameters
#' by sampling from a binomial with a sample size
#' given in the study.
simulate.beta <- function(phi.func, big_T, max_T){
  infected_times <- runif(n=PHI_PARAMS$N_LONG_INFECT,
                          min=big_T,
                          max=max_T)
  recent <- rbinom(n=PHI_PARAMS$N_LONG_INFECT, size=1,
                   p=phi.func(infected_times))
  beta <- sum(recent) / PHI_PARAMS$N_LONG_INFECT
  beta_var <- beta * (1 - beta) / PHI_PARAMS$N_LONG_INFECT
  return(list(est=beta, var=beta_var))
}

#' Simulate infection duration and recency indicators
#' from the study.
#'
simulate.study <- function(phi.func){

  # probability of sampling in each cohort
  p_samp <- (PHI_PARAMS$mean_samp - 1) / (PHI_PARAMS$max_samp - 1)

  # simulate the number of samples per individual in each cohort
  coh_samples <- mapply(rbinom, n=PHI_PARAMS$coh_n,
                        size=PHI_PARAMS$mean_samp,
                        p=p_samp, SIMPLIFY=FALSE)

  # add one to each of the samples in the list
  coh_samples <- lapply(coh_samples, function(x) x + 1)

  # this loops through each cohort to sample the infection durations
  infection_durations <- list()
  for(i in 1:length(PHI_PARAMS$coh)){

    # this is where we'll put all of the infection durations
    # for this cohort
    durations <- c()

    # this gets an array for each person by looping over each person
    # as a num. of samples, and pulling that many observations of a
    # multinomial with the correct bucket probabilities
    people <- lapply(coh_samples[[i]], rmultinom, size=1,
                     prob=PHI_PARAMS$coh_buckets[[i]])

    # this loops through each person in the cohort
    for(person in people){

      # this loops through each person's sample
      # there are ncol(person) number of samples
      for(sample in 1:ncol(person)){

        # this gets which bin the infection duration is in for person/sample
        # and maps that to the duration (lower, upper)
        duration_bin <- PHI_PARAMS$duration_years[[which(person[, sample] == 1)]]

        # this samples the duration from the uniform distribution
        # with lower and upper specified by the duration bin
        duration <- runif(n=1, min=duration_bin[1], max=duration_bin[2])

        # append the current duration to the vector of infection
        # durations, since the person and sample index don't matter
        # at this point
        durations <- append(durations, duration)
      }
    }
    infection_durations[[i]] <- durations
  }

  # generate the recency indicator for the infection duration
  # based on the true phi() function
  indicators <- lapply(infection_durations, function(x) rbinom(n=length(x),
                                                               size=1, prob=phi.func(x)))
  return(list(
    recent=unlist(indicators),
    durations=unlist(infection_durations)
  ))
}

#' Fits a cubic model to indicators and durations
#' and returns the model object.
#'
fit.cubic <- function(recent, durations){
  durations2 <- durations**2
  durations3 <- durations**3

  # fit a logistic regression with cubic polynomials
  model <- glm(recent ~ durations + durations2 + durations3,
               family=binomial(link="logit"))
  return(model)
}

#' Create phi probability predictions
#' based on a model object and d durations
#'
phi.hat <- function(d, model){
  dmat <- cbind(
    rep(1, length(d)),
    d, d**2, d**3
  )
  lin <- dmat %*% matrix(coef(model))
  return(exp(lin) / (1 + exp(lin)))
}

#' Estimate omega_T from a model object
#' Can instead estimate mu hat if you pass in a sufficiently
#' long follow_T parameter (like 9)
integrate.phi <- function(model, max_T){

  # generate a sequence of duration times, starting at 0 and
  # up to max t, just to get the phi_hat
  dt <- 0.001

  ts <- seq(0, max_T, by=dt)
  time <- ts
  ts <- cbind(rep(1, length(ts)),
              ts,
              ts**2,
              ts**3)

  # get the linear predictor
  lin_predictor <- ts %*% matrix(coef(model))

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

  vc <- vcov(model)

  # the variance of omega_t
  variance <- c(t(grad) %*% vc %*% grad)

  return(list(est=estimate, var=variance, ts=time, phi=phi))
}

#' Function that simulates all the properties of
#' the recency assay, and returns all of the ones that we're interested in
#'
#' @export
assay.properties.sim <- function(phi.func, max_T, big_T){
  # SIMULATIONS -----------------------------

  # simulate one phi function and
  # get the estimate and variance
  study <- simulate.study(phi.func)

  model <- fit.cubic(recent=study$recent,
                     durations=study$durations)

  # get mu and omega
  mu_sim <- integrate.phi(model, max_T=max_T)
  omega_sim <- integrate.phi(model, max_T=big_T)

  # simulate beta and get the estimate and variance of them
  beta_sim <- simulate.beta(big_T=big_T, max_T=max_T, phi.func=phi.func)

  result <- list(estimators=list(
                  mu_est=mu_sim$est,
                  mu_var=mu_sim$var,
                  omega_est=omega_sim$est,
                  omega_var=omega_sim$var,
                  beta_est=beta_sim$est,
                  beta_var=beta_sim$var
                 ),
                 curve=data.table(
                   ts=mu_sim$ts,
                   phi=mu_sim$phi
                 )
                )
  return(result)
}

#' Function that simulates all the properties of
#' the recency assay, and returns all of the ones that we're interested in
#'
assay.properties.nsim <- function(n_sims, phi.func, max_T, big_T){
  replicates <- replicate(n_sims, assay.properties.sim(phi.func=phi.func,
                                                       max_T=max_T, big_T=big_T))
  result <- do.call(cbind, replicates["estimators", ]) %>% t %>% data.table
  columns <- colnames(result)
  cols <- lapply(columns, function(x) unlist(result[[x]]))

  df <- do.call(cbind, cols) %>% data.table
  names(df) <- columns
  df[, sim := .I]

  est <- melt(df, id.vars="sim", measure.vars=patterns("_est"), value.name="estimate")
  vars <- melt(df, id.vars="sim", measure.vars=patterns("_var"), value.name="variance")

  est <- est %>% separate(col="variable", into=c("quantity", "stub"), sep="_")
  vars <- vars %>% separate(col="variable", into=c("quantity", "stub"), sep="_")

  est <- est[, .(sim, quantity, estimate)]
  vars <- vars[, .(sim, quantity, variance)]

  df <- merge(est, vars, by=c("sim", "quantity"))
  df[, lower := estimate - qnorm(0.975) * sqrt(variance)]
  df[, upper := estimate + qnorm(0.975) * sqrt(variance)]

  true_mu <- true_mdri(phi.func, follow_T=max_T)
  true_omega <- true_mdri(phi.func, follow_T=big_T)

  df[, truth := NULL]
  df[quantity == "mu", truth := true_mu]
  df[quantity == "omega", truth := true_omega]

  phis <- do.call(cbind, replicates["curve", ])
  n <- ncol(phis)
  phis <- phis[, c(1, seq(2, n, 2)), with=F]
  n <- ncol(phis)
  namecols <- paste0("phi_", 1:(n-1))
  setnames(phis, c("time", namecols))
  phis <- melt(phis, id.vars=c("time"))
  phis <- phis %>% separate(col='variable', into=c("stub", "sim"), sep="_")

  return(list(estimators=df, phi=phis))
}

true_mdri <- function(phi.function, follow_T){
  ts <- seq(0, follow_T, 1e-3)
  return(sum(phi.function(ts) * 1e-3))
}
