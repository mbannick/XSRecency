# -------------------------------------------------------------------------
# EXTERNAL STUDY DATA SIMULATION TO ESTIMATE MU, OMEGA, BETA
# -------------------------------------------------------------------------

source("./external-study-data.R")

#' This simulates the false recency rate using
#' tail probabilities from the phi parameters
#' by sampling from a binomial with a sample size
#' given in the study.
#' 
simulate.beta <- function(){
  false <- rbinom(n=1, size=PHI_PARAMS$N_LONG_INFECT,
                  prob=PHI_PARAMS$TAIL_PROB)
  beta <- false / PHI_PARAMS$N_LONG_INFECT
  beta_var <- beta * (1 - beta) / PHI_PARAMS$N_LONG_INFECT
  return(list(est=beta, var=beta_var))
}

#' Simulate infection duration and recency indicators
#' from the study.
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
integrate.phi <- function(model, follow_T=PHI_PARAMS$FOLLOW_T){
  
  # generate a sequence of duration times, starting at 0 and 
  # up to max t, just to get the phi_hat
  dt <- 0.001
  
  ts <- seq(0, follow_T, by=dt)
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
  
  return(list(est=estimate, var=variance))
}

#' Function that simulates all the properties of
#' the recency assay, and returns all of the ones that we're interested in
assay.properties.sim <- function(phi.func){
  
  # SIMULATIONS -----------------------------
  
  # simulate one phi function and
  # get the estimate and variance
  study <- simulate.study(phi.func)
  
  model <- fit.cubic(recent=study$recent,
                     durations=study$durations)
  
  # get mu and omega
  mu_sim <- integrate.phi(model, follow_T=9.9)
  omega_sim <- integrate.phi(model, follow_T=PHI_PARAMS$FOLLOW_T)
  
  # simulate beta and get the estimate and variance of them
  beta_sim <- simulate.beta()
  
  result <- list(mu=mu_sim,
                 omega=omega_sim,
                 beta=beta_sim)
  return(result)
}
