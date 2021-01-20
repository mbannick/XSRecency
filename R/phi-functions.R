# CHARACTERISTIC FUNCTIONS FOR PHI

#' Gets the Gamma alpha and beta parameters
#' for a given mean window period and shadow
#' @export
get.gamma.params <- function(window, shadow){
  alpha <- window ** 2 / (2 * shadow - window ** 2)
  beta <- window / (2 * shadow - window ** 2)
  return(c(alpha, beta))
}

#' Characteristic function for phi, with a given
#' MDRI and FRR. If FRR = 0, then MDRI is the "window" period.
#'
#' @export
phi.character.1 <- function(t, window=142, frr=0.015, shadow=150){
  params <- get.gamma.params(window=window/356.25, shadow=shadow/365.25)
  phit <- (1 - pgamma(t, shape = params[1], rate = params[2]))*(1-frr)+frr
  return(phit)
}

#' Characteristic function for phi, with
#' increasing phi(t) in the tail
#' such that long-infected individuals have an increasing
#' probability of testing recent
#'
#' @export
phi.character.2 <- function(t, window=142, frr=0.015, shadow=150){
  params <- get.gamma.params(window=window/356.25, shadow=shadow/365.25)
  phit <- (1-pgamma(t, shape=params[1], rate=params[2]))*(1 - frr) + frr + dnorm(t-7, mean=0, sd=1) / 8
  return(phit)
}

#' Non-constant FRR true rate based on the phi function and the tail probability of 0.015
#' @export
true.frr <- function(phi.function, window=142, frr=0.015, shadow=150){
  ts <- seq(PHI_PARAMS$FRR_MIN, PHI_PARAMS$FRR_MAX, 1e-3)
  return(mean(phi.function(ts, window, frr, shadow)))
}

#' When using the shadow period, need to numerically calculate the MDRI
#' the MDRI argument here is really mean window period, just calling it
#' MDRI for convenience
#' @export
true.mdri <- function(phi.function, window=142, frr=0.015, shadow=150){
  if(frr == 0){
    # that means we're calculating the window rather than MDRI
    return(window / 365.25)
  } else {
    ts <- seq(0, PHI_PARAMS$FRR_MIN, 1e-3)
    return(sum(phi.function(ts, window, frr, shadow) * 1e-3))
  }
}
