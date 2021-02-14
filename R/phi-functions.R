# CHARACTERISTIC FUNCTIONS FOR PHI

DELTA <- 1e-3

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
  ts <- seq(PHI_PARAMS$FRR_MIN, PHI_PARAMS$FRR_MAX, DELTA)
  return(mean(phi.function(ts, window, frr, shadow)))
}

#' When using the shadow period, need to numerically calculate the MDRI
#' the MDRI argument here is really mean window period, just calling it
#' MDRI for convenience
#' @export
true.mdri <- function(phi.function, window=142, frr=0.015, shadow=150){
  ts <- seq(0, PHI_PARAMS$FRR_MIN, DELTA)
  return(sum(phi.function(ts, window, frr, shadow) * DELTA))
}

#' @export
#' Numerically calculate the true window period.
true.window <- function(phi.function, window=142, frr=0.015, shadow=150){
  ts <- seq(0, PHI_PARAMS$FRR_MAX, DELTA)
  return(sum(phi.function(ts, window, frr, shadow) * DELTA))
}
#' @export
true.shadow.snap <- function(phi.function, ...){
  ts <- seq(0, PHI_PARAMS$FRR_MAX, DELTA)
  window <- true.window(phi.function, ...)
  return(sum(phi.function(ts, ...) * ts / window * DELTA))
}
#' @export
true.shadow.adj <- function(phi.function, rho, ...){
  ts <- seq(0, PHI_PARAMS$FRR_MIN, DELTA)
  omega <- true.mdri(phi.function, ...)
  beta <- true.frr(phi.function, ...)
  return(sum(ts * (phi.function(ts, ...) - beta)/(omega - beta) * DELTA))
}
#' @export
snap.bias.linear <- function(phi.func, rho, ...){
  shadow <- true.shadow.snap(phi.func, ...)
  return(shadow * rho)
}
#' @export
adj.bias.linear <- function(phi.func, rho, ...){
  shadow <- true.shadow.adj(phi.func, ...)
  return(shadow * rho)
}

#' Expected bias for the snapshot estimator
#' @export
snap.bias <- function(phi.func, inc.func, inc.0, rho, ...){
  window <- true.window(phi.func, ...)
  us <- seq(0, PHI_PARAMS$FRR_MAX, DELTA)

  integrand <- (phi.func(us, ...) / window) * inc.function(t=-us,
                                                           lambda_0=inc.0,
                                                           rho=rho)
  integral <- sum(integrand * DELTA)
  integral <- integral - inc.0
  return(integral)
}

# Expected bias for the snapshot estimator
#' @export
adj.bias <- function(phi.func, inc.func, inc.0, rho, ...){
  mdri <- true.mdri(phi.func, ...)
  frr <- true.frr(phi.func, ...)
  us <- seq(0, PHI_PARAMS$FRR_MIN, DELTA)

  integrand <- (phi.func(us, ...) - frr) / (mdri - PHI_PARAMS$FRR_MIN * frr) * inc.function(t=-us,
                                                                       lambda_0=inc.0,
                                                                       rho=rho)
  integral <- sum(integrand * DELTA)
  integral <- integral - inc.0
  return(integral)
}
