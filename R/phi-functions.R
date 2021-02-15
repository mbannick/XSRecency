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

#' Non-constant FRR true rate based on the phi function and the tail probability of 0.015
#' @export
true.frr <- function(phi.function, bigT, tau){
  ts <- seq(bigT, tau, DELTA)
  return(mean(phi.function(ts)))
}

#' When using the shadow period, need to numerically calculate the MDRI
#' the MDRI argument here is really mean window period, just calling it
#' MDRI for convenience
#' @export
true.window.mdri <- function(phi.function, maxT){
  ts <- seq(0, maxT, DELTA)
  return(sum(phi.function(ts) * DELTA))
}

#' @export
true.shadow.snap <- function(phi.function, tau){
  ts <- seq(0, tau, DELTA)
  window <- true.window(phi.function)
  return(sum(phi.function(ts) * ts / window * DELTA))
}

#' @export
true.shadow.adj <- function(phi.function, bigT, tau, rho){
  ts <- seq(0, bigT, DELTA)
  omega <- true.mdri(phi.function)
  beta <- true.frr(phi.function, bigT=bigT, tau=tau)
  return(sum(ts * (phi.function(ts) - beta)/(omega - beta) * DELTA))
}

#' Expected bias for the snapshot estimator
#' @export
snap.bias <- function(phi.func, inc.func, inc.0, tau, rho){
  window <- true.window.mdri(phi.func, maxT=tau)
  us <- seq(0, tau, DELTA)

  integrand <- (phi.func(us) / window) * inc.function(t=-us,lambda_0=inc.0, rho=rho)
  integral <- sum(integrand * DELTA)
  integral <- integral - inc.0
  return(integral)
}

# Expected bias for the snapshot estimator
#' @export
adj.bias <- function(phi.func, inc.func, inc.0, bigT, tau, rho){
  mdri <- true.window.mdri(phi.func, maxT=bigT)
  frr <- true.frr(phi.func, bigT=bigT, tau=tau)
  us <- seq(0, bigT, DELTA)

  integrand <- (phi.func(us) - frr) / (mdri - bigT * frr) * inc.function(t=-us,
                                                                       lambda_0=inc.0,
                                                                       rho=rho)
  integral <- sum(integrand * DELTA)
  integral <- integral - inc.0
  return(integral)
}
