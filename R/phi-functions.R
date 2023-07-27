# CHARACTERISTIC FUNCTIONS FOR PHI

DELTA <- 1e-3

# Gets Gamma alpha and beta parameters
# for a given mean window period and shadow
# to approximate a phi function.
#
# @param window Window period
# @param shadow Shadow period
get.gamma.params <- function(window, shadow){
  alpha <- window / (2 * shadow - window)
  beta <- 1 / (2 * shadow - window)
  return(c(alpha, beta))
}

OMEGA <- 0.33
BIGT <- 2
SHADOW <- 0.2
FRR <- 0.02

loss <- function(X){
  a <- X[1]
  b <- X[2]

  phi <- function(t) 1-pgamma(t, shape=a, rate=b) + FRR
  phit <- function(t) t * (1-pgamma(t, shape=a, rate=b) + FRR)

  intphi <- integrate(phi, lower=0, upper=BIGT)$value
  intphit <- integrate(phit, lower=0, upper=BIGT)$value

  l1 <- OMEGA - intphi
  l2 <- SHADOW - (intphit - FRR * BIGT**2/2)/(OMEGA - BIGT * FRR)

  return(sum(c(l1, l2)**2))
}

vals <- optim(c(1, 0.01), fn=loss, lower=0, upper=Inf, method="L-BFGS-B")

phi <- function(t){
  value <- (t <= BIGT)* (1 - pgamma(t, shape=vals$par[1], rate=vals$par[2]) + FRR) +
  (t > BIGT)*FRR
  return(value)
}

ts <- seq(0, 4, by=0.02)
plot(phi(ts) ~ ts, type='l')

#' Approximate test-recent function with given window and shadow period
#'
#' Generate a test-recent function that has a desired window and shadow period
#' using a gamma distribution parameters.
#'
#' @param window Window period (in days)
#' @param shadow Shadow period (in days)
#'
#' @export
#' @returns A test-recent function of t
#'
#' @examples
#'
#' phi <- getTestRecentFunc(window=200, shadow=191)
getTestRecentFunc <- function(window, shadow){
  params <- get.gamma.params(window=window/365.25, shadow=shadow/365.25)
  phi.func <- function(t) 1-pgamma(t, shape = params[1], rate = params[2])

  return(phi.func)
}

# Non-constant FRR true rate based on the phi function and the tail probability of 0.015
true.frr <- function(phi.function, bigT, tau){
  # ts <- seq(bigT, tau, DELTA)
  # return(mean(phi.function(ts)))
  val <- integrate(phi.function, bigT, tau)$value/(tau-bigT)
  return(val)
}

#' Integrate a phi function to calculate the MDRI or window period
#'
#' @param phi.function A function of t (years) which has range 0 to at least \code{maxT}
#' @param maxT The t to integrate to
#'
#' @export
integratePhi <- function(phi.function, maxT){
  # ts <- seq(0, maxT, DELTA)
  # return(sum(phi.function(ts) * DELTA))
  val <- integrate(phi.function, 0, maxT)$value
  return(val)
}

true.shadow.snap <- function(phi.function, tau){
  ts <- seq(0, tau, DELTA)
  window <- true.window.mdri(phi.function, maxT=tau)
  print(window)
  return(sum(phi.function(ts) * ts / window * DELTA))
}

true.shadow.adj <- function(phi.function, bigT, tau, rho){
  ts <- seq(0, bigT, DELTA)
  omega <- true.window.mdri(phi.function)
  beta <- true.frr(phi.function, bigT=bigT, tau=tau)
  return(sum(ts * (phi.function(ts) - beta)/(omega - beta) * DELTA))
}

# Expected bias for the snapshot estimator
snap.bias <- function(phi.func, inc.func, inc.0, tau, rho){
  window <- true.window.mdri(phi.func, maxT=tau)
  us <- seq(0, tau, DELTA)

  integrand <- (phi.func(us) / window) * inc.function(t=-us,lambda_0=inc.0, rho=rho)
  integral <- sum(integrand * DELTA)
  integral <- integral - inc.0
  return(integral)
}

# Expected bias for the snapshot estimator
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
