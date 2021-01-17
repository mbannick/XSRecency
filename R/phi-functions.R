# CHARACTERISTIC FUNCTIONS FOR PHI

#' Characteristic function for phi, with a given
#' MDRI and FRR. If FRR = 0, then MDRI is the "window" period.
#'
#' @export
phi.character.1 <- function(t, mdri=142, frr=0.015){
  alpha <- 1
  x <- (mdri/365.25-frr)/(1-frr)
  beta <- alpha/x
  phit <- (1 - pgamma(t, shape = alpha, rate = beta))*(1-frr)+frr
  return(phit)
}

#' Characteristic function for phi, with
#' increasing phi(t) in the tail
#' such that long-infected individuals have an increasing
#' probability of testing recent
#'
#' @export
phi.character.2 <- function(t, mdri=142, frr=0.015){
  alpha <- 1
  x <- (mdri/365.25-frr)/(1-frr)
  beta <- alpha/x
  phit <- (1-pgamma(t, shape=alpha, rate=beta))*(1 - frr) + frr + dnorm(t-7, mean=0, sd=1) / 8
  return(phit)
}

#' Non-constant FRR true rate based on the phi function and the tail probability of 0.015
#' @export
true.frr <- function(phi.function, mdri=142, frr=0.015){
  ts <- seq(PHI_PARAMS$FRR_MIN, PHI_PARAMS$FRR_MAX, 1e-3)
  return(mean(phi.function(ts, mdri, frr)))
}

#' Calculate the shadow period for a given phi function.
#' @export
calculate.shadow <- function(phi.function, tau, mdri=142, frr=0.015){
  ts <- seq(0, tau, 1e-3)
  lomega <- sum(ts * 1e-3 * (phi.function(ts, mdri, frr) - frr) / (mdri/365.25 - frr))
  return(lomega)
}
