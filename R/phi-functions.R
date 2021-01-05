# CHARACTERISTIC FUNCTIONS FOR PHI

#' Characteristic function for phi, with a given
#' MDRI and FRR. If FRR = 0, then MDRI is the "window" period.
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
phi.character.2 <- function(t, mdri=142, frr=0.015){
  alpha <- 1
  x <- (mdri/365.25-frr)/(1-frr)
  beta <- alpha/x
  phit <- (1 - pgamma(t, shape = alpha, rate = beta))*(1-frr)+frr + 
    0.01*(t-2)*(t>2)
  return(phit)
}
