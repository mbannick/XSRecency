#' Snapshot estimator
#'
#' @export
snapshot.estimate <- function(n_r, n_n, mu){
  val <- n_r / (mu * n_n)
  return(val)
}

#' Adjusted estimator from Kassanjee
#'
#' @export
adjusted.estimate <- function(n_r, n_n, n_p, omega, beta, big_T){
  o_b <- omega - beta * big_T
  est <- n_r - beta * n_p
  val <- est / (o_b * n_n)
  return(val)
}

#' Adjusted variance computation, is the same for adjusted
#' and snapshot estimator, just pass in beta_sim=list(est=0, var=0).
#'
#' @export
variance <- function(n_n, n_r, n_p, n, omega, omega_var, beta, beta_var, big_T){

  variance <- (
    n_r * (n_p - n_r) / (n_p * (n_r - n_p * beta)**2) +
      n / (n_p * n_n) +
      beta_var * n_p * (n - n_p) / (n * (n_r - n_p * beta)**2) +
      omega_var / (omega - beta * big_T)**2 +
      beta_var * (
        (n_p * omega - n_r * big_T) /
          ((n_r - n_p * beta) * (omega - beta * big_T))
      ) ** 2
  )
  return(variance)
}

#' Convert log variance to variance
#'
#' @export
var.log.to.var <- function(estimate, variance) (estimate ** 2) * variance

#' Get snapshot estimate and variance.
#'
#' @export
get.snapshot <- function(n_r, n_n, n_p, n, mu, mu_var){
  est <- snapshot.estimate(n_r=n_r, n_n=n_n, mu=mu)
  logvar <- variance(n_n=n_n, n_r=n_r, n_p=n_p, n=n,
                     omega=mu, omega_var=mu_var,
                     beta=0, beta_var=0, big_T=0)
  estvar <- var.log.to.var(est, logvar)
  return(list(est=est, var=estvar))
}

#' Get adjusted estimate and variance.
#'
#' @export
get.adjusted <- function(n_r, n_n, n_p, n, omega, omega_var,
                         beta, beta_var, big_T){
  est <- adjusted.estimate(n_r=n_r, n_n=n_n, n_p=n_p,
                           omega=omega, beta=beta, big_T=big_T)
  logvar <- variance(n_n=n_n, n_r=n_r, n_p=n_p, n=n,
                     omega=omega, omega_var=omega_var,
                     beta=beta, beta_var=beta_var, big_T=big_T)
  estvar <- var.log.to.var(est, logvar)
  return(list(est=est, var=estvar))
}
