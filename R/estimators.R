#' Snapshot estimator
#'
#' @export
snapshot.estimate <- function(n_r, n_n, mu_sim){
  val <- n_r / (mu_sim$est * n_n)
  return(val)
}

#' Adjusted estimator from Kassanjee
#'
#' @export
adjusted.estimate <- function(n_r, n_n, n_p, omega_sim, beta_sim, big_T){
  o_b <- omega_sim$est - beta_sim$est * big_T
  est <- n_r - beta_sim$est * n_p
  val <- est / (o_b * n_n)
  return(val)
}

#' Adjusted variance computation, is the same for adjusted
#' and snapshot estimator, just pass in beta_sim=list(est=0, var=0).
#'
#' @export
variance <- function(n_n, n_r, n_p, n, omega_sim, beta_sim, big_T){

  beta_t <- beta_sim$est
  beta_var <- beta_sim$var
  omega_t <- omega_sim$est
  omega_var <- omega_sim$var

  variance <- (
    n_r * (n_p - n_r) / (n_p * (n_r - n_p * beta_t)**2) +
      n / (n_p * n_n) +
      beta_var * n_p * (n - n_p) / (n * (n_r - n_p * beta_t)**2) +
      omega_var / (omega_t - beta_t * big_T)**2 +
      beta_var * (
        (n_p * omega_t - n_r * big_T) /
          ((n_r - n_p * beta_t) * (omega_t - beta_t * big_T))
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
get.snapshot <- function(n_r, n_n, n_p, n, mu_sim){
  beta_null <- list(est=0, var=0)
  est <- snapshot.estimate(n_r=n_r, n_n=n_n, mu_sim=mu_sim)
  logvar <- variance(n_n=n_n, n_r=n_r, n_p=n_p, n=n,
                     omega_sim=mu_sim, beta_sim=beta_null, big_T=0)
  estvar <- var.log.to.var(est, logvar)
  return(list(est=est, var=estvar))
}

#' Get adjusted estimate and variance.
#'
#' @export
get.adjusted <- function(n_r, n_n, n_p, n, omega_sim, beta_sim, big_T){
  est <- adjusted.estimate(n_r=n_r, n_n=n_n, n_p=n_p,
                           omega_sim=omega_sim, beta_sim=beta_sim, big_T=big_T)
  logvar <- variance(n_n=n_n, n_r=n_r, n_p=n_p, n=n,
                     omega_sim=omega_sim, beta_sim=beta_sim, big_T=big_T)
  estvar <- var.log.to.var(est, logvar)
  return(list(est=est, var=estvar))
}
