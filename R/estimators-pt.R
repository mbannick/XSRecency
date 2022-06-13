#' Adjusted estimator  accounting for prior test results.
#'
#' @inheritParams get.adjusted.pt
#' @export
adjusted.estimate.pt <- function(n_r_pt, n_n, n_p,
                                 omega, beta, big_T,
                                 num_beta, den_omega, den_beta){

  numerator <- n_r_pt - beta * num_beta
  denominator <- n_n * (omega - beta * big_T + (den_omega + beta * den_beta)/n_p)

  val <- numerator / denominator
  return(val)
}

#' Adjusted variance computation, is the same for adjusted
#' and snapshot estimator, just pass in beta_sim=list(est=0, var=0).
#'
#' @export
variance <- function(n_n, n_r, n_p, n, omega, omega_var, beta, beta_var, big_T, q=1){

  # Slight modification for n_{p,test} to give the number of positives tested
  # for recency. Takes the place of n_p in this calculation.
  n_pt <- n_p * q

  variance <- (
    n_r * (n_pt - n_r) / (n_pt * (n_r - n_pt * beta)**2) +
      n / (n_p * n_n) +
      beta_var * n_pt * (n - n_pt) / (n * (n_r - n_pt * beta)**2) +
      omega_var / (omega - beta * big_T)**2 +
      beta_var * (
        (n_pt * omega - n_r * big_T) /
          ((n_r - n_pt * beta) * (omega - beta * big_T))
      ) ** 2
  )
  return(variance)
}

#' Convert log variance to variance
#'
#' @export
var.log.to.var <- function(estimate, variance) (estimate ** 2) * variance

#' Adjusted estimator and variance (Kassanjee et al. 2012)
#'
#'   \deqn{
#'     \hat{\lambda} = \frac{N_{rec}/q - N_{pos} \hat{\beta}_{T^*}}{N_{neg} (\hat{\Omega}_{T^*} - \hat{\beta}_{T^*} T^*)}
#'   }
#'   where \eqn{\hat{\Omega}_{T^*}} is an estimate of the mean duration of recent infection (MDRI),
#'   \eqn{\hat{\beta}_{T^*}} is an estimate of the false recency rate, and \eqn{T^*} is the time cutoff for being
#'   a recent infected or not.
#'
#'   \code{q} is an adjustment for the number of positives that are given recency tests.
#'
#' @export
#' @param n_r_pt Number of recent positives
#' @param n_n Number of negatives
#' @param n_p Number of positives
#' @param n Total number of observations (\eqn{n_n + n_p = n})
#' @param omega Mean duration of recent infection (MDRI) (in years, not days)
#' @param omega_var Variance of the estimator for MDRI (or 0 if MDRI known)
#' @param beta False recency rate (FRR)
#' @param beta_var Variance of the estimator for FRR (or 0 if FRR known)
#' @param big_T The \eqn{T^*} in the equation above
#' @param q The fraction of positives that were given recency tests, defaults to 1
#' @return Returns a list of the estimate and the variance.
#' @examples
#' get.adjusted(n_r=2, n_n=50, n_p=10, n=60,
#'              omega=0.36, omega_var=0, beta=0.02, beta_var=0, big_T=2, q=1)
#' get.adjusted(n_r=c(2, 3), n_n=c(50, 48), n_p=c(10, 12), n=c(60, 60),
#'              omega=0.36, omega_var=0, beta=0.02, beta_var=0, big_T=2)
get.adjusted.pt <- function(n_r_pt, n_n, n_p, n, omega, omega_var,
                         beta, beta_var, big_T, q=1,
                         num_beta, den_omega, den_beta){

  est <- adjusted.estimate.pt(n_r_pt=n_r_pt, n_n=n_n, n_p=n_p,
                              omega=omega, beta=beta, big_T=big_T,
                              num_beta=num_beta,
                              den_omega=den_omega, den_beta=den_beta)
  # logvar <- variance.pt(n_n=n_n, n_r=n_r, n_p=n_p, n=n,
  #                       omega=omega, omega_var=omega_var,
  #                       beta=beta, beta_var=beta_var,
  #                       big_T=big_T, q=q)
  # estvar <- var.log.to.var(est, logvar)
  return(list(est=est, var=NA))
}
