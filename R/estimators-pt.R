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

variance.pt <- function(
    n_n, n_r_pt, n_r, n_p, n,
    omega, omega_var,
    beta, beta_var, big_T,
    r_TA, r_TAprime, r_TAstar,
    omega_TA, omega_TA_var, omega_TAstar, omega_TA2,
    mu_TA, var_TA, mu_TB, var_TB,
    p_A, p_B,
    incidence=NULL
  ){

  p <- n_p / n
  pr <- n_r_pt / n_p # This is P_R^*
  pr_old <- n_r / n_p # This is P_R
  obt <- omega - beta * big_T

  # If incidence is not provided,
  # use the relationship with pr_old instead.
  # This is theoretically more robust to assumption violations.
  if(is.null(incidence)){
    lamp <- (pr_old - beta) / obt
  } else {
    lamp <- incidence * (1-p) / p
    # replace pr_old with the "assumption"-based pr.
    pr_old <- beta + obt * lamp
  }

  EW1 <- n * p * (pr - beta * (1-p_B))
  EW2 <- n * p
  EW3 <- obt
  EW4 <- n * p * p_A * (mu_TA  - omega_TA)
  EW5 <- n * p * beta * p_B * mu_TB

  VW1 <- n * p * (
    pr * (1-pr) + (1-p) * (pr - (1-p_B) * beta)**2 +
    beta_var * (1 - p_B) * (1-p * (1-p_B) + n * p * (1-p_B)) +
    beta * p_B * (beta * (1-p_B) - 2 * (pr - lamp*(
      omega + beta * (mu_TB - big_T)
    )))
  )
  VW2 <- n * p * (1-p)
  VW3 <- omega_var + beta_var * big_T**2
  VW4 <- n * p * p_A * (
    var_TA + omega_TA_var + mu_TA**2 + omega_TA**2 -
    p * p_A * (mu_TA - omega_TA)**2 +
    r_TA + p_A * r_TAprime * (n * p - p) - 2 * omega_TAstar
  )
  VW5 <- n * p * p_B * (
    beta_var * mu_TB**2 * p_B * n * p +
    (beta_var + beta**2) * (var_TB + mu_TB**2 * (1-p_B * p))
  )

  C12 <- n * p * (1-p) * (pr - (1-p_B) * beta)
  C13 <- n * p * beta_var * big_T * (1-p_B)
  C14 <- n * p * p_A * (
    (mu_TA - omega_TA) * (pr_old - p * pr - beta * (1 - p + p * p_B)) +
    lamp * (mu_TA**2 + var_TA + omega_TA**2 + omega_TA_var - 2*omega_TAstar)
  )
  C15 <- n * p * p_B * (
    beta * (
      lamp * (
        mu_TB * (omega - beta * big_T) +
        beta * (var_TB + mu_TB**2)
      ) - p * pr * mu_TB
    ) +
    p * (beta_var + beta**2) * mu_TB * (1 - p_B) -
    beta_var * n * p * (1 - p_B) * mu_TB
  )
  C23 <- rep(0, length(n))
  C24 <- n * p * (1-p) * p_A * (mu_TA - omega_TA)
  C25 <- n * p * (1-p) * beta * p_B * mu_TB
  C34 <- -n * p * p_A * r_TAstar
  C35 <- -n * p * big_T * p_B * mu_TB * beta_var
  C45 <- -n * p^2 * p_A * p_B * beta * mu_TB * (mu_TA - omega_TA)

  components_est <- list(
    EW1=EW1, EW2=EW2, EW3=EW3, EW4=EW4, EW5=EW5,
    VW1=VW1, VW2=VW2, VW3=VW3, VW4=VW4, VW5=VW5,
    C12=C12, C13=C13, C14=C14, C15=C15,
    C23=C23, C24=C24, C25=C25,
    C34=C34, C35=C35,
    C45=C45
  )

  varfunc <- function(ew1, ew2, ew3, ew4, ew5,
                      vw1, vw2, vw3, vw4, vw5,
                      c12, c13, c14, c15,
                      c23, c24, c25,
                      c34, c35,
                      c45,
                      N){

    d <- (ew4 + ew5) / ew2

    delta_gW <- c(
      1/ew1,
      1/(N-ew2) + d/(ew3 * ew2 + ew4 + ew5),
      -1/(ew3+d),
      -1/(ew2 * (ew3 + d)),
      -1/(ew2 * (ew3 + d))
    )

    cov_W <- matrix(
      c(vw1, c12, c13, c14, c15,
        c12, vw2, c23, c24, c25,
        c13, c23, vw3, c34, c35,
        c14, c24, c34, vw4, c45,
        c15, c25, c35, c45, vw5),
      byrow=T,
      nrow=5
    )

    logvar <- delta_gW %*% cov_W %*% delta_gW
    return(logvar)
  }

  logvars <- mapply(
    FUN=varfunc,
    ew1=EW1,
    ew2=EW2,
    ew3=EW3,
    ew4=EW4,
    ew5=EW5,
    vw1=VW1,
    vw2=VW2,
    vw3=VW3,
    vw4=VW4,
    vw5=VW5,
    c12=C12,
    c13=C13,
    c14=C14,
    c15=C15,
    c23=C23,
    c24=C24,
    c25=C25,
    c34=C34,
    c35=C35,
    c45=C45,
    N=n
  )

  return(list(
    logvars=logvars,
    components_est=components_est
  ))
}

#' Convert log variance to variance
#'
#' @export
var.log.to.var <- function(estimate, variance) (estimate ** 2) * variance

#' Enhanced estimator and variance (Bannick and Gao 2023+)
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
#' get.adjusted.pt(n_r_pt=3, n_r=2, n_n=50, n_p=10, n=60,
#'                 omega=0.36, omega_var=0, beta=0.02, beta_var=0, big_T=2, q=1)
#' get.adjusted.pt(n_r=c(2, 3), n_n=c(50, 48), n_p=c(10, 12), n=c(60, 60),
#'                 omega=0.36, omega_var=0, beta=0.02, beta_var=0, big_T=2)
get.adjusted.pt <- function(n_p, n, ptdf,
                            beta, beta_var, big_T,
                            phidat, use_geese, formula, family, plot_phi=TRUE, ...){

  # Summarize data inputs from the prior testing data
  # and an estimate of the phi function
  summdat <- summarize.pt.generator(bigT=big_T, use_geese=use_geese,
                                    formula=formula,
                                    family=family, plot_phi=plot_phi)
  args <- summdat(ptdf=ptdf, n=n, n_p=n_p, phidat=phidat)

  # Set up additional arguments
  args[["beta"]] <- beta
  args[["beta_var"]] <- beta_var
  args[["big_T"]] <- big_T

  # Run estimation function
  funcresult <- R.utils::doCall(get.adjusted.pt.internal, args=args)

  result <- list()
  result[["est"]] <- funcresult$est
  # Use the variance that is robust
  result[["var"]] <- funcresult$var[[1]]

  result[["omega"]] <- args$omega
  result[["omega_var"]] <- args$omega_var
  result[["n_r"]] <- args$n_r
  result[["n_n"]] <- args$n_n
  result[["n_r_pt"]] <- args$n_r_pt
  result[["n_p"]] <- args$n_p
  result[["n"]] <- args$n
  result[["q_eff"]] <- args$q_eff

  return(result)
}

get.adjusted.pt.internal <- function(n_r_pt, n_r, n_n, n_p, n, omega, omega_var,
                                     beta, beta_var, big_T, q=1,
                                     num_beta, den_omega, den_beta,
                                     r_TA, r_TAprime, r_TAstar,
                                     omega_TA, omega_TA_var, omega_TAstar, omega_TA2,
                                     mu_TA, var_TA, mu_TB, var_TB,
                                     p_A, p_B){

  est <- adjusted.estimate.pt(n_r_pt=n_r_pt, n_n=n_n, n_p=n_p,
                              omega=omega, beta=beta, big_T=big_T,
                              num_beta=num_beta,
                              den_omega=den_omega, den_beta=den_beta)

  var <- list()
  components_est <- list()
  for(i in 1:2){
    if(i == 1) inc.i <- NULL
    if(i == 2) inc.i <- est
    logvar <- variance.pt(n_n=n_n, n_r_pt=n_r_pt, n_r=n_r, n_p=n_p, n=n,
                          omega=omega, omega_var=omega_var,
                          beta=beta, beta_var=beta_var,
                          r_TA=r_TA, r_TAprime=r_TAprime, r_TAstar=r_TAstar,
                          omega_TA=omega_TA, omega_TA_var=omega_TA_var, omega_TA2=omega_TA2,
                          omega_TAstar=omega_TAstar,
                          mu_TA=mu_TA, var_TA=var_TA,
                          mu_TB=mu_TB, var_TB=var_TB,
                          p_A=p_A, p_B=p_B,
                          big_T=big_T,
                          incidence=inc.i)
    estvar <- var.log.to.var(est, logvar$logvar)
    var[[i]] <- estvar
    components_est[[i]] <- logvar$components_est
  }
  return(list(est=est, var=var, components_est=components_est))
}
